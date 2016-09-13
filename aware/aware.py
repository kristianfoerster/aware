from __future__ import print_function, division, absolute_import
import aware
import aware.models
import collections
import munch
import numpy as np
import pandas as pd
from . import util

import matplotlib
import matplotlib.pyplot as plt

class Aware(object):
    _results_time_series_columns = '''
        temp
        precip
        snowfall
        rainfall
        swe
        snowmelt
        icemelt
        melt
        snow_outflow
        ice_outflow
        glacier_outflow
        runoff
        sms
        et_pot
        et
        baseflow
        direct_runoff
    '''.strip().split()

    @classmethod
    def run_aware(cls, config):
        aw = cls()
        aw.config = config
        aw.initialize()
        return aw.run()

    def __init__(self):
        self.config = aware.Config()

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, c):
        assert isinstance(c, aware.Config)
        c = munch.munchify(c)
        self._config = c

    def initialize(self):
        dtm = util.read_gdal_file(self.config.dtm_file)
        glaciers = aware.util.read_gdal_file(self.config.glaciers_file)
        catchments = aware.util.read_gdal_file(self.config.catchments_file, fill_value=0)

        self.input_grids = munch.Munch(dtm=dtm, glaciers=glaciers, catchments=catchments)
        self.catchment_ids = np.unique(catchments[catchments > 0])

        default_params = self.config.params.default
        catchment_params = {cid: default_params.copy() for cid in self.catchment_ids}
        if 'catchments' in self.config.params:
            for cid in self.config.params.catchments.keys():
                catchment_params[cid] = default_params.items() + self.config.params.catchments[cid].items()
        self.config.params.catchments = munch.munchify(catchment_params)

        self.catchments = munch.Munch()
        for cid in self.catchment_ids:
            params = self.config.params.catchments[cid]

            soil = aware.models.SoilModel()

            groundwater = aware.models.GroundwaterModel()
            groundwater.K = params.gw_k

            center_lon_lat = util.center_lon_lat(self.config.dtm_file) # TODO calculate this for each catchment
            evapotranspiration = aware.models.EvapotranspirationModel(center_lon_lat[1], J=params.et_j)
            evapotranspiration.delta_daylength = 0
            soil.SMSC  = params.soil_capacity
            soil.DRf   = params.psi
            soil.Rf    = params.percolation_r
            soil.PercF = params.percolation_f

            melt = aware.models.MeltModel()
            melt.melt_temp = params.melt_temp
            melt.trans_temp = params.trans_temp
            melt.temp_trans_range = params.temp_trans_range
            melt.rain_corr = params.rain_corr
            melt.snow_corr = params.snow_corr

            self.catchments[cid] = munch.Munch()
            self.catchments[cid].soil = soil
            self.catchments[cid].groundwater = groundwater
            self.catchments[cid].evapotranspiration = evapotranspiration
            self.catchments[cid].melt = melt
            self.catchments[cid].pixels = (self.input_grids.catchments == cid)

        self.meteo = aware.Meteo(
            self.config,
            dtm,
            {cid: self.catchments[cid].pixels for cid in self.catchment_ids}
        )

    def run(self):
        start_date = pd.Timestamp(self.config.start_date).to_period('M').to_timestamp('M')
        end_date = pd.Timestamp(self.config.end_date).to_period('M').to_timestamp('M')
        assert start_date in self.meteo.dates
        assert end_date in self.meteo.dates

        dates = pd.date_range(start=start_date, end=end_date, freq='M')

        sms = np.zeros(self.input_grids.dtm.shape)
        swe = np.zeros(self.input_grids.dtm.shape)
        icewe = np.zeros(self.input_grids.dtm.shape)
        glaciers = self.input_grids.glaciers
        glaciers_pos = glaciers > 0
        icewe[glaciers_pos] = 1e10

        rts_catchments = collections.OrderedDict()

        for cid in self.catchment_ids:
            params = self.config.params.catchments[cid]
            catchment = self.catchments[cid]
            cpx = self.catchments[cid].pixels

            rts = pd.DataFrame(index=dates, columns=self._results_time_series_columns, dtype=float)

            sms[cpx] += params.sms_init
            gw_storage = params.gw_storage_init

            for date in dates:
                print(date)

                temp, precip = self.meteo.get_meteo(date) # TODO get meteo only for catchment pixels

                cswe, snowmelt, snow_outflow, snowfall, rainfall, melt_avail = catchment.melt.melt(
                    swe[cpx],
                    precip[cpx],
                    temp[cpx],
                    params.ddf_snow,
                    glacier_fraction=None
                )
                swe[cpx] = cswe

                snow_outflow_unglacierized = snow_outflow * (1. - glaciers[cpx])
                snow_outflow_glacierized = np.zeros(self.input_grids.dtm.shape)
                snow_outflow_glacierized = snow_outflow * glaciers[cpx]

                ice_melt_factor = np.minimum(glaciers[cpx], glaciers[cpx] * melt_avail * params.ddf_ice / params.ddf_snow)
                cicewe, icemelt, ice_outflow, _, _, _ = catchment.melt.melt(
                    icewe[cpx],
                    precip[cpx] * 0.0,
                    temp[cpx],
                    params.ddf_ice,
                    glacier_fraction=ice_melt_factor
                )

                glacier_outflow = ice_outflow.mean() + snow_outflow_glacierized.mean()

                if self.config.enable_soil_model:
                    et_pot = catchment.evapotranspiration.monthly_evapotranspiration(
                        date,
                        temp[cpx],
                        n_etp_summer=params.factor_etp_summer
                    )
                    et_pot *= (1. - glaciers[cpx])
            
                    csms, perc, runoff_d, et_act = catchment.soil.soil_water_balance(
                        snow_outflow_unglacierized,
                        et_pot,
                        sms[cpx]
                    )
                    sms[cpx] = csms
                    percolation = perc.mean()
                    direct_runoff = runoff_d.mean()

                    baseflow, gw_storage = catchment.groundwater.groundwater_model(gw_storage, percolation)
                    runoff = (baseflow + direct_runoff) + glacier_outflow
                else:
                    melt_only = snow_outflow_glacierized.mean() + glacier_outflow
                    baseflow, gw_storage = catchment.groundwater.groundwater_model(gw_storage, params.gw_n * melt_only)
                    runoff = baseflow + (1 - params.gw_n) * melt_only

                    et_pot = np.zeros(cswe.shape) * np.nan
                    et_act = np.zeros(cswe.shape) * np.nan
                    direct_runoff = np.nan

                rts_cur = rts.loc[date]
                rts_cur.temp = temp[cpx].mean()
                rts_cur.precip = precip[cpx].mean()
                rts_cur.snowfall = snowfall.mean()
                rts_cur.rainfall = rainfall.mean()
                rts_cur.swe = swe[cpx].mean()
                rts_cur.snowmelt = snowmelt.mean()
                rts_cur.icemelt = icemelt.mean()
                rts_cur.melt = rts_cur.snowmelt + rts_cur.icemelt
                rts_cur.snow_outflow = snow_outflow.mean()
                rts_cur.ice_outflow = ice_outflow.mean()
                rts_cur.glacier_outflow = glacier_outflow
                rts_cur.runoff = runoff
                rts_cur.sms = sms[cpx].mean()
                rts_cur.et_pot = et_pot.mean()
                rts_cur.et = et_act.mean()
                rts_cur.baseflow = baseflow
                rts_cur.direct_runoff = direct_runoff

            rts_catchments[cid] = rts

        results = munch.Munch()
        results.ts = pd.Panel(rts_catchments)

        return results
