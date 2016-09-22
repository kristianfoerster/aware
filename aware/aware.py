from __future__ import print_function, division, absolute_import
import aware
import aware.models
import collections
import munch
import numpy as np
import pandas as pd
from . import util
from . import awarestatevariable


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
        dtm, prj_settings = util.read_gdal_file(self.config.dtm_file, return_prj_settings=True)
        glaciers = aware.util.read_gdal_file(self.config.glaciers_file)
        catchments = aware.util.read_gdal_file(self.config.catchments_file, fill_value=0)

        self.input_grids = munch.Munch(dtm=dtm, glaciers=glaciers, catchments=catchments)
        self.catchment_ids = np.unique(catchments[catchments > 0])

        # define statevariables (glacier reserved for future work)
        self.state_swe          = awarestatevariable.AwareStateVariable('Snow water equivalent SWE', prj_settings, None, None, catchments>0)
        self.state_icewe        = awarestatevariable.AwareStateVariable('Ice water equivalent',      prj_settings, None, None, catchments>0)
        self.state_glacierarea  = awarestatevariable.AwareStateVariable('Glaciated area',            prj_settings, None, None, catchments>0)
        self.state_soilmoisture = awarestatevariable.AwareStateVariable('Soil moisture storage',     prj_settings, None, None, catchments>0)
        self.state_groundwater  = awarestatevariable.AwareStateVariable('Groundwater storage',       prj_settings, None, None, catchments>0)
      
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
            if center_lon_lat[0] < - 180 or center_lon_lat[0] > 180 or center_lon_lat[1] < -90 or center_lon_lat[1] > 90:
                print('Error: gdal might not have correctly detected latitude and longitude!')
                print('Try to use values from config file...')
                if 'latitude' in params:
                    latitude = params['latitude']
                    print('Use value from config file: %f' % latitude)
                else:
                    latitude = 45
                    print('Failed: No valid values found! Using DEFAULT latitude: %f' % latitude)
                    print('Please check configuration since this parameter is crucial for ET computations!')
            else:
                latitude = center_lon_lat[1]

            evapotranspiration = aware.models.EvapotranspirationModel(latitude, J=params.et_j)
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

            cpx = self.catchments[cid].pixels
            self.state_soilmoisture.set_state(params.sms_init, cpx)
            self.state_swe.set_state(0., cpx)
            self.state_icewe.set_state(0., cpx)
            self.state_glacierarea.set_state(glaciers)
            self.state_groundwater.set_state(params.gw_storage_init, cpx)

            # this is a crude assumption which needs to be reconsidered as soon
            # as the glacier model is coupled!
            glaciers_pos = cpx & (glaciers > 0)
            self.state_icewe.set_state(1e10, glaciers_pos)

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

        rts_catchments = collections.OrderedDict()
        rts = pd.DataFrame(index=dates, columns=self._results_time_series_columns, dtype=float)
        for cid in self.catchment_ids:
            rts_catchments[cid] = rts.copy()

        for date in dates:
            print(date)

            temp, precip = self.meteo.get_meteo(date)
            glaciers = self.state_glacierarea.get_state() # couple glacier model here!

            for cid in self.catchment_ids:
                params = self.config.params.catchments[cid]
                catchment = self.catchments[cid]
                cpx = self.catchments[cid].pixels

                cswe, snowmelt, snow_outflow, snowfall, rainfall, melt_avail = catchment.melt.melt(
                    self.state_swe.get_state(cpx),
                    precip[cpx],
                    temp[cpx],
                    params.ddf_snow,
                    glacier_fraction=None
                )
                self.state_swe.set_state(cswe,cpx)

                snow_outflow_unglacierized = snow_outflow * (1. - glaciers[cpx])
                snow_outflow_glacierized = np.zeros(self.input_grids.dtm.shape)
                snow_outflow_glacierized = snow_outflow * glaciers[cpx]

                ice_melt_factor = np.minimum(glaciers[cpx], glaciers[cpx] * melt_avail * params.ddf_ice / params.ddf_snow)
                cicewe, icemelt, ice_outflow, _, _, _ = catchment.melt.melt(
                    self.state_icewe.get_state(cpx),
                    precip[cpx] * 0.0,
                    temp[cpx],
                    params.ddf_ice,
                    glacier_fraction=ice_melt_factor
                )
                self.state_icewe.set_state(cicewe,cpx)
                
                glacier_outflow = ice_outflow + snow_outflow_glacierized

                # get groundwater state
                gw_storage = self.state_groundwater.get_state(cpx)

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
                        self.state_soilmoisture.get_state(cpx)
                    )
                    # sms[cpx] = csms
                    self.state_soilmoisture.set_state(csms,cpx)
                    percolation = perc.mean()
                    direct_runoff = runoff_d.mean()

                    baseflow, gw_storage = catchment.groundwater.groundwater_model(gw_storage, perc)
                    runoff = (baseflow + direct_runoff) + glacier_outflow
                else:
                    melt_only = snow_outflow_glacierized + glacier_outflow
                    baseflow, gw_storage = catchment.groundwater.groundwater_model(gw_storage, params.gw_n * melt_only)
                    runoff = baseflow + (1 - params.gw_n) * melt_only

                    et_pot = np.zeros(cswe.shape) * np.nan
                    et_act = np.zeros(cswe.shape) * np.nan
                    direct_runoff = np.nan
                
                self.state_groundwater.set_state(gw_storage,cpx)                
                
                rts_cur = rts_catchments[cid].loc[date]
                rts_cur.temp = temp[cpx].mean()
                rts_cur.precip = precip[cpx].mean()
                rts_cur.snowfall = snowfall.mean()
                rts_cur.rainfall = rainfall.mean()
                rts_cur.swe = self.state_swe.get_state(cpx).mean()
                rts_cur.snowmelt = snowmelt.mean()
                rts_cur.icemelt = icemelt.mean()
                rts_cur.melt = rts_cur.snowmelt + rts_cur.icemelt
                rts_cur.snow_outflow = snow_outflow.mean()
                rts_cur.ice_outflow = ice_outflow.mean()
                rts_cur.glacier_outflow = glacier_outflow.mean()
                rts_cur.runoff = runoff.mean() # AVERAGE INCLUDING ALL TRIBUTARIES!!!!!!!!!!!!!!!!!!!
                rts_cur.sms = self.state_soilmoisture.get_state(cpx).mean()
                rts_cur.et_pot = et_pot.mean()
                rts_cur.et = et_act.mean()
                rts_cur.baseflow = baseflow.mean()
                rts_cur.direct_runoff = direct_runoff
                # rts_cur.icewe = self.state_icewe.get_state(cpx).mean() # activate if required

        results = munch.Munch()
        results.ts = pd.Panel(rts_catchments)

        return results
