from __future__ import print_function, division, absolute_import
import calendar
import munch
import netCDF4
import numpy as np
import pandas as pd
from . import util
import os

class Meteo(object):
    def __init__(self, config, dtm, catchment_pixels):
        self.config = config
        self.dtm = dtm
        self.catchment_pixels = catchment_pixels

        histalp = munch.Munch()
        self.histalp = histalp
        npz = np.load(config.histalp_mapping_file)
        histalp.mapping_x = npz['mapping_x']
        histalp.mapping_y = npz['mapping_y']

        # read histalp DTM
        # XXX fix this (flip u/d errors, ...)
        # histalp.surface = util.read_gdal_file(config.histalp_dtm_file)
        # @todo: rename histalp to reference_grid or similar
        nc = netCDF4.Dataset(config.histalp_dtm_file, 'r')
        histalp.surface = nc['HSURF'][:, :]

        if config.meteo_type == 'histalp':
            nc_temp = netCDF4.Dataset(config.histalp_temp_file, 'r')
            nc_precip = netCDF4.Dataset(config.histalp_precip_file, 'r')
            self.temp = nc_temp.variables['T_2M']
            self.precip = nc_precip.variables['TOT_PREC']
            histalp.time_temp = util.num2date(nc_temp.variables['time'])
            histalp.time_precip = util.num2date(nc_precip.variables['time'])
            self.dates = sorted(list(set(histalp.time_temp) & set(histalp.time_precip)))
        else:
            npz = np.load(config.meteo_mapping_file)
            self.mapping_x = npz['mapping_x']
            self.mapping_y = npz['mapping_y']

            npz = np.load(config.meteo_bias_file)
            self.temp_slope = npz['temp_slope']
            self.temp_intercept = npz['temp_intercept']
            self.precip_slope = npz['precip_slope']
            self.precip_intercept = npz['precip_intercept']
            self.temp_bias = npz['temp_bias']

            if config.meteo_type == 'cfs':
                self.nc = netCDF4.Dataset(config.meteo_file, 'r')
                self.temp = self.nc.variables['TMP_2maboveground']
                self.precip = self.nc.variables['PRATE_surface']
                time_var = self.nc.variables['time']
                self.dates = util.num2date(time_var)

                # set hour to 0 for all dates
                for date_num, date in enumerate(self.dates):
                    date_new = pd.datetime(date.year, date.month, date.day, 0)
                    self.dates[date_num] = date_new            
            elif config.meteo_type == 'glosea5':
                # in contrast to other data sources, each time step is stored separately
                # therefore, config.meteo_file is the file name without datetime and
                # variable information, e.g.:
                # Amon_GloSea5_horizlResImpact_S19960425_r1i1p1
                self.nc = config.meteo_file
                self.ncpath = config.meteo_path
                self.temp_file_identifier = 'tas'
                self.precip_file_identifier = 'pr'
                # at this stage, only start and end dates are known
                start_date = pd.Timestamp(self.config.start_date).to_period('M').to_timestamp('M')
                end_date   = pd.Timestamp(self.config.end_date).to_period('M').to_timestamp('M')               
                self.dates = pd.date_range(start=start_date, end=end_date, freq='M')
            else:
                raise ValueError('Meteo type %s not supported' % config.meteo_type)

    def _get_meteo_glosea5(self, date):
        date = pd.Timestamp(date)
        yyyy = date.year
        mm = date.month
        string = self.nc
        filename_t = 'tas_%s_%4i%02i-%4i%02i.nc' % (string, yyyy, mm, yyyy, mm)
        filename_p = 'pr_%s_%4i%02i-%4i%02i.nc' % (string, yyyy, mm, yyyy, mm)
        path_t = os.path.join(self.ncpath, filename_t)
        path_p = os.path.join(self.ncpath, filename_p)
        
        # read temperature
        nct = netCDF4.Dataset(path_t, 'r')
        temp_var = nct.variables['tas']
        temp = temp_var[0, :, :]
        nct.close()
        
        # read precipitation
        ncp = netCDF4.Dataset(path_p, 'r')
        prec_var = ncp.variables['pr']
        prec = prec_var[:, :] # precipitation files do not have a time dimension
        ncp.close()
        return self._meteo_mapping(date, temp, prec, temp_const_bias = True)

    def _get_meteo_histalp(self, date):
        date = pd.Timestamp(date)

        histalp_temp_pos = np.where(pd.to_datetime(self.histalp.time_temp) == date)[0][0]
        histalp_precip_pos = np.where(pd.to_datetime(self.histalp.time_precip) == date)[0][0]
        temp_histalp = self.temp[histalp_temp_pos, :, :]
        precip_histalp = self.precip[histalp_precip_pos, :, :]

        return temp_histalp, precip_histalp

    def _get_meteo_cfs(self, date):
        date = pd.Timestamp(date)

        pos = np.where(pd.to_datetime(self.dates) == date)[0][0]

        temp_cfs = self.temp[pos, :, :]
        precip_cfs = self.precip[pos, :, :]
        return self._meteo_mapping(date, temp_cfs, precip_cfs)


    def _meteo_mapping(self, date, temp, precip, temp_const_bias = False):
        mx = self.mapping_x
        my = self.mapping_y
        temp_slope = self.temp_slope[date.month - 1, :, :]
        temp_intercept = self.temp_intercept[date.month - 1, :, :]
        precip_slope = self.precip_slope[date.month - 1, :, :]
        precip_intercept = self.precip_intercept[date.month - 1, :, :]
        temp_bias = self.temp_bias[date.month - 1, :, :]
        num_days_per_month = calendar.monthrange(date.year, date.month)[1]

        if temp_const_bias:
            # print('Bias added.')
            temp_reference = temp[my, mx] - temp_bias - 273.15
        else:
            temp_reference = temp_intercept + temp_slope * temp[my, mx] - 273.15
        
        precip_reference = precip_intercept + precip_slope * precip[my, mx] * num_days_per_month * 86400.
        precip_reference = precip_reference.clip(0.)

        return temp_reference, precip_reference

    def _histalp_to_model_grid(self, date, temp_histalp, precip_histalp):
        temp_lapse_rates = np.zeros(self.dtm.shape) * np.nan
        precip_lapse_rates = np.zeros(self.dtm.shape) * np.nan
        for cid, cpx in self.catchment_pixels.items():
            temp_lapse_rates[cpx] = self.config.params.catchments[cid].temp_lapse_rates[date.month - 1]
            precip_lapse_rates[cpx] = 0.01 * self.config.params.catchments[cid].precip_lapse_rates[date.month - 1]

        mx = self.histalp.mapping_x
        my = self.histalp.mapping_y

        temp = temp_histalp[my, mx] + 273.15 - (self.histalp.surface[my, mx] - self.dtm) * temp_lapse_rates
        precip = precip_histalp[my, mx] - (self.histalp.surface[my, mx] - self.dtm) * precip_histalp[my, mx] * precip_lapse_rates
        precip = precip.clip(0.)

        return temp, precip

    def get_meteo(self, date):
        if self.config.meteo_type == 'histalp':
            temp, precip = self._get_meteo_histalp(date)
        elif self.config.meteo_type == 'cfs':
            temp, precip = self._get_meteo_cfs(date)
        elif self.config.meteo_type == 'glosea5':
            temp, precip = self._get_meteo_glosea5(date)

        return self._histalp_to_model_grid(date, temp, precip)
