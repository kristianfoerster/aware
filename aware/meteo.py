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

        reference = munch.Munch()
        self.reference = reference
        npz = np.load(config.reference_mapping_file)
        reference.mapping_x = npz['mapping_x']
        reference.mapping_y = npz['mapping_y']

        # read histalp DTM
        # XXX fix this (flip u/d errors, ...)
        # histalp.surface = util.read_gdal_file(config.histalp_dtm_file)
        # @todo: rename histalp to reference_grid or similar
        nc = netCDF4.Dataset(config.reference_dtm_file, 'r')
        reference.surface = nc['HSURF'][:, :]

        # get reference climatology
        #np.savez(outfile, mean_temp=month_clim_T, std_temp=month_stdv_T, \
        #        mean_prec=month_clim_P, std_prec=month_stdv_P)
        npz = np.load(config.meteo_ref_climatology_file)
        self.ref_avg_temp = npz['mean_temp']
        self.ref_std_temp = npz['std_temp']
        self.ref_avg_prec = npz['mean_prec']
        self.ref_std_prec = npz['std_prec']
        
        self.cf = False

        if config.meteo_type == 'histalp':
            nc_temp = netCDF4.Dataset(config.histalp_temp_file, 'r')
            nc_precip = netCDF4.Dataset(config.histalp_precip_file, 'r')
            self.temp = nc_temp.variables['T_2M']
            self.precip = nc_precip.variables['TOT_PREC']
            reference.time_temp = util.num2date(nc_temp.variables['time'])
            reference.time_precip = util.num2date(nc_precip.variables['time'])
            self.dates = sorted(list(set(reference.time_temp) & set(reference.time_precip)))
            
            if 'meteo_climatological_forecast' in config:
                self.cf = config.meteo_climatological_forecast
                if self.cf:
                    print('Climatological forecast mode activated.')
        else:
            npz = np.load(config.meteo_mapping_file)
            self.mapping_x = npz['mapping_x']
            self.mapping_y = npz['mapping_y']
            
            # get forecast model climatology
            npz = np.load(config.meteo_mod_climatology_file)
            self.mod_avg_temp = npz['mean_temp']
            self.mod_std_temp = npz['std_temp']
            self.mod_avg_prec = npz['mean_prec']
            self.mod_std_prec = npz['std_prec']
            
            # old bias correction:
            #npz = np.load(config.meteo_bias_file)
            #self.temp_slope = npz['temp_slope']
            #self.temp_intercept = npz['temp_intercept']
            #self.precip_slope = npz['precip_slope']
            #self.precip_intercept = npz['precip_intercept']
            #self.temp_bias = npz['temp_bias']
            
            # check for optimum scale parameters
            #meteo_optimum_scale_active = Tue
            #meteo_optimum_scale_window = 1
            #meteo_optimum_scale_shift  = [0,0]
            if 'meteo_optimum_scale_active' in config:
                self.optimum_scale_activated = config.meteo_optimum_scale_active
                self.optimum_scale_window    = config.meteo_optimum_scale_window
                self.meteo_optimum_scale_shift = config.meteo_optimum_scale_shift
                print('Optimum scale analysis activated with window size = %i and shift %i in x and %i in y direction.' % \
                (self.optimum_scale_window, self.meteo_optimum_scale_shift[0],self.meteo_optimum_scale_shift[1]))
            else:
                self.optimum_scale_activated = False
            
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
        if len(prec_var.shape) > 2:
            # this means that a time dimension is defined:
            prec = prec_var[0,:, :] # precipitation files have a time dimension in seasonal experiment
        else:
            prec = prec_var[:, :] # precipitation files do not have a time dimension in horizontalResImpact experiment
        ncp.close()
        #return self._meteo_mapping(date, temp, prec, temp_const_bias = True)
        if self.optimum_scale_activated:
            return self._meteo_mapping_anom_scale(date, temp, prec, self.optimum_scale_window, self.meteo_optimum_scale_shift)
        else:
            return self._meteo_mapping_anom(date, temp, prec)

    def _get_meteo_histalp(self, date):
        date = pd.Timestamp(date)
        
        if self.cf:
            return self.ref_avg_temp[date.month-1], self.ref_avg_prec[date.month-1]

        histalp_temp_pos = np.where(pd.to_datetime(self.reference.time_temp) == date)[0][0]
        histalp_precip_pos = np.where(pd.to_datetime(self.reference.time_precip) == date)[0][0]
        temp_histalp = self.temp[histalp_temp_pos, :, :]
        precip_histalp = self.precip[histalp_precip_pos, :, :]

        return temp_histalp, precip_histalp

    def _get_meteo_cfs(self, date):
        date = pd.Timestamp(date)

        pos = np.where(pd.to_datetime(self.dates) == date)[0][0]

        temp_cfs = self.temp[pos, :, :]
        precip_cfs = self.precip[pos, :, :]
        #return self._meteo_mapping(date, temp_cfs, precip_cfs)
        #return self._meteo_mapping_anom(date, temp_cfs, precip_cfs)
        if self.optimum_scale_activated:
            return self._meteo_mapping_anom_scale(date, temp_cfs, precip_cfs, self.optimum_scale_window, self.meteo_optimum_scale_shift)
        else:
            return self._meteo_mapping_anom(date, temp_cfs, precip_cfs)


    # new!
    def _meteo_mapping_anom(self, date, temp, prec):
        # forecast model to reference
        mx = self.mapping_x
        my = self.mapping_y

        mi = date.month - 1
        num_days_per_month = calendar.monthrange(date.year, date.month)[1]
        
        # transform standardized anomalies
        temp_reference = (temp[my, mx] - self.mod_avg_temp[mi, my, mx]) / self.mod_std_temp[mi, my, mx] * self.ref_std_temp[mi,:,:] + self.ref_avg_temp[mi,:,:]
        prec_reference = (prec[my, mx] * num_days_per_month * 86400. - self.mod_avg_prec[mi, my, mx]) / self.mod_std_prec[mi, my, mx] * self.ref_std_prec[mi,:,:] + self.ref_avg_prec[mi,:,:]

        return temp_reference, prec_reference

    def _meteo_mapping_anom_scale(self, date, temp, prec, window_size=1, shift=[0,0]):

        temp_reference = np.zeros(self.ref_avg_temp[0,:,:].shape)
        prec_reference = np.zeros(self.ref_avg_prec[0,:,:].shape)

        start = -int(window_size/2)
        stop  = start + window_size

        num_dataframes = 0
        for si in range(start,stop):
            for sj in range(start,stop):

                # forecast model to reference
                # @todo: check boundaries!!!
                mx = self.mapping_x + si + shift[0]
                my = self.mapping_y + sj + shift[1]
        
                mi = date.month - 1
                num_days_per_month = calendar.monthrange(date.year, date.month)[1]
                
                # transform standardized anomalies
                temp_reference += (temp[my, mx] - self.mod_avg_temp[mi, my, mx]) / self.mod_std_temp[mi, my, mx] * self.ref_std_temp[mi,:,:] + self.ref_avg_temp[mi,:,:]
                prec_reference += (prec[my, mx] * num_days_per_month * 86400. - self.mod_avg_prec[mi, my, mx]) / self.mod_std_prec[mi, my, mx] * self.ref_std_prec[mi,:,:] + self.ref_avg_prec[mi,:,:]
                num_dataframes += 1
        
        temp_reference /= num_dataframes
        prec_reference /= num_dataframes        
        
        return temp_reference, prec_reference
        

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

    def _reference_to_model_grid(self, date, temp_reference, precip_reference):
        temp_lapse_rates = np.zeros(self.dtm.shape) * np.nan
        precip_lapse_rates = np.zeros(self.dtm.shape) * np.nan
        for cid, cpx in self.catchment_pixels.items():
            temp_lapse_rates[cpx] = self.config.params.catchments[cid].temp_lapse_rates[date.month - 1]
            precip_lapse_rates[cpx] = 0.01 * self.config.params.catchments[cid].precip_lapse_rates[date.month - 1]

        mx = self.reference.mapping_x
        my = self.reference.mapping_y

        temp = temp_reference[my, mx] + 273.15 - (self.reference.surface[my, mx] - self.dtm) * temp_lapse_rates
        precip = precip_reference[my, mx] - (self.reference.surface[my, mx] - self.dtm) * precip_reference[my, mx] * precip_lapse_rates
        precip = precip.clip(0.)

        return temp, precip

    def get_meteo(self, date):
        if self.config.meteo_type == 'histalp':
            temp, precip = self._get_meteo_histalp(date)
        elif self.config.meteo_type == 'cfs':
            temp, precip = self._get_meteo_cfs(date)
        elif self.config.meteo_type == 'glosea5':
            temp, precip = self._get_meteo_glosea5(date)

        return self._reference_to_model_grid(date, temp, precip)
