from __future__ import print_function, division, absolute_import
import aware
import aware.models
import collections
import munch
import numpy as np
import pandas as pd
from . import util
from . import awarestatevariable
import copy
import os
import dateutil
import pickle

class Aware(object):
    '''This class represents the Alpine WAter balance and Runoff Estimation model
    (AWARE) and its core functionality.
    '''    
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
        # empty variable reserved for intermediate results (required for hotstart capability)
        self.state_swe           = None
        self.state_icewe         = None
        self.state_glacierarea   = None
        self.state_soilmoisture  = None
        self.state_groundwater   = None

    @property
    def config(self):
        return self._config

    @config.setter
    def config(self, c):
        assert isinstance(c, aware.Config)
        c = munch.munchify(c)
        self._config = c

    def initialize(self):
        '''Calling this function prior to the application of the model is mandatory.
        It prepares the model structure (computation of sub-basins) and sets all
        relevant parameters needed for simulations.
        '''
        dtm, prj_settings = util.read_gdal_file(self.config.dtm_file, return_prj_settings=True)
        self.glaciers_init = aware.util.read_gdal_file(self.config.glaciers_file)
        catchments = aware.util.read_gdal_file(self.config.catchments_file, fill_value=0)
        self.hydrographictree = pickle.load(open(self.config.hydrographictree_file,'rb'))

        self.input_grids = munch.Munch(dtm=dtm, glaciers=self.glaciers_init, catchments=catchments)
        self.catchment_ids = np.unique(catchments[catchments > 0])

        # define statevariables (glacier reserved for future work)
        self.state_swe          = awarestatevariable.AwareStateVariable('swe', 'Snow water equivalent SWE', prj_settings, None, None, catchments>0)
        self.state_icewe        = awarestatevariable.AwareStateVariable('iwe', 'Ice water equivalent',      prj_settings, None, None, catchments>0)
        self.state_glacierarea  = awarestatevariable.AwareStateVariable('glc', 'Glaciated area',            prj_settings, None, None, catchments>0)
        self.state_soilmoisture = awarestatevariable.AwareStateVariable('sms', 'Soil moisture storage',     prj_settings, None, None, catchments>0)
        self.state_groundwater  = awarestatevariable.AwareStateVariable('gws', 'Groundwater storage',       prj_settings, None, None, catchments>0)

        # remember projection settings
        # self.prj_settings = prj_settings

        # sub-catchments
        self.catchments = munch.Munch()
        for cid in self.catchment_ids:
            # catchment parameters
            self.catchments[cid] = munch.Munch()

            # hydrographic tree and structure of the basin
            for subtree_id, ti in enumerate(self.hydrographictree):
                tree = ti.get_tree_by_id(cid)
                if tree is not None:
                    break
            # name of catchment
            self.catchments[cid].name = tree.name

            # get tributaries
            ids, areas = tree.get_upstream_areas()
            self.catchments[cid].upstream_ids = np.array(ids)
            self.catchments[cid].upstream_areas = np.array(areas)
            # relative contribution of tributaries
            self.catchments[cid].upstream_areas = self.catchments[cid].upstream_areas / tree.area
            # own contribution
            self.catchments[cid].area = 1. - np.sum(self.catchments[cid].upstream_areas)
            # computation order
            self.catchments[cid].priority = tree.order
            # absolute area
            self.catchments[cid].abs_area = tree.area
            
            self.catchments[cid].next = tree.downstream_node
                                       
            # downstream path to outlet (for assigment of parameters)
            self.catchments[cid].downstream_path = self.hydrographictree[subtree_id].get_downstream_path(cid)
            
            a = np.array([self.catchments[cid].area])
            a = np.append(a, self.catchments[cid].upstream_areas)
            #print(tree.area)
            #print('%i\t' % (cid), (a*tree.area).astype(int), end='')
            #print('\tSumme=%.3f' % np.sum(a))


        # rearange computation order of catchments 
        # dictionary including ids and priorities
        p = {id: self.catchments[id].priority for id in self.catchments.keys()}
        # ordered by priority
        self.computation_order = sorted(p, key=lambda f: p[f])
        # reverse order
        reverse_order = self.computation_order[::-1]

      
        default_params = self.config.params.default
        catchment_params = {cid: default_params.copy() for cid in self.catchment_ids}
        # old method: assign catchment parameters if available. 
        #if 'catchments' in self.config.params:
        #    for cid in self.config.params.catchments.keys():
        #        catchment_params[cid] = default_params.items() + self.config.params.catchments[cid].items()

        # new method: get information from downstream areas as well
        if 'catchments' in self.config.params:
            for cid in reverse_order:
                # assign downstream parameters if relevant                
                for dwid in self.catchments[cid].downstream_path[::-1]:
                    if dwid in self.config.params.catchments.keys():
                        print('[%02i] apply parameters from area (%02i)...' % (cid,dwid))
                        # catchment_params[cid] += self.config.params.catchments[dwid].items()
                        for dii in default_params.keys():
                            if dii in self.config.params.catchments[dwid].keys():
                                param_value = self.config.params.catchments[dwid][dii]

                                if dii not in catchment_params[cid] or catchment_params[cid][dii] != param_value:
                                    catchment_params[cid][dii] = param_value
                                    print('    ... overwrite %s' % dii)

        self.config.params.catchments = munch.munchify(catchment_params)

        
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
        self.write_dates = list()
        if 'write_dates' in self.config:
            for di in self.config.write_dates:
                self.write_dates.append(pd.Timestamp(di).to_period('M').to_timestamp('M'))
        self.is_initialized = True
        self.is_ready       = False

    def reset_state_vars(self):
        '''Resets state variables as prescribed in the config file.
        
        '''       
        if not self.is_initialized:
            return False
        for cid in self.catchment_ids:
            params = self.config.params.catchments[cid]
            cpx = self.catchments[cid].pixels
            self.state_soilmoisture.set_state(params.sms_init, cpx)
            self.state_swe.set_state(0., cpx)
            self.state_icewe.set_state(0., cpx)
            self.state_glacierarea.set_state(self.glaciers_init)
            self.state_groundwater.set_state(params.gw_storage_init, cpx)

            # this is a crude assumption which needs to be reconsidered as soon
            # as the glacier model is coupled!
            glaciers_pos = cpx & (self.glaciers_init > 0)
            self.state_icewe.set_state(1e10, glaciers_pos)
        self.is_ready = True
        self.timestamp = pd.Timestamp(self.config.start_date).to_period('M').to_timestamp('M')
        return True

    def run(self, hotstart=False):
        '''Performs an AWARE simulation and returns the results (time series).
        
        Returns
        ----
        pandas data frame including the results
        '''
        if not self.is_initialized and not self.is_ready:
            print('Error: Model has not been initialized or prepared with initial states.')
            return
        if hotstart:
            if self.is_ready:
                print('AWARE hotstart: Resuming last run ...')
            else:
                print('Waring: cannot resume run in hotstart mode. Using default initialisation!')
        else:
            self.reset_state_vars()
        
        start_date = pd.Timestamp(self.config.start_date).to_period('M').to_timestamp('M')
        end_date = pd.Timestamp(self.config.end_date).to_period('M').to_timestamp('M')
        assert start_date in self.meteo.dates
        assert end_date in self.meteo.dates

        dates = pd.date_range(start=start_date, end=end_date, freq='M')

        rts_catchments = collections.OrderedDict() # results including upstream areas
        rts_catchments_sub_mean = collections.OrderedDict() # results sub-catchment only
        rts = pd.DataFrame(index=dates, columns=self._results_time_series_columns, dtype=float)
        for cid in self.catchment_ids:
            rts_catchments[cid] = rts.copy()
            rts_catchments_sub_mean[cid] = rts.copy()

        for date in dates:
            print(date)

            temp, precip = self.meteo.get_meteo(date)
            glaciers = self.state_glacierarea.get_state() # couple glacier model here!

            for cid in self.computation_order:
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
                    #percolation = perc.mean()
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
                
                
                rts_cur = rts_catchments_sub_mean[cid].loc[date]
                
                # calculate averages for sub-catchment without tributaries
                rts_cur.temp     = temp[cpx].mean()
                rts_cur.precip   = precip[cpx].mean()
                rts_cur.snowfall = snowfall.mean()
                rts_cur.rainfall = rainfall.mean()
                rts_cur.swe      = self.state_swe.get_state(cpx).mean()
                rts_cur.snowmelt = snowmelt.mean()
                rts_cur.icemelt  = icemelt.mean()
                rts_cur.melt     = rts_cur.snowmelt + rts_cur.icemelt
                rts_cur.snow_outflow = snow_outflow.mean()
                rts_cur.ice_outflow = ice_outflow.mean()
                rts_cur.glacier_outflow = glacier_outflow.mean()
                rts_cur.runoff   = runoff.mean()
                rts_cur.sms      = self.state_soilmoisture.get_state(cpx).mean()
                rts_cur.et_pot   = et_pot.mean()
                rts_cur.et       = et_act.mean()
                rts_cur.baseflow = baseflow.mean()
                rts_cur.direct_runoff = direct_runoff
                # rts_cur.icewe = self.state_icewe.get_state(cpx).mean() # activate if required
                
                # prepare averages
                rts_catchments[cid].loc[date] = self.catchments[cid].area * rts_catchments_sub_mean[cid].loc[date]                
                
                # add results of tributaries
                for ii in range(0,len(self.catchments[cid].upstream_ids)):
                    sub_id = self.catchments[cid].upstream_ids[ii]
                    sub_n  = self.catchments[cid].upstream_areas[ii]
                    
                    # tributaries
                    rts_cur_sub = rts_catchments_sub_mean[sub_id].loc[date]
                    # set results to total upstream area
                    rts_cur = rts_catchments[cid].loc[date]

                    rts_cur.temp     += sub_n * rts_cur_sub.temp
                    rts_cur.precip   += sub_n * rts_cur_sub.precip
                    rts_cur.snowfall += sub_n * rts_cur_sub.snowfall
                    rts_cur.rainfall += sub_n * rts_cur_sub.rainfall
                    rts_cur.swe      += sub_n * rts_cur_sub.swe
                    rts_cur.snowmelt += sub_n * rts_cur_sub.snowmelt
                    rts_cur.icemelt  += sub_n * rts_cur_sub.icemelt
                    rts_cur.melt     += sub_n * rts_cur_sub.melt
                    rts_cur.snow_outflow += sub_n * rts_cur_sub.snow_outflow
                    rts_cur.ice_outflow += sub_n * rts_cur_sub.ice_outflow
                    rts_cur.glacier_outflow += sub_n * rts_cur_sub.glacier_outflow
                    rts_cur.runoff   += sub_n * rts_cur_sub.runoff
                    rts_cur.sms      += sub_n * rts_cur_sub.sms
                    rts_cur.et_pot   += sub_n * rts_cur_sub.et_pot
                    rts_cur.et       += sub_n * rts_cur_sub.et
                    rts_cur.baseflow += sub_n * rts_cur_sub.baseflow
                    rts_cur.direct_runoff += sub_n * rts_cur_sub.direct_runoff

            # remember timestamp
            self.timestamp = date
            
            # write stats if required
            if date in self.write_dates:
                self.write_states(add_timestamp = True, verbose=True)

        results = munch.Munch()
        results.ts = pd.Panel(rts_catchments)

        return results
    
    def set_sim_period(self, t1, t2):
        '''Defines the simulation period.
        
        Parameters
        ----
        t1: start time (string, e.g. '2015-10')
        t2: termination time (string, e.g. '2016-07')
        '''       
        self.config.start_date = t1
        self.config.end_date = t2
        self.timestamp = pd.Timestamp(t1).to_period('M').to_timestamp('M')
        
    def copy_state_vars_from(self, other_aware):
        '''Transfers the state variables of another AWARE object to the current
        object.
        
        Parameters
        ----
        other_aware: an instance of an AWARE object whose states will be loaded
        '''       
        self.state_swe           = copy.deepcopy(other_aware.state_swe)
        self.state_icewe         = copy.deepcopy(other_aware.state_icewe)
        self.state_glacierarea   = copy.deepcopy(other_aware.state_glacierarea)
        self.state_soilmoisture  = copy.deepcopy(other_aware.state_soilmoisture)
        self.state_groundwater   = copy.deepcopy(other_aware.state_groundwater)
        self.is_ready = True

    def get_working_directory(self):
        '''Generates the path to the working directory ('out_dir' in the config
        file) taking into account operating system specific ways to concatenate
        sub-directories.
        
        '''       
        
        # check if variable and folder exists
        dir = './'
        if isinstance(self.config['out_dir'],str):
            if os.path.isdir(self.config['out_dir']):
                dir = self.config['out_dir']
        return dir

    def write_states(self, add_timestamp = False, verbose=False):
        '''Exports all relevant state variables to files. The function is
        also capable of handling timestamps included in the file name. In this
        way, system states will be remebered for arbitrary time steps even if
        the model is interrupted.
        
        Parameters
        ----
        add_timestamp: adds the timestamp mode to the file input functions
        verbose: prints details to the screen
        prev_month: reads grids labeled for the previous month (relevant for resuming runs)
    
        '''       
        dir = self.get_working_directory()

        dir = self.get_working_directory()
        if add_timestamp:
            timestamp = self.timestamp
        else:
            timestamp = None
        self.state_swe.export_state(dir, timestamp=timestamp, verbose=verbose)
        self.state_icewe.export_state(dir, timestamp=timestamp, verbose=verbose)
        self.state_glacierarea.export_state(dir, timestamp=timestamp, verbose=verbose)
        self.state_soilmoisture.export_state(dir, timestamp=timestamp, verbose=verbose)
        self.state_groundwater.export_state(dir, timestamp=timestamp, verbose=verbose)

    def load_states(self, add_timestamp = False, verbose=False, prev_month=False):
        '''Loads all relevant state variables needed for AWARE runs from disk.
        This is important to run the model in hotstart mode. The function is
        also capable of handling timestamps included in the file name. In this
        way, system states defined for specified time can be loaded efficiently.
        
        Parameters
        ----
        add_timestamp: adds the timestamp mode to the file input functions
        verbose: prints details to the screen
        prev_month: reads grids labeled for the previous month (relevant for resuming runs)
    
        '''       
        dir = self.get_working_directory()
        if add_timestamp:
            if prev_month:
                timestamp = self.timestamp - dateutil.relativedelta.relativedelta(months=1)
            else:
                timestamp = self.timestamp
            if verbose:
                ts_start = '%4i-%02i' % (self.timestamp.year, self.timestamp.month)
                ts_prev  = '%4i-%02i' % (timestamp.year, timestamp.month)
                if prev_month:
                    print('Simulation start: %s, loading states from previous month %s' % (ts_start, ts_prev))
                else:
                    print('Simulation start: %s, loading states from previous current %s' % (ts_start, ts_prev))
        else:
            timestamp = None

        # $kf 2016-12-14: improved error handling for import functions
        n_errors = 0
        list_failed = []
        
        try:
            self.state_swe.import_state(dir, timestamp=timestamp, verbose=verbose, basin=self.input_grids.catchments)
        except:
            n_errors += 1
            list_failed.append('SWE')
            
        try:
            self.state_icewe.import_state(dir, timestamp=timestamp, verbose=verbose, basin=self.input_grids.catchments)
        except:
            n_errors += 1
            list_failed.append('SWE ice')
            
        try:
            self.state_glacierarea.import_state(dir, timestamp=timestamp, verbose=verbose, basin=self.input_grids.catchments)
        except:
            n_errors += 1
            list_failed.append('Glacier area')
            
        try:
            self.state_soilmoisture.import_state(dir, timestamp=timestamp, verbose=verbose, basin=self.input_grids.catchments)
        except:
            n_errors += 1
            list_failed.append('Soil moisture')
            
        try:
            self.state_groundwater.import_state(dir, timestamp=timestamp, verbose=verbose, basin=self.input_grids.catchments)            
        except:
            n_errors += 1
            list_failed.append('Groundwater storage')
            
        if n_errors < 5:
            self.is_ready = True
        else:
            print('Error: Loading state variables failed.')
        if n_errors > 0 and self.is_ready:
            print('At least one state variable could not be loaded from file.')
            print('Please check the availability of the following state variables:')
            for var in list_failed:
                print(var)

    def get_catchment_info(self):
        df = pd.DataFrame(columns=['name','area','next'], index=self.catchment_ids)
        df['name'] = df['name'].astype(str)
        for ii in self.catchment_ids:
            df.loc[ii]['name'] = self.catchments[ii].name
            df.loc[ii]['area'] = self.catchments[ii].abs_area
            df.loc[ii]['next'] = self.catchments[ii].next
        return df
