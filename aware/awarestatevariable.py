# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:47:07 2016

@author: kristianf
"""

from . import util
import numpy as np
import os

class AwareStateVariable:
    '''AwareStateVariable is class including methods and attributes to manage
    state variables required for hydrological simulations using the AWARE model.
    '''
    
    def __init__(self, name, long_name, prj_settings=None, init_grid=None, init_single_val=None, index=None):
        '''Constructor method for AwareStateVariable.
        
        Parameters
        ----
        name: short name used for naming of files
        long_name: a more detailed label of the object
        prj_settings*: dictionary including the projection settings (rasterio format)
        init_grid*: inital state provided as grid
        init_single_val*: value that will be uniformly applied to the grid
        index*: array index representing a subset of the grid which will be initialized

        '''
        if prj_settings is None:
            raise 'Insufficent information provided! Either projection or x/y dimension needs to be specified.'
        self.name = name
        self.description = long_name
        self.projection_settings = prj_settings
        cols = prj_settings['width']
        rows = prj_settings['height']
        if init_grid is not None:
            self.state = init_grid
        else:
            self.state = np.zeros([rows,cols]) * np.nan
            if init_single_val is not None:
                if index is None:
                    self.state = init_single_val
                else:
                    self.state[index] = init_single_val
    
    def get_state(self, index=None):
        '''Returns the state variable's grid.
        
        Parameters
        ----
        index*: optional argument to request a subset of the grid values
    
        Returns
        ----    
        system state as grid
        '''
        if index is not None:
            return self.state[index]
        else:
            return self.state
    
    def set_state(self, data, index=None):
        '''Passes a grid to the system state grid and overwrites it.
        
        Parameters
        ----
        data: grid representing the new system state
        index*: optional subset definition
    
        '''       
        if index is not None:
            self.state[index] = data
        else:
            self.state = data
    
    def import_state(self, filename, timestamp=None, verbose=False, basin=None):
        '''Read system state from file.
        
        Parameters
        ----
        filename: file that includes the system state grid to be read
        timestamp*: reads grid only that includes this time step in the file name (optional)
        verbose*: optional argument forcing the method to display additional information
    
        '''
        if os.path.isdir(filename):
            # this is a directory, pass automatically generated file name to path
            infile = os.path.join(filename, self.get_state_filename(timestamp))
        else:
            infile = filename
        data, prj_settings = util.read_gdal_file(infile, return_prj_settings=True)
        # force zero assignement for negative values
        data[data<0] = 0
        self.state = data
        self.projection_settings = prj_settings
        if basin is not None:
            i_data_basin = basin>0
            i_invalid = np.isnan(data) & i_data_basin
            n_invalid = i_invalid.sum()
            if n_invalid > 0:
                print('%i invalid cells detected. Set to zero!' % n_invalid)
                data[i_invalid] = 0
        if verbose:
            print('State variable %s updated from file %s.' % (self.name, infile))
                
    def export_state(self, filename, timestamp=None, verbose=False):
        '''Exports system state to a local file.
        
        Parameters
        ----
        filename: target file
        timestamp*: adds a timestamp to the file name (optional)
        verbose*: optional argument forcing the method to display additional information
    
        '''
        if os.path.isdir(filename):
            # this is a directory, pass automatically generated file name to path
            outfile = os.path.join(filename, self.get_state_filename(timestamp))
        else:
            outfile = filename
        if self.projection_settings is not None:
            data_type = self.projection_settings['dtype']
            util.write_gdal_file(outfile, self.state.astype(data_type), self.projection_settings)
            if verbose:
                print('State variable %s written to %s.' % (self.name, outfile))
        else:
            print('No projection settings defined!')

    def get_state_filename(self, timestamp=None):
        '''Generates a file name for i/o of system states taking model time
        into consideration.
        
        Parameters
        ----
        timestamp*: reads grid only that includes this time step in the file name (optional)

        Returns
        ----
        file name string
    
        '''
        extension = '.tif'
        if timestamp is not None:
            extension = '_%4i%02i.tif' % (timestamp.year, timestamp.month)
        
        filename = '%s%s' % (self.name, extension)
        return filename
        