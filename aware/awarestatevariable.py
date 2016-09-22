# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:47:07 2016

@author: kristianf
"""

from . import util
import numpy as np

class AwareStateVariable:
    
    def __init__(self, name, prj_settings=None, init_grid=None, init_single_val=None, index=None):
        if prj_settings is None:
            raise 'Insufficent information provided! Either projection or x/y dimension needs to be specified.'
        self.name = name
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
        if index is not None:
            return self.state[index]
        else:
            return self.state
    
    def set_state(self, data, index=None):
        if index is not None:
            self.state[index] = data
        else:
            self.state = data
    
    def import_state(self, filename):
        data, prj_settings = util.read_gdal_file(self.config.dtm_file, return_prj_settings=True)
        self.state = data
        self.projection_settings = prj_settings
    
    def export_state(self, filename, verbose=False):
        if self.projection_settings is not None:
            data_type = self.projection_settings['dtype']
            util.write_gdal_file(filename, self.state.astype(data_type), self.projection_settings)
            if verbose:
                print('Variable %s written to %s.' % (self.name, filename))
        else:
            print('No projection settings defined!')
    