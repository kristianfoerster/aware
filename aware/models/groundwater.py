from __future__ import print_function, division, absolute_import
import numpy as np

class GroundwaterModel(object):
    def __init__(self):
        self.K = 0.1

    def groundwater_model(self, GWStorage, gw_inflow):
        GWStorage = GWStorage + gw_inflow
        baseflow = GWStorage / self.K
        if baseflow > GWStorage:
            baseflow = GWStorage
        GWStorage = GWStorage - baseflow
        return baseflow, GWStorage
