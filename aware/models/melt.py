from __future__ import print_function, division, absolute_import
import numpy as np

class MeltModel(object):
    def melt(self, swe, prec, temp, ddf, glacier_fraction=None):
        melt = np.zeros(swe.shape)
        rain = np.zeros(swe.shape)
        snow = np.zeros(swe.shape)
        melt_avail_factor = np.zeros(swe.shape)

        # rain or snow
        pos = temp > self.trans_temp + 0.5 * self.temp_trans_range
        rain[pos] = prec[pos]
        snow[pos] = 0

        pos = (temp <= self.trans_temp + 0.5 * self.temp_trans_range) & (temp > self.trans_temp - 0.5 * self.temp_trans_range)
        t1 = self.trans_temp - 0.5 * self.temp_trans_range
        t2 = self.trans_temp + 0.5 * self.temp_trans_range
        rain_fraction = (temp[pos] - t1) / (t2 - t1)
        rain[pos] = rain_fraction * prec[pos]
        snow[pos] = (1 - rain_fraction) * prec[pos]

        pos = temp <= self.trans_temp - 0.5 * self.temp_trans_range
        rain[pos] = 0
        snow[pos] = prec[pos]

        # precipitation correction
        rain *= self.rain_corr
        snow *= self.snow_corr

        # potential melt
        pot_melt = np.maximum(0, ddf * (temp - self.melt_temp) )

        if glacier_fraction is not None:
            pot_melt *= glacier_fraction

        # water balance
        pos = swe + snow >= pot_melt
        pos2 = swe + snow < pot_melt # don't use ~pos here because of possible nan values in snow or pot_melt
        swe[pos] += snow[pos] - pot_melt[pos]
        melt[pos] = pot_melt[pos]
        melt[pos2] = swe[pos2] + snow[pos2]
        swe[pos2] = 0

        # outflow is melt plus rain
        outflow = rain + melt
        
        # if there is some energy for ice melt availbale, return a melt fraction
        pos = pot_melt > 0 | ~np.isnan(pot_melt)
        melt_avail_factor[pos] = 1 - melt[pos] / pot_melt[pos]

        # return swe and melt
        return swe, melt, outflow, snow, rain, melt_avail_factor
