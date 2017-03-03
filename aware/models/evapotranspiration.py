from __future__ import print_function, division, absolute_import
import numpy as np
import calendar

class EvapotranspirationModel:
    def __init__(self,lat, J=None, ts_temp=None):
        if J==None and ts_temp==None:
            raise('Either J or ts_temp must be specified!')
        # empirical correction to account for shadows
        self.delta_daylength = 0 #-5 # Zimmermann (2015)
        # winter time humidity for Ivanov model
        self.wintertime_hum = 85

        # day of mont 15 for each month
        jdays_mi = np.array([15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349])

        # maximum day length for each month
        self.s0_mi = np.maximum(self.max_sunshine_duration(jdays_mi, lat)+self.delta_daylength, 0)

        if J:
            self.J = J
        else:
            if ts_temp is None:
                raise 'Parameter J not specified. Please provide at least a temperature series for estimating J!'
            self.J = EvapotranspirationModel.compute_parameter_j(ts_temp)

        # parameter alpha
        self.alpha = (0.0675 * self.J**3 - 7.71 * self.J**2 + 1792 * self.J +49239) * 1.E-5


    def max_sunshine_duration(self, doy, lat):
        ceta=0.0172*doy-1.39
        s0 = 12.3 + np.sin(ceta) * (4.3 + 1.0*(lat - 51) / 6.0)
        return s0

    def monthly_evapotranspiration(self,date,temp, n_etp_summer=1.):
        # daylength of month
        D = self.s0_mi[date.month-1]
        d = calendar.monthrange(date.year, date.month)

        ET_pot = np.zeros(temp.shape)

        pos_tpos =  temp > 273.15
        # calculation according to Thornthwaite (1948)
        ET_pot[pos_tpos] = n_etp_summer * 0.533 * d[1] * D / 12 * (10 * np.maximum(temp[pos_tpos] - 273.15,0) / self.J )**self.alpha
        
        # calculation according to Ivanov (Wendling & Mueller, 1984)
        ET_pot[~pos_tpos] = np.maximum(0, 0.0011 * (25 + (temp[~pos_tpos] - 273.15))**2 * (100-self.wintertime_hum) )

        return ET_pot

    @staticmethod
    def compute_parameter_j(temperature_time_series):
        # long-term temperature mean
        longterm_t = temperature_time_series.groupby(temperature_time_series.index.month).mean() - 273.15

        # parameter J of Thornthwaite model
        J = 0
        for i in range(0, 12):
            if np.isscalar(longterm_t.values[i]):
                val_i = longterm_t.values[i]
            else:
                val_i = longterm_t.values[i][0]
            ti = np.maximum(val_i, 1.)
            J += (ti / 5)**1.514
        return J

