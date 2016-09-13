from __future__ import print_function, division, absolute_import
import numpy as np
import pandas as pd
import calendar
import matplotlib.pyplot as plt

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
            # long-term temperature mean
            self.longterm_t = ts_temp.groupby(ts_temp.index.month).mean() - 273.15

            # parameter J of Thornthwaite model
            self.J = 0
            for i in range(0, 12):
                ti = np.maximum(self.longterm_t[i+1], 0)
                self.J += (ti / 5)**1.514

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

if __name__ == "__main__":
    data = pd.read_csv('/Users/kristianf/Downloads/rinn_monthly.csv',
                       # 'input/meteo_gepatschdamm_monthly.csv',
                       index_col = 'date', parse_dates=True)
    # data = data[data.index.year>1986]+273.15
    data = data[data.index.year<2008]+273.15
    test = evapotranspiration(47., ts_temp=data['temp'])
    et = np.zeros(len(data))
    for i in range(0,len(data)):
        et[i] = test.monthly_evapotranspiration(data.index[i], data.temp[i])

    fig = plt.figure()
    plt.plot(et)
    plt.xlabel('time (month since simulations started)')
    plt.ylabel('ET [mm/mon]')
    plt.grid(True)
    plt.show()
