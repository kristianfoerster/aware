from __future__ import print_function, division, absolute_import
import numpy as np

class SoilModel(object):
    """The SoilModel class represents a minimal soil model. It is based on
    previous work dedicated to the Thornthwaite water balance model (see, e.g.,
    McCabe and Markstrom, 2007, Raleigh and others).
    """

    def __init__(self):
        # soil-moisture storage capacity (mm)
        self.SMSC  = 100

        # direct runoff factor
        self.DRf   = 0.15

        # percolation factor
        self.Rf    = 0.1

        # percolation fraction of soil moisture
        self.PercF = 0.1

    def soil_water_balance(self,inflow0,ET_pot,SMS):
        """soil_water_balance computes the water balance of a single soil layer
        for each grid cell. The solution of the water balance equation is derived
        for monthly time steps taking into consideration varous fluxes such as
        runoff components and evapotranspiration.
        
        Parameters
        ----
        inflow0: boundary condition - liquid water inflow (rainfall, snowmelt), grid
        ET_pot: boundary condtion - potential evapotranspiration, grid
        SMS: soil moistore storage (state variable), previous time step, grid
        
        returns
        ----
        SMS: soil moistore storage (state variable), new time step, grid
        Perc: water percolating into deeper layers (i.e. groundwater), grid
        Runoff_D: direct runoff, grid
        ET_act: actual (real) evaoptranspiration, grid
        """
        
        # apply dimension of input data to initialization of intermediate results
        Runoff_D = np.zeros(inflow0.shape)    # direct runoff
        STW      = np.zeros(inflow0.shape)    # withdrawal
        ET_act   = np.zeros(inflow0.shape)    # actual evapotranspiration
        #PET_AET  = np.zeros(inflow0.shape)    # evapotranspiration deficit # unused
        Surplus  = np.zeros(inflow0.shape)    # evapotranspiration
        Perc     = np.zeros(inflow0.shape)    # percolation
        inflow   = np.zeros(inflow0.shape)    # inflow WITHOUT infiltration excess
        SMS_1    = SMS.copy()                 # soil moisture of last ts

        # The 2015 version of the soil model was based on two loops iterating
        # over all rows and columns. Since this type of computation is slower
        # than vectorized numpy operations, a revised numerical scheme is
        # applied here
        #for i in range(0, inflow.shape[0]):
        #    for j in range(0, inflow.shape[1]):

        # direct runoff
        Runoff_D = inflow0 * self.DRf
        inflow = inflow0 * (1 - self.DRf)

        # two cases for withdrawal (STW):
        # inflow < ET_pot

        # if inflow[i,j] < ET_pot[i,j]:
        selection = inflow < ET_pot
        STW[selection] = SMS[selection] * (1 - np.exp((inflow[selection]-ET_pot[selection])/self.SMSC))
        
        #if SMS[i,j] < STW[i,j]:
        selection2 = SMS[selection] < STW[selection]
        STW_corr = STW[selection]
        SMS_corr = SMS[selection]
        STW_corr[selection2] = SMS_corr[selection2]
        SMS_corr[selection2] = 0
        # else: # if SMS[i,j] < STW[i,j]:
        SMS_corr[~selection2] = SMS_corr[~selection2] - STW_corr[~selection2]
        # update grids using the modified subsets
        STW[selection] = STW_corr
        SMS[selection] = SMS_corr

        if np.min(STW) < 0:
            print('Warning: Negative soil moisture storage withdrawal!')
        
        ET_act[selection] = inflow[selection] + STW[selection]
        # P_total + STW <= ET_pot(i,1)
        # if inflow[i,j] + STW[i,j] <= ET_pot[i,j]:
        selection2 = inflow[selection] + STW[selection] <= ET_pot[selection]
        # do nothing for cells which fulfill selection2
        ETP_selection = ET_pot[selection]
        ETR_selection = ET_act[selection]
        SMS_corr[~selection2] = SMS_corr[~selection2] + np.abs(ETR_selection[~selection2] - ETP_selection[~selection2])
        ETR_selection[~selection2] = ETP_selection[~selection2]
        # update
        SMS[selection] = SMS_corr
        ET_pot[selection] = ETP_selection
        ET_act[selection] = ETR_selection
        
        # else # inflow[i,j] < ET_pot[i,j]:
        ET_act[~selection] = ET_pot[~selection]
        SMS[~selection] = SMS[~selection] + inflow[~selection] - ET_act[~selection]
        # in case that the capacity is exceeded
        #if SMS[i,j] > self.SMSC:
        selection2 = SMS[~selection] > self.SMSC
        SMS_corr = SMS[~selection]
        Surplus_corr = Surplus[~selection]
        Surplus_corr[selection2] = SMS_corr[selection2] - self.SMSC
        SMS_corr[selection2] = self.SMSC
        # update
        SMS[~selection] = SMS_corr
        Surplus[~selection] = Surplus_corr
                
        Perc = self.PercF * SMS
        SMS = SMS - Perc      

        Perc += Surplus * self.Rf
        #Surplus = Surplus * (1-self.Rf)
        Runoff_D += Surplus * (1 - self.Rf)

        # delta in soil moisture storage
        dS = SMS_1 - SMS
        # mass balance
        mb = dS + inflow0 - ET_act - Perc - Runoff_D
        if (abs(mb) > 1.E-3).any():
            print('mass balance error')
            return

        return SMS, Perc, Runoff_D, ET_act
