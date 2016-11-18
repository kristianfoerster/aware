# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:33:07 2016

@author: kristianf
"""

import netCDF4
import numpy as np
import datetime
import calendar
import pandas as pd

###############################################################################
# useful functions
def read_nc_file(cfsfile, var, no_time_dim=False):
    nc_cfs = netCDF4.Dataset(cfsfile, 'r')
    nc_var = nc_cfs.variables[var]
    if no_time_dim:
        variable = nc_var[:,:][:]
    else:
        variable = nc_var[0,:,:][:] # return first dataset only
    nc_cfs.close()
    return variable
def read_time_dim(array, time_fmt):
    # spacing = time_fmt.units.split()[0] # assumed to be 1 month
    t1      = datetime.datetime.strptime(time_fmt.units.split()[2], '%Y-%m-%d')
    tr      = pd.date_range(t1, periods=len(array), freq='M')
    return tr
def get_time_index(list_ts, date):
    di = 0
    for dii in list_ts:
        date_i = dii.date()
        if  date_i == date:
            break
        di += 1   
    return di
###############################################################################
# settings

# select meteo input
meteo_inp_type = 'cfs'
#meteo_inp_type = 'glosea'
#meteo_inp_type = 'histalp'

# files input
path_histalp_prec = '/Users/kristianf/projekte/MUSICALS/Daten/histalp/HISTALP_precipitation_all_abs_1801-2010_inn.nc'
path_histalp_temp = '/Users/kristianf/projekte/MUSICALS/Daten/histalp/HISTALP_temperature_1780-2008_inn.nc'

if meteo_inp_type == 'cfs':
    path_input = '/Users/kristianf/projekte/MUSICALS/CFS/reforecasts/'
    model_mapping_file = '/Users/kristianf/projekte/MUSICALS/AWARE/awaretestdata/data/cfs_mapping.npz'
    year1 = 1992 # start glosea
    year2 = 2008 # last year histalp
if meteo_inp_type == 'glosea':
    year1 = 1992 # start glosea
    year2 = 2008 # last year histalp
    path_input = '/Users/kristianf/projekte/MUSICALS/Daten/UKMO_seasonal/'
    model_mapping_file = '/Users/kristianf/projekte/MUSICALS/AWARE/awaretestdata/data/glosea_mapping.npz'
if meteo_inp_type == 'histalp':
    year1 = 1992 # start glosea
    year2 = 2008 # last year histalp


basedir = '/Users/kristianf/projekte/MUSICALS/AWARE/awaretestdata/data/'
outfile = '%sclimatology_%s_%i-%i.npz' % (basedir, meteo_inp_type, year1, year2)

###############################################################################
# step 1: read reforecast data
dim_i = 0
dim_j = 0
n     = 0

list_file = list()
list_days_per_month = list()
list_month = list()
skip_months = []

# prepare hindcast files
if  meteo_inp_type == 'cfs':
    field_temp = 'TMP_2maboveground'
    field_prec = 'PRATE_surface'
    n_real = 4
    
    for yy in range(year1,year2+1):
        di = datetime.date(yy,1,1)
        last_month = 0
        while di.year == yy:
            if di.day <= 25: #if di.month != last_month:
                for hh in range(0,24,6):
                    file = '%sflxf.01.%i%02i%02i%02i.nc' % (path_input,di.year,di.month,di.day,hh)
                    list_file.append(file)
                    list_month.append(di.month)
                    list_days_per_month.append(calendar.monthrange(di.year, di.month)[1])
                #last_month = di.month
            # add 5 days
            if di.month == 2 and di.day == 25 and calendar.isleap(yy):
                # no shift of fixed scheme in leap years
                di += datetime.timedelta(6)
            else:
                di += datetime.timedelta(5)
    first_file = list_file[0]

elif meteo_inp_type == 'glosea':
    field_temp = 'tas'
    field_prec = 'pr'
    n_real = 8
    #skip_months = [1,2,3,4,5,6,7,8,9,10,12]
    # months to be evaluated:
    #m = [1,2,11,12] # [5,11]
    m_start = [11] #[5,11]
    all_months = np.arange(1,13)
    #skip_months = list()
    #for mii in all_months:
    #    if mii not in m:
    #        skip_months.append(mii)
    for ii in range(year1,year2+1):
        for jj in m_start:
            for ll in range(0,4):
                month = jj + ll
                year = ii
                if month > 12:
                    month -= 12
                    year += 1
                for kk in range(0,n_real):
                    string0 = 'Amon_GloSea5_seasonal_S%4i%02i25_r%ii1p1_%4i%02i-%4i%02i.nc' % \
                        (ii,(jj-1),(kk+1),year,month,year,month)
                    list_file.append(string0)
                    list_days_per_month.append(calendar.monthrange(ii, jj)[1])
                    list_month.append(month)
                    string1 = 'Amon_GloSea5_seasonal_S%4i%02i01_r%ii1p1_%4i%02i-%4i%02i.nc' % \
                        (ii,jj,(kk+1),year,month,year,month)
                    list_file.append(string1)
                    list_days_per_month.append(calendar.monthrange(year, month)[1])
                    list_month.append(month)
                    if ll > 0:
                        string2 = 'Amon_GloSea5_seasonal_S%4i%02i09_r%ii1p1_%4i%02i-%4i%02i.nc' % \
                            (ii,jj,(kk+1),year,month,year,month)
                        list_file.append(string2)
                        list_days_per_month.append(calendar.monthrange(year, month)[1])
                        list_month.append(month)
        first_file = '%stas_%s' % (path_input, list_file[0])

if meteo_inp_type == 'histalp':
    nc_histalp_t = netCDF4.Dataset(path_histalp_temp,'r')
    nc_histalp_p = netCDF4.Dataset(path_histalp_prec,'r')
    h_time_t     = nc_histalp_t.variables['time']
    h_time_p     = nc_histalp_p.variables['time']
    vtemp_ha     = nc_histalp_t.variables['T_2M']
    vprecip_ha   = nc_histalp_p.variables['TOT_PREC']
            
    list_time_t = read_time_dim(h_time_t[:], h_time_t)
    list_time_p = read_time_dim(h_time_p[:], h_time_p)
    startdate   = datetime.date(year1,1,31) # np.datetime64('%i%02i%02i' % (year1,1,31))
    
    i1_t = get_time_index(list_time_t, startdate)
    i1_p = get_time_index(list_time_p, startdate)
    
    i2_t = i1_t + (year2-year1+1) * 12 # (12-len(skip_months)) - 1
    i2_p = i1_p + (year2-year1+1) * 12 # (12-len(skip_months)) - 1
    
    dim_i = vtemp_ha[0,:,:].shape[0]
    dim_j = vtemp_ha[0,:,:].shape[1]
    
    n_ts = (year2-year1+1) * 12   
    
    stack_T = np.zeros([n_ts,dim_i,dim_j])
    stack_P = np.zeros([n_ts,dim_i,dim_j])
    
    stack_T[:,:,:] = vtemp_ha[i1_t:i2_t,:,:]
    stack_P[:,:,:] = vprecip_ha[i1_p:i2_p,:,:]
    
    list_month = list_time_t[i1_t:i2_t].month

else:
    # get dimension from first file
    data = read_nc_file(first_file,field_temp)
    dim_i = data.shape[0]
    dim_j = data.shape[1]
    
    n_ts = len(list_file)    
    
    # stack of data
    stack_T = np.zeros([n_ts,dim_i,dim_j])
    stack_P = np.zeros([n_ts,dim_i,dim_j])
    
    t = 0
    for fi in list_file:
        #print(fi)
        try:
            if  meteo_inp_type == 'cfs': 
                temp_grid = read_nc_file(fi,field_temp)      
                prec_grid = read_nc_file(fi,field_prec)
                print(fi)
            elif meteo_inp_type == 'glosea':
                path_tfile = '%stas_%s' % (path_input, fi)
                path_pfile = '%spr_%s' % (path_input, fi)
                print(path_tfile)
                temp_grid = read_nc_file(path_tfile,field_temp)
                print(path_pfile)
                prec_grid = read_nc_file(path_pfile,field_prec)
                
    
            stack_T[t,:,:] = temp_grid
            stack_P[t,:,:] = prec_grid * list_days_per_month[t] * 86400
        except:
            print('As the file %s was not found, it will not be included.' % fi)
            # workaround to handle no data
            if t > 0:
                stack_T[t,:,:] = stack_T[t-1,:,:]
                stack_P[t,:,:] = stack_P[t-1,:,:]
                print('Last time step has been used instead.')
            else:
                stack_T[t,:,:] = 0
                stack_P[t,:,:] = 0
                print('Error: no data available!')
        t += 1

###############################################################################
# compute climatology

month_clim_T = np.zeros([12,dim_i,dim_j])
month_clim_P = np.zeros([12,dim_i,dim_j])
month_stdv_T = np.zeros([12,dim_i,dim_j])
month_stdv_P = np.zeros([12,dim_i,dim_j])
n_months = np.zeros([12])

# climatology
for mi in range(0,12):
    for i in range(0,n_ts):
        m = list_month[i]
        if m == mi+1:
            month_clim_T[m-1] += stack_T[i,:,:]
            month_clim_P[m-1] += stack_P[i,:,:]
            n_months[m-1]+=1

for mi in range(0,12):
    if n_months[mi] > 0:
        month_clim_T[mi,:,:] /= n_months[mi]
        month_clim_P[mi,:,:] /= n_months[mi]
    else:
        month_clim_T[mi,:,:] = np.nan
        month_clim_P[mi,:,:] = np.nan

# std_dev.
for mi in range(0,12):
    for i in range(0,n_ts):
        m = list_month[i]
        if m == mi+1:
            month_stdv_T[m-1] += (month_clim_T[mi,:,:] - stack_T[i,:,:])**2
            month_stdv_P[m-1] += (month_clim_P[mi,:,:] - stack_P[i,:,:])**2

for mi in range(0,12):
    month_stdv_T[mi,:,:] = np.sqrt(month_stdv_T[mi,:,:] / (n_months[mi] - 1))
    month_stdv_P[mi,:,:] = np.sqrt(month_stdv_P[mi,:,:] / (n_months[mi] - 1))

np.savez(outfile, mean_temp=month_clim_T, std_temp=month_stdv_T, \
                    mean_prec=month_clim_P, std_prec=month_stdv_P)
print('Results saved to file %s.' % outfile)
