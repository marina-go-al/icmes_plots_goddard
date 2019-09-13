#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 09:08:04 2018

@author: Marina
"""

import numpy as np
import pandas as pd
import math
#import os
#from read_ST import read_STEREO_par, read_STEREO

#from readfiles_wi import read_wi_mag, read_YYYYpar, read_wi_swep
def beta_average(average,dfB):
    
    change_day = 0
    first_yyyyB = dfB.year.iloc[0]
    last_yyyyB = dfB.year.iloc[-1]
    dfB.index = range(len(dfB.doy))
    if first_yyyyB != last_yyyyB:
        #index_lastyyyyB = (dfB.year == last_yyyyB)
        max_doy = max(dfB.doy)
        #dfB['doy_equiv'] = dfB.doy
        a = dfB.doy.values
        c = dfB.year.values
        b = a
        for index in range(0,len(a)):
            if int(c[index])==int(last_yyyyB):
                b[index] = a[index]+max_doy    #La columna dfB.doy ya esta mal y, por tanto, el max_doy tambien
        dfB['doy_equiv'] = b
        #dfB.doy_equiv[index_lastyyyyB] = dfB.doy[index_lastyyyyB] + max_doy
        num_days = max(dfB.doy_equiv) - min(dfB.doy_equiv) + 1
        change_day = 1
    else:
        num_days = max(dfB.doy) - min(dfB.doy) + 1
        
    ndt = num_days*24*average
    n_min = np.zeros(average)
    dt_doy = np.zeros(ndt)
    dt_day = np.zeros(ndt)
    dt_hour = np.zeros(ndt)
    dt_min = np.zeros(ndt)
    
    interval_min = 60/average
    step_min = interval_min/2
    k = 0
    
    for i in range(0,average):
        n_min[i] = step_min + i*interval_min
    for m in range(0,num_days):
        dt_day[k:(m+1)*24*average] = dfB.day.iloc[0]+m
        if change_day != 1:
            dt_doy[k:(m+1)*24*average] = min(dfB.doy)+m
        else:
            dt_doy[k:(m+1)*24*average] = min(dfB.doy_equiv)+m
        k0 = k
        k = (m+1)*24*average
        hr = 0
        for l in range(k0,k,average):
            dt_hour[l:l+average] = hr
            dt_min[l:l+average] = n_min
            hr = hr+1
    """
    max_dayB = max(dfB.day)
    filter_dayB = dt_day >= 30
    filter_dayB2 = dt_day[filter_dayB] > max_dayB
    filter_dayB3 = filter_dayB[filter_dayB2]
    dt_day[filter_dayB3] = dt_day[filter_dayB3]-max_dayB
    """
    max_dayB = max(dfB.day)
    filter_dayB = [i for i in range(0,len(dt_day)) if dt_day[i]>=30]
    aux_dt_day = [dt_day[i] for i in filter_dayB]
    filter_dayB2 = [i for i in range(0,len(aux_dt_day)) if aux_dt_day[i] > max_dayB]
    aux_filter_day = [filter_dayB[i] for i in filter_dayB2]
    if len(aux_filter_day) >= 1:
        for i in aux_filter_day:
            dt_day[i] = dt_day[i]-max_dayB
    
    dt_B = np.zeros(ndt)
    dt_Bx = np.zeros(ndt)
    dt_By = np.zeros(ndt)
    dt_Bz = np.zeros(ndt)
    dt_Np = np.zeros(ndt)
    dt_Vth = np.zeros(ndt)
    dt_Vp = np.zeros(ndt)
    
    for index in range(0,ndt):
        dfB_filter = dfB
        dfB_filter = dfB_filter[(dfB_filter.day == dt_day[index]) & (dfB_filter.hour == dt_hour[index]) & (dfB.minute >= dt_min[index]-step_min) & (dfB.minute <= dt_min[index]+step_min)]
        
        if dfB_filter.empty == True:
            dt_B[index] = 0
            dt_Bx[index] = 0
            dt_By[index] = 0
            dt_Bz[index] = 0
            dt_Np[index] = 0
            dt_Vth[index] = 0
            dt_Vp[index] = 0
        else:
            dt_B[index] = dfB_filter.B.mean()
            dt_Bx[index] = dfB_filter.Bx.mean()
            dt_By[index] = dfB_filter.By.mean()
            dt_Bz[index] = dfB_filter.Bz.mean()
            dt_Np[index] = dfB_filter.Np.mean()
            dt_Vth[index] = dfB_filter.Vth.mean() 
            dt_Vp[index] = dfB_filter.Vp.mean() 
        
    dt_ddoy = dt_doy+(dt_hour+(dt_min/60))/24
        
    filter_dtB = dt_B != 0
    dt_B = dt_B[filter_dtB]
    dt_Bx = dt_Bx[filter_dtB]
    dt_By = dt_By[filter_dtB]
    dt_Bz = dt_Bz[filter_dtB]
    dt_Np = dt_Np[filter_dtB]
    dt_Vth = dt_Vth[filter_dtB]
    dt_Vp = dt_Vp[filter_dtB]
    dt_ddoy = dt_ddoy[filter_dtB]
    
    kBz = 1.3806503e-23
    Mp = 1.672621777e-27
    mu_0 = math.pi*4e-7
    
    dt_Tp = Mp*(1000*dt_Vth)**2/(2*kBz)
    P_plasma = (100**3)*dt_Np*kBz*dt_Tp
    P_mag = (dt_B*1e-9)**2/(2*mu_0)
    beta = P_plasma/P_mag
    
    d = {'prom_B':dt_B, 'prom_Bx':dt_Bx, 'prom_By':dt_By, 'prom_Bz':dt_Bz,
         'prom_Np':dt_Np, 'prom_Tp':dt_Tp, 'prom_beta':beta, 'prom_ddoy':dt_ddoy,
         'prom_Vp':dt_Vp}
    df = pd.DataFrame(data=d)
    return df
'''
rootdir = '/Users/Marina/Desktop/ZPython/'
year = '1995' # Selection of the year to be analyzed
filepath1 = rootdir + 'data/' + year + '/' + year + 'par.txt'
dfpar = read_YYYYpar(filepath1)
filepath2 = rootdir + '/data/' + year + '/' + dfpar.save_name[0] + '_WI_H0_MFI.txt'
dfB = read_wi_mag(filepath2)
'''
'''
rootdir = os.path.dirname(os.path.realpath(__file__))

year = '2013' # Selection of the year to be analyzed
spacecraft = 'A'

filepath1 = rootdir + '/data/' + year + '/' + year + spacecraft + 'par.xlsx'
dfpar = read_STEREO_par(filepath1)
filepath2 = rootdir + '/data/' + year + '/' + dfpar.save_name[0] + '.txt'
df = read_STEREO(filepath2)
#filepath3 = rootdir + '/data/' + year + '/' + dfpar.save_name[0] + '_WI_K0_SWE.txt'
#dfP  =  read_wi_swep(filepath3)

average = 4

df = beta_average(average,df)
'''