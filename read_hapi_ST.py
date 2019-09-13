#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 16:00:43 2018

@author: mgonza30
"""
import pandas as pd
from datetime import datetime,timedelta
import copy
from hapiclient.hapi import hapi
import os
import math

def read_STEREO_par(filepath,spacecraft):
    '''
    Function that reads Lan Jian's excel catalog and creates a DataFrame with
    the following columns (Output):
        
        01) ICME_st -> Format: YYYY DOY M/D H:MM
                Notes: ICME Start
                
        02) MO_st -> Format: YYYY DOY M/D H:MM
                Notes: Magnetic Obstacle Start
                
        03) ICME_ed -> Format: YYYY DOY M/D H:MM
                Notes: ICME End
                
        04) ICME_datestart -> Format: YYYY-MM-DDTHH:MM:SS
                Notes: ICME Start with another format
    
        05) ICME_doystart -> Format: DOY
                Notes: ICME Doy Start
    
        06) ICME_dateend   -> Format: YYYY-MM-DDTHH:MM:SS
                Notes: ICME End with another format
                
        07) ICME_doyend -> Format: DOY
                Notes: ICME Doy End
                
        08) ddoy_Ist -> Format: DDOY
                Notes: ICME Decimal Doy Start
                
        09) ddoy_Ied -> Format: DDOY
                Notes: ICME Decimal Doy End
                
        10) ddoy_st -> Format: DDOY
                Notes: MO Decimal Doy Start
                
        11) ddoy_ed -> Format: DDOY
                Notes: MO Decimal Doy End
                
        12) save_name -> Format: YYYYMMDD_DOY_STA or _STB
                Notes: Name for the plot files that will be saved (.jpg)
    
        13) download_start -> Format: YYYY-MM-DDT00:00:00
                Notes: Start date for the data download
    
        14) download_end -> Format: YYYY-MM-DDT23:59:30
                Notes: End date for the data download
    
    '''
    df = pd.read_excel(filepath,header = None,names  =  ["ICME_st", "MO_st", "ICME_ed"]) #Reading Excel to create DataFrame
    date_start = [] #Empty lists
    doy_start = []
    date_end = []
    doy_end = []
    aux_ddoy_Ist = []
    aux_ddoy_Ied = []
    aux_ddoy_st = []
    aux_stfile = []
    aux_edfile = []
    aux_save_name = []
    
    for index in range(0,df.ICME_st.count()):
        
        aux_ICMEst = df.ICME_st[index].split()
        aux_start_str = str(aux_ICMEst[0] + '-' + aux_ICMEst[2] + '-' + aux_ICMEst[3]) #datetime str type YYYY-MM/DD-HH:MM
        aux_date_start = datetime.strptime(aux_start_str,"%Y-%m/%d-%H:%M")
        date_start.append(datetime.strftime(aux_date_start,"%Y-%m-%dT%H:%M:00"))
        doy_start.append(aux_ICMEst[1])
        hour_Ist = aux_date_start.timetuple().tm_hour
        minute_Ist = aux_date_start.timetuple().tm_min
        sec_Ist = aux_date_start.timetuple().tm_sec
        day_fraction_Ist = float(hour_Ist*3600+minute_Ist*60+sec_Ist)/float(24*3600)
        aux_ddoy_Ist.append(float(aux_ICMEst[1])+day_fraction_Ist) #decimal doy ICME start calculation
        aux_date = datetime.strftime(aux_date_start,"%Y-%m-%d")
        aux_date = datetime.strptime(aux_date,"%Y-%m-%d") - timedelta(days=2)
        aux_stfile.append(datetime.strftime(aux_date,"%Y-%m-%dT00:00:00"))
        
        aux_ICMEed = df.ICME_ed[index].split()
        aux_end_str = str(aux_ICMEed[0] + '-' + aux_ICMEed[2] + '-' + aux_ICMEed[3]) #datetime str type YYYY-MM/DD-HH:MM
        aux_date_end = datetime.strptime(aux_end_str,"%Y-%m/%d-%H:%M")
        date_end.append(datetime.strftime(aux_date_end,"%Y-%m-%dT%H:%M:00"))
        doy_end.append(aux_ICMEed[1])
        hour_Ied = aux_date_end.timetuple().tm_hour
        minute_Ied = aux_date_end.timetuple().tm_min
        sec_Ied = aux_date_end.timetuple().tm_sec
        day_fraction_Ied = float(hour_Ied*3600+minute_Ied*60+sec_Ied)/float(24*3600)
        aux_ddoy_Ied.append(float(aux_ICMEed[1])+day_fraction_Ied) #decimal doy ICME start calculation
        aux_date = datetime.strftime(aux_date_end,"%Y-%m-%d")
        aux_date = datetime.strptime(aux_date,"%Y-%m-%d") + timedelta(days=2)
        aux_edfile.append(datetime.strftime(aux_date,"%Y-%m-%dT23:59:30"))
        
        aux_MOst = df.MO_st[index].split()
        aux_start_str_MO = str(aux_MOst[0] + '-' + aux_MOst[2] + '-' + aux_MOst[3]) #datetime str type YYYY-MM/DD-HH:MM
        aux_date_start_MO = datetime.strptime(aux_start_str_MO,"%Y-%m/%d-%H:%M")
        hour_st = aux_date_start_MO.timetuple().tm_hour
        minute_st = aux_date_start_MO.timetuple().tm_min
        sec_st = aux_date_start_MO.timetuple().tm_sec
        day_fraction_st = float(hour_st*3600+minute_st*60+sec_st)/float(24*3600)
        aux_ddoy_st.append(float(aux_MOst[1])+day_fraction_st) #decimal doy ICME start calculation
        
        year_Ist = aux_date_start.timetuple().tm_year
        month_Ist = aux_date_start.timetuple().tm_mon
        if month_Ist < 10:
            month_Ist = '0'+str(month_Ist)
        else:
            month_Ist = str(month_Ist)
        day_Ist = aux_date_start.timetuple().tm_mday
        if day_Ist < 10:
            day_Ist = '0'+str(day_Ist)
        else:
            day_Ist = str(day_Ist)
        aux_doy_str = int(aux_ICMEst[1])
        if aux_doy_str < 10:
            aux_doy_str = '00'+str(aux_doy_str)
        elif aux_doy_str < 100:
            aux_doy_str = '0'+str(aux_doy_str)
        else:
            aux_doy_str = str(aux_doy_str)
        save_name = str(year_Ist) + month_Ist + day_Ist + '_' + aux_doy_str + '_ST' + spacecraft
        aux_save_name.append(save_name)
        
    df['ICME_datestart']  = copy.deepcopy(date_start) #new columns to DataFrame
    df['ICME_doystart']   = copy.deepcopy(doy_start)
    df['ICME_dateend']    = copy.deepcopy(date_end)
    df['ICME_doyend']     = copy.deepcopy(doy_end)
    df['ddoy_Ist']        = copy.deepcopy(aux_ddoy_Ist) 
    df['ddoy_Ied']        = copy.deepcopy(aux_ddoy_Ied)
    df['ddoy_st']         = copy.deepcopy(aux_ddoy_st)
    df['ddoy_ed']         = copy.deepcopy(aux_ddoy_Ied)
    df['save_name']       = copy.deepcopy(aux_save_name)
    df['download_start']  = copy.deepcopy(aux_stfile)
    df['download_end']    = copy.deepcopy(aux_edfile)
    return df

def read_STEREO(download_start,download_end,spacecraft):
    server     = 'https://cdaweb.gsfc.nasa.gov/hapi'
    dataset    = 'ST' + spacecraft + '_L2_MAGPLASMA_1M'
    start      = download_start
    stop       = download_end
    parameters = 'BTOTAL,BFIELDRTN,Np,Vp,Tp,Vth,Beta,Total_Pressure,Cone_Angle,Magnetic_Pressure,Dynamic_Pressure,Vp_RTN,HEE,R'
    #parameters = 'BTOTAL,BFIELDRTN,Np,Vp,Tp,Vth,Beta,Total_Pressure,Cone_Angle,Magnetic_Pressure,Dynamic_Pressure,Vp_RTN,HEE,R'
    opts       = {'logging': True, 'use_cache': True}    
    data,meta = hapi(server, dataset, parameters, start, stop, **opts)
    
    aux_date_time = [str(data['Time'][i]) for i in range(0,len(data['Time']))]
    date_time = [datetime.strptime(aux_date_time[i],"b'%Y-%m-%dT%H:%M:%S.%fZ'") for i in range(0,len(aux_date_time))]
    #save_time = [datetime.strftime(date_time[0],"%Y%m%d")]
    #save_time = save_time[0] #+ '_' + str(dfpar.ICME_doystart[index])
    
    HEEQx = data['HEE'][:,0]
    HEEQy = data['HEE'][:,1]
    HEEQz = data['HEE'][:,2]
    Angle = [round(math.atan2(HEEQy[index],HEEQx[index])*180/math.pi,1) for index in range(len(HEEQx))]
    Latitude = [round(HEEQz[index]*180/math.pi,1) for index in range(len(HEEQz))]
    
    d = {'date_time':date_time,'B':data['BTOTAL'],'Bx':data['BFIELDRTN'][:,0],
         'By':data['BFIELDRTN'][:,1],'Bz':data['BFIELDRTN'][:,2],'Np':data['Np'],
         'Vp':data['Vp'],'Tp':data['Tp'],'Vth':data['Vth'],'beta':data['Beta'],
         'Ptot':data['Total_Pressure'],'Vsw':data['Vp_RTN'],'Angle':Angle,
         'Latitude':Latitude,'HEEx':HEEQx,'HEEy':HEEQy,'HEEz':HEEQz,'Radi':data['R']}
         #'Latitude':Latitude,'Radi':data['R']}
    df = pd.DataFrame(data=d)
    
    aux_datetitle = []
    aux_xticks = []
    aux_doy = []
    aux_ddoy = []
    aux_year = []
    aux_day = []
    aux_hour = []
    aux_minute = []
    
    for index in range(0,len(date_time)):
        aux_datetitle.append(date_time[index].strftime('%b %-d, %Y %H:%M:%S'))
        aux_xticks.append(str(date_time[index].strftime('%b %-d')))
        year = date_time[index].timetuple().tm_year #year, day, hour, minute, sec as floats
        day = date_time[index].timetuple().tm_mday
        hour = date_time[index].timetuple().tm_hour
        minute = date_time[index].timetuple().tm_min
        sec = date_time[index].timetuple().tm_sec
        day_fraction = float(hour*3600+minute*60+sec)/float(24*3600)
        aux_doy.append(date_time[index].timetuple().tm_yday) #doy calculation
        aux_ddoy.append(float(date_time[index].timetuple().tm_yday)+day_fraction) #decimal doy calculation
        aux_year.append(float(year))
        aux_day.append(float(day))
        aux_hour.append(float(hour))
        aux_minute.append(float(minute))

    aux_xticks = [(str(round(aux_ddoy[i],1)) + '\n' + str(aux_xticks[i])) for i in range(0,len(aux_doy))] #date format for xticks axis
    df['doy']  =  copy.deepcopy(aux_doy) #new columns to DataFrame
    df['ddoy']  =  copy.deepcopy(aux_ddoy)
    df['datetime']  =  copy.deepcopy(aux_datetitle)
    df['minute'] = copy.deepcopy(aux_minute)
    df['hour'] = copy.deepcopy(aux_hour)
    df['day'] = copy.deepcopy(aux_day)
    df['year'] = copy.deepcopy(aux_year)
    df['xplot'] = copy.deepcopy(aux_xticks)
    
    filter_vals = (df.B > 1e-10) & (df.Vth < 90) & (df.Vth > 0)
    df = df[filter_vals]
    
    return df

'''
rootdir = os.path.dirname(os.path.realpath(__file__))

year = '2013' # Selection of the year to be analyzed
spacecraft = 'A'

filepath = rootdir + '/data/' + year + '/' + year + spacecraft + 'par.xlsx'
dfpar = read_STEREO_par(filepath,spacecraft)
df1 = read_STEREO(dfpar.download_start[0],dfpar.download_end[0],spacecraft)
#print(df.dtypes)
    
    #f = open(rootdir + '/data_down/' + 'ST' + spacecraft + '_' + save_time + '.txt',"w+")
    #df.to_csv(f,index=None,sep=' ',mode='w+')
'''