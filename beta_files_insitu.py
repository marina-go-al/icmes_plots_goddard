#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 16:02:09 2018

@author: mgonza30
"""

import numpy as np
from st_average import beta_average
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from astropy.convolution import convolve, Box1DKernel
import os
from read_hapi_ST import read_STEREO_par, read_STEREO
import pandas as pd
from datetime import datetime
import PIL

rootdir = os.path.dirname(os.path.realpath(__file__))

year = '2014' # Selection of the year to be analyzed
spacecraft = 'B'

filepath = rootdir + '/data/' + year + '/' + year + spacecraft + 'insitu.xlsx'
dfpar = read_STEREO_par(filepath,spacecraft)
nelem = dfpar.ICME_datestart.count()

average_vec = [1,2,3,4,6,10,20]

box_w = 0.95
labels_fontsize = 11

save_name = [word[:-4] for word in dfpar.save_name.tolist()]

for index in range(0,nelem): 
    count = 0
    for average_beta in average_vec:
        df = read_STEREO(dfpar.download_start[index],dfpar.download_end[index],spacecraft)
        dfbeta = beta_average(average_beta,df)
        
        filter_beta = (dfbeta.prom_ddoy > dfpar.ddoy_st[index]) & (dfbeta.prom_ddoy < dfpar.ddoy_ed[index]) # Only work with values INSIDE X range
        #filter_beta = (dfbeta.prom_ddoy > dfpar.ddoy_Ist[index]) & (dfbeta.prom_ddoy < dfpar.ddoy_Ied[index])
        dfbeta  =  dfbeta[filter_beta]
        dfbeta.to_csv(rootdir + '/list_STEREO/ST' + spacecraft + '/' + year + '/' + save_name[index] + '/' + dfpar.save_name[index] + '_dtfit_' + str(average_beta) + '#hr.csv', index = False)
        count += 1
        print('File:',index+1,'->',count*100/len(average_vec),'% completed.')