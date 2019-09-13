#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 16:02:09 2018

@author: mgonza30
"""

import numpy as np
#from read_hapi import read_hapiwi_mag, read_YYYY, read_YYYYpar, read_hapiwi_swep, read_hapiPADswe
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
#import matplotlib.gridspec as gridspec
#from matplotlib.ticker import AutoMinorLocator
#rootdir = '/Volumes/Marina1000/MiPython/Copy4/'

rootdir = os.path.dirname(os.path.realpath(__file__))

year = '2013' # Selection of the year to be analyzed
spacecraft = 'B' #Selection of the spacecraft: A or B

filepath = rootdir + '/data/' + year + '/' + year + spacecraft + 'par.xlsx' #Input file directory (Lan's catalog excel)
dfpar = read_STEREO_par(filepath,spacecraft) #Input file to DataFrame
nelem = dfpar.ICME_datestart.count() #Number of events to be analyzed

average_beta = 8 #Number of averages per hour for beta plot
average_odog = 2 #Number of averages per hour for hodograms plot

box_w = 0.95 #Graphics box size
labels_fontsize = 11 #Fontsize

#for index in range(23,24):
for index in range(14,15):
    df = read_STEREO(dfpar.download_start[index],dfpar.download_end[index],spacecraft) #Getting all the data for each event
    df_Odog = df

    dfbeta = beta_average(average_beta,df) #Calculating the averages
    dfp_Odog = beta_average(average_odog,df)

    fig1 = plt.figure(1) #Organizing the graphs' structure and disposition
    gs1 = gridspec.GridSpec(6,3)
    gs1.update(hspace = 0, bottom = 0.27)
    ax1 = plt.subplot(gs1[0,:])
    ax2 = plt.subplot(gs1[1,:])
    ax3 = plt.subplot(gs1[2,:])
    ax4 = plt.subplot(gs1[3,:])
    ax5 = plt.subplot(gs1[4,:])
    ax6 = plt.subplot(gs1[5,:])

    gs2 = gridspec.GridSpec(1,3)
    gs2.update(wspace = 0.4,top = 0.18,bottom=0.05,right=0.873)
    ax8a = plt.subplot(gs2[0,0])
    ax8b = plt.subplot(gs2[0,1])
    ax8c = plt.subplot(gs2[0,2])
    
    #################--- PLOT 1: B ---###################
    
    xrange = [dfpar.ddoy_Ist[index]-1,dfpar.ddoy_Ied[index]+1] # X axis range definition
    
    filter_vals1 = (df.ddoy > xrange[0]) & (df.ddoy < xrange[1]) # Only work with values INSIDE X range
    df  =  df[filter_vals1]      
    
    if df.empty == False:
    
        if xrange[0] < min(df.ddoy):
            xrange[0] = min(df.ddoy)
        elif xrange[1] > max(df.ddoy):
            xrange[1] = max(df.ddoy)
    
        df_plot = df.values # We create a np.matrix from the DataFrame to do the plots
    
        title = fig1.text(0.5,0.93,'STEREO ' + spacecraft + '\n IMPACT, PLASTIC Combined Data \n'+ # Title
                     str(df_plot[0,-6]) + '  to  ' + str(df_plot[-1,-6]) + '\n'
                    'RTN Coordinate System. HEEQ Longitude: ' + str(df_plot[0,12]) + 'º \n Time Resolution: 1m',horizontalalignment='center',
                     verticalalignment='center',fontsize=12)
        
        yrange = [0,max(df_plot[:,1])+6] # Y range definition
        
        ddoy = df_plot[:,-7] # Keep variables to be plotted
        B = df_plot[:,1]
        
        B_plot, =ax1.plot(ddoy, B,'b.',markersize = 1)#linewidth = 0.75) #Plot B
        
        ICME_start_plot1, = ax1.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ist[index]],
                 yrange,'k-',linewidth=0.75) # ICME_start plot
        CME_start_plot1, = ax1.plot([dfpar.ddoy_st[index],dfpar.ddoy_st[index]],
                 yrange,'g--',linewidth=1) #CME_start plot
        CME_end_plot1, = ax1.plot([dfpar.ddoy_ed[index],dfpar.ddoy_ed[index]],
                 yrange,'g--',linewidth=1) # CME_end plot
        
        ICME_dur_plot, = ax1.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ied[index]],
                 [yrange[1]-0.7,yrange[1]-0.7],'k-',linewidth=2) # ICME duration plpt
        ICME_dur_plot.set_label('ICME')
        MO_dur_plot, = ax1.plot([dfpar.ddoy_st[index],dfpar.ddoy_ed[index]],
                 [yrange[1]-1.4,yrange[1]-1.4],'m-',linewidth=2) # MO duration plot
        MO_dur_plot.set_label('MO')
        
        box1 = ax1.get_position()
        ax1.set_position([box1.x0, box1.y0, box1.width * box_w, box1.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        #ax1.tick_params(right = True, top = True)
        ax1.yaxis.set_major_locator(ticker.MultipleLocator(10)) # Properties for major and minor locators
        ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        
        if dfpar.ddoy_Ist[index] != dfpar.ddoy_st[index]:
            ax1.text(dfpar.ddoy_Ist[index],yrange[1]-3.2,str(round(dfpar.ddoy_Ist[index],2)))
        ax1.text(dfpar.ddoy_st[index],yrange[1]-3.5,str(round(dfpar.ddoy_st[index],2)))
        ax1.text(dfpar.ddoy_ed[index],yrange[1]-3.5,str(round(dfpar.ddoy_ed[index],2)))
        
        ax1.set_xlim(xrange)
        ax1.set_ylim(yrange)
        ax1.set_xticks([])
        ax1.set_ylabel(r'$B$ (nT)',fontsize = labels_fontsize)
        
        ###############--- PLOT 2: Bx, By, Bz ---################
        Bx = df_plot[:,2]
        By = df_plot[:,3]
        Bz = df_plot[:,4]
        
        aux_yrange_bi = max([max(abs(df.Bx)),max(abs(df.By)),max(abs(df.Bz))])
        yrange_bi = [-aux_yrange_bi,aux_yrange_bi]
        Bx_plot, = ax2.plot(ddoy,Bx,'r.',markersize = 1)#linewidth = 0.75)
        Bx_plot.set_label(r'$B_R$')
        By_plot, = ax2.plot(ddoy,By,'b.',markersize = 1)#linewidth = 0.75)
        By_plot.set_label(r'$B_T$')
        Bz_plot, = ax2.plot(ddoy,Bz,'g.',markersize = 1)#linewidth = 0.75)
        Bz_plot.set_label(r'$B_Z$')
        
        ICME_start_plot2, = ax2.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ist[index]],
                 yrange_bi,'k-',linewidth=0.75) # ICME_start plot
        CME_start_plot2, = ax2.plot([dfpar.ddoy_st[index],dfpar.ddoy_st[index]],
                 yrange_bi,'g--',linewidth=1) #CME_start plot
        CME_end_plot2, = ax2.plot([dfpar.ddoy_ed[index],dfpar.ddoy_ed[index]],
                 yrange_bi,'g--',linewidth=1) # CME_end plot
        Horizontal_line, = ax2.plot(xrange,[0,0],'k--',linewidth=0.75) # CME_end plot
        
        box2 = ax2.get_position()
        ax2.set_position([box2.x0, box2.y0, box2.width * box_w, box2.height])
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        ax2.set_xlim(xrange)
        ax2.set_ylim(yrange_bi) 
        ax2.set_xticks([])
        ax2.set_ylabel('(nT)',fontsize = labels_fontsize)
        ax2.yaxis.set_major_locator(ticker.MultipleLocator(5)) # Properties for major and minor locators
        ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        
        ###############--- PLOT 3: Np ---################
    
        Np = df_plot[:,5]
        
        aux_yrange_Np = df.Np.median()
        yrange_Np = [0,5*aux_yrange_Np]
        
        Np_plot, =ax3.plot(ddoy, Np,'b.',markersize = 1)#linewidth = 0.75) #Plot Np
        
        ICME_start_plot3, = ax3.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ist[index]],
                 yrange_Np,'k-',linewidth=0.75) # ICME_start plot
        CME_start_plot3, = ax3.plot([dfpar.ddoy_st[index],dfpar.ddoy_st[index]],
                 yrange_Np,'g--',linewidth=1) #CME_start plot
        CME_end_plot3, = ax3.plot([dfpar.ddoy_ed[index],dfpar.ddoy_ed[index]],
                 yrange_Np,'g--',linewidth=1) # CME_end plot
        
        box3 = ax3.get_position()
        ax3.set_position([box3.x0, box3.y0, box3.width * box_w, box3.height])
        #ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        ax3.set_xlim(xrange)
        ax3.set_ylim(yrange_Np)
        ax3.set_xticks([])
        
        ax3.yaxis.set_major_locator(ticker.MultipleLocator(10)) # Properties for major and minor locators
        ax3.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax3.set_ylabel(r'$N_p$(*/cc)',fontsize = labels_fontsize)
        
        ###############--- PLOT 4: Vth ---################
        Vth = df_plot[:,8]
        
        aux_yrange_Vth = df.Vth.median()
        yrange_Vth = [0,4*aux_yrange_Vth]
        
        Vth_plot, = ax4.plot(ddoy, Vth,'b.',markersize = 1) #Plot Vth
        
        ICME_start_plot4, = ax4.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ist[index]],
                 yrange_Vth,'k-',linewidth=0.75) # ICME_start plot
        CME_start_plot4, = ax4.plot([dfpar.ddoy_st[index],dfpar.ddoy_st[index]],
                 yrange_Vth,'g--',linewidth=1) #CME_start plot
        CME_end_plot4, = ax4.plot([dfpar.ddoy_ed[index],dfpar.ddoy_ed[index]],
                 yrange_Vth,'g--',linewidth=1) # CME_end plot
        
        box4 = ax4.get_position()
        ax4.set_position([box4.x0, box4.y0, box4.width * box_w, box4.height])
        #ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        ax4.set_xlim(xrange)
        ax4.set_ylim(yrange_Vth)
        ax4.set_xticks([])
        
        ax4.yaxis.set_major_locator(ticker.MultipleLocator(25)) # Properties for major and minor locators
        ax4.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax4.set_ylabel(r'$V_{th}$(km/s)',fontsize = labels_fontsize) 
        
        ###############--- PLOT 5: Beta ---################
        
        #beta = df_plot[:,7]
        
        filter_vals2 = (dfbeta.prom_ddoy > xrange[0]) & (dfbeta.prom_ddoy < xrange[1]) # Only work with values INSIDE X range
        dfbeta  =  dfbeta[filter_vals2]
        
        dfbeta_plot = dfbeta.values
        ddoy_beta = dfbeta_plot[:,-2] # Keep variables to be plotted
        beta_prom = dfbeta_plot[:,-3]
        
        yrange_beta = [0,1]
        
        #Beta_plot, = ax5.plot(ddoy, beta,'b.',markersize = 1) #Plot beta
        Beta_prom_plot, = ax5.plot(ddoy_beta, beta_prom,'g.',markersize = 2) #Plot beta
        
        ICME_start_plot5, = ax5.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ist[index]],
                 yrange_beta,'k-',linewidth=0.75) # ICME_start plot
        CME_start_plot5, = ax5.plot([dfpar.ddoy_st[index],dfpar.ddoy_st[index]],
                 yrange_beta,'g--',linewidth=1) #CME_start plot
        CME_end_plot5, = ax5.plot([dfpar.ddoy_ed[index],dfpar.ddoy_ed[index]],
                 yrange_beta,'g--',linewidth=1) # CME_end plot
        
        box5 = ax5.get_position()
        ax5.set_position([box5.x0, box5.y0, box5.width * box_w, box5.height])
        #ax5.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        ax5.set_xlim(xrange)
        ax5.set_ylim(yrange_beta)
        ax5.set_xticks([])
        
        ax5.yaxis.set_major_locator(ticker.MultipleLocator(0.25)) # Properties for major and minor locators
        ax5.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax5.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax5.set_ylabel(r'$\beta_{proton}$',fontsize = labels_fontsize)
            
        ###############--- PLOT 6: Vsw ---################
        Vsw = df_plot[:,11]
        
        aux_Vsw = df.Vsw
        aux_Vsw_filter = (df.Vsw < 1000) & (df.Vsw > 0)
        aux_Vsw = aux_Vsw[aux_Vsw_filter]
        if len(aux_Vsw) > 1:
            minVsw = min(aux_Vsw)
            maxVsw = max(aux_Vsw)
            yrange_Vsw = [minVsw-50,maxVsw+50]
        else:
            yrange_Vsw = [0,700]
        
        Vsw_plot, = ax6.plot(ddoy, Vsw,'b.',markersize = 1) #Plot Vth
        
        ICME_start_plot6, = ax6.plot([dfpar.ddoy_Ist[index],dfpar.ddoy_Ist[index]],
                 yrange_Vsw,'k-',linewidth=0.75) # ICME_start plot
        CME_start_plot6, = ax6.plot([dfpar.ddoy_st[index],dfpar.ddoy_st[index]],
                 yrange_Vsw,'g--',linewidth=1) #CME_start plot
        CME_end_plot6, = ax6.plot([dfpar.ddoy_ed[index],dfpar.ddoy_ed[index]],
                 yrange_Vsw,'g--',linewidth=1) # CME_end plot
        
        box6 = ax6.get_position()
        ax6.set_position([box6.x0, box6.y0, box6.width * box_w, box6.height])
        #ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        ax6.set_xlim(xrange)
        ax6.set_ylim(yrange_Vsw)
        ax6.set_xlabel(year, fontsize = 'large')
        
        ax6.yaxis.set_major_locator(ticker.MultipleLocator(100)) # Properties for major and minor locators
        ax6.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))
        ax6.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax6.set_ylabel(r'$V_{SW}$(km/s)',fontsize = labels_fontsize)
        
        x_min, x_max = ax6.get_xlim()
        ticks = [(tick - x_min)/(x_max-x_min) for tick in ax6.get_xticks()]
        x_value = ax6.get_xticks()
        x_pos = []
        aux = df.ddoy.values
        for m in range(1,len(x_value)-1):
            count = np.arange(len(df.ddoy))
            filter_aux = abs(aux-x_value[m]) < 0.01
            count = count[filter_aux]
            if len(count) > 0:
                x_pos.append(int(max(count)))
                #x_pos.append(int((max(count)+min(count))/2))
        
        if len(x_pos) == len(x_value[1:-1]):
            labels = df.xplot.values
            x_labels = [labels[i] for i in x_pos]
            ax6.set_xticklabels(x_labels)
            ax6.set_xticks(x_value[1:-1])
        
        ###############--- PLOT 8: Odograms ---################
        
        filter_vals = (df_Odog.ddoy > dfpar.ddoy_st[index]) & (df_Odog.ddoy < dfpar.ddoy_ed[index])
        df_Odog = df_Odog[filter_vals]
        
        df_Odog_plot = df_Odog.values
        
        filter_vals = (dfp_Odog.prom_ddoy > dfpar.ddoy_st[index]) & (dfp_Odog.prom_ddoy < dfpar.ddoy_ed[index])
        dfp_Odog = dfp_Odog[filter_vals]
        
        dfp_Odog_plot = dfp_Odog.values
        
        #yrange = [0,max(dfB_plot[:,2])+6] # Y range definition
    
        Bx_Odog = df_Odog_plot[:,2] # Keep variables to be plotted
        By_Odog = df_Odog_plot[:,3]
        Bz_Odog = df_Odog_plot[:,4]
        
        Bxprom_Odog = dfp_Odog_plot[:,1]
        Byprom_Odog = dfp_Odog_plot[:,2]
        Bzprom_Odog = dfp_Odog_plot[:,3]
    
        Odog_a_plot1, = ax8a.plot(Bx_Odog, By_Odog,'b.',markersize=0.8)
        Odog_a_plot2, = ax8a.plot(Bxprom_Odog, Byprom_Odog,'m-',linewidth=0.8)
        Odog_a_plot3, = ax8a.plot(Bxprom_Odog[0], Byprom_Odog[0],'ro',markersize=6)
        
        box8a = ax8a.get_position()
        ax8a.set_position([box8a.x0, box8a.y0, box8a.width * box_w, box8a.height])
        #ax8a.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        aux_Bt_range = max(abs(max(By_Odog)),abs(min(By_Odog)))
        Bt_range = [-aux_Bt_range,aux_Bt_range]
        aux_Br_range = max(abs(max(Bx_Odog)),abs(min(Bx_Odog)))
        Br_range = [-aux_Br_range,aux_Br_range]
        aux_Bn_range = max(abs(max(Bz_Odog)),abs(min(Bz_Odog)))
        Bn_range = [-aux_Bn_range,aux_Bn_range]
        
        ax8a.yaxis.set_major_locator(ticker.MultipleLocator(5)) # Properties for major and minor locators
        ax8a.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax8a.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax8a.set_ylim(Bt_range)
        ax8a.set_xlim(Br_range)
        ax8a.set_ylabel(r'$B_T$',fontsize = labels_fontsize)
        ax8a.set_xlabel(r'$B_R$',fontsize = labels_fontsize)
        ax8a.set_title("(nT)")
        
        Odog_b_plot1, =ax8b.plot(Bx_Odog, Bz_Odog,'b.',markersize=0.8)
        Odog_b_plot2, = ax8b.plot(Bxprom_Odog, Bzprom_Odog,'m-',linewidth=0.8)
        Odog_b_plot3, = ax8b.plot(Bxprom_Odog[0], Bzprom_Odog[0],'ro',markersize=6)
        
        box8b = ax8b.get_position()
        ax8b.set_position([box8b.x0, box8b.y0, box8b.width * box_w, box8b.height])
        #ax8b.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        ax8b.yaxis.set_major_locator(ticker.MultipleLocator(5)) # Properties for major and minor locators
        ax8b.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax8b.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax8b.set_ylim(Bn_range)
        ax8b.set_xlim(Br_range)
        ax8b.set_ylabel(r'$B_N$',fontsize = labels_fontsize)
        ax8b.set_xlabel(r'$B_R$',fontsize = labels_fontsize)
        ax8b.set_title("(nT)")
        
        Odog_c_plot1, =ax8c.plot(Bz_Odog, By_Odog,'b.',markersize=0.8)
        Odog_c_plot2, = ax8c.plot(Bzprom_Odog, Byprom_Odog,'m-',linewidth=0.8)
        Odog_c_plot3, = ax8c.plot(Bzprom_Odog[0], Byprom_Odog[0],'ro',markersize=6)
        
        box8c = ax8c.get_position()
        ax8c.set_position([box8c.x0, box8c.y0, box8c.width * box_w, box8c.height])
        #ax8c.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # Legend
        
        ax8c.yaxis.set_major_locator(ticker.MultipleLocator(5)) # Properties for major and minor locators
        ax8c.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax8c.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax8c.set_ylim(Bt_range)
        ax8c.set_xlim(Bn_range)
        ax8c.set_ylabel(r'$B_T$',fontsize = labels_fontsize)
        ax8c.set_xlabel(r'$B_N$',fontsize = labels_fontsize)
        ax8c.set_title("(nT)")
        
        # Saving files
        fig1.set_size_inches(w=9,h=12)
        '''
        name_a = rootdir + '/grafs/' + str(dfpar.save_name[index]) + '.png'
        fig1save = fig1.savefig(name_a)
        plt.close(fig1)
        '''
