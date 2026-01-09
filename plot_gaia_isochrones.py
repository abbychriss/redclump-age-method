#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:20:56 2024

@author: abbychriss
"""

import numpy as np
from numpy import median
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import *
from astropy.table import Table
from astropy.io import ascii
import matplotlib.cm as cm
from mpl_point_clicker import clicker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import random
from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance

rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = True

cluster_keep = ['King_5','Melotte_66','Melotte_71','NGC_1245','NGC_2420','NGC_2477',
                'NGC_2506','NGC_2509','NGC_2627','NGC_2682','NGC_6208','NGC_6791',
                'NGC_7789','Ruprecht_68']

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.where(np.isin(np.array(cluster_info['Cluster']),cluster_keep))]
selected_cluster_info.rename_column('Input file name', 'Gaia file')

plx = selected_cluster_info['Plx']

#open table of cluster parameters (cluster name and E(B-V)) from Kharchenko N. V.,
#Piskunov A. E., Schilbach E., Röser S., Scholz R. D., 2013, A&A, 558, A53
cluster_parameters = Table.read('cluster_parameters.dat', format='ascii')
cluster_parameters = cluster_parameters.group_by('Cluster')

#Isochrone function
# for convenience, get the list of ages and metallicities (the latter all zero
# in this case.)
a,feh=np.loadtxt('/Users/abbychriss/Desktop/WSU/Isochrones/pops.solar11',usecols=(0,1),unpack=True)

rgb_errors = [0.03,0.01,0.01,0.02,0.01,0.01,0.04,0.01,0.03,0.015,0.025,0.095,0.025,0.015]

#Marigo 2007 [Fe/H]=0.00:
def fetchiso(iiso):
    if iiso > 65:
        print('Error! iiso must be 65 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('/Users/abbychriss/Desktop/WSU/Isochrones/hess_gaia_ev6res0mixss.out','r') as filehandle:
        filecontents = filehandle.readlines()
        #print(filecontents)

    teff=[];logl=[];gmag=[];bmag=[];rmag=[];nst=[];indx=[]
    nHA = 0
    popcount = -1
    for line in filecontents:
        row = line.split()
        if row[0] == '#(':     # end old block, make a plot start new one
            popcount = popcount + 1
            ireadthis = 0
            if popcount == iiso:     # store this one. Otherwise, spin on.
                ireadthis = 1
        elif row[0] == '#':
            # do nothing
            fiddle = 0
        elif ireadthis == 1:
            nHA = nHA + 1
            # collect data
            teff.append(float(row[3]))
            logl.append(float(row[2]))
            nst.append(float(row[4]))
            gmag.append(float(row[5]))
            bmag.append(float(row[6]))
            rmag.append(float(row[7]))
            indx.append(int(row[0]))
            #print(popcount,nHA,row[0],row[1],row[2])
        else:
            # do nothing
            fiddle = 0
    # Teff, log L/Lo, number of stars in bin
    teff = np.array(teff)
    logl = np.array(logl)
    nst  = np.array(nst)
    # Gaia BP, Gaia RP, integer array of bin number
    bmag = np.array(bmag)
    gmag = np.array(gmag)
    rmag = np.array(rmag)
    indx = np.array(indx)
    # nHA is just an integer - the number of point in the isochrone
    return teff,logl,nst,gmag,bmag,rmag,indx,nHA

#Basti 5.0.1 [Fe/H]=+0.06
abasti,fehbasti=np.loadtxt('Isochrones/pops.solarbasti',usecols=(0,1),unpack=True,skiprows=1)
def fetchbasti(iiso):
    if iiso > 46:
        print('Error! iiso must be 46 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('Isochrones/hess_gaia_ev5res0mixss.out','r') as filehandle:
        filecontents = filehandle.readlines()

    teff=[];logl=[];gmag=[];bmag=[];rmag=[];nst=[];indx=[]
    nHA = 0
    popcount = -1
    for line in filecontents:
        row = line.split()
        if row[0] == '#(':     # end old block, make a plot start new one
            popcount = popcount + 1
            ireadthis = 0
            if popcount == iiso:     # store this one. Otherwise, spin on.
                ireadthis = 1
        elif row[0] == '#':
            # do nothing
            fiddle = 0
        elif ireadthis == 1:
            nHA = nHA + 1
            # collect data
            teff.append(float(row[3]))
            logl.append(float(row[2]))
            nst.append(float(row[4]))
            gmag.append(float(row[5]))
            bmag.append(float(row[6]))
            rmag.append(float(row[7]))
            indx.append(int(row[0]))
            #print(popcount,nHA,row[0],row[1],row[2])
        else:
            # do nothing
            fiddle = 0
    # Teff, log L/Lo, number of stars in bin
    teff = np.array(teff)
    logl = np.array(logl)
    nst  = np.array(nst)
    # Gaia BP, Gaia RP, integer array of bin number
    bmag = np.array(bmag)
    rmag = np.array(rmag)
    indx = np.array(indx)
    # nHA is just an integer - the number of point in the isochrone
    return teff,logl,nst,gmag,bmag,rmag,indx,nHA

#BASTI 6 [Fe/H]=+0.06
abasti6,fehbasti6=np.loadtxt('Isochrones/pops_solarbasti6.tex',usecols=(0,1),unpack=True,skiprows=1)

def fetchbasti6(iiso):
    if iiso > 46:
        print('Error! iiso must be 46 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('Isochrones/hess_gaia_ev7res0mixss.out','r') as filehandle:
        filecontents = filehandle.readlines()

    teff=[];logl=[];gmag=[];bmag=[];rmag=[];nst=[];indx=[]
    nHA = 0
    popcount = -1
    for line in filecontents:
        row = line.split()
        if row[0] == '#(':     # end old block, make a plot start new one
            popcount = popcount + 1
            ireadthis = 0
            if popcount == iiso:     # store this one. Otherwise, spin on.
                ireadthis = 1
        elif row[0] == '#':
            # do nothing
            fiddle = 0
        elif ireadthis == 1:
            nHA = nHA + 1
            # collect data
            teff.append(float(row[3]))
            logl.append(float(row[2]))
            nst.append(float(row[4]))
            gmag.append(float(row[5]))
            bmag.append(float(row[6]))
            rmag.append(float(row[7]))
            indx.append(int(row[0]))
            #print(popcount,nHA,row[0],row[1],row[2])
        else:
            # do nothing
            fiddle = 0
    # Teff, log L/Lo, number of stars in bin
    teff = np.array(teff)
    logl = np.array(logl)
    nst  = np.array(nst)
    # Gaia BP, Gaia RP, integer array of bin number
    bmag = np.array(bmag)
    rmag = np.array(rmag)
    indx = np.array(indx)
    # nHA is just an integer - the number of points in the isochrone
    return teff,logl,nst,gmag,bmag,rmag,indx,nHA

#PARSEC+COLIBRI [Fe/H]=0.00
aparsec,fehparsec=np.loadtxt('Isochrones/pops_solarparsec.tex',usecols=(0,1),unpack=True,skiprows=1)
def fetchparsec(iiso):
    if iiso > 35:
        print('Error! iiso must be 35 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('Isochrones/hess_gaia_ev8res0mixss.out','r') as filehandle:
        filecontents = filehandle.readlines()
        #print(filecontents)

    teff=[];logl=[];gmag=[];bmag=[];rmag=[];nst=[];indx=[]
    nHA = 0
    popcount = -1
    for line in filecontents:
        row = line.split()
        if row[0] == '#(':     # end old block, make a plot start new one
            popcount = popcount + 1
            ireadthis = 0
            if popcount == iiso:     # store this one. Otherwise, spin on.
                ireadthis = 1
        elif row[0] == '#':
            # do nothing
            fiddle = 0
        elif ireadthis == 1:
            nHA = nHA + 1
            # collect data
            teff.append(float(row[3]))
            logl.append(float(row[2]))
            nst.append(float(row[4]))
            gmag.append(float(row[5]))
            bmag.append(float(row[6]))
            rmag.append(float(row[7]))
            indx.append(int(row[0]))
            #print(popcount,nHA,row[0],row[1],row[2])
        else:
            # do nothing
            fiddle = 0
    # Teff, log L/Lo, number of stars in bin
    teff = np.array(teff)
    logl = np.array(logl)
    nst  = np.array(nst)
    # Gaia BP, Gaia RP, integer array of bin number
    bmag = np.array(bmag)
    rmag = np.array(rmag)
    indx = np.array(indx)
    # nHA is just an integer - the number of points in the isochrone
    return teff,logl,nst,gmag,bmag,rmag,indx,nHA

prob_cutoff=[0.4, #rup68
      0.4, #ngc7789
      0.95, #ngc6791
      0.1, #ngc6208
      0.95, #ngc2682
      0.4, #ngc2627
      0.25, #ngc2509
      0.95, #ngc2506
      0.6, #ngc2477
      0.8, #ngc2420
      0.3, #ngc1245
      0.3, #mel71
      0.75, #mel66
      0.35] #king_5

g_bounds = [[7,-3], #rup68
      [9,-4], #ngc7789
      [8.5,-2.5], #ngc6791
      [10,-3], #ngc6208
      [11,-3], #ngc2682
      [10,-3], #ngc2627
      [7,-2], #ngc2509
      [9,-2.5], #ngc2506
      [10,-3.75], #ngc2477
      [9,-4], #ngc2420
      [8,-4], #ngc1245
      [9,-4], #mel71
      [8,-4], #mel66
      [8,-4]] #king_5

bprp_bounds = [[0,2], #rup68
      [0,2.5], #ngc7789
      [0,2.5], #ngc6791
      [-1,3], #ngc6208
      [-0.25,3], #ngc2682
      [-0.,2.5], #ngc2627
      [0,2], #ngc2509
      [-0.5,2.5], #ngc2506
      [-0.25,2.5], #ngc2477
      [-0.25,2.25], #ngc2420
      [-0.25,2.0], #ngc1245
      [0,2.25], #mel71
      [-0.2,2.25], #mel66
      [0,2.5]] #king_5

iiso_list_basti =[25, #rup68
      25, #ngc7789
      36, #ngc6791
      26, #ngc6208
      31, #ngc2682
      25, #ngc2627
      24, #ngc2509
      26, #ngc2506
      24, #ngc2477
      27, #ngc2420
      22, #ngc1245
      24, #mel71
      29, #mel66
      23] #king_5

iiso_list_parsec =[17, #rup68
      17, #ngc7789
      26, #ngc6791
      18, #ngc6208
      22, #ngc2682
      17, #ngc2627
      16, #ngc2509
      18, #ngc2506
      16, #ngc2477
      18, #ngc2420
      14, #ngc1245
      15, #mel71
      20, #mel66
      14] #king_5

iso_bounds_basti = [[0,336,380,446,453,-1], #rup68
      [0,336,380,446,453,-1], #ngc7789
      [0,435,455,540,550,700], #ngc6791
      [0,340,380,452,459,-1], #ngc6208
      [0,400,434,453,465,-1], #ngc2682
      [0,350,380,452,459,-1], #ngc2627
      [0,350,380,452,459,-1], #ngc2509
      [0,350,380,452,459,-1], #ngc2506
      [0,350,380,452,459,-1], #ngc2477
      [0,347,380,440,452,-1], #ngc2420
      [0,350,380,452,459,-1], #ngc1245
      [0,350,380,452,459,-1], #mel71
      [0,360,360,410,435,-1], #mel66
      [0,350,380,452,459,-1]] #king_5

iso_bounds_parsec = [[0,121,140,165,175,250], #rup68
      [0,121,140,165,175,250], #ngc7789
      [0,123,139,150,166,230], #ngc6791
      [0,113,131,160,175,240], #ngc6208
      [0,132,160,160,160,210], #ngc2682
      [0,113,131,160,175,240], #ngc2627
      [0,113,131,160,175,240], #ngc2509
      [0,113,131,160,175,240], #ngc2506
      [0,113,131,150,175,240], #ngc2477
      [0,113,131,158,162,225], #ngc2420
      [0,113,131,160,175,240], #ngc1245
      [0,113,131,160,175,240], #mel71
      [0,105,120,150,159,222], #mel66
      [0,113,131,160,175,240]] #king_5

med_rc_bprp_list = []
med_rc_g_list = []
med_rgb_bprp_list = []
d_bprp_list = []
e_bprps = []
total_errors = []
fe_h = []
rc_numbers = []
rgb_numbers = []
n_rc80_list = []
for i in range(len(selected_cluster_info['Cluster'])):
    
    cluster=selected_cluster_info['Cluster'][i]
    
    #CALCULATE REDENNING AND DISTANCE CORRECTIONS
    #M = m + 5 - 5 log d - A
    for k in range(len(cluster_parameters['Cluster'])):
        #apply reddening and distance modulus to get absolute magnitudes for each cluster
        if cluster_parameters['Cluster'][k] == cluster:
            g_corr = 5 - (5 * np.log10(1000 / plx[i])) - (2.740 * float(cluster_parameters['E(B-V)'][k]))
            E_BP_RP = (1.339 * float(cluster_parameters['E(B-V)'][k]))
            e_bprps.append(E_BP_RP)
            fe_h.append(cluster_parameters['[Fe/H]'][k])

    cluster_out = Table.read('/Users/abbychriss/Desktop/WSU/'+cluster+'_out2.dat', format='ascii')
    bprp=cluster_out['BP-RP']
    g=cluster_out['G']
    prob=cluster_out['Prob']
    
    #PROBABILITY FILTER
    bprp_filt = []
    g_filt = []
    prob_filt = []
    for j in range(len(bprp)):
        if prob[j] >= prob_cutoff[i]:
            bprp_filt.append(bprp[j])
            g_filt.append(g[j])
            prob_filt.append(prob[j])
    
    #APPLY REDENNING CORRECTIONS
    corr_bprp_filt = np.array([x - E_BP_RP for x in bprp_filt])
    corr_g_filt = np.array([y + g_corr for y in g_filt])
    corr_colors = np.zeros(len(corr_bprp_filt), dtype={'names':('corr_bprp_filt', 'corr_g_filt'),
                              'formats':('f8', 'f8')})
    corr_colors['corr_bprp_filt'] = corr_bprp_filt
    corr_colors['corr_g_filt'] = corr_g_filt
    corr_colors = np.sort(corr_colors, order='corr_g_filt')
    
    bprp_filt_80 = []
    g_filt_80 = []
    prob_filt_80 = []
    for j in range(len(bprp)):
        if prob[j] >= 0.6:
            bprp_filt_80.append(bprp[j])
            g_filt_80.append(g[j])
            prob_filt_80.append(prob[j])
 
    corr_bprp_filt_80 = np.array([x - E_BP_RP for x in bprp_filt_80])
    corr_g_filt_80 = np.array([y + g_corr for y in g_filt_80])
    
    #GET ISOCHRONES
    iiso_basti = iiso_list_basti[i]
    teff,logl,nst,gmag_basti,bmag_basti,rmag_basti,indx,nHA=fetchbasti6(iiso_basti)
    br_basti = bmag_basti - rmag_basti
    
    iiso_parsec = iiso_list_parsec[i]
    teff,logl,nst,gmag_parsec,bmag_parsec,rmag_parsec,indx,nHA=fetchparsec(iiso_parsec)
    br_parsec = bmag_parsec - rmag_parsec
    
    #PLOT CMD
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
    im = scatter(corr_colors['corr_bprp_filt'],corr_colors['corr_g_filt'], s=5, c=prob_filt, cmap=cm.plasma_r)#, alpha=0.75)
    
    #PLOT ISOCHRONES
    plot(br_basti[iso_bounds_basti[i][0]:iso_bounds_basti[i][1]],
         gmag_basti[iso_bounds_basti[i][0]:iso_bounds_basti[i][1]],
         label=str(abasti6[iiso_basti])+' Gyr',color='#55cbf2') #(BaSTI (2018) [Fe/H]=+0.06)
    plot(br_basti[iso_bounds_basti[i][2]:iso_bounds_basti[i][3]],
         gmag_basti[iso_bounds_basti[i][2]:iso_bounds_basti[i][3]],color='#55cbf2')
    plot(br_basti[iso_bounds_basti[i][4]:iso_bounds_basti[i][5]],
         gmag_basti[iso_bounds_basti[i][4]:iso_bounds_basti[i][5]],color='#55cbf2')
    
    plot(br_parsec[iso_bounds_parsec[i][0]:iso_bounds_parsec[i][1]],
         gmag_parsec[iso_bounds_parsec[i][0]:iso_bounds_parsec[i][1]],
         label=str(aparsec[iiso_parsec])+' Gyr',color='#61e861')#, alpha=0.5) #(PARSEC (2012) [Fe/H]=0.00)
    plot(br_parsec[iso_bounds_parsec[i][2]:iso_bounds_parsec[i][3]],
         gmag_parsec[iso_bounds_parsec[i][2]:iso_bounds_parsec[i][3]],color='#61e861')
    plot(br_parsec[iso_bounds_parsec[i][4]:iso_bounds_parsec[i][5]],
         gmag_parsec[iso_bounds_parsec[i][4]:iso_bounds_parsec[i][5]],color='#61e861')
    
    #ylim([g_bounds[i][0],g_bounds[i][1]-1])
    #xlim(bprp_bounds[i])
    ylim(7,-2.5)
    xlim(0,2.5)
    #title('Color Magnitude Diagram: ' + cluster.replace('_',' '))
    xlabel(r'$BP-RP$ (mag)',fontsize=50)
    ylabel(r'$G$ (mag)',fontsize=50)
    yticks(fontsize=50)
    xticks(fontsize=50)
    legend(frameon=False,fontsize=40,loc=2)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(im, cax=cax)
    cax.tick_params(labelsize=50)
    cb.set_label(label='Prob', size=50)
    savefig('/Users/abbychriss/Desktop/'+selected_cluster_info['Cluster'][i]+'_iso_cmd.pdf', bbox_inches='tight', format='pdf')
    close()
    #show()
