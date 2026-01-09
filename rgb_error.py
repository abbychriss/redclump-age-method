#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 18:41:52 2023

@author: abbychriss
"""

import numpy as np
from numpy import median
import math
import matplotlib.pyplot as plt
from pylab import *
from astropy.table import Table
import matplotlib.cm as cm
from mpl_point_clicker import clicker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import random
from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance

cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)
#open table of cluster parameters (cluster name and E(B-V)) from Kharchenko N. V.,
#Piskunov A. E., Schilbach E., Röser S., Scholz R. D., 2013, A&A, 558, A53
cluster_parameters = Table.read('cluster_parameters.dat', format='ascii')
cluster_parameters = cluster_parameters.group_by('Cluster')
plx = [0.320,0.358,0.472,2.290,0.850,0.518,0.211,0.848,1.144,0.322,0.560,
       0.360,0.267,0.686,0.951,0.379,0.330,0.619,0.285,0.432,0.187,0.344,
       0.282,0.380,0.429]


#Padova isochrone function:

# for convenience, get the list of ages and metallicities (the latter all zero
# in this case.)
a,feh=np.loadtxt('/Users/abbychriss/Desktop/WSU/Isochrones/pops.solar11',usecols=(0,1),unpack=True)
#print(len(a),len(feh))

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

#New BASTI isochrone function:
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

e_bprps = []
rgb_err = []
slopes = []
for i in range(len(cluster_info['Cluster'])):
    cluster=cluster_info['Cluster'][i]
    #M = m + 5 - 5 log d - A
    for k in range(len(cluster_parameters['Cluster'])):
        if cluster_parameters['Cluster'][k] == cluster:
            g_corr = 5 - (5 * np.log10(1000 / plx[i])) - (2.740 * float(cluster_parameters['E(B-V)'][k]))
            E_BP_RP = (1.339 * float(cluster_parameters['E(B-V)'][k]))
            e_bprps.append(E_BP_RP)

    cluster_out = Table.read('/Users/abbychriss/Desktop/WSU/'+cluster+'_out.dat', format='ascii')
    bprp=cluster_out['BP-RP']
    g=cluster_out['G']
    prob=cluster_out['Prob']

    
    bprp_filt = []
    g_filt = []
    prob_filt = []
    for j in range(len(bprp)):
        #apply reddening and distance modulus to get absolute magnitudes for each cluster
        """if cluster=='NGC_7789' or cluster=='NGC_2660' or cluster=='Melotte_71' or cluster=='NGC_2236' or cluster=='NGC_2420' or cluster=='NGC_2506' or cluster=='NGC_6791' or cluster=='NGC_2682' or cluster=='Melotte_66':
            if prob[j] >= 0.7:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                prob_filt.append(prob[j])
        elif cluster=='NGC_1245':
            if prob[j] >= 0.8: 
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                prob_filt.append(prob[j])
        elif cluster=='NGC_2627':
            if prob[j] >= 0.96: 
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                prob_filt.append(prob[j])
        elif cluster=='NGC_2447' or cluster=='NGC_6940':
            if prob[j] >= 0.75: 
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                prob_filt.append(prob[j])
        elif cluster=='Collinder_110' or cluster=='Czernik_37':
            if prob[j] >= 0.7:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                prob_filt.append(prob[j])
        else:"""
        if prob[j] >= 0.7:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                prob_filt.append(prob[j])
    """
    #plot Padova isochrones
    iiso = 42
    teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchiso(iiso)
    br = bmag - rmag
    plot(br,gmag,label='Age: '+str(a[iiso])+' Gyr (Padova)',color='#14B4D7',linewidth=1.5)
    """
    
    #plot BASTI isochrones
    cluster_iiso = [21,18,21,21,17,20,39,22,25,17,21,20,20,
                    19,14,22,16,15,18,19,24,18,18,17,21]
    iiso = cluster_iiso[i]
    teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchbasti(iiso)
    br = bmag - rmag
    for j in np.linspace(-0.25,0.15,40):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
        plot(br + j,gmag,label='Age '+str(abasti[iiso])+' Gyr (BASTI)',color='#3AD714')
        im = scatter([x - E_BP_RP for x in bprp_filt],[y + g_corr for y in g_filt], s=5, c=prob_filt, cmap=cm.plasma_r)
        ylim([9.5,-5.5])
        xlim([-0.5,3.0])
        title('Color Magnitude Diagram: ' + cluster.replace('_',' '))
        xlabel('BP-RP (mag)')
        ylabel('G (mag)')
        legend()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        colorbar(im, cax=cax, label='Prob')
        show()
    
    rgb_errors = [0.03, 0.06, 0.01, 0.01, 0.035, 0.015, 0.01, 0.02, 0.01, 0.035, 0.01, 0.04, 0.01, 0.03, 
                  0.00, 0.015, 0.06, 0.11, 0.025, 0.095, 0.025, 0.015, 0.00, 0.00, 0.025]

    
    """klicker = clicker(ax, ["event"], markers=["x"])
    
    print(klicker.get_positions()["event"])"""
    
    #find slope of best fit isochrone for each cluster:
    #for 1.25 Gyr use stars 128-148
    #for 7.5 Gyr use stars 74 to 84
    m = (gmag[148] - gmag[128]) / (br[148] - br[128])
    slopes.append(m)

#points on the RGB in absolute magnitude with dust corrections
rgb_point = [[ 1.33726076, -0.10013136],
[ 1.38989049, -1.30411958],
[ 1.34431378, -0.19511843],
[1.14394512, 0.80268595],
[ 1.09043202, -0.31284084],
[1.25291812, 0.32484399],
[1.25793968, 2.24338   ],
[ 1.32497082, -0.18838196],
[1.16868885, 1.0061351 ],
[1.03965532, 0.19560634],
[1.06581738, 1.08407123],
[ 1.18826963, -0.08953871],
[1.17943447, 0.01113173],
[ 1.27856351, -0.01219992],
[ 1.11493376, -0.60364311],
[1.18673861, 0.54146989],
[ 1.0011218,  -0.60838674],
[ 1.05289967, -0.17524034],
[1.14161332, 0.21812206],
[1.19702673, 0.46630966],
[1.20523775, 0.57977034],
[ 1.17693509, -0.36573527],
[ 1.05388534, -0.14505678],
[ 1.18691454, -1.13441465],
[1.23746425, 0.27312082]]
