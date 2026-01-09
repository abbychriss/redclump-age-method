#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 18:22:29 2024

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

g_bounds_zoom = [[7,-3], #rup68
      [2.25,-2.25], #ngc7789
      [3.75,-1.5], #ngc6791
      [10,-3], #ngc6208
      [3.25,-2], #ngc2682
      [10,-3], #ngc2627
      [7,-2], #ngc2509
      [3,-1.75], #ngc2506
      [10,-3.75], #ngc2477
      [2.75,-2], #ngc2420
      [8,-4], #ngc1245
      [9,-4], #mel71
      [3,-2.5], #mel66
      [8,-4]] #king_5

bprp_bounds_zoom = [[0,2], #rup68
      [0.75,1.8], #ngc7789
      [1,2], #ngc6791
      [-1,3], #ngc6208
      [0.75,1.8], #ngc2682
      [-0.,2.5], #ngc2627
      [0,2], #ngc2509
      [0.75,1.5], #ngc2506
      [-0.25,2.5], #ngc2477
      [0.6,1.75], #ngc2420
      [-0.25,2.0], #ngc1245
      [0,2.25], #mel71
      [0.5,1.75], #mel66
      [0,2.5]] #king_5

iiso_list=[26, #rup68
      24, #ngc7789
      26, #ngc6791
      22, #ngc6208
      23, #ngc2682
      18, #ngc2627
      18, #ngc2509
      23, #ngc2506
      19, #ngc2477
      22, #ngc2420
      19, #ngc1245
      26, #mel71
      28, #mel66
      23] #king_5

iso_bounds = [[67,120], #rup68
      [70,126], #ngc7789
      [68,130], #ngc6791
      [68,130], #ngc6208
      [70,128], #ngc2682
      [90,160], #ngc2627
      [90,140], #ngc2509
      [70,140], #ngc2506
      [90,140], #ngc2477
      [65,135], #ngc2420
      [90,140], #ngc1245
      [60,130], #mel71
      [55,115], #mel66
      [70,130]] #king_5

iso_shift = [[0,-1.6], #rup68
      [0.01,-1.5], #ngc7789
      [0,0], #ngc6791
      [-.04,-1], #ngc6208
      [-.12,-.28], #ngc2682
      [-.2,-.2], #ngc2627
      [-.15,-.7], #ngc2509
      [-.14,-.9], #ngc2506
      [0.013,-.94], #ngc2477
      [-.1,-.8], #ngc2420
      [-.1,-.9], #ngc1245
      [0.03,-2.7], #mel71
      [-.12,-1.15], #mel66
      [-0.06,-1.98]] #king_5

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
    
    #GET ISOCHRONE
    iiso = iiso_list[i]
    teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchparsec(iiso)
    br = bmag - rmag
    
    #PLOT CMD
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 15))
    im = scatter(corr_colors['corr_bprp_filt'],corr_colors['corr_g_filt'], s=10, c=prob_filt, cmap=cm.plasma_r)#, alpha=0.75)
    
    #PLOT ISOCHRONE SNIPPET
    plot(np.array(br[iso_bounds[i][0]:iso_bounds[i][1]])+iso_shift[i][0],
         np.array(gmag[iso_bounds[i][0]:iso_bounds[i][1]])+iso_shift[i][1],
         #label=str(aparsec[iiso])+' Gyr',
         color='#48e8bb')#, alpha=0.5)
    
    #Absolute magnitudes of bp-rp and g limits of red clump
    rc_bprp_lim = [[0.95,1.24],[1.05,1.2],[1.21,1.413],[1.05,1.20],[1.1,1.22],[0.9,1.15],
                [0.96,1.02],[1.00,1.15],[0.95,1.18],[1.12,1.18],[0.95,1.1],[0.97,1.15],
                [1.075,1.1865],[1.00,1.11]]
    
    rc_g_lim = [[0.4,-0.2],[0.68,-0.04],[0.66,-0.12],[0.5,0.1],[0.4,0.2],[0.4,-0.3],
             [0.5,0.2],[0.3,-0.1],[0.8,-0.5],[0.21,-0.06],[0.6,-0.4],[0.72,0.08],
             [0.2,-0.16],[0.4,-0.4]]
    
    rc_rgb_bprp_lim = [[0.95,1.86],[1.0,1.95],[1.14,2.31],[1.0,1.4],[1.0,2.0],[0.85,1.4],
                       [0.9,1.4],[0.95,2.0],[0.94,2.0],[1.0,1.7],[0.95,1.3],[0.9,1.9],
                       [0.924,1.8],[0.95,1.42]]
    
    rc_rgb_g_lim = [[1.5,-2.7],[1.95,-3],[3.75,-1.35],[1.2,-1.0],[2.984,-2],[2.0,-2.3],
                    [1.5,-1.4],[2.35,-2.4],[1.5,-3.0],[2.285,-2.0],[1.0,-1.5],[1.0,-2.0],
                    [3.1,-3.0],[0.4,-2.0]]
    
    #FIND MEDIAN LOCATION OF RED CLUMP
    rc_bprp = []
    rc_g = []
    nrgb_rc=0
    #collect all the stars in the red clump by taking the stars within the bp-rp and g limits above
    for j in range(len(bprp_filt)):
        if rc_bprp_lim[i][0] <= corr_colors['corr_bprp_filt'][j] <= rc_bprp_lim[i][1] and rc_g_lim[i][1] <= corr_colors['corr_g_filt'][j] <= rc_g_lim[i][0]:
            rc_bprp.append(corr_colors['corr_bprp_filt'][j])
            rc_g.append(corr_colors['corr_g_filt'][j])
        
        if rc_rgb_bprp_lim[i][0] <= corr_colors['corr_bprp_filt'][j] <= rc_rgb_bprp_lim[i][1] and rc_rgb_g_lim[i][1] <= corr_colors['corr_g_filt'][j] <= rc_rgb_g_lim[i][0]:
            nrgb_rc+=1
    rc_numbers.append(len(rc_bprp))
    n_rgb = nrgb_rc - len(rc_bprp)
    rgb_numbers.append(n_rgb)
    
    #use astropy statistics biweight_location to find the central location of the sample
    #if the median absolute deviation is zero, the biweight location is simply the median
    med_rc_bprp = round(biweight_location(rc_bprp),6)
    med_rc_g = round(biweight_location(rc_g),6)
    med_rc_bprp_list.append(med_rc_bprp)
    med_rc_g_list.append(med_rc_g)
    
    #CALCUATE NUMBER OF RC STARS WITH PROB>=80
    n_rc80=0
    #collect all the stars in the red clump by taking the stars within the bp-rp and g limits above
    for j in range(len(bprp_filt_80)):
        if rc_bprp_lim[i][0] <= corr_bprp_filt_80[j] <= rc_bprp_lim[i][1] and rc_g_lim[i][1] <= corr_g_filt_80[j] <= rc_g_lim[i][0]:
            n_rc80+=1
    n_rc80_list.append(n_rc80)
    
    scatter(med_rc_bprp,med_rc_g,s=200,color='white',edgecolor='black',marker='*')

    #FIND COLOR OF RGB AT LEVEL OF RC
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx],idx
    
    iso_rgb_bprp = np.array(br[iso_bounds[i][0]:iso_bounds[i][1]])+iso_shift[i][0] #bp-rp of RGB of isochrone, shifted to fit cluster
    iso_rgb_g = np.array(gmag[iso_bounds[i][0]:iso_bounds[i][1]])+iso_shift[i][1] #gmag of RGB of isochrone, shifted to fit cluster
    
    rgb_g,idx = find_nearest(iso_rgb_g, med_rc_g)
    med_rgb_bprp = iso_rgb_bprp[idx]
    med_rgb_bprp_list.append(round(med_rgb_bprp,6))
    scatter(med_rgb_bprp,med_rc_g,s=200,color='white',edgecolor='black',marker='*')
    
    d_bprp = med_rgb_bprp - med_rc_bprp
    d_bprp_list.append(round(d_bprp,6))
    
    #plot tangent line to RGB at median RGB location
    iso_slope = (iso_rgb_g[idx]-iso_rgb_g[idx-1])/(iso_rgb_bprp[idx]-iso_rgb_bprp[idx-1])
    #plt.plot(np.linspace(1,1.7,80), iso_slope*(np.linspace(1,1.7,80)- med_rgb_bprp) + med_rc_g,color='orange',alpha=0.5)
    
    #DRAW BOX AROUND RED CLUMP
    plot(np.linspace(rc_bprp_lim[i][0], rc_bprp_lim[i][1],80), np.full(80, rc_g_lim[i][0]), color='r',lw=1)
    plot(np.linspace(rc_bprp_lim[i][0], rc_bprp_lim[i][1],80), np.full(80, rc_g_lim[i][1]), color='r',lw=1)
    plot(np.full(80, rc_bprp_lim[i][0]), np.linspace(rc_g_lim[i][0],rc_g_lim[i][1],80), color='r',lw=1)
    plot(np.full(80, rc_bprp_lim[i][1]), np.linspace(rc_g_lim[i][1], rc_g_lim[i][0], 80), color='r',lw=1)

    #ylim([g_bounds_zoom[i][0],g_bounds_zoom[i][1]])
    #xlim([bprp_bounds_zoom[i][0],bprp_bounds_zoom[i][1]])
    ylim(3,-1.5)
    xlim(0.5,2)
    #title('Color Magnitude Diagram: ' + cluster.replace('_',' '))
    xlabel(r'$BP-RP$ (mag)',fontsize=50)
    ylabel(r'$G$ (mag)',fontsize=50)
    yticks(fontsize=50)
    xticks(fontsize=50)
    #legend(frameon=False,fontsize=25)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(im, cax=cax)
    cax.tick_params(labelsize=50)
    cb.set_label(label='Prob', size=50)
    fig.set_size_inches(12,15)
    savefig('/Users/abbychriss/Desktop/'+selected_cluster_info['Cluster'][i]+'_cmd.pdf', bbox_inches='tight', format='pdf')
    close()
    #show()
    
    #find total error for d_bprp:
    
    #log error due to G mag is approximately -3.5, so G mag error is 10^-3.5 = 0.000316
    #fluctuation in G mag contributes 0.000316/iso_slope - Divide by slope to quantify how much a change in G affects the color
    
    #determine variance of red clump bp-rp and g using astropy statistics biweight_midvariance
    
    d_bprp_error = math.sqrt((rgb_errors[i]**2/n_rgb) + (biweight_midvariance(rc_bprp)/len(rc_bprp)) + ((0.000316/iso_slope)**2/nrgb_rc) + ((10**(-3))**2/nrgb_rc) + ((10**(-3.25))**2/nrgb_rc)) 
    d_bprp_error_new = math.sqrt(d_bprp_error**2 + 0.004**2)
    total_errors.append(round(d_bprp_error_new,6))
    #total error is the sqrt of the red clump standard error^2, photometric g mag error/sqrt(N red clump stars)^2, photometric BP error/sqrt(N red clump stars)^2, photometric RP error/sqrt(N red clump stars)^2

print('Median BP-RP red clump: ',med_rc_bprp_list)
print('Median G red clump: ',med_rc_g_list)
print('Median BP-RP red giant branch: ',med_rgb_bprp_list)
print('d(BP-RP): ',d_bprp_list)
print('Total d(BP-RP) error:', total_errors)
print('nRC: ', [(selected_cluster_info['Cluster'][i], rc_numbers[i]) for i in range(len(rc_numbers))])
print('nRGB: ', [(selected_cluster_info['Cluster'][i], rgb_numbers[i]) for i in range(len(rgb_numbers))])

#open list of cluster parameters from Cantat-Gaudin, 2020 to get ages
table1_2020 = Table.read('/Users/abbychriss/Desktop/WSU/table1_2020.dat.txt', format='ascii')
headers = ['Cluster', 'RAdeg', 'DEdeg', 'GLON', 'GLAT', 'r50', 'nbstars07', 'pmRA*', 'e_pmRA*', 'pmDE', 'e_pmDE', 'plx', 'e_plx', 'Flag', 'AgeNN', 'AVNN', 'DMNN', 'DistPc', 'X', 'Y', 'Z', 'Rgc']
for i in range(len(headers)):
    table1_2020.rename_column('col'+str(i+1), headers[i])
table1_2020 = table1_2020.group_by('Cluster')

cluster_ages = []
for i in range(len(table1_2020['Cluster'])):
    if np.flip(table1_2020['Cluster'])[i] in selected_cluster_info['Cluster']:
        cluster_ages.append(np.flip(table1_2020['AgeNN'])[i])
cluster_ages = np.array(cluster_ages,dtype=float)

age_log_err = []
for age in cluster_ages:
    if age <= 9.25:
        age_log_err.append(0.175)
    else:
        age_log_err.append(0.15)
        
age_err = [((age_log_err[i] * np.log(10) * 10**cluster_ages[i])/10**9) for i in range(len(age_log_err))]

final_data = np.zeros(14, dtype={'names':('Cluster', 'd(BP-RP)', 'error d(BP-RP)', 'log age', 'error log age','nRC80', '[Fe/H]', 'E(BP-RP)'),
                          'formats':('U20', 'f8', 'f8', 'f8', 'f8', 'int', 'f8', 'f8')})
final_data['Cluster'] = selected_cluster_info['Cluster']
final_data['d(BP-RP)'] = d_bprp_list
final_data['error d(BP-RP)'] = total_errors
final_data['log age'] = cluster_ages
final_data['error log age'] = age_log_err
final_data['nRC80'] = n_rc80_list
final_data['E(BP-RP)'] = e_bprps
final_data['[Fe/H]'] = fe_h
final_data = np.sort(final_data, order='log age')
ascii.write(final_data, 'final_data_d_bprp.dat', overwrite=True)
final_table = Table.read('final_data_d_bprp.dat',format='ascii')
print(final_table)

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6))
scatter(final_data['log age'],final_data['d(BP-RP)'],color='#8748B4')
#title('Age Calibration Fit')
xlabel('Age (log(years))',fontsize=20)
ylabel(r'$d_{BP-RP}$ (mag)',fontsize=20)
yticks(fontsize=20)
xticks(fontsize=20)
xlim([8.1,10.3])
errorbar(final_data['log age'],final_data['d(BP-RP)'], xerr=age_log_err, yerr=total_errors, fmt="o",capsize=4,capthick=1,color='#8748B4')#,alpha=0.3) #put in error bars
#savefig('/Users/abbychriss/Desktop/age_calibration_plot.pdf', bbox_inches='tight', format='pdf')
#close()
show()

#Error bar implementation:
#Include error due to Gaia photometry, reddening parameter?, and distance modulus (for the level of red clump)
#Include uncertainty in choice of clump starts through bootstrap analysis
#Include error due to age estimates in Cantat-Gaudin, 2020 - "the uncertainty on the determination of log t ranges
#from 0.15 to 0.25 for young clusters and from 0.1 to 0.2 for old clusters"
#we are in the intermediate regime, so we should use 0.175 (log age) for intermediate age error, 0.15 for old
#old clusters are any cluster older than 9 log years??

print('E(BP-RP) =', [str(e_bprps[i])+ ' ' + cluster_info['Cluster'][i] for i in range(len(e_bprps))])

#do the K-S test from scipy
from scipy import stats
x = d_bprp_list
#remove the clusters older than 9.5 log years (mel 66, ngc 6791, ngc 2682)
for i in sorted([2,4,12], reverse=True):
    del x[i]
    del e_bprps[i]
    del fe_h[i]
print(stats.kstest(x, stats.norm.cdf))

