#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 22:28:52 2024

@author: abbychriss
"""

import numpy as np
from numpy import median
import math
import matplotlib as mpl
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

rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = True

cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)
#open table of cluster parameters (cluster name and E(B-V)) from Kharchenko N. V.,
#Piskunov A. E., Schilbach E., Röser S., Scholz R. D., 2013, A&A, 558, A53
cluster_parameters = Table.read('cluster_parameters.dat', format='ascii')
cluster_parameters = cluster_parameters.group_by('Cluster')
plx = [0.320,0.358,0.472,2.290,0.850,0.518,0.211,0.848,1.144,0.322,0.560,
       0.360,0.267,0.686,0.951,0.379,0.330,0.619,0.285,0.432,0.187,0.344,
       0.282,0.380,0.429]

n_rgb = [6,5,73,4,6,18,222,5,32,8,4,9,43,48,1,20,6,3,13,6,104,9,1,1,40]

rgb_errors = [0.03, 0.06, 0.01, 0.01, 0.035, 0.015, 0.01, 0.02, 0.01, 0.035, 0.01, 0.04, 0.01, 0.03, 
              0.00, 0.015, 0.06, 0.11, 0.025, 0.095, 0.025, 0.015, 0.00, 0.00, 0.025] #in magnitude

e_bprps = []
rc_numbers = []

#point slope form of a line: y - y1 = m(x-x1)
#function that finds the point along the red giant branch at y = med_g_rc
def f(x1,y1,y):
    x = ((y - y1)/-8.367) + x1
    return x
#x = med_rgb_bprp, y = med_rc_g, y1 = rgb_g, x1 = rgb_bprp

for i in range(len(cluster_info['Cluster'])):
    cluster=cluster_info['Cluster'][i]
    #M = m + 5 - 5 log d - A
    for k in range(len(cluster_parameters['Cluster'])):
        #apply reddening and distance modulus to get absolute magnitudes for each cluster
        if cluster_parameters['Cluster'][k] == cluster:
            g_corr = 5 - (5 * np.log10(1000 / plx[i])) - (2.740 * float(cluster_parameters['E(B-V)'][k]))
            E_BP_RP = (1.339 * float(cluster_parameters['E(B-V)'][k]))
            e_bprps.append(E_BP_RP)

    cluster_out = Table.read('/Users/abbychriss/Desktop/WSU/'+cluster+'_out2.dat', format='ascii')
    bprp=cluster_out['BP-RP']
    g=cluster_out['G']
    total_probs2=cluster_out['Prob']
    
    bprp_filt = []
    g_filt = []
    total_probs2_filt = []
    for j in range(len(bprp)): #filter out stars of probability less than ___
        if cluster_info['Cluster'][i] == 'Ruprecht_68' or cluster_info['Cluster'][i] == 'Czernik_37' or cluster_info['Cluster'][i] == 'FSR_1252' or cluster_info['Cluster'][i] == 'King_5' or cluster_info['Cluster'][i] == 'NGC_1907' or cluster_info['Cluster'][i] == 'NGC_2509':
            if total_probs2[j] >= 0.3:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
        if cluster_info['Cluster'][i] == 'NGC_7789' or cluster_info['Cluster'][i] == 'NGC_2477' or cluster_info['Cluster'][i] == 'NGC_6208' or cluster_info['Cluster'][i] == 'NGC_6940' or cluster_info['Cluster'][i] == 'Melotte_66':
            if total_probs2[j] >= 0.7:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
        if cluster_info['Cluster'][i] == 'NGC_6791':
            if total_probs2[j] >= 0.9:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
        else:
            if total_probs2[j] >= 0.6:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
    corr_bprp_filt = np.array([x - E_BP_RP for x in bprp_filt])
    corr_g_filt = np.array([y + g_corr for y in g_filt])
    corr_colors = np.zeros(len(corr_bprp_filt), dtype={'names':('corr_bprp_filt', 'corr_g_filt'),
                              'formats':('f8', 'f8')})
    corr_colors['corr_bprp_filt'] = corr_bprp_filt
    corr_colors['corr_g_filt'] = corr_g_filt
    corr_colors = np.sort(corr_colors, order='corr_g_filt')

    #Apparent magnitudes of bp-rp and g limits of red clump
    bprp_lim = [[1.4265,1.66],[1.80,1.968],[1.3564,1.5858],[1.1312,1.1951],
                [1.2372,1.33425],[1.481,1.62],[1.42,1.56],[1.38172,1.44891],
                [1.2258,1.25855],[1.52,1.672],[1.2,1.2755],[1.0166,1.2765],
                [1.06,1.21185],[1.38,1.53],[1.045,1.07752],[1.1525,1.17],
                [1.5526,1.8224],[1.4285,1.7286],[1.31,1.428],[1.1283,1.2864],
                [1.24497,1.34523],[1.9,2.0],[1.43,1.5],[2.626,2.8683],[1.53,1.9]]
    
    g_lim = [[12.0,14.3],[12.744,14.0],[12.08,13.135],[8.6,8.9],[10.336,11.434],
             [12.2,13.0],[13.75,14.4],[11.1232,11.4155],[10.14,10.25],[13.65,13.84],
             [11.45,11.8117],[12.5,13.214],[12.7,13.15],[11.23,12.426],[9.8716,10.1058],
             [12.26,12.4],[12.8258,13.7037],[11.218,12.305],[13.015,13.949],[12.043,12.7218],
             [13.8645,14.2108],[13.4,14.4],[13.41,13.5661],[13.8,14.4],[11.86,13.726]]

    #Find location of red clump
    rc_bprp = []
    rc_g = []
    #collect all the stars in the red clump by taking the stars within the bp-rp and g limits above
    for j in range(len(bprp_filt)):
        if bprp_lim[i][0] <= bprp_filt[j] <= bprp_lim[i][1] and g_lim[i][0] <= g_filt[j] <= g_lim[i][1]:
            rc_bprp.append(bprp_filt[j])
            rc_g.append(g_filt[j])
    rc_numbers.append(len(rc_bprp))
    
    #use astropy statistics biweight_location to find the central location of the sample
    #if the median absolute deviation is zero, the biweight location is simply the median
    med_rc_bprp = round(biweight_location(rc_bprp) - E_BP_RP,4) #correct for reddening
    med_rc_g = round(biweight_location(rc_g) + g_corr,4) #correct for distance and reddening

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
    
    #use line euqation from above to find median rgb location at height of rc
    med_rgb_bprp = f(rgb_point[i][0],rgb_point[i][1],med_rc_g)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    
    #plot cmd
    im = scatter(corr_colors['corr_bprp_filt'],corr_colors['corr_g_filt'], s=5, c=total_probs2_filt, cmap=cm.plasma_r)
    #plot rgb and rc median locations with big red stars
    scatter(med_rgb_bprp, med_rc_g, marker="*", color='#FF2222', s=20)
    scatter(med_rc_bprp, med_rc_g, marker="*", color='#FF2222', s=20)
    #plot rgb line segment (should intersect rgb star)
    x_rgb = np.linspace(0.8,(med_rgb_bprp+.3),80)
    y_rgb = -8.367*(x_rgb - med_rgb_bprp) + med_rc_g
    plot(x_rgb,y_rgb,color='red',alpha=0.4)
    
    import os
    os.chdir('/Users/abbychriss/Desktop')

    ylim([4,-4])
    xlim([0,2.5])
    title('Color Magnitude Diagram: ' + cluster.replace('_',' '))
    xlabel(r'$BP-RP$ (mag)',fontsize=25)
    ylabel(r'$G$ (mag)',fontsize=25)
    yticks(fontsize=25)
    xticks(fontsize=25)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = colorbar(im, cax=cax)
    cax.tick_params(labelsize=25)
    cb.set_label(label='Prob', size=25)
    #show()
    savefig(cluster_info['Cluster'][i]+'_cmd.pdf', bbox_inches='tight', format='pdf')
    close()
    
    """klicker = clicker(ax, ["event"], markers=["x"])
    
    print(klicker.get_positions()["event"])"""