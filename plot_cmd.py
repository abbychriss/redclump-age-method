#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 15:29:01 2023

@author: abbychriss
"""

import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import math

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from sklearn.ensemble import RandomForestClassifier
import scipy.stats
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

#filtered_stars_plx_t = Table.read('/Users/abbychriss/Desktop/WSU/filtered_stars_plx_t.dat', format='ascii')#, names=['Cluster','proba','RAdeg','DEdeg','pmRA*','pmDE','Plx','BP-RP','Gmag'])

filtered_stars_plx_t_2 = Table.read('/Users/abbychriss/Desktop/WSU/filtered_stars_plx_t_2.dat', format='ascii')#, names=['Cluster','proba','RAdeg','DEdeg','pmRA*','pmDE','Plx','BP-RP','Gmag'])
filtered_stars_plx_t_2 = filtered_stars_plx_t_2.group_by('Cluster')

#open list of probabilities
proba = np.array(filtered_stars_plx_t_2['proba'], dtype=float)

#convert list of numerical strings to floats for BP-RP and Gmag
bp_rp = filtered_stars_plx_t_2['BP-RP']
bp_rp_f = []
for i in range(len(bp_rp)):
    if bp_rp[i] == 'BP-RP':
        pass
    else:
        bp_rp_f.append(float(bp_rp[i]))
        
gmag = filtered_stars_plx_t_2['Gmag']
gmag_f = []
for i in range(len(gmag)):
    if gmag[i] != 'Gmag':
        gmag_f.append(float(gmag[i]))

clusters = np.array(filtered_stars_plx_t_2['Cluster'],dtype=str)

#create list of indexes of the last star in each cluster
stop_indexes = []
for i in range(len(clusters)):
    if i == 0:
        stop_indexes.append(i)
    elif i == 84966:
        stop_indexes.append(i)
    elif clusters[i] != clusters[i+1]:
        stop_indexes.append(i+1)

#test with first cluster
"""
br = []
gm = []
prob = []
for i in range(stop_indexes[0], stop_indexes[1]):
    br.append(bp_rp_f[i])
    gm.append(gmag_f[i])
    prob.append(proba[i])
scatter(br, gm, marker='.', c=prob, cmap=cm.plasma_r)
xlabel(r'BP - RP (mag)')
ylabel(r'G (mag)')
xlim([-1.0,6.0])
ylim([20.0,-10.0])
colorbar(extendfrac=None, drawedges=False,label='Probability')
title(clusters[0].replace('_',' '))
show()
"""

for j in range(len(stop_indexes)): #iterate from 0 to 632
    br = []
    gm = []
    prob = []
    if j != 632:
        for k in range(stop_indexes[j],stop_indexes[j+1]-1):
            br.append(bp_rp_f[k])
            gm.append(gmag_f[k])
            prob.append(proba[k])
        scatter(br, gm, marker='.', c=prob, cmap=cm.plasma_r)
        title(clusters[stop_indexes[j]].replace('_',' '))
        xlabel(r'BP - RP (mag)')
        ylabel(r'G (mag)')
        xlim([-1.0,5])
        ylim([22.0,2.0])
        colorbar(extendfrac=None, drawedges=False,label='Probability')
        show()
        
#make array of cluster names that are rich, have sizable clumps, and RGBs
RC_RGB_clusters = np.array(['ASCC_90', 'ESO_518_03', 'IC_2714', 'IC_4651', 'IC_4756', 
                   'LP_145', 'LP_1994', 'LP_2117', 'LP_658', 'NGC_1817', 'NGC_2099',
                   'NGC_2354',' NGC_2360', 'NGC_2423', 'NGC_2437', 'NGC_2447',
                   'NGC_2477', 'NGC_2548', 'NGC_2627', 'NGC_3532', 'NGC_5822',
                   'NGC_6134', 'NGC_6152', 'NGC_6208', 'NGC_6939', 'NGC_6940',
                   'NGC_6991', 'NGC_752', 'Ruprecht_111', 'UBC_284', 'UPK_549'])

RC_cluster_indexes = []
for row in range(len(filtered_stars_plx_t_2['Cluster'])):
    if filtered_stars_plx_t_2['Cluster'][row] in RC_RGB_clusters:
        RC_cluster_indexes.append(row)

RC_candidates = Table()
RC_candidates['Cluster'] = [np.array(filtered_stars_plx_t_2['Cluster'])[i] for i in RC_cluster_indexes]
RC_candidates['proba'] = [np.array(filtered_stars_plx_t_2['proba'])[i] for i in RC_cluster_indexes]
RC_candidates['RAdeg'] = [np.array(filtered_stars_plx_t_2['RAdeg'])[i] for i in RC_cluster_indexes]
RC_candidates['DEdeg'] = [np.array(filtered_stars_plx_t_2['DEdeg'])[i] for i in RC_cluster_indexes]
RC_candidates['pmRA*'] = [np.array(filtered_stars_plx_t_2['pmRA*'])[i] for i in RC_cluster_indexes]
RC_candidates['pmDE'] =[np.array(filtered_stars_plx_t_2['pmDE'])[i] for i in RC_cluster_indexes]
RC_candidates['Plx'] = [np.array(filtered_stars_plx_t_2['Plx'])[i] for i in RC_cluster_indexes]
RC_candidates['BP-RP'] = [np.array(filtered_stars_plx_t_2['BP-RP'])[i] for i in RC_cluster_indexes]
RC_candidates['Gmag'] = [np.array(filtered_stars_plx_t_2['Gmag'])[i] for i in RC_cluster_indexes]
ascii.write(RC_candidates, 'RC_candidates.dat', overwrite=True)

RC_candidates = Table.read('/Users/abbychriss/Desktop/WSU/RC_candidates.dat', format='ascii')

print('\ndone :)')
    
