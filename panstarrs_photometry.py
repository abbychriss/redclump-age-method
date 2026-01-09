#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 12:40:56 2024

@author: abbychriss
"""

import astropy.units as u
from astropy.io import fits
from astropy.table import Table, vstack
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt

cluster_names = np.array(['Collinder_110', 'King_5', 'Melotte_71', 'NGC_1245', 'NGC_1907', 'NGC_2420', 'NGC_2506', 'NGC_2509', 'NGC_2627', 'NGC_2682', 'NGC_6791', 'NGC_6939', 'NGC_6940', 'NGC_752', 'NGC_7789'])

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.where(np.isin(np.array(cluster_info['Cluster']),cluster_names))]
selected_cluster_info.rename_column('Input file name', 'Gaia file')

panstarrs_file = np.array([selected_cluster_info['Gaia file'][i].replace('_cone.csv', '_pan_{j}.csv') for i in range(len(selected_cluster_info['Cluster']))])
selected_cluster_info.add_column(panstarrs_file, name='PanSTARRS file', index=1)


for i in range(len(selected_cluster_info)):
    cluster = selected_cluster_info['Cluster'][i]
    if cluster == 'NGC_6791':
        k=19
    if cluster == 'NGC_7789':
        k=9
    if cluster=='NGC_2682':
        k=7
    if cluster == 'NGC_2506':
        k=6    
    if cluster == 'Collinder_110' or cluster == 'NGC_6940':
        k=5
    if cluster == 'NGC_6939' or cluster=='NGC_1245':
        k=4   
    if cluster == 'Melotte_71' or cluster=='NGC_1245' or cluster=='NGC_2420':
        k=3
    if cluster == 'King_5' or cluster=='NGC_752' or cluster=='NGC_1907' or cluster=='NGC_2509' or cluster=='NGC_2627':
        k=2
    panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][i].format(j=l),format='csv') for l in np.arange(k)+1 ] )
    g_i = panstarrs_info['gMeanPSFMag'] - panstarrs_info['iMeanPSFMag']
    g_mag = panstarrs_info['gMeanPSFMag']
    
    #panstarrs clusters
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    im = plt.scatter(g_i,g_mag, s=4)#, c=darkblue)
    plt.xlim([-1.0,4.0])
    if cluster=='NGC_6791':
        plt.ylim([20,12.0])
    else:
        plt.ylim([18,12.0])   
    plt.title('Color Magnitude Diagram (Pan-STARRS): '+selected_cluster_info['Cluster'][i].replace('_',' '))
    plt.xlabel('g-i (mag)')
    plt.ylabel('g (mag)')
    plt.show()