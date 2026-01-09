#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:20:33 2024

@author: abbychriss
"""

import requests
import numpy as np
from astropy.table import Table

cluster_names = np.array(['Collinder_110', 'King_5', 'Melotte_71', 'NGC_1245', 'NGC_1907', 'NGC_2420', 'NGC_2506', 'NGC_2509', 'NGC_2627', 'NGC_2682', 'NGC_6791', 'NGC_6939', 'NGC_6940', 'NGC_752', 'NGC_7789'])

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.isin(np.array(cluster_info['Cluster']), cluster_names)]

"""for i in range(len(cluster_names)):
    requests.post('https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean/crossmatch/upload')
    file=selected_cluster_info['Cluster'][i]+'_cross_ref.csv'
    header=True
    resolve=True
    radius=0.0008
    ra_name=ra
    dec_name=dec
    target_name=target"""

file_names = ['NGC_7789_cross_ref.csv',
'NGC_752_cross_ref.csv',
'NGC_6940_cross_ref.csv',
'NGC_6939_cross_ref.csv',
'NGC_6791_cross_ref.csv',
'NGC_2682_cross_ref.csv',
'NGC_2627_cross_ref.csv',
'NGC_2509_cross_ref.csv',
'NGC_2506_cross_ref.csv',
'NGC_2420_cross_ref.csv',
'NGC_1907_cross_ref.csv',
'NGC_1245_cross_ref.csv',
'Melotte_71_cross_ref.csv',
'King_5_cross_ref.csv',
'Collinder_110_cross_ref.csv']

for i in range(len(cluster_names)):
    print("curl -F radius=0.000833 -F"+" 'file=@{}'".format(file_names[i])+" \ -X POST https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean/crossmatch/upload.csv >" + ' {}.'.format(file_names[i]).replace('csv.','ps1.csv'),'\n')