#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 20:30:41 2024

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

cluster_names = np.array(['Collinder_110', 'King_5', 'Melotte_71', 'NGC_1245', 'NGC_1907', 'NGC_2420', 'NGC_2506', 'NGC_2509', 'NGC_2627', 'NGC_2682', 'NGC_6791', 'NGC_6939', 'NGC_6940', 'NGC_752', 'NGC_7789'])

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.isin(np.array(cluster_info['Cluster']), cluster_names)]

for i in range(len(cluster_names)):
    csv_file = Table.read(selected_cluster_info['Input file name'][i], format='csv')
    csv_file.rename_column('RA','ra')
    csv_file.rename_column('DEC','dec')
    csv_file.rename_column('Gaia ID', 'target')
    csv_file.keep_columns(['target','ra','dec'])
    csv_file.write('/Users/abbychriss/Desktop/WSU/'+selected_cluster_info['Cluster'][i]+'_cross_ref.csv',overwrite=True)
    