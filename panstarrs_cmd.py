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

selected_cluster_info = cluster_info[np.isin(np.array(cluster_info['Cluster']), cluster_names)]
selected_cluster_info.rename_column('Input file name', 'Gaia file')

panstarrs_file = np.array([selected_cluster_info['Gaia file'][i].replace('_cone.csv', '_pan_{j}.csv') for i in range(len(selected_cluster_info['Cluster']))])
selected_cluster_info.add_column(panstarrs_file, name='PanSTARRS file', index=1)

print(selected_cluster_info)

for i in range(len(selected_cluster_info)):
    try: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(19)+1 ] )
    except: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(7)+1 ] )
    except: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(6)+1 ] )
    except: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(5)+1 ] )
    except: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(4)+1 ] )
    except: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(3)+1 ] )
    except: panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][-1].format(j=l),format='csv') for l in np.arange(2)+1 ] )
        
    g_i = col110_panstarrs['gMeanPSFMag'] - col110_panstarrs['iMeanPSFMag']
    g_mag = col110_panstarrs['gMeanPSFMag']
    
    print(col110_panstarrs)
    
    #panstarrs clusters
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    im = plt.scatter(g_i,g_mag, s=4)#, c=darkblue)
    plt.xlim([-1.0,4.0])
    plt.ylim([17.5,12.0])   
    plt.title('Color Magnitude Diagram (Pan-STARRS): Collinder 110')#+ selected_cluster_info['Cluster'][i].replace('_',' '))
    plt.xlabel('g-i (mag)')
    plt.ylabel('g (mag)')
    plt.show()