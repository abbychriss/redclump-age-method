#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:41:29 2023

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

#import members.dat file from  T. Cantat-Gaudin et al., 2018, "A Gaia DR2 view of the open cluster population in the Milky Way"
#downloaded from https://cdsarc.cds.unistra.fr/ftp/J/A+A/618/A93/
#members_2018 = Table.read('/Users/abbychriss/Desktop/WSU/members_2018.dat', format= 'csv')
#print(members)

#import table1 from same paper
#table1_2018 = Table.read('/Users/abbychriss/Desktop/WSU/table1_2018.fits')
#for explanations of headers for 2018 data see ReadMe: https://cdsarc.cds.unistra.fr/ftp/J/A+A/618/A93/ReadMe

#import nodup (table of individual stars) from T. Cantat-Gaudin et al., 2020, "Painting a portrait of the Galactic disc with its stellar clusters"
nodup_2020 = Table.read('/Users/abbychriss/Desktop/WSU/nodup_2020.dat', format = 'csv')
print(nodup_2020)

#import table1 from same paper (table of cluster parameters)
table1_2020 = Table.read('/Users/abbychriss/Desktop/WSU/table1_2020.dat.txt', format='ascii')

headers = ['Cluster', 'RAdeg', 'DEdeg', 'GLON', 'GLAT', 'r50', 'nbstars07', 'pmRA*', 'e_pmRA*', 'pmDE', 'e_pmDE', 'plx', 'e_plx', 'Flag', 'AgeNN', 'AVNN', 'DMNN', 'DistPc', 'X', 'Y', 'Z', 'Rgc']
#for explanations on headers for 2020 data, see ReadMe: https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A+A/640/A1?format=html&tex=true

for i in range(len(headers)):
    table1_2020.rename_column('col'+str(i+1), headers[i])

#failed AGB stars: - need turn off masses, most massive star remaining in cluster ~1.25 M_s to 5 M_s translates to few hundred mil yrs to 2 billion
#take between 8.3 and 9.3 AgeNN

#we begin with the 2020 data set to narrow down our list of clusters because it provides decently accurate estimates of ages

#convert strings in column of ages to floats
cluster_names = np.array(table1_2020['Cluster'])
age_s = np.array(table1_2020['AgeNN'])
#with np.printoptions(threshold=np.inf):
    #print(age_s)
age_f = []
ageless = 0
for row in range(len(age_s)):
    if age_s[row] == '---':
        #print('Cluster ' + cluster_names[row] + ' does not have an age.')
        age_f.append(0.)
        ageless += 1
    else:
        age_f.append(float(age_s[row]))
print('\nNumber of ageless stars = ' + str(ageless) + '\n')

age_f = np.array(age_f)

#from Cantat-Gaudin, 2020: the uncertainty on the determination of log t ranges
#from 0.15 to 0.25 for young clusters and from 0.1 to 0.2 for old clusters.

#convert strings in plx column to floats
plx_s = np.array(table1_2020['plx'])

plx_f = []
for row in range(len(plx_s)):
    plx_f.append(float(plx_s[row]))

plx_f = np.array(plx_f)

#filter clusters so that plx > 0.25 (d < 4000 pc) AND 8.3 <= AgeNN <= 9.3 (200 Myrs to 2 Gyrs)

cluster_indexes_plx_t_2 = []
for i in range(len(table1_2020['Cluster'])):
    if 8.3 <= age_f[i] <= 9.3 and plx_f[i] > 0.25:
        cluster_indexes_plx_t_2.append(i)
        
#how many clusters, examine all robustly, maybe lists are ok for quick look
#figure out if they did or did not use colors and magnitudes
#if they did not use then that would be good news
#did not use photometry to determine age or membership

#Write ASCII table for clusters filtered with respect to plx > 0.25 (d < 4 kpc)
#and age 8.3 <= log t <= 9.3 (200 Myrs to 2 Gyrs) stored in file called 'clusterplx_ages.dat'
filtered_clusters_plx_t_2 = Table()
filtered_clusters_plx_t_2['Cluster'] = [np.array(cluster_names[i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['AgeNN'] = [np.array(age_f[i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['nbstars07'] = [np.array(table1_2020['nbstars07'][i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['RAdeg'] = [np.array(table1_2020['RAdeg'][i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['DEdeg'] = [np.array(table1_2020['DEdeg'][i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['pmRA*'] = [np.array(table1_2020['pmRA*'][i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['pmDE'] = [np.array(table1_2020['pmDE'][i]) for i in cluster_indexes_plx_t_2]
filtered_clusters_plx_t_2['plx'] = [np.array(table1_2020['plx'][i]) for i in cluster_indexes_plx_t_2]
ascii.write(filtered_clusters_plx_t_2, 'clusterplx_ages_2.dat', overwrite=True)

filtered_clusters_plx_t_2 = Table.read('/Users/abbychriss/Desktop/WSU/clusterplx_ages_2.dat', format='ascii')

print(filtered_clusters_plx_t_2['Cluster'])

print(filtered_clusters_plx_t_2['Cluster','AgeNN','plx'][np.where(filtered_clusters_plx_t_2['Cluster']=='NGC_7789')[0]])

print('\ndone :)')





