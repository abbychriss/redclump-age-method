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
#print(nodup_2020)

#import table1 from same paper
table1_2020 = Table.read('/Users/abbychriss/Desktop/WSU/table1_2020.dat.txt', format='ascii')

headers = ['Cluster', 'RAdeg', 'DEdeg', 'GLON', 'GLAT', 'r50', 'nbstars07', 'pmRA*', 'e_pmRA*', 'pmDE', 'e_pmDE', 'plx', 'e_plx', 'Flag', 'AgeNN', 'AVNN', 'DMNN', 'DistPc', 'X', 'Y', 'Z', 'Rgc']
#for explanations on headers for 2020 data, see ReadMe: https://cdsarc.cds.unistra.fr/viz-bin/ReadMe/J/A+A/640/A1?format=html&tex=true

for i in range(len(headers)):
    table1_2020.rename_column('col'+str(i+1), headers[i])
    
print("\n",table1_2020)

print("\n",table1_2020['Cluster','AgeNN'])

#failed AGB stars: - need turn off masses, most massive star remaining in cluster ~1.25 M_s to 5 M_s translates to few hundred mil yrs to 2 billion
# take between 8.3 and 9.3 AgeNN

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
        print('Cluster ' + cluster_names[row] + ' does not have an age.')
        age_f.append(0.)
        ageless += 1
    else:
        age_f.append(float(age_s[row]))
print('Number of ageless stars = ' + str(ageless))

age_f = np.array(age_f)

#filter clusters to include only those in our age range of interest
#from 8.3 to 9.3 in log years
new_ages = []
cluster_indeces_t = []
for row in range(len(age_f)):
    if 8.3 <= age_f[row] <= 9.3:
        new_ages.append(age_f[row])
        cluster_indeces_t.append(row)
        
#make new lists for certain parameters for filtered clusters      
new_clusters_t = [cluster_names[i] for i in cluster_indeces_t]
new_nbstars07_t = [np.array(table1_2020['nbstars07'])[i] for i in cluster_indeces_t]
new_RAdeg_t = [np.array(table1_2020['RAdeg'])[i] for i in cluster_indeces_t]
new_DEdeg_t = [np.array(table1_2020['DEdeg'])[i] for i in cluster_indeces_t]
new_pmRA_t = [np.array(table1_2020['pmRA*'])[i] for i in cluster_indeces_t]
new_pmDE_t = [np.array(table1_2020['pmDE'])[i] for i in cluster_indeces_t]
new_plx_t = [np.array(table1_2020['plx'])[i] for i in cluster_indeces_t]

#Write ASCII table stored in file called 'clusterages.dat'
filtered_clusters_t = Table()
filtered_clusters_t['Cluster'] = np.array(new_clusters_t, dtype=str)
filtered_clusters_t['AgeNN'] = np.array(new_ages, dtype=float)
filtered_clusters_t['nbstars07'] = np.array(new_nbstars07_t)
filtered_clusters_t['RAdeg'] = np.array(new_RAdeg_t)
filtered_clusters_t['DEdeg'] = np.array(new_DEdeg_t)
filtered_clusters_t['pmRA*'] = np.array(new_pmRA_t)
filtered_clusters_t['pmDE'] = np.array(new_pmDE_t)
filtered_clusters_t['plx'] = np.array(new_plx_t)
ascii.write(filtered_clusters_t, 'clusterages.dat', overwrite=True)

#print(Table.read('/Users/abbychriss/Desktop/WSU/clusterages.dat', format='ascii'))

#from Cantat-Gaudin, 2020: the uncertainty on the determination of log t ranges from 0.15 to 0.25 for young clusters and from 0.1 to 0.2 for old clusters.

#how many clusters, examine all robustly, maybe lists are ok for quick look - figure out if they did or did not use colors and magnitudes, if they did not use then that would be good news

#now eliminate all clusters whose distance is > 2kpc = 2000 pc
#convert strings in plx column to floats
plx_s = np.array(table1_2020['plx'])

plx_f = []
for row in range(len(plx_s)):
    plx_f.append(float(plx_s[row]))

plx_f = np.array(plx_f)

#filter clusters to include only those in our age range of interest: plx > 0.5 
new_plx = []
cluster_indeces_plx = []
for row in range(len(plx_f)):
    if plx_f[row] > 0.5:
        new_plx.append(plx_f[row])
        cluster_indeces_plx.append(row)

#make new lists for certain parameters for clusters filtered with respect to plx > 0.5 (equiv, dist < 2 kpc)   
new_ages_plx = [age_f[i] for i in cluster_indeces_plx]
new_clusters_plx = [cluster_names[i] for i in cluster_indeces_plx]
new_nbstars07_plx = [np.array(table1_2020['nbstars07'])[i] for i in cluster_indeces_plx]
new_RAdeg_plx = [np.array(table1_2020['RAdeg'])[i] for i in cluster_indeces_plx]
new_DEdeg_plx = [np.array(table1_2020['DEdeg'])[i] for i in cluster_indeces_plx]
new_pmRA_plx = [np.array(table1_2020['pmRA*'])[i] for i in cluster_indeces_plx]
new_pmDE_plx = [np.array(table1_2020['pmDE'])[i] for i in cluster_indeces_plx]
new_plx_plx = [np.array(table1_2020['plx'])[i] for i in cluster_indeces_plx]

#Write ASCII table for clusters with average distance < 2kpc stored in file called 'clusterplx.dat'
filtered_clusters_plx = Table()
filtered_clusters_plx['Cluster'] = np.array(new_clusters_plx, dtype=str)
filtered_clusters_plx['AgeNN'] = np.array(new_ages_plx, dtype=float)
filtered_clusters_plx['nbstars07'] = np.array(new_nbstars07_plx)
filtered_clusters_plx['RAdeg'] = np.array(new_RAdeg_plx)
filtered_clusters_plx['DEdeg'] = np.array(new_DEdeg_plx)
filtered_clusters_plx['pmRA*'] = np.array(new_pmRA_plx)
filtered_clusters_plx['pmDE'] = np.array(new_pmDE_plx)
filtered_clusters_plx['plx'] = np.array(new_plx_plx)
ascii.write(filtered_clusters_plx, 'clusterplx.dat', overwrite=True)

print(Table.read('/Users/abbychriss/Desktop/WSU/clusterplx.dat', format='ascii'))

print('\ndone :)')





