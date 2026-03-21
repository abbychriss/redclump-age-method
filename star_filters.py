#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 16:34:08 2023

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

#import table of filtered clusters (8.3 <= log yrs <= 9.3)
filtered_clusters_t = Table.read('/Users/abbychriss/Desktop/WSU/clusterages.dat', format='ascii')

filtered_clusters_plx_t = Table.read('/Users/abbychriss/Desktop/WSU/clusterplx_ages.dat', format='ascii')

filtered_clusters_plx_t_2 = Table.read('/Users/abbychriss/Desktop/WSU/clusterplx_ages_2.dat', format='ascii')

#import table of members

#members_2020 = Table.read('/Users/abbychriss/Desktop/WSU/nodup_2020.dat', format='csv', guess=False, delimiter= '\s') #names = ['RAdeg', 'DEdeg', 'GaiaDR2', 'GLON', 'GLAT', 'Plx', 'e_Plx', 'pmRA*', 'e_pmRA*', 'pmDE', 'e_pmDE', 'RADEcor', 'RAPlxcor', 'RApmRAcor', 'RApmDEcor', 'DEPlxcor', 'DEpmRAcor', 'DEpmDEcor', 'PlxpmRAcor', 'PlxpmDEcor', 'pmRApmDEcor', 'o_Gmag', 'Gmag', 'BP-RP', 'proba', 'Cluster', 'Teff50'])

#members_2020 = np.genfromtxt('/Users/abbychriss/Desktop/WSU/nodup_2020.dat',
                     #delimiter = (21,22,20,21,23,23,23,21,22,22,22,22,25,15,15,15,15,15,15,15,15,15,15,5,11,14,4,23,14), 
                     #usecols = range(1,29) )#, skip_footer = 1)
                     
members_2020 = ascii.read('/Users/abbychriss/Desktop/WSU/nodup_2020.dat', format='fixed_width_no_header',
                names = ('RAdeg', 'DEdeg', 'GaiaDR2', 'GLON', 'GLAT', 'Plx', 'e_Plx', 'pmRA*', 'e_pmRA*', 'pmDE', 'e_pmDE', 'RV', 'e_RV', 'RADEcor', 'RAPlxcor', 'RApmRAcor', 'RApmDEcor', 'DEPlxcor', 'DEpmRAcor', 'DEpmDEcor', 'PlxpmRAcor', 'PlxpmDEcor', 'pmRApmDEcor', 'o_Gmag', 'Gmag', 'BP-RP', 'proba', 'Cluster', 'Teff50'),
                col_starts=(0, 21, 43, 63, 84, 107, 130, 153, 174, 196, 218, 240, 262, 287, 302, 317, 332, 347, 362, 377, 392, 407, 422, 437, 442, 453, 467, 471, 494))

print(members_2020.info)

print(members_2020)

"""
#use 2018 members table because the 2020 table was not entered properly
members_2018 = Table.read('/Users/abbychriss/Desktop/WSU/members_2018.dat', format= 'ascii')

names=['RAdeg', 'DEdeg', 'Source', 'GLON', 'GLAT', 'plx', 'e_plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'RADEcor', 'RAPlxcor', 'RApmRAcor', 'RApmDEcor', 'DEPlxcor', 'DEpmRAcor', 'DEpmDEcor', 'PlxpmRAcor', 'PlxpmDEcor', 'pmRApmDEcor', 'o_Gmag', 'Gmag', 'BP-RP', 'PMemb', 'Cluster']
for i in range(len(names)):
    members_2018.rename_column('col'+str(i+1), names[i]) """

#want to write a program that runs through members and keeps only the rows that contain stars from clusters in filtered_clusters_t and filtered_clusters_plx_t
#member_clusters_2018 = members_2018['Cluster']
member_clusters_2020 = members_2020['Cluster']
print(member_clusters_2020)

"""
#count how many stars remain if we restrict age to 8.3 <= log t <= 9.3
member_indexes_t = []
for row in range(len(member_clusters)):
    if member_clusters[row] in filtered_clusters_t['Cluster']:
        member_indexes_t.append(row)
print('\nStars remaining from age filter= ' + str(len(member_indeces_t)))

#count how many stars remain if we restrict age to 8.3 <= log t <= 9.3 AND plx > 0.5
member_indexes_plx_t = []
for row in range(len(member_clusters_2020)):
    if member_clusters_2020[row] in filtered_clusters_plx_t['Cluster']:
        member_indexes_plx_t.append(row)
print('\nStars remaining from age and parallax filter= ' + str(len(member_indexes_plx_t)))
"""

#count how many stars remain if we restrict age to 8.3 <= log t <= 9.3 and plx > 0.25
member_indexes_plx_t_2 = []
for row in range(len(member_clusters_2020)):
    if member_clusters_2020[row] in filtered_clusters_plx_t_2['Cluster']:
        member_indexes_plx_t_2.append(row)
print('\nStars remaining from age and parallax filter= ' + str(len(member_indexes_plx_t_2)))

"""
#use the age and plx filter to create a new table containing stars from filtered_clusters_plx_t

#make new lists for certain parameters for stars filtered for 8.3 <= log t <= 9.3 and plx > 0.5
star_GaiaDR2_plx_t = [np.array(members_2020['GaiaDR2'][i]) for i in member_indexes_plx_t]
star_clusters_plx_t = [np.array(members_2020['Cluster'][i]) for i in member_indexes_plx_t]
#star_GLON_plx_t = [np.array(members_2020['GLON'][i]) for i in member_indexes_plx_t]
#star_GLAT_plx_t = [np.array(members_2020['GLAT'][i]) for i in member_indexes_plx_t]
star_proba_plx_t = [np.array(members_2020['proba'][i]) for i in member_indexes_plx_t]
star_RAdeg_plx_t = [np.array(members_2020['RAdeg'][i]) for i in member_indexes_plx_t]
star_DEdeg_plx_t = [np.array(members_2020['DEdeg'][i]) for i in member_indexes_plx_t]
star_pmRA_plx_t = [np.array(members_2020['pmRA*'][i]) for i in member_indexes_plx_t]
#star_e_pmRA_plx_t = [np.array(members_2020['e_pmRA'])[i] for i in member_indexes_plx_t]
star_pmDE_plx_t = [np.array(members_2020['pmDE'][i]) for i in member_indexes_plx_t]
#star_e_pmDE_plx_t = [np.array(members_2020['e_pmDE'])[i] for i in member_indexes_plx_t]
star_plx_plx_t = [np.array(members_2020['Plx'][i]) for i in member_indexes_plx_t]
#star_e_plx_plx_t = [np.array(members_2020['e_plx'][i]) for i in member_indexes_plx_t]
star_BP_RP_plx_t = [np.array(members_2020['BP-RP'][i]) for i in member_indexes_plx_t]
star_Gmag_plx_t = [np.array(members_2020['Gmag'][i]) for i in member_indexes_plx_t]

#Write ASCII table for clusters with average distance < 2kpc and between 8.3 and 9.3 log yrs
#stored in file called 'filtered_stars_plx_t.dat'
filtered_stars_plx_t = Table()
filtered_stars_plx_t['GaiaDR2'] = np.array(star_GaiaDR2_plx_t)
filtered_stars_plx_t['Cluster'] = np.array(star_clusters_plx_t, dtype=str)
#filtered_stars_plx_t['GLON'] = np.array(star_GLON_plx_t)
#filtered_stars_plx_t['GLAT'] = np.array(star_GLAT_plx_t)
filtered_stars_plx_t['proba'] = np.array(star_proba_plx_t)
filtered_stars_plx_t['RAdeg'] = np.array(star_RAdeg_plx_t)
filtered_stars_plx_t['DEdeg'] = np.array(star_DEdeg_plx_t)
filtered_stars_plx_t['pmRA*'] = np.array(star_pmRA_plx_t)
#filtered_stars_plx_t['e_pmRA'] = np.array(star_e_pmRA_plx_t)
filtered_stars_plx_t['pmDE'] = np.array(star_pmDE_plx_t)
#filtered_stars_plx_t['e_pmDE'] = np.array(star_e_pmDE_plx_t)
filtered_stars_plx_t['Plx'] = np.array(star_plx_plx_t)
#filtered_stars_plx_t['e_plx'] = np.array(star_e_plx_plx_t)
filtered_stars_plx_t['BP-RP'] = np.array(star_BP_RP_plx_t)
filtered_stars_plx_t['Gmag'] = np.array(star_Gmag_plx_t)
ascii.write(filtered_stars_plx_t, 'filtered_stars_plx_t_2.dat', overwrite=True)

filtered_stars_plx_t = Table.read('/Users/abbychriss/Desktop/WSU/filtered_stars_plx_t.dat', format='ascii')
"""

#Write ASCII table for clusters with average distance < 4kpc and between 8.3 and 9.3 log yrs
#stored in file called 'filtered_stars_plx_t_2.dat'

star_GaiaDR2_plx_t_2 = [np.array(members_2020['GaiaDR2'][i]) for i in member_indexes_plx_t_2]
star_clusters_plx_t_2 = [np.array(members_2020['Cluster'][i]) for i in member_indexes_plx_t_2]
#star_GLON_plx_t_2 = [np.array(members_2020['GLON'][i]) for i in member_indexes_plx_t_2]
#star_GLAT_plx_t_2 = [np.array(members_2020['GLAT'][i]) for i in member_indexes_plx_t_2]
star_proba_plx_t_2 = [np.array(members_2020['proba'][i]) for i in member_indexes_plx_t_2]
star_RAdeg_plx_t_2 = [np.array(members_2020['RAdeg'][i]) for i in member_indexes_plx_t_2]
star_DEdeg_plx_t_2 = [np.array(members_2020['DEdeg'][i]) for i in member_indexes_plx_t_2]
star_pmRA_plx_t_2 = [np.array(members_2020['pmRA*'][i]) for i in member_indexes_plx_t_2]
#star_e_pmRA_plx_t_2 = [np.array(members_2020['e_pmRA'])[i] for i in member_indexes_plx_t_2]
star_pmDE_plx_t_2 = [np.array(members_2020['pmDE'][i]) for i in member_indexes_plx_t_2]
#star_e_pmDE_plx_t_2 = [np.array(members_2020['e_pmDE'])[i] for i in member_indexes_plx_t_2]
star_plx_plx_t_2 = [np.array(members_2020['Plx'][i]) for i in member_indexes_plx_t_2]
#star_e_plx_plx_t_2 = [np.array(members_2020['e_plx'][i]) for i in member_indexes_plx_t_2]
star_BP_RP_plx_t_2 = [np.array(members_2020['BP-RP'][i]) for i in member_indexes_plx_t_2]
star_Gmag_plx_t_2 = [np.array(members_2020['Gmag'][i]) for i in member_indexes_plx_t_2]

filtered_stars_plx_t_2 = Table()
filtered_stars_plx_t_2['GaiaDR2'] = star_GaiaDR2_plx_t_2
filtered_stars_plx_t_2['Cluster'] = star_clusters_plx_t_2
#filtered_stars_plx_t_2['GLON'] = star_GLON_plx_t_2
#filtered_stars_plx_t_2['GLAT'] = star_GLAT_plx_t_2 
filtered_stars_plx_t_2['proba'] = star_proba_plx_t_2
filtered_stars_plx_t_2['RAdeg'] = star_RAdeg_plx_t_2
filtered_stars_plx_t_2['DEdeg'] = star_DEdeg_plx_t_2
filtered_stars_plx_t_2['pmRA*'] = star_pmRA_plx_t_2
#filtered_stars_plx_t_2['e_pmRA'] = star_e_pmRA_plx_t_2
filtered_stars_plx_t_2['pmDE'] = star_pmDE_plx_t_2
#filtered_stars_plx_t_2['e_pmDE'] = star_e_pmDE_plx_t_2
filtered_stars_plx_t_2['Plx'] = star_plx_plx_t_2
#filtered_stars_plx_t_2['e_plx'] = star_e_plx_plx_t_2
filtered_stars_plx_t_2['BP-RP'] = star_BP_RP_plx_t_2
filtered_stars_plx_t_2['Gmag'] = star_Gmag_plx_t_2
ascii.write(filtered_stars_plx_t_2, '/Users/abbychriss/Desktop/WSU/filtered_stars_plx_t_2.dat', overwrite=True)

filtered_stars_plx_t_2 = Table.read('/Users/abbychriss/Desktop/WSU/filtered_stars_plx_t_2.dat', format='ascii')

print(filtered_clusters_plx_t_2.group_by('Cluster'))

print(filtered_stars_plx_t_2['Cluster','Plx'][np.where(filtered_stars_plx_t_2['Cluster']=='Melotte_25')[0]])

print('\ndone :)')

#y axis = G,