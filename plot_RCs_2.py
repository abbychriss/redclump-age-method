#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 15:46:20 2023

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

#Isochrone function:

# for convenience, get the list of ages and metallicities (the latter all zero
# in this case.)
a,feh=np.loadtxt('/Users/abbychriss/Desktop/WSU/pops.solar11',usecols=(0,1),unpack=True)
#print(len(a),len(feh))

def fetchiso(iiso):
    if iiso > 65:
        print('Error! iiso must be 65 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('/Users/abbychriss/Desktop/WSU/hess_gaia_ev6res0mixss.out','r') as filehandle:
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


#Magnitude correcting parameters and functions:
    
#open table of cluster parameters (cluster name and E(B-V)) from Kharchenko N. V.,
#Piskunov A. E., Schilbach E., Röser S., Scholz R. D., 2013, A&A, 558, A53
cluster_parameters_np = ascii.read('/Users/abbychriss/Desktop/WSU/cluster_parameters_Kharchenko_2013.txt',
                                   format='fixed_width_no_header', 
                                   names=('MWSC', 'Cluster', 'Type', 'n_Type',
                                          'RAhour', 'DEdeg', 'GLON', 'GLAT',
                                          'r0', 'r1', 'r2', 'pmRA', 'pmDE',
                                          'e_pm', 'RV', 'e_RV', 'o_RV', 'N1sr0',
                                          'N1sr1', 'N1sr2', 'd', 'E(B-V)', 'MOD',
                                          'E(J-Ks)', 'E(J-H)', 'dH', 'logt',
                                          'e_logt', 'Nt', 'rc', 'e_rc', 'rt',
                                          'e_rt', 'k', 'e_k', 'Src', 'SType',
                                          '[Fe/H]', 'e_[Fe/H]', 'o_[Fe/H]'), 
                                   col_starts=(0, 5, 23, 24, 27, 35, 43, 51, 60,
                                               67, 74, 81, 88, 96, 101, 111, 120,
                                               126, 132, 139, 144, 152, 158, 166,
                                               173, 179, 187, 194, 201, 207, 215,
                                               222, 230, 237, 246, 252, 257, 263,
                                               271, 278))


cluster_parameters = Table()
cluster_parameters['Cluster'] = cluster_parameters_np['Cluster']
cluster_parameters['d'] = cluster_parameters_np['d']
cluster_parameters['E(B-V)'] = cluster_parameters_np['E(B-V)']
cluster_parameters = ascii.write(cluster_parameters, 'cluster_parameters.dat', overwrite=True)
cluster_parameters = Table.read('cluster_parameters.dat', format='ascii')

#open table of stars in clusters selected for age (8.3 <= log t <= 9.3) and distance (< 4 kpc) from Cantat-Gaudin 2020
filtered_stars_plx_t_2 = Table.read('/Users/abbychriss/Desktop/WSU/filtered_stars_plx_t_2.dat', format='ascii')#, names=['Cluster','proba','RAdeg','DEdeg','pmRA*','pmDE','Plx','BP-RP','Gmag'])
filtered_stars_plx_t_2 = filtered_stars_plx_t_2.group_by('Cluster')

RC_clusters = ['IC_2714', 'IC_4651', 'IC_4756', 'NGC_2099', 'NGC_2447', 'NGC_2477',
               'NGC_2627', 'NGC_6208', 'NGC_6939', 'NGC_6940', 'NGC_752', 'BH_211',
               'Berkeley_9', 'COIN-Gaia_6', 'Collinder_110', 'Collinder_277',
               'Czernik_37', 'FSR_0893', 'FSR_0942', 'FSR_1252', 'FSR_1663',
               'King_23', 'King_5', 'LP_1540', 'LP_947', 'Melotte_71', 'NGC_1245',
               'NGC_1907', 'NGC_2236', 'NGC_2420', 'NGC_2432', 'NGC_2489', 'NGC_2506',
               'NGC_2509', 'NGC_2627', 'NGC_2660', 'NGC_3496', 'NGC_3532', 'NGC_5381',
               'NGC_559', 'NGC_7245', 'NGC_7789', 'Pismis_18', 'Ruprecht_101', 'Ruprecht_112',
               'Ruprecht_68', 'Skiff_J0058+68.4', 'Trumpler_23', 'UBC_200']

RC_star_indexes = []
for row in range(len(filtered_stars_plx_t_2['Cluster'])):
    if filtered_stars_plx_t_2['Cluster'][row] in RC_clusters:
        RC_star_indexes.append(row)

RC_candidates = np.zeros(22775, dtype={'names':('Cluster', 'proba', 'RAdeg', 'DEdeg', 'pmRA*', 'pmDE', 'Plx', 'BP-RP', 'Gmag'),
                          'formats':('U18', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
RC_candidates['Cluster'] = [np.array(filtered_stars_plx_t_2['Cluster'])[i] for i in RC_star_indexes]
RC_candidates['proba'] = [np.array(filtered_stars_plx_t_2['proba'])[i] for i in RC_star_indexes]
RC_candidates['RAdeg'] = [np.array(filtered_stars_plx_t_2['RAdeg'])[i] for i in RC_star_indexes]
RC_candidates['DEdeg'] = [np.array(filtered_stars_plx_t_2['DEdeg'])[i] for i in RC_star_indexes]
RC_candidates['pmRA*'] = [np.array(filtered_stars_plx_t_2['pmRA*'])[i] for i in RC_star_indexes]
RC_candidates['pmDE'] =[np.array(filtered_stars_plx_t_2['pmDE'])[i] for i in RC_star_indexes]
RC_candidates['Plx'] = [np.array(filtered_stars_plx_t_2['Plx'])[i] for i in RC_star_indexes]
RC_candidates['BP-RP'] = [np.array(filtered_stars_plx_t_2['BP-RP'])[i] for i in RC_star_indexes]
RC_candidates['Gmag'] = [np.array(filtered_stars_plx_t_2['Gmag'])[i] for i in RC_star_indexes]

print(RC_candidates)

#ESO 518 03 also known as UBC 327 and MWSC 2462
#LP 145 ([LP2019] 145 in simbad) also known as UBC 551 and UFMG 5
#LP 1994 ([LP2019] 1994 in simbad) also known as UBC 288
#LP 2117 ([LP2019] 2117 in simbad) also known as UBC 111
#LP 658 ([LP2019] 658 in simbad) also known as UBC 211
#NGC 2360 also known as [NGC Cl] Melotte 64, [KC2019] Theia 1473, [KPS2012] MWSC 1165, [C] 0715-155, [NAME] Caroline's Cluster, [KPR2004b] 132
#NGC 6991 also known as [KC2019] Theia 1231 and [KPR2004b] 494
#UBC 284 has no other names in simbad
#UPK 549 also known as UBC 256 and [LP2019] 2236

def check_cluster(name):
    if name in cluster_parameters['Cluster']:
        return (name + ' is in cluster_parameters at index ' + str(np.where(cluster_parameters['Cluster']==name)[0][0]))
        return True
    else:
        return False

for cluster in RC_clusters:
    if check_cluster(cluster) == False:
        RC_clusters.remove(cluster)

#we will have to do without those eight clusters, because they don't seem to be in cluster_parameters under any of their aliases

#covert color magnitudes from BP-RP to B-V to factor in reddening
#Casagrande L., VandenBerg D. A., 2018, Monthly Notices of the Royal Astronomical Society: Letters, 479, L102
#give the conversions from E(B-V) to A and E(BP-RP)
#use distance modulus formula m - M = 5 log d -5 + A ==> M = m + 5 - 5 log d - A
E_BP_RP = []
Gmag_correction = []
new_clusters = []
for row in range(len(cluster_parameters['Cluster'])):
    if cluster_parameters['Cluster'][row] in RC_clusters:
        Gmag_correction.append(5 - (5 * np.log10(float(cluster_parameters['d'][row]))) - (2.740 * float(cluster_parameters['E(B-V)'][row])))
        E_BP_RP.append(1.339 * float(cluster_parameters['E(B-V)'][row]))
        new_clusters.append(cluster_parameters['Cluster'][row])

#create a structured numpy array from the above lists of magnitude corrections 
corrections = np.zeros(43, dtype={'names':('Cluster', 'E(BP-RP)', 'Gmag Correction'),
                          'formats':('U18', 'f8', 'f8')})
corrections['Cluster'] = new_clusters
corrections['E(BP-RP)'] = E_BP_RP
corrections['Gmag Correction'] = Gmag_correction
print(corrections)

#now let's plot some color magnitude diagrams and overlay some isochrones.
 
RC_star_indexes = []
for row in range(len(filtered_stars_plx_t_2)):
    if filtered_stars_plx_t_2['Cluster'][row] in corrections['Cluster']:
        RC_star_indexes.append(row)

#make new table with new stars with only the indexes that are in the isochrone star indexes list
RC_stars = np.zeros(22373, dtype={'names':('Cluster', 'Plx', 'proba', 'RAdeg', 'DEdeg', 'BP-RP', 'E(BP-RP)', 'Gmag', 'Gmag Correction'),
                          'formats':('U18', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
RC_stars['Cluster'] = [np.array(filtered_stars_plx_t_2['Cluster'][i]) for i in RC_star_indexes]
RC_stars['Plx'] = [np.array(filtered_stars_plx_t_2['Plx'][i]) for i in RC_star_indexes]
RC_stars['proba'] = [np.array(filtered_stars_plx_t_2['proba'][i]) for i in RC_star_indexes]
RC_stars['RAdeg'] = [np.array(filtered_stars_plx_t_2['RAdeg'][i]) for i in RC_star_indexes]
RC_stars['DEdeg'] = [np.array(filtered_stars_plx_t_2['DEdeg'][i]) for i in RC_star_indexes]
RC_stars['BP-RP'] = [np.array(filtered_stars_plx_t_2['BP-RP'][i]) for i in RC_star_indexes]
RC_stars['Gmag'] = [np.array(filtered_stars_plx_t_2['Gmag'][i]) for i in RC_star_indexes]

RC_stars_e_bp_rp = []
RC_stars_gmag_corr = []
# for every star in RC_stars, in new columns called E(BP-RP) and Gmag correction,
#add the E(BP-RP) and Gmag correction associated with the cluster from the corrections array
for i in range(len(RC_stars['Cluster'])):
    for j in range(len(corrections['Cluster'])):
        if corrections['Cluster'][j]==RC_stars['Cluster'][i]:
            RC_stars_e_bp_rp.append(corrections['E(BP-RP)'][j])
            RC_stars_gmag_corr.append(corrections['Gmag Correction'][j])
            
RC_stars['E(BP-RP)'] = np.array(RC_stars_e_bp_rp)
RC_stars['Gmag Correction'] = np.array(RC_stars_gmag_corr)

#create list of indexes of the first star in each cluster
RC_stop_indexes = []
for i in range(len(RC_stars['Cluster'])):
    if i == 0:
        RC_stop_indexes.append(i)
    elif i == 22372:
        RC_stop_indexes.append(i)
    elif RC_stars['Cluster'][i] != RC_stars['Cluster'][i+1]:
        RC_stop_indexes.append(i+1)

"""#single cluster test plot: Alessi 1
br = []
gm = []
prob = []
for i in range(stop_indexes[0], stop_indexes[1]-1):
    br.append(bp_rp_f[i] - 0.1339)
    gm.append(gmag_f[i] - 0.274 -9.375)
    prob.append(proba[i])
scatter(br, gm, marker='.', c=prob, cmap=cm.plasma_r)
xlabel(r'BP - RP (mag)')
ylabel(r'G (mag)')
xlim([-1.0,6.0])
ylim([20.0,-10.0])
colorbar(extendfrac=None, drawedges=False,label='Probability')
title(isochrone_stars['Cluster'][0].replace('_',' '))
show()"""

#plot CMDs of the clusters in our filtered list (age and distance) over isochrones, subtracting E(BP-RP)
#from BP-RP and Gmag correction from Gmag to plot corrected values
for j in range(len(RC_stop_indexes)): #iterate from 0 to 307
    br = []
    gm = []
    prob = []
    if j != 43:
        for k in range(RC_stop_indexes[j],RC_stop_indexes[j+1]-1):
            br.append(RC_stars['BP-RP'][k] - RC_stars['E(BP-RP)'][k]) #subtract E(BP-RP)
            gm.append(RC_stars['Gmag'][k] + RC_stars['Gmag Correction'][k]) #add Gmag correction = 5 - 5log(d) - A_G since then we are subtracting A_G
            prob.append(RC_stars['proba'][k])
        scatter(br, gm, marker='.', c=prob, cmap=cm.plasma_r)
        title(RC_stars['Cluster'][RC_stop_indexes[j]].replace('_',' '))

        iiso = 35
        teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchiso(iiso)
        br = bmag - rmag
        plot(br,gmag,label='Age '+str(a[iiso])+' Gyr')
        
        iiso = 25
        teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchiso(iiso)
        br = bmag - rmag
        plot(br,gmag,label='Age '+str(a[iiso])+' Gyr')
        
        iiso = 45
        teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchiso(iiso)
        br = bmag - rmag
        plot(br,gmag,label='Age '+str(a[iiso])+' Gyr')
        
        #make list of rich clusters and with at least 3 clump stars
        #also throw out clusters without a red giant branch, should have about a dozen left
        
        xlim([-0.5,2.0])
        ylim([5.0,-5.0])
        legend(frameon=False)
        xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
        ylabel(r'G (mag)')
        colorbar(extendfrac=None, drawedges=False,label='Probability')
        show()
#king model instead of gaussian for RA-DEC
#Ruprecht 68,  Pismis 18, NGC 7789, NGC 752, NGC 6940, NGC 6939, NGC 6208, NGC 2660, NGC 2627, NGC 2509, NGC 2506, NGC 2477, NGC 2447, NGC 2420, NGC 2236, NGC 1907, NGC 1245, Melotte 71 (analysis), King 5, IC 4756, IC 4651, IC 2714, FSR 1252, Czernik 37, Collinder 110      