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
a,feh=np.loadtxt('/Users/abbychriss/Desktop/WSU/Isochrones/pops.solar11',usecols=(0,1),unpack=True)
#print(len(a),len(feh))

def fetchiso(iiso):
    if iiso > 65:
        print('Error! iiso must be 65 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('/Users/abbychriss/Desktop/WSU/Isochrones/hess_gaia_ev6res0mixss.out','r') as filehandle:
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
cluster_parameters['[Fe/H]'] = cluster_parameters_np['[Fe/H]']
cluster_parameters = ascii.write(cluster_parameters, 'cluster_parameters.dat', overwrite=True)
cluster_parameters = Table.read('cluster_parameters.dat', format='ascii')

#open table of star parameters for stars in clusters with red clumps and RGBs
#in age (between 200 Myrs and 2 Gyrs) and distance (less than 2 kpc) ranges
RC_candidates = Table.read('/Users/abbychriss/Desktop/WSU/RC_candidates.dat', format='ascii')
        
RC_clusters = ['ASCC_90', 'ESO_518_03', 'IC_2714', 'IC_4651', 'IC_4756', 
                   'LP_145', 'LP_1994', 'LP_2117', 'LP_658', 'NGC_1817', 'NGC_2099',
                   'NGC_2354','NGC_2360', 'NGC_2423', 'NGC_2437', 'NGC_2447',
                   'NGC_2477', 'NGC_2548', 'NGC_2627', 'NGC_3532', 'NGC_5822',
                   'NGC_6134', 'NGC_6152', 'NGC_6208', 'NGC_6939', 'NGC_6940',
                   'NGC_6991', 'NGC_752', 'Ruprecht_111', 'UBC_284', 'UPK_549']

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
corrections = np.zeros(23, dtype={'names':('Cluster', 'E(BP-RP)', 'Gmag Correction'),
                          'formats':('U18', 'f8', 'f8')})
corrections['Cluster'] = new_clusters
corrections['E(BP-RP)'] = E_BP_RP
corrections['Gmag Correction'] = Gmag_correction
print(corrections)

#now let's plot some color magnitude diagrams and overlay some isochrones.
#we will need the table of member stars from Cantat-Gaudin 2020, filtered for age (8.3 to 9.3 log yrs) and distance (< 2kpc)
#also filtered out those that don't have members in cluster_parameters
isochrone_stars = Table.read('/Users/abbychriss/Desktop/WSU/isochrone_stars.dat', format='ascii')
 
RC_star_indexes = []
for row in range(len(isochrone_stars)):
    if isochrone_stars['Cluster'][row] in corrections['Cluster']:
        RC_star_indexes.append(row)

#make new table with new stars with only the indexes that are in the isochrone star indexes list
RC_stars = np.zeros(15693, dtype={'names':('Cluster', 'Plx', 'proba', 'RAdeg', 'DEdeg', 'BP-RP', 'E(BP-RP)', 'Gmag', 'Gmag Correction'),
                          'formats':('U18', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')})
RC_stars['Cluster'] = [np.array(isochrone_stars['Cluster'][i]) for i in RC_star_indexes]
RC_stars['Plx'] = [np.array(isochrone_stars['Plx'][i]) for i in RC_star_indexes]
RC_stars['proba'] = [np.array(isochrone_stars['proba'][i]) for i in RC_star_indexes]
RC_stars['RAdeg'] = [np.array(isochrone_stars['RAdeg'][i]) for i in RC_star_indexes]
RC_stars['DEdeg'] = [np.array(isochrone_stars['DEdeg'][i]) for i in RC_star_indexes]
RC_stars['BP-RP'] = [np.array(isochrone_stars['BP-RP'][i]) for i in RC_star_indexes]
RC_stars['Gmag'] = [np.array(isochrone_stars['Gmag'][i]) for i in RC_star_indexes]

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
    elif i == 15692:
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
    if j != 23:
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

#Yes: NGC 752, NGC 6940, NGC 6939, NGC 6208, NGC 2627, NGC 2447, NGC 2099, IC 4756, IC 4651, IC 2714
#Do your own membership study : NGC 2477

#NGC 7789

