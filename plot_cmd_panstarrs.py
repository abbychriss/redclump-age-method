#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 19:33:28 2024

@author: abbychriss
"""

import numpy as np
from numpy import median
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.table import Table, vstack
from astropy.io import ascii
import matplotlib.cm as cm
from mpl_point_clicker import clicker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
import random
from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['text.usetex'] = True

cluster_names = np.array(['Collinder_110', 'King_5', 'Melotte_71', 'NGC_1245', 'NGC_1907', 'NGC_2420', 'NGC_2506', 'NGC_2509', 'NGC_2627', 'NGC_2682', 'NGC_6791', 'NGC_6939', 'NGC_6940', 'NGC_752', 'NGC_7789'])

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.where(np.isin(np.array(cluster_info['Cluster']),cluster_names))]
selected_cluster_info.rename_column('Input file name', 'Gaia file')

panstarrs_file = np.array([selected_cluster_info['Gaia file'][i].replace('_cone.csv', '_pan_{j}.csv') for i in range(len(selected_cluster_info['Cluster']))])
selected_cluster_info.add_column(panstarrs_file, name='PanSTARRS file', index=1)

plx = selected_cluster_info['Plx']

e_bprps = []
N_80_list = []
fe_h= []
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
        
    #load panstarrs data cross-matched to Gaia data of stars in clusters with greater than 20% membership probability
    panstarrs_info = vstack( [Table.read(selected_cluster_info['PanSTARRS file'][i].format(j=l),format='csv') for l in np.arange(k)+1 ] )
    panstarrs_gi = panstarrs_info['gMeanPSFMag'] - panstarrs_info['iMeanPSFMag']
    panstarrs_g = panstarrs_info['gMeanPSFMag']
    panstarrs_ra = panstarrs_info['ra']
    panstarrs_dec = panstarrs_info['dec']
    
    #load gaia data for stars in each cluster
    import csv
    ID=[];ra=[];dec=[];plx=[];eplx=[];pmra=[];epmra=[];pmdec=[];epmdec=[];g=[];bp=[];rp=[];bprp=[];rv=[];teff=[];ag=[];ebprp=[];vflag=[]
    with open(selected_cluster_info['Gaia file'][i],'r',newline='') as csvfile:
         cr = csv.reader(csvfile)
         ict = 0
         for row in cr:
              if ict == 0:
                   columnnames = row
                   ict = ict + 1
              else:
                   ict = ict + 1
                   ID.append(row[0])
                   try:
                        ra.append(float(row[1]))  # The downside of csv format is
                   except:                        # having to cast the floats
                        ra.append(-999.)          # instead of simply reading them
                   try:                           # in as floats to start with.
                        dec.append(float(row[2]))
                   except:
                        dec.append(-999.)
                   try:
                        plx.append(float(row[3]))
                   except:
                        plx.append(-999.)
                   try:
                        eplx.append(float(row[4]))
                   except:
                        eplx.append(-999.)
                   try:
                        pmra.append(float(row[5]))
                   except:
                        pmra.append(-999.)
                   try:
                        epmra.append(float(row[6]))
                   except:
                        epmra.append(-999.)
                   try:
                        pmdec.append(float(row[7]))
                   except:
                        pmdec.append(-999.)
                   try:
                        epmdec.append(float(row[8]))
                   except:
                        epmdec.append(-999.)
                   try:
                        g.append(float(row[9]))
                   except:
                        g.append(-999.)
                   try:
                        bprp.append(float(row[12]))
                   except:
                        bprp.append(-999.)
    
    ra = np.array(ra); dec=np.array(dec); pmra=np.array(pmra); pmdec=np.array(pmdec)
    g = np.array(g); bprp=np.array(bprp); plx=np.array(plx); ag=np.array(ag); ebprp=np.array(ebprp)
    
    colors = ['']*len(ID)
    
    if cluster == 'NGC_7789':
        for alpha in range(len(ra)):
            if ra[alpha] > 300.00:
                ra[alpha] = ra[alpha] - 360.0000
    
    #define colors for plots:
    green = '#B7F4C7'
    purple = '#C698EF'
    pink = '#F4B7DC'
    orange = '#ffb65c'
    darkblue = '#1C0075'
    lightblue = '#AAD9EB'
    
    #COMMENTING OUT ALL OF THE "OLD PROBABILITY" METHOD EQUATIONS, VARIABLES, ETC
    #Put each equation for each category and each cluster in a structured array
    #Every equation has a 'c' coefficient: that will act as the constant field star background
    #Prob = (f(x) - c) / (f(x))   [ old "or" logic]
    #x will be dpm, dpos, dplx in each equation  
    ppmc=[];ppmf=[]
    pposc=[];pposf=[]
    pplxc=[];pplxf=[]
    total_probs2 = []
    for j in range(len(ra)):
        
        # --------
        # PROPER MOTION
        # --------
      
        dpm = math.sqrt( (pmra[j] - selected_cluster_info['PMra'][i])**2 + (pmdec[j] - selected_cluster_info['PMdec'][i])**2 )
        
        """ppm_eqns = [(12429.746 * np.exp(-dpm*11.399)) / (111.261 + 12429.746 * np.exp(-dpm*11.399)), #Ruprecht 68
        (7889.704 * np.exp(-dpm*10.716))/ (317.237 + 7889.704 * np.exp(-dpm*10.716)), #Pismis 18
        (41242.811 * np.exp(-dpm*7.103)) / (1225.596 + 41242.811 * np.exp(-dpm*7.103)), #NGC 7789
        (820.675* np.exp(-dpm*3.205)) / (15.690 + 820.675* np.exp(-dpm*3.205)), #NGC 752
        (7707.473 * np.exp(-dpm*6.728)) / (588.992 + 7707.473 * np.exp(-dpm*6.728)), #NGC 6940
        (13557.738 * np.exp(-dpm*9.377)) / (99.831 + 13557.738 * np.exp(-dpm*9.377)), #NGC 6939
        (47524.902 * np.exp(-dpm*8.087)) / (130.203 + 47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        (5469.882 * np.exp(-dpm*7.669)) / (682.752 + 5469.882 * np.exp(-dpm*7.669)), #NGC 6208
        (6386.179 * np.exp(-dpm*4.596)) / (14.800+ 6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        (12625.640 * np.exp(-dpm*10.175)) / (155.176 + 12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (3115.683 * np.exp(-dpm*1.959)) / (53.852 + 3115.683 * np.exp(-dpm*1.959)), #NGC 2627
        (17280.461 * np.exp(-dpm*17.437)) / (159.432 + 17280.461 * np.exp(-dpm*17.437)), #NGC 2509
        (33119.650 * np.exp(-dpm*8.828)) / (117.545 + 33119.650 * np.exp(-dpm*8.828)), #NGC 2506
        (21370.314 * np.exp(-dpm*5.392)) / (854.059 + 21370.314 * np.exp(-dpm*5.392)), #NGC 2477
        (8745.948 * np.exp(-dpm*6.120)) / (1355.076 + 8745.948 * np.exp(-dpm*6.120)), #NGC 2447
        (12374.103 * np.exp(-dpm*9.113)) / (120.148 + 12374.103 * np.exp(-dpm*9.113)), #NGC 2420
        (10819.833 * np.exp(-dpm*9.048)) / (428.479 + 10819.833 * np.exp(-dpm*9.048)), #NGC 2236
        (3789.570 * np.exp(-dpm*6.141)) / (123.164 + 3789.570 * np.exp(-dpm*6.141)), #NGC 1907
        (13338.112 * np.exp(-dpm*7.465)) / (389.433 + 13338.112 * np.exp(-dpm*7.465)), #NGC 1245
        (16787.800 * np.exp(-dpm*9.424)) / (136.353 + 16787.800 * np.exp(-dpm*9.424)), #Melotte 71
        (23473.741 * np.exp(-dpm*9.355)) / (139.946 + 23473.741 * np.exp(-dpm*9.355)),#Melotte 66
        (7973.347 * np.exp(-dpm*9.869)) / (609.280 + 7973.347 * np.exp(-dpm*9.869)), #King 5
        (3992.037 * np.exp(-dpm*9.735)) / (143.573 + 3992.037 * np.exp(-dpm*9.735)), #FSR 1252
        (5841.210 * np.exp(-dpm*10.104)) / (321.770 + 5841.210 * np.exp(-dpm*10.104)), #Czernik 37
        (28809.904 * np.exp(-dpm*10.447)) / (955.309 + 28809.904 * np.exp(-dpm*10.447))] #Collinder 110"""
        
        phi_cpm = [#(12429.746 * np.exp(-dpm*11.399)) , #Ruprecht 68
        #(7889.704 * np.exp(-dpm*10.716) ), #Pismis 18
        (41242.811 * np.exp(-dpm*7.103)) , #NGC 7789
        (820.675* np.exp(-dpm*3.205)) , #NGC 752
        (7707.473 * np.exp(-dpm*6.728)) , #NGC 6940
        (13557.738 * np.exp(-dpm*9.377)) , #NGC 6939
        (47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        #(5469.882 * np.exp(-dpm*7.669)) , #NGC 6208
        (6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        #(12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (3115.683 * np.exp(-dpm*1.959)) , #NGC 2627
        (17280.461 * np.exp(-dpm*17.437)) , #NGC 2509
        (33119.650 * np.exp(-dpm*8.828)) , #NGC 2506
        #(21370.314 * np.exp(-dpm*5.392)) , #NGC 2477
        #(8745.948 * np.exp(-dpm*6.120)) , #NGC 2447
        (12374.103 * np.exp(-dpm*9.113)) , #NGC 2420
        #(10819.833 * np.exp(-dpm*9.048)) , #NGC 2236
        (3789.570 * np.exp(-dpm*6.141)) , #NGC 1907
        (13338.112 * np.exp(-dpm*7.465)) , #NGC 1245
        (16787.800 * np.exp(-dpm*9.424)) , #Melotte 71
        #(23473.741 * np.exp(-dpm*9.355)) , #Melotte 66
        (7973.347 * np.exp(-dpm*9.869)) , #King 5
        #(3992.037 * np.exp(-dpm*9.735)) , #FSR 1252
        #(5841.210 * np.exp(-dpm*10.104)) , #Czernik 37
        (28809.904 * np.exp(-dpm*10.447)) ] #Collinder 110 

        phi_fpm = [#111.261,
                   #317.237,
                   1225.596,
                   15.690,
                   588.992,
                   99.831, 
                   130.203,
                   #682.752,
                   14.800,
                   #155.176,
                   53.852,
                   159.432, 
                   117.545,
                   #854.059,
                   #1355.076,
                   120.148,
                   #428.479, 
                   123.164, 
                   389.433,
                   136.353, 
                   #139.946, 
                   609.280, 
                   #143.573, 
                   #321.770, 
                   955.309]
        
        #ppm.append(ppm_eqns[i]) # add a proper motion probability instance
        #ppmc.append(Anorm_pmC*phi_cpm[i]) # normalized cluster probability instance
        #ppmf.append(Anorm_pmF*phi_fpm[i]) # normalized field probability instance
        ppmc.append(phi_cpm[i])
        ppmf.append(phi_fpm[i])

        # --------
        # POSITION
        # --------
        cosdec = math.cos(selected_cluster_info['Dec'][i]*math.pi/180.0) #use dec
        if cluster == 'NGC_7789':
            dpos = math.sqrt( ((ra[j] - (selected_cluster_info['Ra'][i]-360.000))*cosdec)**2 + (dec[j] - selected_cluster_info['Dec'][i])**2 )
        else: dpos = math.sqrt( ((ra[j] - selected_cluster_info['Ra'][i])*cosdec)**2 + (dec[j] - selected_cluster_info['Dec'][i])**2 )
        
        """ppos_eqns = [(30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))) / (73332.225 + 30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))), #Ruprecht 68
                     (170521.883 * np.exp(-((dpos-0.002)**2/(2*0.016**2)))) / (261896.000 + 170521.883 * np.exp(-((dpos-0.002)**2/(2*0.016**2)))), #Pismis 18
                     (55926.096 * np.exp(-((dpos)**2/(2*0.090**2)))) / (66608.971 + 55926.096 * np.exp(-((dpos)**2/(2*0.090**2)))), #NGC 7789
                     (2715.388 * np.exp(-dpos*0.048)) / (2039.455 + 2715.388 * np.exp(-dpos*0.048)), #NGC 752
                     (12646.468 * np.exp(-((dpos)**2/(2*0.280**2)))) / (93792.077 + 12646.468 * np.exp(-((dpos)**2/(2*0.280**2)))), #NGC 6940
                     (51920.447 * np.exp(-dpos*22.830)) / (32503.010 + 51920.447 * np.exp(-dpos*22.830)), #NGC 6939
                     (479976.011 * np.exp(-dpos*21.117)) / (53557.792 + 479976.011 * np.exp(-dpos*21.117)), #NGC 6791
                     (25208.794 * np.exp(-((dpos+0.091)**2/(2*0.1**2)))) / (106234.735 + 25208.794 * np.exp(-((dpos+0.091)**2/(2*0.1**2)))), #NGC 6208
                     (20857.247 * np.exp(-dpos*14.331)) / (6243.826 + 20857.247 * np.exp(-dpos*14.331)), #NGC 2682
                     (414912.100 * np.exp(-dpos*59.109)) / (139757.004 + 414912.100 * np.exp(-dpos*59.109)), #NGC 2660
                     (23000.000 * np.exp(-((dpos-0.022)**2/(2*0.039**2)))) / (58828.279 + 23000.000 * np.exp(-((dpos-0.022)**2/(2*0.039**2)))), #NGC 2627
                     (62734.146 * np.exp(-dpos*39.211)) / (52748.495 + 62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (163447.774 * np.exp((-dpos*17.868))) / (28867.075 + 163447.774 * np.exp((-dpos*17.868))), #NGC 2506
                     (84765.416 * np.exp(-dpos*8.043)) / (72807.519 + 84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     (16000.000 * np.exp(-((dpos-0.030)**2/(2*0.025**2)))) / (85553.016 + 16000.000 * np.exp(-((dpos-0.030)**2/(2*0.025**2)))), #NGC 2447
                     (139661.988 * np.exp(-dpos*38.210)) / (14699.250 + 139661.988 * np.exp(-dpos*38.210)), #NGC 2420
                     (62730.170 * np.exp(-dpos*20.570)) / (41666.921 + 62730.170 * np.exp(-dpos*20.570)), #NGC 2236
                     (58995.329 * np.exp(-dpos*39.285)) / (44185.987 + 58995.329 * np.exp(-dpos*39.285)), #NGC 1907
                     (77263.108 * np.exp(-dpos*17.225)) / (27484.123 + 77263.108 * np.exp(-dpos*17.225)), #NGC 1245
                     (76703.328 * np.exp(-dpos*18.636)) / (52110.055 + 76703.328 * np.exp(-dpos*18.636)), #Melotte 71
                     (118498.343 * np.exp(-dpos*15.582)) / (23179.392 + 118498.343 * np.exp(-dpos*15.582)), #Melotte 66
                     (70172.124 * np.exp(-dpos*27.637)) / (34381.117 + 70172.124 * np.exp(-dpos*27.637)), #King 5
                     (27350.619*np.exp(-((dpos-0.000)**2/(2*0.030**2)))) / (65091.819 + 27350.619*np.exp(-((dpos-0.000)**2/(2*0.030**2)))), #FSR 1252
                     (117912.880*np.exp(-((dpos-0.004)**2/(2*0.019**2)))) / (143655.820 + 117912.880*np.exp(-((dpos-0.004)**2/(2*0.019**2)))), #Czernik 37
                     (30646.548 * np.exp(-dpos*6.359)) / (37600.457 + 30646.548 * np.exp(-dpos*6.359)) ] #Collidner 110"""


        phi_cpos = [#(30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))), #Ruprecht 68
                     #(170521.883 * np.exp(-((dpos-0.002)**2/(2*0.016**2)))), #Pismis 18
                     (55926.096 * np.exp(-((dpos)**2/(2*0.090**2)))), #NGC 7789
                     (2715.388 * np.exp(-dpos*0.048)), #NGC 752
                     (12646.468 * np.exp(-((dpos)**2/(2*0.280**2)))), #NGC 6940
                     (51920.447 * np.exp(-dpos*22.830)), #NGC 6939
                     (479976.011 * np.exp(-dpos*21.117)), #NGC 6791
                     #(25208.794 * np.exp(-((dpos+0.091)**2/(2*0.1**2)))), #NGC 6208
                     (20857.247 * np.exp(-dpos*14.331)), #NGC 2682
                     #(414912.100 * np.exp(-dpos*59.109)), #NGC 2660
                     (23000.000 * np.exp(-((dpos-0.022)**2/(2*0.039**2)))), #NGC 2627
                     (62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (163447.774 * np.exp((-dpos*17.868))), #NGC 2506
                     #(84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     #(16000.000 * np.exp(-((dpos-0.030)**2/(2*0.025**2)))), #NGC 2447
                     (139661.988 * np.exp(-dpos*38.210)), #NGC 2420
                     #(62730.170 * np.exp(-dpos*20.570)), #NGC 2236
                     (58995.329 * np.exp(-dpos*39.285)), #NGC 1907
                     (77263.108 * np.exp(-dpos*17.225)), #NGC 1245
                     (76703.328 * np.exp(-dpos*18.636)), #Melotte 71
                     #(118498.343 * np.exp(-dpos*15.582)), #Melotte 66
                     (70172.124 * np.exp(-dpos*27.637)), #King 5
                     #(27350.619*np.exp(-((dpos-0.000)**2/(2*0.030**2)))), #FSR 1252
                     #(117912.880*np.exp(-((dpos-0.004)**2/(2*0.019**2)))), #Czernik 37
                     (30646.548 * np.exp(-dpos*6.359))] #Collinder 110 
        
        phi_fpos = [#73332.225, #Ruprecht 68
                    #261896.000, #Pismis 18
                    66608.971, #NGC 7789
                    2039.455, #NGC 752
                    93792.077, #NGC 6940
                    32503.010, #NGC 6939
                    53557.792, #NGC 6791
                    #106234.735, #NGC 6208
                    6456.678, #NGC 2682
                    #139757.004, #NGC 2660
                    58828.279, #NGC 2627
                    52748.495, #NGC 2509
                    28867.075, #NGC 2506
                    #72807.519, #NGC 2477
                    #85553.016, #NGC 2447
                    14699.250, #NGC 2420
                    #41666.921, #NGC 2236
                    44185.987, #NGC 1907
                    27484.123, #NGC 1245
                    52110.055, #Melotte 71
                    #23179.392, #Melotte 66
                    34381.117, #King 5
                    #65091.820, #FSR 1252
                    #143655.820, #Czernik 37
                    37600.457] #Collinder 110
        
        #ppos.append(ppos_eqns[i]) #distance from cluster center
        #pposc.append(Anorm_posC*phi_cpos[i])
        #pposf.append(Anorm_posF*phi_fpos[i])
        pposc.append(phi_cpos[i])
        pposf.append(phi_fpos[i])
        
        # --------
        # PARALLAX
        # --------
        
        dplx = plx[j] #- selected_cluster_info['Plx'][i] 
        
        """pplx_eqns = [(75.000*np.exp(-((dplx-.320)**2/(2*0.052**2)))) / (75.000*np.exp(-((dplx-.320)**2/(2*0.052**2))) + 200.000*np.exp(-((dplx-0.256)**2/(2*0.232**2))) + 23.672), #Ruprecht 68
        (43.500*np.exp(-((dplx-0.358)**2/(2*0.040**2)))) / (43.500*np.exp(-((dplx-0.358)**2/(2*0.040**2))) + 140.000*np.exp(-((dplx-0.339)**2/(2*0.200**2))) + 21.057), #Pismis 18
        (540.000*np.exp(-((dplx-0.472)**2/(2*0.020**2)))) / (540.000*np.exp(-((dplx-0.472)**2/(2*0.020**2))) + 663.369*np.exp(-((dplx-0.316)**2/(2*0.234**2))) + 83.308), #NGC 7789
        (50.000*np.exp(-((dplx-2.293)**2/(2*0.050**2)))) / (50.000*np.exp(-((dplx-2.293)**2/(2*0.050**2))) + 493.452*np.exp(-((dplx-0.500)**2/(2*0.900**2))) + 31.763), #NGC 752
        (925.842*np.exp(-((dplx-0.850)**2/(2*0.273**2)))) / (925.842*np.exp(-((dplx-0.850)**2/(2*0.273**2))) + 4900.910*np.exp(-((dplx-0.263)**2/(2*0.230**2))) + 267.896), #NGC 6940
        (146.938*np.exp(-((dplx-0.518)**2/(2*0.026**2)))) / (146.938*np.exp(-((dplx-0.518)**2/(2*0.026**2))) + 93.331*np.exp(-((dplx-0.350)**2/(2*0.266**2))) + 13.438), #NGC 6939
        (247.012*np.exp(-((dplx-0.211)**2/(2*0.073**2)))) / (247.012*np.exp(-((dplx-0.211)**2/(2*0.073**2))) + 163.732*np.exp(-((dplx-0.273)**2/(2*0.275**2))) + 19.346), #NGC 6791
        (120.000*np.exp(-((dplx-0.855)**2/(2*0.020**2)))) / (120.000*np.exp(-((dplx-0.855)**2/(2*0.020**2))) + 2628.768*np.exp(-((dplx)**2/(2*0.486**2))) + 193.809), #NGC 6208
        (130.000*np.exp(-((dplx-1.144)**2/(2*0.050**2)))) / (130.000*np.exp(-((dplx-1.144)**2/(2*0.050**2))) + 48.367*np.exp(-((dplx-0.563)**2/(2*0.622**2))) + 10.000), #NGC 2682
        (57.559*np.exp(-((dplx-0.307)**2/(2*0.070**2)))) / (57.559*np.exp(-((dplx-0.307)**2/(2*0.070**2))) + 80.000*np.exp(-((dplx-0.267)**2/(2*0.300**2))) + 5.618), #NGC 2660
        (130.000*np.exp(-((dplx-0.545)**2/(2*0.050**2)))) / (130.000*np.exp(-((dplx-0.545)**2/(2*0.050**2))) + 247.682*np.exp(-((dplx-0.240)**2/(2*0.169**2))) + 30.000), #NGC 2627
        (65.731*np.exp(-((dplx-0.360)**2/(2*0.018**2)))) / (65.731*np.exp(-((dplx-0.360)**2/(2*0.018**2))) + 94.471*np.exp(-((dplx-0.297)**2/(2*0.250**2))) + 9.676), #NGC 2509
        (186.345*np.exp(-((dplx-0.267)**2/(2*0.051**2)))) / (186.345*np.exp(-((dplx-0.267)**2/(2*0.051**2))) + 125.332*np.exp(-((dplx-0.313)**2/(2*0.228**2))) + 19.544), #NGC 2506
        (570.000*np.exp(-((dplx-0.686)**2/(2*0.041**2)))) / (570.000*np.exp(-((dplx-0.686)**2/(2*0.041**2))) + 818.705*np.exp(-((dplx-0.245)**2/(2*0.245**2))) + 100.748),#NGC 2477
        (699.896*np.exp(-((dplx-0.800)**2/(2*0.260**2)))) / (699.896*np.exp(-((dplx-0.800)**2/(2*0.260**2))) + 3391.611*np.exp(-((dplx-0.258)**2/(2*0.192**2))) + 145.156), #NGC 2447
        (65.812*np.exp(-((dplx-0.379)**2/(2*0.033**2)))) / (65.812*np.exp(-((dplx-0.379)**2/(2*0.033**2))) + 38.696*np.exp(-((dplx-0.369)**2/(2*0.216**2))) + 5.554), #NGC 2420
        (70.304*np.exp(-((dplx-0.301)**2/(2*0.123**2)))) / (70.304*np.exp(-((dplx-0.301)**2/(2*0.123**2))) + 61.964*np.exp(-((dplx-0.351)**2/(2*0.358**2))) + 9.219), #NGC 2236
        (30.000*np.exp(-((dplx-0.623)**2/(2*0.053**2)))) / (30.000*np.exp(-((dplx-0.623)**2/(2*0.053**2))) + 86.597*np.exp(-((dplx-0.372)**2/(2*0.285**2))) + 14.261), #NGC 1907
        (112.810*np.exp(-((dplx-0.285)**2/(2*0.070**2)))) / (112.810*np.exp(-((dplx-0.285)**2/(2*0.070**2))) + 130.000*np.exp(-((dplx-0.370)**2/(2*0.281**2))) + 20.000), #NGC 1245
        (87.411*np.exp(-((dplx-0.432)**2/(2*0.020**2)))) / (87.411*np.exp(-((dplx-0.432)**2/(2*0.020**2))) + 157.641*np.exp(-((dplx-0.317)**2/(2*0.235**2))) + 19.697), #Melotte 71
        (138.563*np.exp(-((dplx-0.188)**2/(2*0.063**2)))) / (138.563*np.exp(-((dplx-0.188)**2/(2*0.063**2))) + 96.834*np.exp(-((dplx-0.337)**2/(2*0.300**2))) + 9.903), #Melotte 66
        (29.727*np.exp(-((dplx-0.336)**2/(2*0.100**2)))) / (29.727*np.exp(-((dplx-0.336)**2/(2*0.100**2))) + 110.000*np.exp(-((dplx-0.372)**2/(2*0.278**2))) + 17.263), #King 5
        (37.282*np.exp(-((dplx-0.285)**2/(2*0.066**2)))) / (37.282*np.exp(-((dplx-0.285)**2/(2*0.066**2))) + (45.188*np.exp(-((dplx-0.324)**2/(2*0.314**2))) + 6.308)), #FSR 1252
        (65.000*np.exp(-((dplx-0.380)**2/(2*0.055**2)))) / (65.000*np.exp(-((dplx-0.380)**2/(2*0.055**2))) + 95.000*np.exp(-((dplx-0.421)**2/(2*0.267**2))) + 15.000), #Czernik 37
        (144.891*np.exp(-((dplx-0.430)**2/(2*0.060**2)))) / (144.891*np.exp(-((dplx-0.430)**2/(2*0.060**2))) + 570.000*np.exp(-((dplx-0.336)**2/(2*0.328**2))) + 100.000)] #Collinder 110
"""
        phi_cpx = [#(75.000*np.exp(-((dplx-.320)**2/(2*0.052**2)))), #Ruprecht 68
                   #(43.500*np.exp(-((dplx-0.358)**2/(2*0.040**2)))), #Pismis 18
                   (540.000*np.exp(-((dplx-0.472)**2/(2*0.020**2)))), #NGC 7789
                   (50.000*np.exp(-((dplx-2.293)**2/(2*0.050**2)))), #NGC 752
                   (925.842*np.exp(-((dplx-0.850)**2/(2*0.273**2)))), #NGC 6940
                   (146.938*np.exp(-((dplx-0.518)**2/(2*0.026**2)))), #NGC 6939
                   (247.012*np.exp(-((dplx-0.211)**2/(2*0.073**2)))), #NGC 6791
                   #(120.000*np.exp(-((dplx-0.855)**2/(2*0.020**2)))), #NGC 6208
                   (130.000*np.exp(-((dplx-1.144)**2/(2*0.050**2)))), #NGC 2682
                   #(57.559*np.exp(-((dplx-0.307)**2/(2*0.070**2)))), #NGC 2660
                   (130.000*np.exp(-((dplx-0.545)**2/(2*0.050**2)))), #NGC 2627
                   (65.731*np.exp(-((dplx-0.360)**2/(2*0.018**2)))), #NGC 2509
                   (186.345*np.exp(-((dplx-0.267)**2/(2*0.051**2)))), #NGC 2506
                   #(570.000*np.exp(-((dplx-0.686)**2/(2*0.041**2)))), #NGC 2477
                   #(699.896*np.exp(-((dplx-0.800)**2/(2*0.260**2)))), #NGC 2447
                   (65.812*np.exp(-((dplx-0.379)**2/(2*0.033**2)))), #NGC 2420
                   #(70.304*np.exp(-((dplx-0.301)**2/(2*0.123**2)))), #NGC 2236
                   (30.000*np.exp(-((dplx-0.623)**2/(2*0.053**2)))), #NGC 1907
                   (112.810*np.exp(-((dplx-0.285)**2/(2*0.070**2)))), #NGC 1245
                   (87.411*np.exp(-((dplx-0.432)**2/(2*0.020**2)))), #Melotte 71
                   #(138.563*np.exp(-((dplx-0.188)**2/(2*0.063**2)))), #Melotte 66
                   (29.727*np.exp(-((dplx-0.336)**2/(2*0.100**2)))), #King 5
                   #(37.282*np.exp(-((dplx-0.285)**2/(2*0.066**2)))), #FSR 1252
                   #(65.000*np.exp(-((dplx-0.380)**2/(2*0.055**2)))), #Czernik 37
                   (144.891*np.exp(-((dplx-0.430)**2/(2*0.060**2))))] #Collinder 110


        phi_fpx = [#(200.000*np.exp(-((dplx-0.256)**2/(2*0.232**2))) + 23.672), #Ruprecht 68
                   #(140.000*np.exp(-((dplx-0.339)**2/(2*0.200**2))) + 21.057), #Pismis 18
                   (663.369*np.exp(-((dplx-0.316)**2/(2*0.234**2))) + 83.308), #NGC 7789
                   (493.452*np.exp(-((dplx-0.500)**2/(2*0.900**2))) + 31.763), #NGC 752
                   (4900.910*np.exp(-((dplx-0.263)**2/(2*0.230**2))) + 267.896), #NGC 6940
                   (93.331*np.exp(-((dplx-0.350)**2/(2*0.266**2))) + 13.438), #NGC 6939
                   (163.732*np.exp(-((dplx-0.273)**2/(2*0.275**2))) + 19.346), #NGC 6791
                   #(2628.768*np.exp(-((dplx)**2/(2*0.486**2))) + 193.809), #NGC 6208
                   (48.367*np.exp(-((dplx-0.563)**2/(2*0.622**2))) + 10.000), #NGC 2682
                   #(80.000*np.exp(-((dplx-0.267)**2/(2*0.300**2))) + 5.618), #NGC 2660
                   (247.682*np.exp(-((dplx-0.240)**2/(2*0.169**2))) + 30.000), #NGC 2627
                   (94.471*np.exp(-((dplx-0.297)**2/(2*0.250**2))) + 9.676), #NGC 2509
                   (125.332*np.exp(-((dplx-0.313)**2/(2*0.228**2))) + 19.544), #NGC 2506
                   #(818.705*np.exp(-((dplx-0.245)**2/(2*0.245**2))) + 100.748), #NGC 2477
                   #(3391.611*np.exp(-((dplx-0.258)**2/(2*0.192**2))) + 145.156), #NGC 2447
                   (38.696*np.exp(-((dplx-0.369)**2/(2*0.216**2))) + 5.554), #NGC 2420
                   #(61.964*np.exp(-((dplx-0.351)**2/(2*0.358**2))) + 9.219), #NGC 2236
                   (86.597*np.exp(-((dplx-0.372)**2/(2*0.285**2))) + 14.261), #NGC 1907
                   (130.000*np.exp(-((dplx-0.370)**2/(2*0.281**2))) + 20.000), #NGC 1245
                   (157.641*np.exp(-((dplx-0.317)**2/(2*0.235**2))) + 19.697), #Melotte 71
                   #(96.834*np.exp(-((dplx-0.337)**2/(2*0.300**2))) + 9.903), #Melotte 66
                   (110.000*np.exp(-((dplx-0.372)**2/(2*0.278**2))) + 17.263), #King 5
                   #(45.188*np.exp(-((dplx-0.324)**2/(2*0.314**2))) + 6.308), #FSR 1252
                   #(95.000*np.exp(-((dplx-0.421)**2/(2*0.267**2))) + 15.000), #Czernik 37
                   (570.000*np.exp(-((dplx-0.336)**2/(2*0.328**2))) + 100.000) ] #Collinder 110
        
        #pplx.append(pplx_eqns[i]) #parallax distribution
        #pplxc.append(Anorm_plxC*phi_cpx[i])
        #pplxf.append(Anorm_plxF*phi_fpx[i])
        pplxc.append(phi_cpx[i])
        pplxf.append(phi_fpx[i])
        
        #calculate total probability from marginals:
        #total_probs.append(1 - (1-ppm[j])*(1-pplx[j])*(1-ppos[j]))
        #phic = n_c*ppmc[j]*pposc[j]*pplxc[j]
        #phif = n_f*ppmf[j]*pposf[j]*pplxf[j]
        if cluster == 'NGC_6208' or cluster == 'NGC_6940' or cluster == 'Czernik_37' or cluster == 'FSR_1252':
            phic = ppmc[j]
            phif = ppmf[j]
        elif  cluster == 'NGC_2447' or cluster == 'Collinder_110':
            phic = ppmc[j]*pplxc[j]
            phif = ppmf[j]*pplxf[j]
        elif  cluster == 'Pismis_18' or cluster == 'King_5':
            phic = ppmc[j]*pposf[j]
            phif = ppmf[j]*pposc[j]
        else:
            phic = ppmc[j]*pposc[j]*pplxc[j]
            phif = ppmf[j]*pposf[j]*pplxf[j]
        total_probs2.append(phic / (phic + phif) )

    """#histogram of membership probabilities
    counts, edges, plot = plt.hist(total_probs2, bins=40, color=orange,rwidth = 1.8)
    plt.title('Membership Probability Distribution: ' + cluster.replace('_',' '))
    plt.ylim([0.0,300.0])
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.xlabel('Probability', fontsize=15)
    plt.ylabel('N', fontsize=15)
    #plt.savefig(cluster + '_prob_dist.pdf', bbox_inches='tight', format='pdf')
    #plt.close()
    plt.show()"""
    panstarrs_prob = []
    panstarrs_g_filt = []
    panstarrs_gi_filt = []
    for j in range(len(ra)):
        if np.round(dec[j],8) in np.round(panstarrs_dec,8):
            l = int(np.where(np.round(panstarrs_dec,8)==np.round(dec[j],8))[0][0])
            panstarrs_prob.append(total_probs2[j])
            panstarrs_g_filt.append(panstarrs_g[l])
            panstarrs_gi_filt.append(panstarrs_gi[l])
            
    #plot CMDs using panstarrs photometry and Gaia membership probabilities
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 10))

    im = plt.scatter(panstarrs_gi_filt,panstarrs_g_filt, s=4, c=panstarrs_prob, cmap=cm.plasma_r)
    
    if cluster=='NGC_6791':
        ylim=[19.5,13.0]
        xlim=[-1.,3.0]
    if cluster=='King_5':
        ylim=[18.0,14.0]
        xlim=[0.5,3.0]
    if cluster=='NGC_6940':
        ylim=[17.5,10.0]
        xlim=[-0.75,2.0]
    if cluster=='NGC_2682' or cluster=='NGC_1907':
        ylim=[17.5,10.0]
        xlim=[-0.75,3.0]
    if cluster=='NGC_2627' or cluster=='NGC_2509':
        ylim=[17.0,11.0]
        xlim=[-1,1.75]
    if cluster=='NGC_2420':
        ylim=[17.0,11.5]
        xlim=[-0.5,1.25]
    if cluster=='NGC_7789':
        ylim=[21.0,8.0]
        xlim=[0,2.5]
    if cluster=='NGC_752':
        ylim=[18.0,8.0]
        xlim=[-0.75,2.5]
    if cluster=='NGC_2506':
        ylim=[16.5,12.0]
        xlim=[-0.75,1.5]
    if cluster=='NGC_1245' or cluster=='Melotte_71':
        ylim=[17.0,12.0]
        xlim=[-1.,2.0]
    if cluster=='Collinder_110':
        ylim=[17.0,12.5]
        xlim=[0.,2.5]
    if cluster=='NGC_6939':
        ylim=[17.0,12.0]
        xlim=[-.5,2.0]
    plt.xlim(xlim)
    plt.ylim(ylim)
    #plt.title('Color Magnitude Diagram (Pan-STARRS): '+cluster.replace('_',' '))
    plt.xlabel('g-i (mag)',fontsize=20)
    plt.ylabel('g (mag)',fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label('Prob', fontsize=20)
    cax.tick_params(labelsize=20)
    plt.savefig('/Users/abbychriss/Desktop/'+cluster+ '_panstarrs_cmd.pdf', bbox_inches='tight', format='pdf')
    plt.close()
    #plt.show()

    cluster_out = Table()
    cluster_out['g-i'] = panstarrs_gi_filt
    cluster_out['g'] = panstarrs_g_filt
    cluster_out['Prob'] = panstarrs_prob
    ascii.write(cluster_out, cluster+'_out_panstarrs.dat', overwrite=True)

    N_80 = 0
    for j in panstarrs_prob:
        if j >= 0.8:
            N_80 +=1
    N_80_list.append(N_80)
print(N_80_list)
print('normal stop')