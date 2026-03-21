#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 15:54:13 2023

@author: abbychriss
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
import math

import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from sklearn.ensemble import RandomForestClassifier
from scipy.optimize import curve_fit
import scipy.stats
import seaborn as sns

#Clusters in study: Ruprecht 68, Pismis 18, NGC 7789, NGC 752, NGC 6940, NGC 6939
#NGC 6791, #NGC 6208, NGC 2682, NGC 2660, NGC 2627, NGC 2509, NGC 2506, NGC 2477,
#NGC 2447, NGC 2420, NGC 2236, NGC 1907, NGC 1245, Melotte 71, Melotte 66, King 5,
#IC 4756, IC 4651, IC 2714, FSR 1252, Czernik 37, Collinder 110

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['text.usetex'] = True

# read cone search results from a file
cluster_info = Table()
cluster_info['Input file name'] = np.array(['/Users/abbychriss/Desktop/WSU/rup68_cone.csv', '/Users/abbychriss/Desktop/WSU/pismis18_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc7789_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc752_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc6940_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc6939_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc6791_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc6208_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2682_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc2660_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2627_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2509_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc2506_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2477_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2447_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc2420_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2236_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc1907_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc1245_cone.csv', '/Users/abbychriss/Desktop/WSU/mel71_cone.csv', '/Users/abbychriss/Desktop/WSU/mel66_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/king5_cone.csv', '/Users/abbychriss/Desktop/WSU/fsr1252_cone.csv', '/Users/abbychriss/Desktop/WSU/czernik37_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/col110_cone.csv'],dtype=str)
cluster_info['Cluster'] = np.array(['Ruprecht_68', 'Pismis_18', 'NGC_7789', 'NGC_752', 'NGC_6940', 'NGC_6939', 'NGC_6791', 'NGC_6208', 'NGC_2682', 'NGC_2660', 'NGC_2627', 'NGC_2509', 'NGC_2506', 'NGC_2477','NGC_2447', 'NGC_2420', 'NGC_2236', 'NGC_1907', 'NGC_1245', 'Melotte_71', 'Melotte_66', 'King_5', 'FSR_1252', 'Czernik_37', 'Collinder_110'],dtype=str)
cluster_info['Ra'] = np.array([131.1470, 204.2270, 359.3340, 029.2230, 308.6260, 307.9170, 290.2210, 252.3360, 132.8460, 130.6670, 129.3090, 120.2010, 120.0100, 118.0460, 116.1410, 114.6020, 097.4160, 082.0330, 048.6910, 114.3830, 111.5730, 048.6820, 110.6590, 268.3200, 099.6770],dtype=float)
cluster_info['Dec'] = np.array([-35.9000,-62.0910,+56.7260,+37.7940,+28.2780,+60.6530,+37.7780,-53.7140,+11.8140,-47.2010,-29.9520,-19.0560,-10.7730,-38.5370,-23.8530,+21.5750,+06.8340,+35.3300,+47.2350,-12.0650,-47.6850,+52.6950,-18.7030,-27.3730,+02.0690],dtype=float)
cluster_info['Radius'] = np.array([0.074,0.033,0.211,0.485,0.25,0.123,0.068,0.162,0.208,0.035,0.087,0.06,0.088,0.15,0.202,0.053,0.07,0.075,0.109,0.075,0.088,0.091,0.045,0.0417,0.157],dtype=float)
cluster_info['PMra'] = np.array([-2.619,-5.658,-0.922,9.8092,-1.954,-1.841,-0.421,-0.966,-10.9737,-2.763,-2.381,-2.708,-2.571,-2.449,-3.5680,-1.190,-0.752,-0.040,0.504,-2.446,-1.474,-0.282,-1.822,0.459,-1.091],dtype=float)
cluster_info['PMdec'] = np.array([5.610,-2.286,-1.933,-11.7637,-9.413,-5.413,-2.269,-1.519,-2.9396,5.165,2.879,0.764,3.912,0.870,5.0434,-2.125,0.010,-3.418,-1.579,4.210,2.745,-1.200,2.990,-0.445,-2.049],dtype=float)
cluster_info['Plx'] = np.array([0.301,0.332,0.453,2.2304,0.947,0.506,0.192,0.830,1.1325,0.308,0.521,0.363,0.292,0.665,0.9603,0.363,0.350,0.611,0.289,0.434,0.183,0.367,0.263,0.406,0.425],dtype=float)
ascii.write(cluster_info, 'RC_cluster_info.dat', overwrite=True)
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

opening_angles = []
N_80_list = []
for i in range(len(cluster_info['Cluster'])):
    import csv
    ID=[];ra=[];dec=[];plx=[];eplx=[];pmra=[];epmra=[];pmdec=[];epmdec=[];g=[];bp=[];rp=[];bprp=[];rv=[];teff=[];ag=[];ebprp=[];vflag=[]
    with open(cluster_info['Input file name'][i],'r',newline='') as csvfile:
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
    
    #define colors for plots:
    green = '#B7F4C7'
    purple = '#C698EF'
    pink = '#F4B7DC'
    orange = '#ffb65c'
    darkblue = '#1C0075'
    lightblue = '#AAD9EB'
    
    # 1. plot proper motions to assess cluster locus
    plt.plot(pmra,pmdec,'b+')
    plt.xlabel(r'$\mu_\alpha$ (mas/yr)',fontsize=15)
    plt.ylabel(r'$\mu_\delta$ (mas/yr)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    close()
    #title('Proper Motions: ' + str(cluster_info['Cluster'][i]).replace('_', ' '))
    #show()
    
    # parameters: center = pmra pmdec, goodradius = radius, definitelybadradius = radius*2
    
    # compute star density as a function of radial distance
    pm_mag = math.sqrt((cluster_info['PMra'][i])**2 + (cluster_info['PMdec'][i])**2)
    if cluster_info['Cluster'][i] == 'Collinder_110' or cluster_info['Cluster'][i] == 'King_5':
        lim = 0.8
    if cluster_info['Cluster'][i] == 'Czernik_37':
        lim = 0.5
    if cluster_info['Cluster'][i] == 'FSR_1252':
        lim = 0.75
    if cluster_info['Cluster'][i] == 'King_5':
        lim = 0.55
    if cluster_info['Cluster'][i] == 'Melotte_66':
        lim = 1.5
    if cluster_info['Cluster'][i] == 'Melotte_71':
        lim = 0.8
    if cluster_info['Cluster'][i] == 'NGC_1245':
        lim = 1.0
    if cluster_info['Cluster'][i] == 'NGC_1907':
        lim = 1.0
    if cluster_info['Cluster'][i] == 'NGC_2236':
        lim = 0.75
    if cluster_info['Cluster'][i] == 'NGC_2420':
        lim = 0.7
    if cluster_info['Cluster'][i] == 'NGC_2447':
        lim = 1.25
    if cluster_info['Cluster'][i] == 'NGC_2477':
        lim = 1.1
    if cluster_info['Cluster'][i] == 'NGC_2506':
        lim = 1.0
    if cluster_info['Cluster'][i] == 'NGC_2509':
        lim = 0.5
    if cluster_info['Cluster'][i] == 'NGC_2627':
        lim = 3.75
    if cluster_info['Cluster'][i] == 'NGC_2660':
        lim = 1.6
    if cluster_info['Cluster'][i] == 'NGC_2682':
        lim = 1.5
    if cluster_info['Cluster'][i] == 'NGC_6208':
        lim = 1.8
    if cluster_info['Cluster'][i] == 'NGC_6791':
        lim = 1.5
    if cluster_info['Cluster'][i] == 'NGC_6939':
        lim = 0.75
    if cluster_info['Cluster'][i] == 'NGC_6940':
        lim = 1.25
    if cluster_info['Cluster'][i] == 'NGC_752':
        lim = 1.6
    if cluster_info['Cluster'][i] == 'NGC_7789':
        lim = 1.5
    if cluster_info['Cluster'][i] == 'Pismis_18':
        lim = 1.4
    if cluster_info['Cluster'][i] == 'Ruprecht_68':
        lim = 0.6
    bsum = np.zeros(12,dtype=float)
    bins = np.linspace(0.0,lim,13) # 13 bin edges
    for j in range(len(ra)):
         d = math.sqrt( (pmra[j] - cluster_info['PMra'][i])**2 + (pmdec[j] - cluster_info['PMdec'][i])**2 )
         ipos = np.searchsorted(bins,d) - 1
         if ipos < 12:
              bsum[ipos] = bsum[ipos] + 1
    areas = math.pi*bins**2 #units: (mas/yr)^2
    density_pm = []
    for k in range(12):
         bsum[k] = bsum[k]/(areas[k+1] - areas[k])
         density_pm.append(bsum[k])
    x = bins[-12:] - lim/24
    bar(x,bsum,width=(lim/10),color=lightblue)
    
    #use scipy.optomize.curve_fit to fit a curve to our data
    
    #define a function where the machine learning algorithm can input coefficients
    def func(x, a, b, c):
        return c + a * np.exp(-x*b)
    bound = ([0,0,min(density_pm)-10], [np.inf,np.inf,min(density_pm)])
    popt, pcov = curve_fit(func, x, density_pm, bounds=bound)
    plt.plot(np.linspace(0.0,lim,80), func(np.linspace(0.0,lim,80), *popt), 'r-', label=r'$%5.1f$exp$\{-%5.1fx\} + %5.1f$' % tuple(popt))
    coeff_dist_pm = 'a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
    print(str(cluster_info['Cluster'][i]) + ' proper motion: ' + str(coeff_dist_pm))    
    plt.xlabel(r'Distance in proper motion space (mas/yr)',fontsize=15)
    plt.ylabel(r'Star density (N yr$^2$ mas$^{-2}$)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.title('Proper Motion: ' + cluster_info['Cluster'][i].replace('_', ' '))
    plt.ylim((0,max(density_pm)+max(density_pm)/12))
    plt.legend(frameon=False,fontsize=15)
    plt.savefig(cluster_info['Cluster']+'_pm_dist.pdf', bbox_inches='tight', format='pdf')
    plt.close()
    #show()
    
    # generate colors for proper motion candidates
    #center around pmra pmdec
    category = []
    for j in range(len(ra)):
         d = math.sqrt( (pmra[j] - cluster_info['PMra'][i])**2 + (pmdec[j] - cluster_info['PMdec'][i])**2 )
         if d < cluster_info['Radius'][i]:
              colors[j] = 'xkcd:black'
              category.append(2)
         elif d > 2*cluster_info['Radius'][i]:
              colors[j] = 'xkcd:light green'
              category.append(0)
         else:
              colors[j] = 'xkcd:light purple'
              category.append(1)
    category = np.array(category)
    
    # 2. plot RA/dec and model cluster/field
    if cluster_info['Cluster'][i] == 'NGC_7789':
        for alpha in range(len(ra)):
            if ra[alpha] > 300.00:
                ra[alpha] = ra[alpha] - 360.0000
    plt.scatter(ra,dec,c=colors,s=2)
    plt.xlabel(r'$\alpha$ ($^{\circ}$)',fontsize=15)
    plt.ylabel(r'$\delta$ ($^{\circ}$)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.title('Sky Positions: ' + str(cluster_info['Cluster'][i]).replace('_', ' '))
    plt.savefig(cluster_info['Cluster']+'_sky_pos.pdf', bbox_inches='tight', format='pdf')
    #show()
    plt.close()
    opening_angles.append((max(dec)-min(dec))/2)

    #center around ra, dec from simbad
    # repeat code from above to translate to star density
    if cluster_info['Cluster'][i] == 'NGC_6208' or cluster_info['Cluster'][i] == 'NGC_2447' :
        lim=1.25
    if cluster_info['Cluster'][i] == 'NGC_1907' or cluster_info['Cluster'][i] == 'NGC_7789' or cluster_info['Cluster'][i] == 'NGC_2682' or cluster_info['Cluster'][i] == 'NGC_6939':
        lim = 2
    if cluster_info['Cluster'][i] == 'NGC_2660' or cluster_info['Cluster'][i] == 'FSR_1252' or cluster_info['Cluster'][i] == 'NGC_6791' or cluster_info['Cluster'][i] == 'Melotte_66' or cluster_info['Cluster'][i] == 'NGC_1245' or cluster_info['Cluster'][i] == 'NGC_2627' or cluster_info['Cluster'][i] == 'Pismis_18' or cluster_info['Cluster'][i] == 'Ruprecht_68' or cluster_info['Cluster'][i] == 'King_5' or cluster_info['Cluster'][i] == 'NGC_2477' or cluster_info['Cluster'][i] == 'Czernik_37' or cluster_info['Cluster'][i] == 'Melotte_71' or cluster_info['Cluster'][i] == 'NGC_2236' or cluster_info['Cluster'][i] == 'NGC_2506' or cluster_info['Cluster'][i] == 'NGC_2509':
        lim = 3
    if cluster_info['Cluster'][i] == 'Collinder_110' or cluster_info['Cluster'][i] == 'NGC_752' or cluster_info['Cluster'][i] == 'NGC_2420' or cluster_info['Cluster'][i] == 'NGC_6940':
        lim = 4
    
    bins=np.linspace(0.0,lim*cluster_info['Radius'][i],13) # 17 bin edges
    bsum = np.zeros(12,dtype=float)
    cosdec = math.cos(cluster_info['Dec'][i]*math.pi/180.0) #use dec
    for j in range(len(ra)):
        if cluster_info['Cluster'][i] == 'NGC_7789':
            d = math.sqrt( ((ra[j] - (cluster_info['Ra'][i] - 360.000))*cosdec)**2 + (dec[j] - cluster_info['Dec'][i])**2 )
        else:
            d = math.sqrt( ((ra[j] - cluster_info['Ra'][i])*cosdec)**2 + (dec[j] - cluster_info['Dec'][i])**2 )
        ipos = np.searchsorted(bins,d) - 1
        if ipos < 12:
             bsum[ipos] = bsum[ipos] + 1
    areas = math.pi*bins**2 #units: ()
    density_center = []
    for k in range(12):
        bsum[k] = bsum[k]/(areas[k+1] - areas[k])
        density_center.append(bsum[k])
    x = bins[-12:] - lim*cluster_info['Radius'][i]/24
    plt.bar(x,bsum,width=(lim*(cluster_info['Radius'][i]/10)),color=pink)
    
    #use scipy.optomize.curve_fit to fit a curve to our data
    
    #define a function where the machine learning algorithm can input coefficients
    def func(x, a, b, c, d):
        if cluster_info['Cluster'][i] == 'NGC_752' or cluster_info['Cluster'][i] == 'NGC_1907' or cluster_info['Cluster'][i] == 'Collinder_110' or cluster_info['Cluster'][i] == 'NGC_6208' or cluster_info['Cluster'][i] == 'NGC_2660' or cluster_info['Cluster'][i] == 'NGC_6940' or cluster_info['Cluster'][i] == 'NGC_2477' or cluster_info['Cluster'][i] == 'NGC_2509' or cluster_info['Cluster'][i] == 'NGC_6939' or cluster_info['Cluster'][i] == 'Melotte_71' or cluster_info['Cluster'][i] == 'King_5' or cluster_info['Cluster'][i] == 'NGC_1245' or cluster_info['Cluster'][i] == 'NGC_2236' or cluster_info['Cluster'][i] == 'NGC_2420' or cluster_info['Cluster'][i] == 'NGC_2506':
            return c + a * np.exp(-x*b)
        elif cluster_info['Cluster'][i] == 'Pismis_18' or cluster_info['Cluster'][i] == 'NGC_2447' or cluster_info['Cluster'][i] == 'FSR_1252' or cluster_info['Cluster'][i] == 'Ruprecht_68':
            return d*np.exp(-(((x-b)/a)**2)) + c
        else:
            return a*np.exp(-((x/b)**2)) + c
        #1d gaussian with mean=0 and sigma ~ cluster radius
    if cluster_info['Cluster'][i] == 'Ruprecht_68':
        bound = ([0.01,0.,min(density_center),20000],[0.05,0.04,80000.,40000])
    elif cluster_info['Cluster'][i] == 'NGC_2477':
        bound = ([0.,0.,40000.,0.],[100000.,np.inf,min(density_center),np.inf])
    elif cluster_info['Cluster'][i] == 'NGC_2447':
        bound = ([0.025,0.03,min(density_center),16000.],[0.08,0.05,88000.,18000])
    elif cluster_info['Cluster'][i] == 'NGC_2509':
        bound = (0,[np.inf,np.inf,min(density_center),np.inf])
    elif cluster_info['Cluster'][i] == 'NGC_752':
        bound = ([0.,0.,min(density_center),0.],[500.,10,min(density_center) + 200,np.inf])
    elif cluster_info['Cluster'][i] == 'Pismis_18':
        bound = ([0.,-0.01,min(density_center),50000.],[0.02,0.01,290000.,200000.])
    elif cluster_info['Cluster'][i] == 'FSR_1252':
        bound = ([0.04,0.013,min(density_center),23000.],[0.05,0.015,66000.,25000.])
    else:
        bound = (0,[np.inf,np.inf,np.inf,np.inf])
    popt, pcov = curve_fit(func, x, density_center, bounds=bound)
    plt.plot(np.linspace(0.0,lim*cluster_info['Radius'][i],80), func(np.linspace(0.0,lim*cluster_info['Radius'][i],80), *popt), 'r-', label = r'$%5.1f$exp$\bigl\{-(\frac{x}{%5.3f})^2\bigl\} + %5.1f$' %tuple([popt[0],popt[1],popt[2]]))
    coeff_dist_cent = 'a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt)
    print(str(cluster_info['Cluster'][i]) + ' position: ' + str(coeff_dist_cent))
    plt.xlabel(r'Distance from cluster center (deg)',fontsize=15)
    plt.ylabel(r'Star density (N/deg$^2$)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.ylim((0,max(density_center)+max(density_center)/16))
    #title('Position: ' + cluster_info['Cluster'][i].replace('_', ' '))
    plt.legend(frameon=False,fontsize=15)
    plt.savefig(cluster_info['Cluster']+'_pos_dist.pdf', bbox_inches='tight', format='pdf')
    plt.close()
    #show()
    
    # 3. and, finally, model parallax
   """ n_plx = []
    if cluster_info['Cluster'][i] == 'IC_4651' or cluster_info['Cluster'][i] == 'NGC_6208':
        x2 = np.linspace(0.5,float(cluster_info['Plx'][i])+1,125)
    else:
        x2 = np.linspace(float(cluster_info['Plx'][i])-1,float(cluster_info['Plx'][i])+1,125)
    for j in range(len(x2)):
        n = 0
        for k in range(len(plx)):
            if j < 124:
                if x2[j] <= plx[k] < x2[j+1]:
                    n += 1
        n_plx.append(n)
            
    comp1=plx[(category==0)]
    comp2=plx[(category==1)]
    comp3=plx[(category==2)]
    if cluster_info['Cluster'][i] == 'IC_4651' or cluster_info['Cluster'][i] == 'NGC_6208':                
        plt.hist([comp1,comp2,comp3],bins=125,range=(0.5,float(cluster_info['Plx'][i])+1),stacked=True,color=[green,purple,darkblue], rwidth=1.5)
    else:
        plt.hist([comp1,comp2,comp3],bins=125,range=(float(cluster_info['Plx'][i])-1,float(cluster_info['Plx'][i])+1),stacked=True,color=[green,purple,darkblue])
    hist_plx = np.histogram(plx, range=(np.min(plx),np.max(plx)), bins=125)"""
    
    #use scipy.optomize.curve_fit to fit a curve to our data
    
    #define a function where the machine learning algorithm can input coefficients
    #have the function try different types of equations for a given data set
    def func(x, a, b, c, d, e, f, g):
        return a*np.exp(-((x-b)**2/(2*c**2))) + d*np.exp(-((x-e)**2/(2*f**2))) + g
        #two 1d gaussians: one with mean = cluster parallax, one representing field star background plus constant pedestal
        #          a     b     c     d   e    f     g
    if cluster_info['Cluster'][i] == 'Collinder_110':
        bound = ([500., 0.25, 0.35, 70, 0.40, 0.01, 100.], [650., 0.32, 0.44, 120, 0.5, 0.05, 130.])
    if cluster_info['Cluster'][i] == 'Czernik_37':
        bound = ([90., 0.42, 0.34, 68, 0.38, 0.098, 0.], [140., 0.6, 24., 0.6, 80, 0.5, 0.1])
    if cluster_info['Cluster'][i] == 'FSR_1252':
        bound = ([40., 0.25, 0.3, 50, 0.27, 0.01, 0.], [60., 0.45, 10., 0.45, 60, 0.3, 0.1])
    if cluster_info['Cluster'][i] == 'King_5':
        bound = ([100, 0.3, 0.4, 10., 0.3, 0.07, .1], [120, 0.7, 20, 1., 80, 0.4, 0.3])
    if cluster_info['Cluster'][i] == 'Melotte_66':
        bound = ([95., 0.2, 0.3, 120, 0.17, 0.01, 5.], [150., 0.5, 0.5, 170, 0.25, 0.1, 30.])
    if cluster_info['Cluster'][i] == 'Melotte_71':
        bound = ([150., 0.28, 0.3, 85, 0.4, 0.01, 18], [180., 0.33, 0.36, 150, 0.45, 0.05, 22])
    if cluster_info['Cluster'][i] == 'NGC_1245':
        bound = ([130., 0.35, 0.25, 80, 0.28, 0.07, 20.], [140., 0.4, 0.4, 150, 0.3, 0.1, 30.])
    if cluster_info['Cluster'][i] == 'NGC_1907':
        bound = ([80., 0.35, 0.4, 40, 0.60, 0.053, 10.], [100., 0.4, 0.45, 50, 0.65, 0.07, 20.])
    if cluster_info['Cluster'][i] == 'NGC_2236':
        bound = ([85., 0.3, 0.4, 75, 0.33, 0.12, 5.], [90., 0.4, 0.5, 85, 0.35, 0.3, 15.])
    if cluster_info['Cluster'][i] == 'NGC_2420':
        bound = ([30., 0.35, 0.1, 50, 0.3, 0.01, 5.], [50., 0.45, 0.37, 80, 0.4, 0.05, 7.])
    if cluster_info['Cluster'][i] == 'NGC_2447':
        bound = ([3000, 0.25, 0.35, 400, 0.9, 0.15, 80.], [4000, 0.29, 0.5, 550, 1.0, 0.25, 200.])
    if cluster_info['Cluster'][i] == 'NGC_2477':
        bound = ([650, 0.2, 0.45, 350, 0.65, 0.05, 40.], [1000, 0.35, 0.5, 650, 0.8, 0.1, 150.])
    if cluster_info['Cluster'][i] == 'NGC_2506':
        bound = ([100., 0.25, 0.2, 180, 0.26, 0.05, 10.], [130., 0.35, 0.35, 220, 0.3, 0.15, 30.])
    if cluster_info['Cluster'][i] == 'NGC_2509':
        bound = ([90., 0.25, 0.29, 115, 0.36, 0.018, 10.], [105., 0.33, 0.4, 150, 0.38, 0.05, 51.])
    if cluster_info['Cluster'][i] == 'NGC_2627':
        bound = ([150, 0.25, 0.35, 100, 0.56, 0.01, 10.], [230, 0.3, 0.45, 170, 0.6, 0.12, 40.])
    if cluster_info['Cluster'][i] == 'NGC_2660':
        bound = ([80., 0.25, 0.3, 55, 0.3, 0.01, 5.], [100., 0.27, 0.4, 70, 0.35, 0.07, 15.])
    if cluster_info['Cluster'][i] == 'NGC_2682':
        bound = ([40., 0.4, 0.3, 130, 1.1, 0.05, 10.], [70., 0.7, 0.75, 160, 1.25, 0.2, 20.])
    if cluster_info['Cluster'][i] == 'NGC_6208':
        bound = ([0., 0., 0.1, 120., 0.845, 0.02, 170], [3000, 0.6, 0.84, 190, 0.86, 0.1, 250])
    if cluster_info['Cluster'][i] == 'NGC_6791':
        bound = ([160.,0.2,0.25,130,0.2,0., 5.], [300., 0.3, 0.4, 250, 0.275, 0.15, 50])
    if cluster_info['Cluster'][i] == 'NGC_6939':
        bound = ([0,0,0,0,0,0,0], [120, 0.35, 0.6, 150, 0.55, 0.05, 15])
    if cluster_info['Cluster'][i] == 'NGC_6940':
        bound = ([4500, 0.25, 0.28, 300, 0.85, 0.3, 200], [5000, 0.3, 0.4, 1000, 0.97, 0.35, 400])
    if cluster_info['Cluster'][i] == 'NGC_752':
        bound = ([350, 0.5, 0.9, 50, 2.25, 0.05, 10.], [500, 2.0, 1.2, 70, 2.4, 0.2, 50.])
    if cluster_info['Cluster'][i] == 'NGC_7789':
        bound = ([620, 0.25, 0.2, 540, 0.46, 0.01, 80], [670, 0.34, 0.35, 550, 0.55, 0.1, 100])
    if cluster_info['Cluster'][i] == 'Pismis_18':
        bound = ([140, 0.32, 0.2, 43.5, 0.358, 0.04, 10], [150, 0.35, 0.45, 100, 0.45, 0.13, 30])
    if cluster_info['Cluster'][i] == 'Ruprecht_68':
        bound = ([200, 0.2, 0.2, 75, 0.32, 0.052, 20], [230, 0.27, 0.35, 80, 0.37, 0.06, 30])
    """popt, pcov = curve_fit(func, x2, n_plx, bounds=bound)
    plt.plot(x2, func(x2, *popt), 'r-', label= r'$%5.1f$exp$\bigl\{\frac{-(x-%5.2f)^2}{2 \cdot %5.3f^2}\bigr\} + %5.1f$exp$\bigl\{\frac{-(x-%5.2f)^2}{2 \cdot %5.3f^2}\bigr\} + %5.1f$' %tuple(popt))
    #title("Parallax: " + str(cluster_info['Cluster'][i].replace('_', ' ')))
    coeff_plx = 'a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f, f=%5.3f, g=%5.3f' % tuple(popt)
    print(str(cluster_info['Cluster'][i]) + ' parallax: ' + str(coeff_plx))
    plt.legend(frameon=False,fontsize=12)
    plt.ylim([0.0,1400])
    plt.xlabel(r'Parallax (mas)',fontsize=15)
    plt.ylabel('N',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.savefig(cluster_info['Cluster']+'_plx_dist.pdf', bbox_inches='tight', format='pdf')
    plt.close()"""
    #show()
    
    # 4. put it all together for a probability. 
     
    #Now that we have coefficients for the curves for each plot in each cluster,
    #we can take the coefficients and write equations for each plot.
    #Put each equation for each category and each cluster in a structured array
    #Every equation has a 'c' coefficient: that will act as the constant field star background
    #Prob = (f(x)) / (f(x) + c)
    #x will be dpm, dpos, dplx in each equation  
    ppm = []
    ppos = []
    pplx = []
    total_probs = []
    for j in range(len(ra)):
        dpm = math.sqrt( (pmra[j] - cluster_info['PMra'][i])**2 + (pmdec[j] - cluster_info['PMdec'][i])**2 )
        
        ppm_eqns = [(12429.746 * np.exp(-dpm*11.399)) / (111.261 + 12429.746 * np.exp(-dpm*11.399)), #Ruprecht 68
        (7889.704 * np.exp(-dpm*10.716)) / (317.237 + 7889.704 * np.exp(-dpm*10.716)), #Pismis 18
        (41242.811 * np.exp(-dpm*7.103)) / (1225.596 + 41242.811 * np.exp(-dpm*7.103)), #NGC 7789
        (814.720 * np.exp(-dpm*3.192)) / (2.456 + 814.720 * np.exp(-dpm*3.192)), #NGC 752
        (7707.473 * np.exp(-dpm*6.728)) / (588.992 + 7707.473 * np.exp(-dpm*6.728)), #NGC 6940
        (13557.738 * np.exp(-dpm*9.377)) / (99.831 + 13557.738 * np.exp(-dpm*9.377)), #NGC 6939
        (47524.902 * np.exp(-dpm*8.087)) / (130.203 + 47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        (5329.684 * np.exp(-dpm*7.083)) / (682.752 + 5329.684 * np.exp(-dpm*7.083)), #NGC 6208
        (6386.179 * np.exp(-dpm*4.596)) / (14.800 + 6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        (12625.640 * np.exp(-dpm*10.175)) / (155.176 + 12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (3115.683 * np.exp(-dpm*1.959)) / (53.852 + 3115.683 * np.exp(-dpm*1.959)),  #NGC 2627
        (17280.461 * np.exp(-dpm*17.437)) / (159.432 + 17280.461 * np.exp(-dpm*17.437)), #NGC 2509
        (33119.650 * np.exp(-dpm*8.828)) / (117.545 + 33119.650 * np.exp(-dpm*8.828)), #NGC 2506
        (21370.314 * np.exp(-dpm*5.392)) / (854.059 + 21370.314 * np.exp(-dpm*5.392)), #NGC 2477
        (8745.948 * np.exp(-dpm*6.120)) / (1355.076 + 8745.948 * np.exp(-dpm*6.120)), #NGC 2447
        (12374.103 * np.exp(-dpm*9.113)) / (120.148 + 12374.103 * np.exp(-dpm*9.113)), #NGC 2420
        (10819.833 * np.exp(-dpm*9.048)) / (428.479 + 10819.833 * np.exp(-dpm*9.048)), #NGC 2236
        (3789.570 * np.exp(-dpm*6.141)) / (123.524 + 3789.570 * np.exp(-dpm*6.141)), #NGC 1907
        (13338.112 * np.exp(-dpm*7.465)) / (389.433 + 13338.112 * np.exp(-dpm*7.465)), #NGC 1245
        (16787.800 * np.exp(-dpm*9.424)) / (136.353 + 16787.800 * np.exp(-dpm*9.424)), #Melotte 71
        (23473.741 * np.exp(-dpm*9.355)) / (139.946 + 23473.741 * np.exp(-dpm*9.355)), #Melotte 66
        (7973.347 * np.exp(-dpm*9.869)) / (609.280 + 7973.347 * np.exp(-dpm*9.869)), #King 5
        (3992.037 * np.exp(-dpm*9.735)) / (143.573 + 3992.037 * np.exp(-dpm*9.735)), #FSR 1252
        (5841.210 * np.exp(-dpm*10.104)) / (321.770 + 5841.210 * np.exp(-dpm*10.104)), #Czernik 37
        (28809.904 * np.exp(-dpm*10.447)) / (955.309 + 28809.904 * np.exp(-dpm*10.447))] #Collinder 110 
        
        ppm.append(ppm_eqns[i]) #proper motion distribution
        
        if cluster_info['Cluster'][i] == 'NGC_7789':
            dpos = math.sqrt( ((ra[j] - (cluster_info['Ra'][i]-360.000))*cosdec)**2 + (dec[j] - cluster_info['Dec'][i])**2 )
        else: dpos = math.sqrt( ((ra[j] - cluster_info['Ra'][i])*cosdec)**2 + (dec[j] - cluster_info['Dec'][i])**2 )

        ppos_eqns = [(28629.163*np.exp(-(((dpos-0.024)/0.044)**2))) / (28629.163*np.exp(-(((dpos-0.024)/0.044)**2)) + 74723.268), #Ruprecht 68
                     (162983.751*np.exp(-(((dpos-0.004)/0.020)**2))) / (162983.751*np.exp(-(((dpos-0.004)/0.020)**2)) + 266119.992), #Pismis 18
                     (55926.119*np.exp(-((dpos/0.127)**2))) / (55926.119*np.exp(-((dpos/0.127)**2)) + 66608.996), #NGC 7789
                     (291.079 * np.exp(-dpos*0.794)) / (4485.024 + 291.079 * np.exp(-dpos*0.794)), #NGC 752
                     (18525.755 * np.exp(-dpos*2.182)) / (90713.287 + 18525.755 * np.exp(-dpos*2.182)), #NGC 6940
                     (51920.463 * np.exp(-dpos*22.830)) / (32503.016 + 51920.463 * np.exp(-dpos*22.830)), #NGC 6939
                     (357387.563 * np.exp(-(((dpos)/0.054)**2))) / ((357387.563 * np.exp(-((dpos)/0.054)**2)) + 79040.411), #NGC 6791
                     (19653.147 * np.exp(-dpos*11.648)) / (104341.075 + 19653.147 * np.exp(-dpos*11.648)), #NGC 6208
                     (15717.902 * np.exp(-(((dpos)/0.081)**2))) / ((15717.902 * np.exp(-((dpos)/0.081)**2)) + 6981.739), #NGC 2682
                     (414911.408 * np.exp(-dpos*59.109)) / (139756.945 + 414911.408 * np.exp(-dpos*59.109)), #NGC 2660
                     (26196.647 * np.exp(-((dpos/0.079)**2))) / (26196.647 * np.exp(-((dpos/0.079)**2)) + 57973.031), #NGC 2627
                     (62734.146 * np.exp(-dpos*39.211)) / (52748.495 + 62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (163447.781 * np.exp(-dpos*17.868)) / (28867.082 + 163447.781 * np.exp(-dpos*17.868)), #NGC 2506
                     (84765.416 * np.exp(-dpos*8.043)) / (72807.519 + 84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     (17316.935 * np.exp(-(((dpos-0.030)/0.025)**2))) / (17316.935 * np.exp(-(((dpos-0.030)/0.025)**2)) + 86163.938), #NGC 2447
                     (139661.322 * np.exp(-dpos*38.210)) / (14699.237 + 139661.322 * np.exp(-dpos*38.210)), #NGC 2420
                     (62730.500 * np.exp(-dpos*20.570)) / (41666.997 + 62730.500 * np.exp(-dpos*20.570)), #NGC 2236
                     (58995.150 * np.exp(-dpos*39.285)) / (44185.898 + 58995.150 * np.exp(-dpos*39.285)), #NGC 1907
                     (77263.047 * np.exp(-dpos*17.225)) / (27484.122 + 77263.047 * np.exp(-dpos*17.225)), #NGC 1245
                     (76702.600 * np.exp(-dpos*18.635)) / (52110.072 + 76702.600 * np.exp(-dpos*18.635)), #Melotte 71
                     (88115.524 * np.exp(-(((dpos)/0.072)**2))) / ((88115.524 * np.exp(-((dpos)/0.072)**2)) + 30223.704), #Melotte 66
                     (70172.166 * np.exp(-dpos*27.637)) / (34381.118 + 70172.166 * np.exp(-dpos*27.637)), #King 5
                     (23813.002 * np.exp(-(((dpos-0.013)/0.040)**2))) / (23813.002 * np.exp(-(((dpos-0.013)/0.040)**2)) + 64227.803), #FSR 1252
                     (123709.357 * np.exp(-((dpos/0.031)**2))) / (123709.357 * np.exp(-((dpos/0.031)**2)) + 143239.265), #Czernik 37
                     (37600.378 * np.exp(-dpos*6.359)) / (37600.378 + 30646.870 * np.exp(-dpos*6.359))] #Collinder 110 
        
        ppos.append(ppos_eqns[i]) #distance from cluster center
        
        dplx = plx[j] - cluster_info['Plx'][i] 
        
        pplx_eqns = [(75.000*np.exp(-((dplx-0.320)/0.052)**2)) / (200.000*np.exp(-((dplx-0.260)/0.323)**2) + 75.000*np.exp(-((dplx-0.320)/0.052)**2) + 25.933), #Ruprecht 68
        (43.500*np.exp(-((dplx-0.358)/0.040)**2)) / (140.000*np.exp(-((dplx-0.340)/0.273)**2) + 43.500*np.exp(-((dplx-0.358)/0.040)**2) + 22.867), #Pismis 18
        (540.000*np.exp(-((dplx-0.472)/0.028)**2)) / (663.369*np.exp(-((dplx-0.316)/0.331)**2) + 540.000*np.exp(-((dplx-0.472)/0.028)**2) + 83.308), #NGC 7789
        (50.000*np.exp(-((dplx-2.290)/0.076)**2)) / (431.931*np.exp(-((dplx-0.500)/1.184)**2) + 50.000*np.exp(-((dplx-2.290)/0.076)**2) + 33.142), #NGC 752
        (952.959*np.exp(-((dplx-0.850)/0.350)**2)) / (4952.871*np.exp(-((dplx-0.265)/0.318)**2) + 952.959*np.exp(-((dplx-0.850)/0.350)**2) + 295.351), #NGC 6940
        (146.938*np.exp(-((dplx-0.518)/0.036)**2)) / (93.331*np.exp(-((dplx-0.350)/0.377)**2) + 146.938*np.exp(-((dplx-0.518)/0.036)**2) + 13.438), #NGC 6939
        (247.012*np.exp(-((dplx-0.211)/0.104)**2)) / (163.732*np.exp(-((dplx-0.273)/0.389)**2) + 247.012*np.exp(-((dplx-0.211)/0.104)**2) + 19.346), #NGC 6791
        (120.000*np.exp(-((dplx-0.848)/0.020)**2)) /  (2609.216*np.exp(-((dplx)/0.704)**2) + 120.000*np.exp(-((dplx-0.848)/0.020)**2) + 191.989), #NGC 6208
        (140.641*np.exp(-((dplx-1.144)/0.057)**2)) / (48.819*np.exp(-((dplx-0.642)/0.750)**2) + 140.641*np.exp(-((dplx-1.144)/0.057)**2) + 12.076), #NGC 2682
        (55.000*np.exp(-((dplx-0.322)/0.051)**2)) / (93.578*np.exp(-((dplx-0.267)/0.328)**2) + 55.000*np.exp(-((dplx-0.322)/0.051)**2) + 10.934), #NGC 2660
        (100.000*np.exp(-((dplx-0.560)/0.023)**2)) / (226.573*np.exp(-((dplx-0.278)/0.350)**2) + 100.000*np.exp(-((dplx-0.560)/0.023)**2) + 20.928), #NGC 2627
        (115.000*np.exp(-((dplx-0.360)/0.018)**2)) / (94.407*np.exp(-((dplx-0.291)/0.315)**2) + 115.000*np.exp(-((dplx-0.360)/0.018)**2) + 12.518), #NGC 2509
        (186.344*np.exp(-((dplx-0.267)/0.072)**2)) / (125.333*np.exp(-((dplx-0.313)/0.322)**2) + 186.344*np.exp(-((dplx-0.267)/0.072)**2) + 19.544), #NGC 2506
        (386.274*np.exp(-((dplx-0.686)/0.050)**2)) / (787.544*np.exp(-((dplx-0.268)/0.450)**2) + 386.274*np.exp(-((dplx-0.686)/0.050)**2) + 57.019), #NGC 2477
        (502.166*np.exp(-((dplx-0.951)/0.250)**2)) / (3211.549*np.exp(-((dplx-0.275)/0.350)**2) + 502.166*np.exp(-((dplx-0.951)/0.250)**2) + 131.877), #NGC 2447
        (65.811*np.exp(-((dplx-0.379)/0.046)**2)) / (38.698*np.exp(-((dplx-0.369)/0.305)**2) + 65.811*np.exp(-((dplx-0.379)/0.046)**2) + 5.554), #NGC 2420
        (70.000*np.exp(-((dplx-0.330)/0.120)**2)) / (70.000*np.exp(-((dplx-0.330)/0.120)**2) + 80.000*np.exp(-((dplx-0.323)/0.478)**2) + 6.541), #NGC 2236
        (40.000*np.exp(-((dplx-0.619)/0.053)**2)) / (40.000*np.exp(-((dplx-0.619)/0.053)**2) + 86.755*np.exp(-((dplx-0.374)/0.406)**2) + 14.150), #NGC 1907
        (115.595*np.exp(-((dplx-0.285)/0.094)**2)) / (130.000*np.exp(-((dplx-0.370)/0.398)**2) + 115.595*np.exp(-((dplx-0.285)/0.094)**2) + 20.000), #NGC 1245
        (87.410*np.exp(-((dplx-0.432)/0.028)**2)) / (157.641*np.exp(-((dplx-0.317)/0.332)**2) + 87.410*np.exp(-((dplx-0.432)/0.028)**2) + 19.697), #Melotte 71
        (136.157*np.exp(-((dplx-0.187)/0.086)**2)) / (99.081*np.exp(-((dplx-0.331)/0.397)**2) + 136.157*np.exp(-((dplx-0.187)/0.086)**2) + 11.957), #Melotte 66
        (60.000*np.exp(-((dplx-0.344)/0.126)**2)) / (100.000*np.exp(-((dplx-0.383)/0.483)**2) + 60.000*np.exp(-((dplx-0.344)/0.126)**2) + 10.000), #King 5
        (50.000*np.exp(-((dplx-0.282)/0.075)**2)) / (43.985*np.exp(-((dplx-0.328)/0.450)**2) + 50.000*np.exp(-((dplx-0.282)/0.075)**2) + 6.339), #FSR 1252
        (68.000*np.exp(-((dplx-0.380)/0.098)**2)) / (68.000*np.exp(-((dplx-0.380)/0.098)**2) + 90.000*np.exp(-((dplx-0.425)/0.361)**2) + 18.000), #Czernik 37
        (99.331*np.exp(-((dplx-0.429)/0.039)**2)) / (634.585*np.exp(-((dplx-0.320)/0.405)**2) + 99.331*np.exp(-((dplx-0.429)/0.039)**2) + 121.966)] #Collinder 110 
    
        pplx.append(pplx_eqns[i]) #parallax distribution
        
        #calculate total probability from marginals:
        total_probs.append(1 - (1-ppm[j])*(1-pplx[j])*(1-ppos[j]))
        
    
    #plot a histograms of the probabilities for stars in each cluster
    hist(ppm, range=(min(ppm),max(ppm)), bins=40, color='#E982C4', rwidth = 1.2)
    #title('Proper motion probability: ' + cluster_info['Cluster'][i].replace('_',' '))
    xlabel('Probability',fontsize=15)
    ylim([0.0,3100.0])
    ylabel('N',fontsize=15)
    close()
    #show()
    
    hist(ppos, range=(min(ppos),max(ppos)), bins=40, color='#E98291', rwidth = 1.2)
    #title('Position probability: ' + cluster_info['Cluster'][i].replace('_',' '))
    ylim([0.0,3100.0])
    xlabel('Probability',fontsize=15)
    ylabel('N',fontsize=15)
    close()
    #show()
    
    hist(pplx, range=(min(pplx),max(pplx)), bins=40, color='#E9DB82', rwidth = 1.2)
    #title('Parallax probability: ' + cluster_info['Cluster'][i].replace('_',' '))
    ylim([0.0,200.0])
    xlabel('Probability', fontsize=15)
    ylabel('N', fontsize=15)
    close()
    #show()
    
    counts, edges, plot = hist(total_probs, range=(min(total_probs),max(total_probs)), bins=40, color=orange,rwidth = 1.2)
    #title('Membership Probability Distribution: ' + cluster_info['Cluster'][i].replace('_',' '))
    ylim([0.0,3100.0])
    yticks(fontsize=15)
    xticks(fontsize=15)
    xlabel('Probability', fontsize=15)
    ylabel('N', fontsize=15)
    savefig(cluster_info['Cluster'][i]+'prob_dist.pdf', bbox_inches='tight', format='pdf')
    close()
    #show()
    
    print(cluster_info['Cluster'][i] + ' n>80: ' + str(counts[0]))
    
    #let's try something a bit...different for the probability based on Griggio and Bedin MNRAS 551, 4702–4709 (2022)
    # add the field pm, pos, plx distributions then divide by sum of total distributions
    better_total_prob_dist = []
    cluster_combined_dist = []
    field_plus_cluster_combined_dist = []
    for j in range(len(ra)):
        dpm = math.sqrt( (pmra[j] - cluster_info['PMra'][i])**2 + (pmdec[j] - cluster_info['PMdec'][i])**2 )
        
        ppm_cluster = [(12429.746 * np.exp(-dpm*11.399)), #Ruprecht 68
        (7889.704 * np.exp(-dpm*10.716)), #Pismis 18
        (41242.811 * np.exp(-dpm*7.103)), #NGC 7789
        (814.720 * np.exp(-dpm*3.192)), #NGC 752
        (7707.473 * np.exp(-dpm*6.728)), #NGC 6940
        (13557.738 * np.exp(-dpm*9.377)), #NGC 6939
        (47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        (5329.684 * np.exp(-dpm*7.083)), #NGC 6208
        (6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        (12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (3115.683 * np.exp(-dpm*1.959)),  #NGC 2627
        (17280.461 * np.exp(-dpm*17.437)), #NGC 2509
        (33119.650 * np.exp(-dpm*8.828)), #NGC 2506
        (21370.314 * np.exp(-dpm*5.392)), #NGC 2477
        (8745.948 * np.exp(-dpm*6.120)), #NGC 2447
        (12374.103 * np.exp(-dpm*9.113)), #NGC 2420
        (10819.833 * np.exp(-dpm*9.048)), #NGC 2236
        (3789.570 * np.exp(-dpm*6.141)), #NGC 1907
        (13338.112 * np.exp(-dpm*7.465)), #NGC 1245
        (16787.800 * np.exp(-dpm*9.424)), #Melotte 71
        (23473.741 * np.exp(-dpm*9.355)), #Melotte 66
        (7973.347 * np.exp(-dpm*9.869)), #King 5
        (3992.037 * np.exp(-dpm*9.735)), #FSR 1252
        (5841.210 * np.exp(-dpm*10.104)), #Czernik 37
        (28809.904 * np.exp(-dpm*10.447))] #Collinder 110 
        
        ppm_total = [(111.261 + 12429.746 * np.exp(-dpm*11.399)), #Ruprecht 68
        (317.237 + 7889.704 * np.exp(-dpm*10.716)), #Pismis 18
        (1225.596 + 41242.811 * np.exp(-dpm*7.103)), #NGC 7789
        (2.456 + 814.720 * np.exp(-dpm*3.192)), #NGC 752
        (588.992 + 7707.473 * np.exp(-dpm*6.728)), #NGC 6940
        (99.831 + 13557.738 * np.exp(-dpm*9.377)), #NGC 6939
        (130.203 + 47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        (682.752 + 5329.684 * np.exp(-dpm*7.083)), #NGC 6208
        (14.800 + 6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        (155.176 + 12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (53.852 + 3115.683 * np.exp(-dpm*1.959)),  #NGC 2627
        (159.432 + 17280.461 * np.exp(-dpm*17.437)), #NGC 2509
        (117.545 + 33119.650 * np.exp(-dpm*8.828)), #NGC 2506
        (854.059 + 21370.314 * np.exp(-dpm*5.392)), #NGC 2477
        (1355.076 + 8745.948 * np.exp(-dpm*6.120)), #NGC 2447
        (120.148 + 12374.103 * np.exp(-dpm*9.113)), #NGC 2420
        (428.479 + 10819.833 * np.exp(-dpm*9.048)), #NGC 2236
        (123.524 + 3789.570 * np.exp(-dpm*6.141)), #NGC 1907
        (389.433 + 13338.112 * np.exp(-dpm*7.465)), #NGC 1245
        (136.353 + 16787.800 * np.exp(-dpm*9.424)), #Melotte 71
        (139.946 + 23473.741 * np.exp(-dpm*9.355)), #Melotte 66
        (609.280 + 7973.347 * np.exp(-dpm*9.869)), #King 5
        (143.573 + 3992.037 * np.exp(-dpm*9.735)), #FSR 1252
        (321.770 + 5841.210 * np.exp(-dpm*10.104)), #Czernik 37
        (955.309 + 28809.904 * np.exp(-dpm*10.447))] #Collinder 110 
        
        if cluster_info['Cluster'][i] == 'NGC_7789':
            dpos = math.sqrt( ((ra[j] - (cluster_info['Ra'][i]-360.000))*cosdec)**2 + (dec[j] - cluster_info['Dec'][i])**2 )
        else: dpos = math.sqrt( ((ra[j] - cluster_info['Ra'][i])*cosdec)**2 + (dec[j] - cluster_info['Dec'][i])**2 )

        ppos_cluster = [(28629.163*np.exp(-(((dpos-0.024)/0.044)**2))), #Ruprecht 68
                     (162983.751*np.exp(-(((dpos-0.004)/0.020)**2))), #Pismis 18
                     (55926.119*np.exp(-((dpos/0.127)**2))), #NGC 7789
                     (291.079 * np.exp(-dpos*0.794)), #NGC 752
                     (18525.755 * np.exp(-dpos*2.182)), #NGC 6940
                     (51920.463 * np.exp(-dpos*22.830)), #NGC 6939
                     (357387.563 * np.exp(-(((dpos)/0.054)**2))), #NGC 6791
                     (19653.147 * np.exp(-dpos*11.648)), #NGC 6208
                     (15717.902 * np.exp(-(((dpos)/0.081)**2))), #NGC 2682
                     (414911.408 * np.exp(-dpos*59.109)), #NGC 2660
                     (26196.647 * np.exp(-((dpos/0.079)**2))), #NGC 2627
                     (62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (163447.781 * np.exp(-dpos*17.868)), #NGC 2506
                     (84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     (17316.935 * np.exp(-(((dpos-0.030)/0.025)**2))), #NGC 2447
                     (139661.322 * np.exp(-dpos*38.210)), #NGC 2420
                     (62730.500 * np.exp(-dpos*20.570)), #NGC 2236
                     (58995.150 * np.exp(-dpos*39.285)), #NGC 1907
                     (77263.047 * np.exp(-dpos*17.225)), #NGC 1245
                     (76702.600 * np.exp(-dpos*18.635)), #Melotte 71
                     (88115.524 * np.exp(-(((dpos)/0.072)**2))), #Melotte 66
                     (70172.166 * np.exp(-dpos*27.637)), #King 5
                     (23813.002 * np.exp(-(((dpos-0.013)/0.040)**2))), #FSR 1252
                     (123709.357 * np.exp(-((dpos/0.031)**2))), #Czernik 37
                     (37600.378 * np.exp(-dpos*6.359))] #Collinder 110 
        
        ppos_total = [(28629.163*np.exp(-(((dpos-0.024)/0.044)**2)) + 74723.268), #Ruprecht 68
                     (162983.751*np.exp(-(((dpos-0.004)/0.020)**2)) + 266119.992), #Pismis 18
                     (55926.119*np.exp(-((dpos/0.127)**2)) + 66608.996), #NGC 7789
                     (4485.024 + 291.079 * np.exp(-dpos*0.794)), #NGC 752
                     (90713.287 + 18525.755 * np.exp(-dpos*2.182)), #NGC 6940
                     (32503.016 + 51920.463 * np.exp(-dpos*22.830)), #NGC 6939
                     ((357387.563 * np.exp(-((dpos)/0.054)**2)) + 79040.411), #NGC 6791
                     (104341.075 + 19653.147 * np.exp(-dpos*11.648)), #NGC 6208
                     ((15717.902 * np.exp(-((dpos)/0.081)**2)) + 6981.739), #NGC 2682
                     (139756.945 + 414911.408 * np.exp(-dpos*59.109)), #NGC 2660
                     (26196.647 * np.exp(-((dpos/0.079)**2)) + 57973.031), #NGC 2627
                     (52748.495 + 62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (28867.082 + 163447.781 * np.exp(-dpos*17.868)), #NGC 2506
                     (72807.519 + 84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     (17316.935 * np.exp(-(((dpos-0.030)/0.025)**2)) + 86163.938), #NGC 2447
                     (14699.237 + 139661.322 * np.exp(-dpos*38.210)), #NGC 2420
                     (41666.997 + 62730.500 * np.exp(-dpos*20.570)), #NGC 2236
                     (44185.898 + 58995.150 * np.exp(-dpos*39.285)), #NGC 1907
                     (27484.122 + 77263.047 * np.exp(-dpos*17.225)), #NGC 1245
                     (52110.072 + 76702.600 * np.exp(-dpos*18.635)), #Melotte 71
                     ((88115.524 * np.exp(-((dpos)/0.072)**2)) + 30223.704), #Melotte 66
                     (34381.118 + 70172.166 * np.exp(-dpos*27.637)), #King 5
                     (23813.002 * np.exp(-(((dpos-0.013)/0.040)**2)) + 64227.803), #FSR 1252
                     (123709.357 * np.exp(-((dpos/0.031)**2)) + 143239.265), #Czernik 37
                     (37600.378 + 30646.870 * np.exp(-dpos*6.359))] #Collinder 110 
        
        dplx = plx[j] - cluster_info['Plx'][i] 
        
        pplx_cluster = [(75.000*np.exp(-((dplx-0.320)/0.052)**2)), #Ruprecht 68
        (43.500*np.exp(-((dplx-0.358)/0.040)**2)), #Pismis 18
        (540.000*np.exp(-((dplx-0.472)/0.028)**2)), #NGC 7789
        (50.000*np.exp(-((dplx-2.290)/0.076)**2)), #NGC 752
        (952.959*np.exp(-((dplx-0.850)/0.350)**2)), #NGC 6940
        (146.938*np.exp(-((dplx-0.518)/0.036)**2)), #NGC 6939
        (247.012*np.exp(-((dplx-0.211)/0.104)**2)), #NGC 6791
        (120.000*np.exp(-((dplx-0.848)/0.020)**2)), #NGC 6208
        (140.641*np.exp(-((dplx-1.144)/0.057)**2)), #NGC 2682
        (55.000*np.exp(-((dplx-0.322)/0.051)**2)), #NGC 2660
        (100.000*np.exp(-((dplx-0.560)/0.023)**2)), #NGC 2627
        (115.000*np.exp(-((dplx-0.360)/0.018)**2)), #NGC 2509
        (186.344*np.exp(-((dplx-0.267)/0.072)**2)), #NGC 2506
        (386.274*np.exp(-((dplx-0.686)/0.050)**2)), #NGC 2477
        (502.166*np.exp(-((dplx-0.951)/0.250)**2)), #NGC 2447
        (65.811*np.exp(-((dplx-0.379)/0.046)**2)), #NGC 2420
        (70.000*np.exp(-((dplx-0.330)/0.120)**2)), #NGC 2236
        (40.000*np.exp(-((dplx-0.619)/0.053)**2)), #NGC 1907
        (115.595*np.exp(-((dplx-0.285)/0.094)**2)), #NGC 1245
        (87.410*np.exp(-((dplx-0.432)/0.028)**2)), #Melotte 71
        (136.157*np.exp(-((dplx-0.187)/0.086)**2)), #Melotte 66
        (60.000*np.exp(-((dplx-0.344)/0.126)**2)), #King 5
        (50.000*np.exp(-((dplx-0.282)/0.075)**2)), #FSR 1252
        (68.000*np.exp(-((dplx-0.380)/0.098)**2)), #Czernik 37
        (99.331*np.exp(-((dplx-0.429)/0.039)**2))] #Collinder 110 
    
        pplx_total = [(200.000*np.exp(-((dplx-0.260)/0.323)**2) + 75.000*np.exp(-((dplx-0.320)/0.052)**2) + 25.933), #Ruprecht 68
        (140.000*np.exp(-((dplx-0.340)/0.273)**2) + 43.500*np.exp(-((dplx-0.358)/0.040)**2) + 22.867), #Pismis 18
        (663.369*np.exp(-((dplx-0.316)/0.331)**2) + 540.000*np.exp(-((dplx-0.472)/0.028)**2) + 83.308), #NGC 7789
        (431.931*np.exp(-((dplx-0.500)/1.184)**2) + 50.000*np.exp(-((dplx-2.290)/0.076)**2) + 33.142), #NGC 752
        (4952.871*np.exp(-((dplx-0.265)/0.318)**2) + 952.959*np.exp(-((dplx-0.850)/0.350)**2) + 295.351), #NGC 6940
        (93.331*np.exp(-((dplx-0.350)/0.377)**2) + 146.938*np.exp(-((dplx-0.518)/0.036)**2) + 13.438), #NGC 6939
        (163.732*np.exp(-((dplx-0.273)/0.389)**2) + 247.012*np.exp(-((dplx-0.211)/0.104)**2) + 19.346), #NGC 6791
        (2609.216*np.exp(-((dplx)/0.704)**2) + 120.000*np.exp(-((dplx-0.848)/0.020)**2) + 191.989), #NGC 6208
        (48.819*np.exp(-((dplx-0.642)/0.750)**2) + 140.641*np.exp(-((dplx-1.144)/0.057)**2) + 12.076), #NGC 2682
        (93.578*np.exp(-((dplx-0.267)/0.328)**2) + 55.000*np.exp(-((dplx-0.322)/0.051)**2) + 10.934), #NGC 2660
        (226.573*np.exp(-((dplx-0.278)/0.350)**2) + 100.000*np.exp(-((dplx-0.560)/0.023)**2) + 20.928), #NGC 2627
        (94.407*np.exp(-((dplx-0.291)/0.315)**2) + 115.000*np.exp(-((dplx-0.360)/0.018)**2) + 12.518), #NGC 2509
        (125.333*np.exp(-((dplx-0.313)/0.322)**2) + 186.344*np.exp(-((dplx-0.267)/0.072)**2) + 19.544), #NGC 2506
        (787.544*np.exp(-((dplx-0.268)/0.450)**2) + 386.274*np.exp(-((dplx-0.686)/0.050)**2) + 57.019), #NGC 2477
        (3211.549*np.exp(-((dplx-0.275)/0.350)**2) + 502.166*np.exp(-((dplx-0.951)/0.250)**2) + 131.877), #NGC 2447
        (38.698*np.exp(-((dplx-0.369)/0.305)**2) + 65.811*np.exp(-((dplx-0.379)/0.046)**2) + 5.554), #NGC 2420
        (70.000*np.exp(-((dplx-0.330)/0.120)**2) + 80.000*np.exp(-((dplx-0.323)/0.478)**2) + 6.541), #NGC 2236
        (40.000*np.exp(-((dplx-0.619)/0.053)**2) + 86.755*np.exp(-((dplx-0.374)/0.406)**2) + 14.150), #NGC 1907
        (130.000*np.exp(-((dplx-0.370)/0.398)**2) + 115.595*np.exp(-((dplx-0.285)/0.094)**2) + 20.000), #NGC 1245
        (157.641*np.exp(-((dplx-0.317)/0.332)**2) + 87.410*np.exp(-((dplx-0.432)/0.028)**2) + 19.697), #Melotte 71
        (99.081*np.exp(-((dplx-0.331)/0.397)**2) + 136.157*np.exp(-((dplx-0.187)/0.086)**2) + 11.957), #Melotte 66
        (100.000*np.exp(-((dplx-0.383)/0.483)**2) + 60.000*np.exp(-((dplx-0.344)/0.126)**2) + 10.000), #King 5
        (43.985*np.exp(-((dplx-0.328)/0.450)**2) + 50.000*np.exp(-((dplx-0.282)/0.075)**2) + 6.339), #FSR 1252
        (68.000*np.exp(-((dplx-0.380)/0.098)**2) + 90.000*np.exp(-((dplx-0.425)/0.361)**2) + 18.000), #Czernik 37
        (634.585*np.exp(-((dplx-0.320)/0.405)**2) + 99.331*np.exp(-((dplx-0.429)/0.039)**2) + 121.966)] #Collinder 110 
        
        cluster_combined_dist.append(ppm_cluster[i] + ppos_cluster[i] + pplx_cluster[i])
        field_plus_cluster_combined_dist.append(ppm_total[i] + ppos_total[i] + pplx_total[i])
    better_total_prob_dist = [cluster_combined_dist[j] / field_plus_cluster_combined_dist[j] for j in range(len(ra))]
    
    #plot the better probability distribution (let's hope it's not wack)
    counts, edges, plot = hist(better_total_prob_dist, range=(min(better_total_prob_dist),max(better_total_prob_dist)), bins=40, color=orange,rwidth = 1.2)
    title('Membership Probability Distribution: ' + cluster_info['Cluster'][i].replace('_',' '))
    ylim([0.0,3100.0])
    yticks(fontsize=15)
    xticks(fontsize=15)
    xlabel('Probability', fontsize=15)
    ylabel('N', fontsize=15)
    #savefig('ngc7789_prob_dist.pdf', bbox_inches='tight', format='pdf')
    #close()
    show()
    
    # plot RA vs Dec, PMRa vs PMDec, CMDs with probabilities from total_probs
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    ax1.scatter(ra, dec, s=5, c=better_total_prob_dist, cmap=cm.plasma)
    ax1.set_title("RA - Dec: " + str(cluster_info['Cluster'][i]).replace('_',' '))
    ax1.set_ylabel(r'Dec (deg)')
    ax1.set_xlabel(r'RA (deg)')
    
    im = ax2.scatter(pmra, pmdec, s=5, c=better_total_prob_dist, cmap=cm.plasma)
    ax2.set_title("pmRA - pmDEC: " + str(cluster_info['Cluster'][i]).replace('_',' '))
    ax2.set_ylabel(r'pmdec (mas yr$^{-1}$)')
    ax2.set_xlabel(r'pmRA (mas yr$^{-1}$)')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, label='Prob')
    plt.show()
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    bprp_filt = []
    g_filt = []
    total_probs_filt = []
    for j in range(len(bprp)):
        if cluster_info['Cluster'][i] == "NGC_7789":
            if better_total_prob_dist[j] >= 0.45:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs_filt.append(better_total_prob_dist[j])
        else:
            if better_total_prob_dist[j] >= 0.1:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs_filt.append(better_total_prob_dist[j])
    im = plt.scatter(bprp_filt,g_filt, s=5, c=total_probs_filt, cmap=cm.plasma_r)
    plt.xlim([-2.0,6.0])
    plt.ylim([22.0,6.0])
    plt.title('Color Magnitude Diagram: ' + str(cluster_info['Cluster'][i]).replace('_',' '))
    plt.xlabel('BP-RP (mag)',fontsize=15)
    plt.ylabel('G (mag)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax, label='Prob')
    plt.show()
    
    N_80 = 0
    for j in total_probs:
        if j >= 0.8:
            N_80 +=1
    N_80_list.append(N_80)
    
    cluster_out = Table()
    cluster_out['BP-RP'] = bprp
    cluster_out['G'] = g
    cluster_out['Prob'] = total_probs
    ascii.write(cluster_out, cluster_info['Cluster'][i]+'_out.dat', overwrite=True)

#lists of opening angles and number of stars with proability > 80 for each cluster
print([str(N_80_list[i]) + ' ' + cluster_info['Cluster'][i] for i in range(len(N_80_list))])
print([str(opening_angles[i]) + ' ' + cluster_info['Cluster'][i] for i in range(len(N_80_list))])

print('normal stop')