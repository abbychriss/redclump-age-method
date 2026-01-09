#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 13:33:01 2024

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

rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = True

#load in names of clusters that we are using in alphabetical order
cluster_names = np.array(['Collinder_110', 'King_5', 'Melotte_66', 'Melotte_71', 'NGC_1245', 'NGC_1907', 'NGC_2420', 'NGC_2477', 'NGC_2506', 'NGC_2509', 'NGC_2627', 'NGC_2682', 'NGC_6208', 'NGC_6791', 'NGC_6939', 'NGC_6940', 'NGC_752', 'NGC_7789', 'Ruprecht_68'])

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.isin(np.array(cluster_info['Cluster']), cluster_names)]
selected_cluster_info.rename_column('Input file name', 'Gaia file')

panstarrs_file = np.array([selected_cluster_info['Gaia file'][i].replace('.csv', '_panstarrs.csv') for i in range(len(selected_cluster_info['Cluster']))])
selected_cluster_info.add_column(panstarrs_file, name='PanSTARRS file', index=1)
for i in range(len(selected_cluster_info['Cluster'])):
    if selected_cluster_info['Cluster'][i] in ['NGC_2660', 'NGC_2477', 'NGC_6208', 'Ruprecht_68', 'Melotte_66']:
        selected_cluster_info['PanSTARRS file'][i] = 'None'

print(selected_cluster_info)

"""
#load in gaia astrometry for stars in each cluster
opening_angles = []
N_80_list = []
for i in range(len(selected_cluster_info['Cluster'])):
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
    
    
    #read csv panstarrs files to get g and i values for stars in clusters that have panstarrs photometry
    if selected_cluster_info['Cluster'][i] not in ['NGC_2660', 'NGC_2477', 'NGC_6208', 'Ruprecht_68', 'Melotte_66']:
        
    
    #define colors for plots:
    green = '#B7F4C7'
    purple = '#C698EF'
    pink = '#F4B7DC'
    orange = '#ffb65c'
    darkblue = '#1C0075'
    lightblue = '#AAD9EB'
    
    # 1. plot proper motions to assess cluster locus
    plot(pmra,pmdec,'b+')
    xlabel(r'$\mu_\alpha$ (mas yr$^{-1}$)')
    ylabel(r'$\mu_\delta$ (mas ys$^{-1}$)')
    title('Proper Motions: ' + str(selected_cluster_info['Cluster'][i]).replace('_', ' '))
    show()
    
    # parameters: center = pmra pmdec, goodradius = radius, definitelybadradius = radius*2
    
    # compute star density as a function of radial distance
    pm_mag = math.sqrt((selected_cluster_info['PMra'][i])**2 + (selected_cluster_info['PMdec'][i])**2)
    if selected_cluster_info['Cluster'][i] == 'Collinder_110' or selected_cluster_info['Cluster'][i] == 'King_5':
        lim = 0.8
    if selected_cluster_info['Cluster'][i] == 'King_5':
        lim = 0.55
    if selected_cluster_info['Cluster'][i] == 'Melotte_66':
        lim = 1.5
    if selected_cluster_info['Cluster'][i] == 'Melotte_71':
        lim = 0.8
    if selected_cluster_info['Cluster'][i] == 'NGC_1245':
        lim = 1.0
    if selected_cluster_info['Cluster'][i] == 'NGC_1907':
        lim = 1.0
    if selected_cluster_info['Cluster'][i] == 'NGC_2236':
        lim = 0.75
    if selected_cluster_info['Cluster'][i] == 'NGC_2420':
        lim = 0.7
    if selected_cluster_info['Cluster'][i] == 'NGC_2447':
        lim = 1.25
    if selected_cluster_info['Cluster'][i] == 'NGC_2477':
        lim = 1.1
    if selected_cluster_info['Cluster'][i] == 'NGC_2506':
        lim = 1.0
    if selected_cluster_info['Cluster'][i] == 'NGC_2509':
        lim = 0.5
    if selected_cluster_info['Cluster'][i] == 'NGC_2627':
        lim = 3.75
    if selected_cluster_info['Cluster'][i] == 'NGC_2660':
        lim = 1.6
    if selected_cluster_info['Cluster'][i] == 'NGC_2682':
        lim = 1.5
    if selected_cluster_info['Cluster'][i] == 'NGC_6208':
        lim = 1.8
    if selected_cluster_info['Cluster'][i] == 'NGC_6791':
        lim = 1.5
    if selected_cluster_info['Cluster'][i] == 'NGC_6939':
        lim = 0.75
    if selected_cluster_info['Cluster'][i] == 'NGC_6940':
        lim = 1.25
    if selected_cluster_info['Cluster'][i] == 'NGC_752':
        lim = 1.6
    if selected_cluster_info['Cluster'][i] == 'NGC_7789':
        lim = 1.5
    if selected_cluster_info['Cluster'][i] == 'Ruprecht_68':
        lim = 0.6
    bsum = np.zeros(12,dtype=float)
    bins = np.linspace(0.0,lim,13) # 13 bin edges
    for j in range(len(ra)):
         d = math.sqrt( (pmra[j] - selected_cluster_info['PMra'][i])**2 + (pmdec[j] - selected_cluster_info['PMdec'][i])**2 )
         ipos = np.searchsorted(bins,d) - 1
         if ipos < 12:
              bsum[ipos] = bsum[ipos] + 1
    areas = math.pi*bins**2
    density_pm = []
    for k in range(12):
         bsum[k] = bsum[k]/(areas[k+1] - areas[k])
         density_pm.append(bsum[k])
    x = bins[-12:] - lim/24
    #plt.bar(x,bsum,width=(lim/10),color=lightblue)
    
    #use scipy.optimize.curve_fit to fit a curve to our data
    
    #define a function where the machine learning algorithm can input coefficients
    def func(x, a, b, c):
        return c + a * np.exp(-x*b)
    bound = ([0,0,min(density_pm)-10], [np.inf,np.inf,min(density_pm)])
    popt, pcov = curve_fit(func, x, density_pm, bounds=bound)
    plt.plot(np.linspace(0.0,lim,80), func(np.linspace(0.0,lim,80), *popt), 'r-', label=r'$%5.1f$exp$\{-%5.1fx\} + %5.1f$' % tuple(popt))
    coeff_dist_pm = 'a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)
    print(str(selected_cluster_info['Cluster'][i]) + ' proper motion: ' + str(coeff_dist_pm))    
    plt.xlabel(r'Distance in proper motion space (mas/yr)',fontsize=15)
    plt.ylabel(r'Star density (N yr$^2$ mas$^{-2}$)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.title('Proper Motion: ' + selected_cluster_info['Cluster'][i].replace('_', ' '))
    plt.ylim((0,max(density_pm)+max(density_pm)/12))
    plt.legend(frameon= False, fontsize = 15)
    #plt.savefig(selected_cluster_info['Cluster']+'_pm_dist.pdf', bbox_inches='tight', format='pdf')
    #plt.close()
    plt.show()
    
    # generate colors for first-pass estimate of proper motion member candidates
    #center around pmra pmdec
    pm_min = 0.0     # new
    pm_max = -1.0    # new
    category = []
    for j in range(len(ra)):
         d = math.sqrt( (pmra[j] - selected_cluster_info['PMra'][i])**2 + (pmdec[j] - selected_cluster_info['PMdec'][i])**2 )
         pm_max = max(d,pm_max)    # new. Find the radius of the farthest star.
         if d < selected_cluster_info['Radius'][i]:
              colors[j] = 'xkcd:black'
              category.append(2)
         elif d > 2*selected_cluster_info['Radius'][i]:
              colors[j] = 'xkcd:light green'
              category.append(0)
         else:
              colors[j] = 'xkcd:light purple'
              category.append(1)
    category = np.array(category)
    
    # 2. plot RA/dec and model cluster/field
    if selected_cluster_info['Cluster'][i] == 'NGC_7789':
        for alpha in range(len(ra)):
            if ra[alpha] > 300.00:
                ra[alpha] = ra[alpha] - 360.0000
    plt.scatter(ra,dec,c=colors,s=2)
    plt.xlabel(r'$\alpha$ ($^{\circ}$)',fontsize=15)
    plt.ylabel(r'$\delta$ ($^{\circ}$)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.title('Sky Positions: ' + str(selected_cluster_info['Cluster'][i]).replace('_', ' '))
    #plt.savefig(selected_cluster_info['Cluster'][i]+'_sky_pos.pdf', bbox_inches='tight', format='pdf')
    #plt.close()
    plt.show()
    opening_angles.append((max(dec)-min(dec))/2)
    
    #center around ra, dec from simbad
    # repeat code from above to translate to star density
    if selected_cluster_info['Cluster'][i] == 'NGC_6208' or selected_cluster_info['Cluster'][i] == 'NGC_2447' :
        lim=1.25
    if selected_cluster_info['Cluster'][i] == 'NGC_1907' or selected_cluster_info['Cluster'][i] == 'NGC_7789' or selected_cluster_info['Cluster'][i] == 'NGC_2682' or selected_cluster_info['Cluster'][i] == 'NGC_6939':
        lim = 2
    if selected_cluster_info['Cluster'][i] == 'NGC_2660' or selected_cluster_info['Cluster'][i] == 'FSR_1252' or selected_cluster_info['Cluster'][i] == 'NGC_6791' or selected_cluster_info['Cluster'][i] == 'Melotte_66' or selected_cluster_info['Cluster'][i] == 'NGC_1245' or selected_cluster_info['Cluster'][i] == 'NGC_2627' or selected_cluster_info['Cluster'][i] == 'Pismis_18' or selected_cluster_info['Cluster'][i] == 'Ruprecht_68' or selected_cluster_info['Cluster'][i] == 'King_5' or selected_cluster_info['Cluster'][i] == 'NGC_2477' or selected_cluster_info['Cluster'][i] == 'Czernik_37' or selected_cluster_info['Cluster'][i] == 'Melotte_71' or selected_cluster_info['Cluster'][i] == 'NGC_2236' or selected_cluster_info['Cluster'][i] == 'NGC_2506' or selected_cluster_info['Cluster'][i] == 'NGC_2509':
        lim = 3
    if selected_cluster_info['Cluster'][i] == 'Collinder_110' or selected_cluster_info['Cluster'][i] == 'NGC_752' or selected_cluster_info['Cluster'][i] == 'NGC_2420' or selected_cluster_info['Cluster'][i] == 'NGC_6940':
        lim = 4
    
    bins=np.linspace(0.0,lim*selected_cluster_info['Radius'][i],13) # 17 bin edges
    bsum = np.zeros(12,dtype=float)
    cosdec = math.cos(selected_cluster_info['Dec'][i]*math.pi/180.0) #use dec
    pos_min = 0.0
    pos_max = -1.0
    for j in range(len(ra)):
        if selected_cluster_info['Cluster'][i] == 'NGC_7789':
            d = math.sqrt( ((ra[j] - (selected_cluster_info['Ra'][i] - 360.000))*cosdec)**2 + (dec[j] - selected_cluster_info['Dec'][i])**2 )
        else:
            d = math.sqrt( ((ra[j] - selected_cluster_info['Ra'][i])*cosdec)**2 + (dec[j] - selected_cluster_info['Dec'][i])**2 )
        pos_max = max(d,pos_max)  # find maximum delta-angle from cluster center
        ipos = np.searchsorted(bins,d) - 1
        if ipos < 12:
             bsum[ipos] = bsum[ipos] + 1
    areas = math.pi*bins**2
    density_center = []
    for k in range(12):
        bsum[k] = bsum[k]/(areas[k+1] - areas[k])
        density_center.append(bsum[k])
    x = bins[-12:] - lim*selected_cluster_info['Radius'][i]/24
    #plt.bar(x,bsum,width=(lim*(selected_cluster_info['Radius'][i]/10)),color=pink)
    
    #use scipy.optimize.curve_fit to fit a curve to our data
    #define a function where the machine learning algorithm can input coefficients
    def func(x, a, b, c, d):
        if selected_cluster_info['Cluster'][i] == 'NGC_752' or selected_cluster_info['Cluster'][i] == 'NGC_2682' or selected_cluster_info['Cluster'][i] == 'NGC_1907' or selected_cluster_info['Cluster'][i] == 'Collinder_110' or selected_cluster_info['Cluster'][i] == 'NGC_2660' or selected_cluster_info['Cluster'][i] == 'NGC_2477' or selected_cluster_info['Cluster'][i] == 'NGC_2509' or selected_cluster_info['Cluster'][i] == 'NGC_6939' or selected_cluster_info['Cluster'][i] == 'Melotte_71' or selected_cluster_info['Cluster'][i] == 'Melotte_66' or selected_cluster_info['Cluster'][i] == 'King_5' or selected_cluster_info['Cluster'][i] == 'NGC_1245' or selected_cluster_info['Cluster'][i] == 'NGC_2236' or selected_cluster_info['Cluster'][i] == 'NGC_2420' or selected_cluster_info['Cluster'][i] == 'NGC_2506' or selected_cluster_info['Cluster'][i] == 'NGC_6791':
            return c + a*np.exp(-x*b)
        else:
            return d*np.exp(-((x-b)**2/(2*a**2))) + c
        #1d gaussian with mean=0 and sigma ~ cluster radius
    if selected_cluster_info['Cluster'][i] == 'NGC_6208':
        bound = ([0.,-0.1,min(density_center),10000],[0.1,0.1,120000.,30000])
    elif selected_cluster_info['Cluster'][i] == 'NGC_2477':
        bound = ([0.,0.,40000.,0.],[100000.,np.inf,min(density_center),np.inf])
    elif selected_cluster_info['Cluster'][i] == 'NGC_6940':
        bound = ([0,0,80000,10000],[3,0.01,np.inf,np.inf])
    elif selected_cluster_info['Cluster'][i] == 'NGC_2447':
        bound = ([0.025,0.03,min(density_center),16000.],[0.08,0.05,88000.,18000])
    elif selected_cluster_info['Cluster'][i] == 'NGC_2509':
        bound = (0,[np.inf,np.inf,min(density_center),np.inf])
        bound = ([0.,-0.1,261895.,70000.],[5,0.1,261896.,200000.])
    elif selected_cluster_info['Cluster'][i] == 'Ruprecht_68':
        #print('density_center = ',density_center)
        bound = ([0.05,0.,min(density_center),30000],[0.06,0.04,80000.,40000])
    elif selected_cluster_info['Cluster'][i] == 'IC_2714':
        bound = ([0,-0.5,6000.,20000.],[1.5,0.5,15000.,30000.])
    elif selected_cluster_info['Cluster'][i] == 'NGC_7789':
        bound = ([0,0,0,40000.],[0.4,0.1,100000.,100000])
    elif selected_cluster_info['Cluster'][i] == 'NGC_2682':
        bound = ([20000,0,3000.,0.],[40000,200,8000.,np.inf])
    elif selected_cluster_info['Cluster'][i] == 'NGC_2627':
        bound = ([0,0,58000.,20000],[0.1,0.06,61000.,23000])
    elif selected_cluster_info['Cluster'][i] == 'Czernik_37':
        bound = ([0,0,140000.,100000],[0.1,0.06,160000.,150000])
    elif selected_cluster_info['Cluster'][i] == 'FSR_1252':
        bound = ([0.03,0,48000.,20000],[0.05,0.06,70000.,40000])
    else:
        bound = (-10,[np.inf,np.inf,np.inf,np.inf])
    popt, pcov = curve_fit(func, x, density_center, bounds=bound)
    plt.plot(np.linspace(0.0,lim*selected_cluster_info['Radius'][i],80), func(np.linspace(0.0,lim*selected_cluster_info['Radius'][i],80), *popt), 'r-', label = r'$%5.2f$exp$\bigl\{-(\frac{(x-%5.2f)^2}{2 \cdot %5.3f^2})\bigl\} + %5.1f$' %tuple(popt))
    coeff_dist_cent = 'a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt)
    print(str(selected_cluster_info['Cluster'][i]) + ' Gaia position: ' + str(coeff_dist_cent))
    plt.xlabel(r'Distance from cluster center (deg)',fontsize=15)
    plt.ylabel(r'Star density (N/deg$^2$)',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    plt.ylim((0,max(density_center)+max(density_center)/16))
    plt.title('Gaia Position: '+selected_cluster_info['Cluster'][i].replace('_', ' '))
    plt.legend(frameon=False,fontsize=15)
    #plt.savefig(selected_cluster_info['Cluster'][i]+'_pos_dist.pdf', bbox_inches='tight', format='pdf')
    #plt.close()
    plt.show()
    
    #repeat position histogram+curve fit for panstarrrs ra,dec
    if selected_cluster_info['Cluster'][i] not in ['NGC_2660', 'NGC_2477', 'NGC_6208', 'Ruprecht_68', 'Melotte_66']:
        for j in range(len(panstarrs_ra)):
            if selected_cluster_info['Cluster'][i] == 'NGC_7789':
                d = math.sqrt( ((panstarrs_ra[j] - (selected_cluster_info['Ra'][i] - 360.000))*cosdec)**2 + (panstarrs_dec[j] - selected_cluster_info['Dec'][i])**2 )
            else:
                d = math.sqrt( ((panstarrs_ra[j] - selected_cluster_info['Ra'][i])*cosdec)**2 + (panstarrs_dec[j] - selected_cluster_info['Dec'][i])**2 )
            pos_max = max(d,pos_max)  # find maximum delta-angle from cluster center
            ipos = np.searchsorted(bins,d) - 1
            if ipos < 12:
                 bsum[ipos] = bsum[ipos] + 1
        areas = math.pi*bins**2
        density_center = []
        for k in range(12):
            bsum[k] = bsum[k]/(areas[k+1] - areas[k])
            density_center.append(bsum[k])
        x = bins[-12:] - lim*selected_cluster_info['Radius'][i]/24
        #plt.bar(x,bsum,width=(lim*(selected_cluster_info['Radius'][i]/10)),color=pink)
        
        #use scipy.optimize.curve_fit to fit a curve to our data
        
        #define a function where the machine learning algorithm can input coefficients
        def func(x, a, b, c, d):
            if selected_cluster_info['Cluster'][i] == 'NGC_752' or selected_cluster_info['Cluster'][i] == 'NGC_2682' or selected_cluster_info['Cluster'][i] == 'NGC_1907' or selected_cluster_info['Cluster'][i] == 'Collinder_110' or selected_cluster_info['Cluster'][i] == 'NGC_2660' or selected_cluster_info['Cluster'][i] == 'NGC_2477' or selected_cluster_info['Cluster'][i] == 'NGC_2509' or selected_cluster_info['Cluster'][i] == 'NGC_6939' or selected_cluster_info['Cluster'][i] == 'Melotte_71' or selected_cluster_info['Cluster'][i] == 'Melotte_66' or selected_cluster_info['Cluster'][i] == 'King_5' or selected_cluster_info['Cluster'][i] == 'NGC_1245' or selected_cluster_info['Cluster'][i] == 'NGC_2236' or selected_cluster_info['Cluster'][i] == 'NGC_2420' or selected_cluster_info['Cluster'][i] == 'NGC_2506' or selected_cluster_info['Cluster'][i] == 'NGC_6791':
                return c + a*np.exp(-x*b)
            else:
                return d*np.exp(-((x-b)**2/(2*a**2))) + c
            #1d gaussian with mean=0 and sigma ~ cluster radius
        if selected_cluster_info['Cluster'][i] == 'NGC_6208':
            bound = ([0.,-0.1,min(density_center),10000],[0.1,0.1,120000.,30000])
        elif selected_cluster_info['Cluster'][i] == 'NGC_2477':
            bound = ([0.,0.,40000.,0.],[100000.,np.inf,min(density_center),np.inf])
        elif selected_cluster_info['Cluster'][i] == 'NGC_6940':
            bound = ([0,0,80000,10000],[3,0.01,np.inf,np.inf])
        elif selected_cluster_info['Cluster'][i] == 'NGC_2447':
            bound = ([0.025,0.03,min(density_center),16000.],[0.08,0.05,88000.,18000])
        elif selected_cluster_info['Cluster'][i] == 'NGC_2509':
            bound = (0,[np.inf,np.inf,min(density_center),np.inf])
            bound = ([0.,-0.1,261895.,70000.],[5,0.1,261896.,200000.])
        elif selected_cluster_info['Cluster'][i] == 'Ruprecht_68':
            #print('density_center = ',density_center)
            bound = ([0.05,0.,min(density_center),30000],[0.06,0.04,80000.,40000])
        elif selected_cluster_info['Cluster'][i] == 'IC_2714':
            bound = ([0,-0.5,6000.,20000.],[1.5,0.5,15000.,30000.])
        elif selected_cluster_info['Cluster'][i] == 'NGC_7789':
            bound = ([0,0,0,40000.],[0.4,0.1,100000.,100000])
        elif selected_cluster_info['Cluster'][i] == 'NGC_2682':
            bound = ([20000,0,3000.,0.],[40000,200,8000.,np.inf])
        elif selected_cluster_info['Cluster'][i] == 'NGC_2627':
            bound = ([0,0,58000.,20000],[0.1,0.06,61000.,23000])
        elif selected_cluster_info['Cluster'][i] == 'Czernik_37':
            bound = ([0,0,140000.,100000],[0.1,0.06,160000.,150000])
        elif selected_cluster_info['Cluster'][i] == 'FSR_1252':
            bound = ([0.03,0,48000.,20000],[0.05,0.06,70000.,40000])
        else:
            bound = (-10,[np.inf,np.inf,np.inf,np.inf])
        popt, pcov = curve_fit(func, x, density_center, bounds=bound)
        plt.plot(np.linspace(0.0,lim*selected_cluster_info['Radius'][i],80), func(np.linspace(0.0,lim*selected_cluster_info['Radius'][i],80), *popt), 'r-', label = r'$%5.2f$exp$\bigl\{-(\frac{(x-%5.2f)^2}{2 \cdot %5.3f^2})\bigl\} + %5.1f$' %tuple(popt))
        coeff_dist_cent = 'a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt)
        print(str(selected_cluster_info['Cluster'][i]) + ' PanSTARRS position: ' + str(coeff_dist_cent))
        plt.xlabel(r'Distance from cluster center (deg)',fontsize=15)
        plt.ylabel(r'Star density (N/deg$^2$)',fontsize=15)
        plt.yticks(fontsize=15)
        plt.xticks(fontsize=15)
        plt.ylim((0,max(density_center)+max(density_center)/16))
        plt.title('PanSTARRS Position: '+selected_cluster_info['Cluster'][i].replace('_', ' '))
        plt.legend(frameon=False,fontsize=15)
        #plt.savefig(selected_cluster_info['Cluster'][i]+'_pos_dist.pdf', bbox_inches='tight', format='pdf')
        #plt.close()
        plt.show()
        
    # 3. and, finally, model parallax
    n_plx = []
    plx_min = np.min(plx)
    plx_min = min(plx_min,-1.0)
    plx_max = np.max(plx)
    plx_max = max(plx_max, 1.5)
    if selected_cluster_info['Cluster'][i] == 'NGC_6208':
        x2 = np.linspace(0.5,float(selected_cluster_info['Plx'][i])+1,125)
    else:
        x2 = np.linspace(float(selected_cluster_info['Plx'][i])-1,float(selected_cluster_info['Plx'][i])+1,125)
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
    if selected_cluster_info['Cluster'][i] == 'NGC_6208':                
        plt.hist([comp1,comp2,comp3],bins=125,range=(0.5,float(selected_cluster_info['Plx'][i])+1),stacked=True,color=[green,purple,darkblue], rwidth=1.5)
    else:
        plt.hist([comp1,comp2,comp3],bins=125,rwidth=1.2,range=(float(selected_cluster_info['Plx'][i])-1,float(selected_cluster_info['Plx'][i])+1),stacked=True,color=[green,purple,darkblue])
    
    #use scipy.optimize.curve_fit to fit a curve to our data
    
    #define a function where the machine learning algorithm can input coefficients
    #have the function try different types of equations for a given data set
    def func(x, a, b, c, d, e, f, g):
        return a*np.exp(-((x-b)**2/(2*c**2))) + d*np.exp(-((x-e)**2/(2*f**2))) + g
        #two 1d gaussians: one with mean = cluster parallax, one representing field star background plus constant pedestal
        #          a     b     c     d    e    f     g
    if selected_cluster_info['Cluster'][i] == 'Collinder_110':
        bound = ([550., 0.25, 0.2, 70, 0.43, 0.04, 70.], [570., 0.4, 0.44, 150, 0.53, 0.06, 100.])
    if selected_cluster_info['Cluster'][i] == 'Czernik_37':
        bound = ([95., 0.35, 0.1, 65, 0.38, 0.04, 15], [130., 0.5, .3, 100, 0.5, 0.06, 25])
    if selected_cluster_info['Cluster'][i] == 'FSR_1252':
        bound = ([40., 0.25, 0.3, 10, 0.27, 0.01, 5], [80., 0.45, 10., 40, 0.35, 0.3, 20])
    if selected_cluster_info['Cluster'][i] == 'King_5':
        bound = ([110, 0.3, 0.2, 20, 0.3, 0.1, 5], [150, 0.7, 0.4, 60, 0.5, 0.23, 30])
    if selected_cluster_info['Cluster'][i] == 'Melotte_66':
        bound = ([95., 0.2, 0.3, 120, 0.17, 0.01, 5.], [150., 0.5, 0.5, 170, 0.25, 0.1, 30.])
    if selected_cluster_info['Cluster'][i] == 'Melotte_71':
        bound = ([100., 0.28, 0.2, 85, 0.4, 0.01, 15], [180., 0.33, 0.36, 150, 0.45, 0.05, 22])
    if selected_cluster_info['Cluster'][i] == 'NGC_1245':
        bound = ([130., 0.35, 0.25, 80, 0.28, 0.07, 20.], [140., 0.4, 0.4, 150, 0.3, 0.1, 30.])
    if selected_cluster_info['Cluster'][i] == 'NGC_1907':
        bound = ([60., 0.35, 0.25, 30, 0.60, 0.053, 5.], [100., 0.4, 0.45, 50, 0.65, 0.07, 20.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2236':
        bound = ([60., 0.3, 0.25, 50, 0.28, 0.05, 2.], [90., 0.4, 0.5, 85, 0.35, 0.3, 15.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2420':
        bound = ([30., 0.35, 0.1, 50, 0.3, 0.01, 5.], [50., 0.45, 0.37, 80, 0.4, 0.05, 7.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2447':
        bound = ([3000, 0.15, 0.1, 200, 0.8, 0.1, 100.], [4500, 0.29, 0.5, 1000, 1.5, 0.4, 150.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2477':
        bound = ([650, 0.2, 0.2, 570, 0.65, 0.03, 100.], [1000, 0.35, 0.5, 650, 0.8, 0.05, 200.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2506':
        bound = ([100., 0.25, 0.2, 180, 0.26, 0.05, 10.], [130., 0.35, 0.35, 220, 0.3, 0.15, 30.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2509':
        bound = ([70., 0.25, 0.25, 60, 0.36, 0.01, 8.], [95., 0.33, 0.4, 150, 0.38, 0.05, 15.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2627':
        bound = ([150, 0.2, 0.1, 130, 0.52, 0.016, 20.], [250, 0.3, 0.35, 160, 0.56, 0.05, 30.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2660':
        bound = ([80., 0.25, 0.3, 55, 0.3, 0.01, 5.], [100., 0.27, 0.4, 70, 0.35, 0.07, 15.])
    if selected_cluster_info['Cluster'][i] == 'NGC_2682':
        bound = ([40., 0.4, 0.3, 130, 1.1, 0.05, 10.], [70., 0.7, 0.75, 160, 1.25, 0.2, 20.])
    if selected_cluster_info['Cluster'][i] == 'NGC_6208':
        bound = ([0., 0., 0.1, 120., 0.845, 0.02, 170], [3000, 0.6, 0.84, 190, 0.86, 0.1, 250])
    if selected_cluster_info['Cluster'][i] == 'NGC_6791':
        bound = ([160.,0.2,0.25,130,0.2,0., 5.], [300., 0.3, 0.4, 250, 0.275, 0.15, 50])
    if selected_cluster_info['Cluster'][i] == 'NGC_6939':
        bound = ([0,0,0,140,0.5,0.025,0], [120, 0.35, 0.6, 180, 0.7, 0.05, 15])
    if selected_cluster_info['Cluster'][i] == 'NGC_6940':
        bound = ([4800, 0.25, 0.23, 600, 0.85, 0.2, 150], [5500, 0.3, 0.28, 1000, 0.9, 0.35, 300])
    if selected_cluster_info['Cluster'][i] == 'NGC_752':
        bound = ([350, 0.5, 0.9, 50, 2.25, 0.05, 10.], [500, 2.0, 1.2, 70, 2.4, 0.2, 50.])
    if selected_cluster_info['Cluster'][i] == 'NGC_7789':
        bound = ([620, 0.25, 0.2, 540, 0.46, 0.01, 80], [670, 0.34, 0.35, 550, 0.55, 0.1, 100])
    if selected_cluster_info['Cluster'][i] == 'Pismis_18':
        bound = ([140, 0.32, 0.2, 43.5, 0.358, 0.04, 10], [150, 0.35, 0.45, 100, 0.45, 0.13, 30])
    if selected_cluster_info['Cluster'][i] == 'Ruprecht_68':
        bound = ([200, 0.2, 0.2, 75, 0.32, 0.052, 20], [230, 0.27, 0.35, 80, 0.37, 0.06, 30])
    popt, pcov = curve_fit(func, x2, n_plx, bounds=bound)
    plt.plot(x2, func(x2, *popt), 'r-', label= r'$%5.1f$exp$\bigl\{\frac{-(x-%5.2f)^2}{2 \cdot %5.3f^2}\bigr\} + %5.1f$exp$\bigl\{\frac{-(x-%5.2f)^2}{2 \cdot %5.3f^2}\bigr\} + %5.1f$' %tuple(popt))
    plt.title("Parallax: " + str(selected_cluster_info['Cluster'][i].replace('_', ' ')))
    coeff_plx = 'a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f, e=%5.3f, f=%5.3f, g=%5.3f' % tuple(popt)
    plt.legend(frameon=False,fontsize=12)
    plt.ylim([0.0,1400])
    #plt.xlim([0.75,1.25])
    plt.xlabel(r'Parallax (mas)',fontsize=15)
    plt.ylabel('N',fontsize=15)
    plt.yticks(fontsize=15)
    plt.xticks(fontsize=15)
    #plt.savefig(selected_cluster_info['Cluster'][i]+'_plx_dist.pdf', bbox_inches='tight', format='pdf')
    #plt.close()
    plt.show()
 
    # 4. put it all together for a probability. 
     
    #Now that we have coefficients for the curves for each plot in each cluster,
    #we can take the coefficients and write equations for each plot.

    # compute N_cluster, N_field, and normalization factors for the cluster and field functional distributions
    # Ensure that the integral over all stars: A* INT_0^max f(v)dv = 1, where A is a constant
    # Ensure this separately for cluster and field (but the dv and range must be the same for both)
    bin_edges = np.linspace(pm_min,pm_max,num=3951,endpoint=True)
    # note large number of points. We need to fairly sample the tiny cluster locus AND the huge field locus in PM space
    bin_centers = bin_edges[0:-1] + 0.5*( bin_edges[1]-bin_edges[0])
    dbin =  0.5*( bin_edges[1]-bin_edges[0])
    pmsumc = 0.0; pmsumf = 0.0
    for j in range(len(bin_centers)):
        dpm = bin_centers[j]
        phi_cpm = [(12429.746 * np.exp(-dpm*11.399)) , #Ruprecht 68
        (7889.704 * np.exp(-dpm*10.716) ), #Pismis 18
        (41242.811 * np.exp(-dpm*7.103)) , #NGC 7789
        (820.675* np.exp(-dpm*3.205)) , #NGC 752
        (7707.473 * np.exp(-dpm*6.728)) , #NGC 6940
        (13557.738 * np.exp(-dpm*9.377)) , #NGC 6939
        (47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        (5469.882 * np.exp(-dpm*7.669)) , #NGC 6208
        (6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        (12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (3115.683 * np.exp(-dpm*1.959)) , #NGC 2627
        (17280.461 * np.exp(-dpm*17.437)) , #NGC 2509
        (33119.650 * np.exp(-dpm*8.828)) , #NGC 2506
        (21370.314 * np.exp(-dpm*5.392)) , #NGC 2477
        (8745.948 * np.exp(-dpm*6.120)) , #NGC 2447
        (12374.103 * np.exp(-dpm*9.113)) , #NGC 2420
        (10819.833 * np.exp(-dpm*9.048)) , #NGC 2236
        (3789.570 * np.exp(-dpm*6.141)) , #NGC 1907
        (13338.112 * np.exp(-dpm*7.465)) , #NGC 1245
        (16787.800 * np.exp(-dpm*9.424)) , #Melotte 71
        (23473.741 * np.exp(-dpm*9.355)) , #Melotte 66
        (7973.347 * np.exp(-dpm*9.869)) , #King 5
        (3992.037 * np.exp(-dpm*9.735)) , #FSR 1252
        (5841.210 * np.exp(-dpm*10.104)) , #Czernik 37
        (28809.904 * np.exp(-dpm*10.447)) ] #Collinder 110 
                
        peak_cpm = [(12429.746) , #Ruprecht 68
        (7889.704), #Pismis 18
        (41242.811) , #NGC 7789
        (820.675) , #NGC 752
        (7707.473 ) , #NGC 6940
        (13557.738) , #NGC 6939
        (47524.902), #NGC 6791
        (5469.882) , #NGC 6208
        (6386.179), #NGC 2682
        (12625.640 ), #NGC 2660
        (3115.683 ) , #NGC 2627
        (17280.461 ) , #NGC 2509
        (33119.650 ) , #NGC 2506
        (21370.314 ) , #NGC 2477
        (8745.948 ) , #NGC 2447
        (12374.103 ) , #NGC 2420
        (10819.833 ) , #NGC 2236
        (3789.570) , #NGC 1907
        (13338.112 ) , #NGC 1245
        (16787.800) , #Melotte 71
        (23473.741), #Melotte 66
        (7973.347) , #King 5
        (3992.037 ) , #FSR 1252
        (5841.210 ) , #Czernik 37
        (28809.904 ) ] #Collinder 110 

        phi_fpm = [111.261, 317.237, 1225.596, 15.690, 588.992, 99.831, 
                   130.203, 682.752, 14.800, 155.176, 53.852, 159.432, 
                   117.545, 854.059, 1355.076, 120.148, 428.479, 123.164, 
                   389.433, 136.353, 139.946, 609.280, 143.573, 321.770, 955.309]

        pmsumc = pmsumc + dbin*phi_cpm[i]
        pmsumf = pmsumf + dbin*phi_fpm[i]
        pmpeakc = peak_cpm[i]
        pmpeakf = phi_fpm[i]
   
    N_cluster = pmsumc  # N_cluster and N_field should be broadly proportional to the actual values
    N_field = pmsumf
    Anorm_pmC = 1/pmsumc
    Anorm_pmF = 1/pmsumf
    Anorm_N = 1.0/(pmsumc + pmsumf)
    n_c = Anorm_N*pmsumc  # normalized number of cluster members
    n_f = Anorm_N*pmsumf  # normalized number of field members
    print('PM: Ncluster, Nfield, Anorms = ',N_cluster,N_field,Anorm_pmC,Anorm_pmF)
    print( 'PM: Number of cluster, field, and total = ',n_c,n_f,n_c+n_f)
    print('(The last value should equal 1)')
    
    # Repeat for position
    bin_edges = np.linspace(pos_min,pos_max,num=151,endpoint=True)
    bin_centers = bin_edges[0:-1] + 0.5*( bin_edges[1]-bin_edges[0])
    dbin =  0.5*( bin_edges[1]-bin_edges[0])
    possumc = 0.0; possumf = 0.0
    for j in range(len(bin_centers)):
        dpos = bin_centers[j]
        phi_cpos = [(30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))), #Ruprecht 68
                     (170521.883 * np.exp(-((dpos-0.002)**2/(2*0.016**2)))), #Pismis 18
                     (55926.096 * np.exp(-((dpos)**2/(2*0.090**2)))), #NGC 7789
                     (2715.388 * np.exp(-dpos*0.048)), #NGC 752
                     (12646.468 * np.exp(-((dpos)**2/(2*0.280**2)))), #NGC 6940
                     (51920.447 * np.exp(-dpos*22.830)), #NGC 6939
                     (479976.011 * np.exp(-dpos*21.117)), #NGC 6791
                     (25208.794 * np.exp(-((dpos+0.091)**2/(2*0.1**2)))), #NGC 6208
                     (20857.247 * np.exp(-dpos*14.331)), #NGC 2682
                     (414912.100 * np.exp(-dpos*59.109)), #NGC 2660
                     (23000.000 * np.exp(-((dpos-0.022)**2/(2*0.039**2)))), #NGC 2627
                     (62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (163447.774 * np.exp((-dpos*17.868))), #NGC 2506
                     (84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     (16000.000 * np.exp(-((dpos-0.030)**2/(2*0.025**2)))), #NGC 2447
                     (139661.988 * np.exp(-dpos*38.210)), #NGC 2420
                     (62730.170 * np.exp(-dpos*20.570)), #NGC 2236
                     (58995.329 * np.exp(-dpos*39.285)), #NGC 1907
                     (77263.108 * np.exp(-dpos*17.225)), #NGC 1245
                     (76703.328 * np.exp(-dpos*18.636)), #Melotte 71
                     (118498.343 * np.exp(-dpos*15.582)), #Melotte 66
                     (70172.124 * np.exp(-dpos*27.637)), #King 5
                     (27350.619*np.exp(-((dpos-0.000)**2/(2*0.030**2)))), #FSR 1252
                     (117912.880*np.exp(-((dpos-0.004)**2/(2*0.019**2)))), #Czernik 37
                     (30646.548 * np.exp(-dpos*6.359))] #Collinder 110 
        
        peak_cpos = [(30000.000), #Ruprecht 68
                     (170521.883), #Pismis 18
                     (55926.096), #NGC 7789
                     (2715.388), #NGC 752
                     (12646.468), #NGC 6940
                     (51920.447), #NGC 6939
                     (479976.011), #NGC 6791
                     (25208.794), #NGC 6208
                     (20857.247), #NGC 2682
                     (414912.100), #NGC 2660
                     (23000.000), #NGC 2627
                     (62734.146), #NGC 2509
                     (163447.774), #NGC 2506
                     (84765.416), #NGC 2477
                     (16000.000), #NGC 2447
                     (139661.988), #NGC 2420
                     (62730.170), #NGC 2236
                     (58995.329), #NGC 1907
                     (77263.108), #NGC 1245
                     (76703.328), #Melotte 71
                     (118498.343), #Melotte 66
                     (70172.124), #King 5
                     (27350.619), #FSR 1252
                     (117912.880), #Czernik 37
                     (30646.548)] #Collinder 110 
        
        phi_fpos = [73332.225, #Ruprecht 68
                    261896.000, #Pismis 18
                    66608.971, #NGC 7789
                    2039.455, #NGC 752
                    93792.077, #NGC 6940
                    32503.010, #NGC 6939
                    53557.792, #NGC 6791
                    106234.735, #NGC 6208
                    6456.678, #NGC 2682
                    139757.004, #NGC 2660
                    58828.279, #NGC 2627
                    52748.495, #NGC 2509
                    28867.075, #NGC 2506
                    72807.519, #NGC 2477
                    85553.016, #NGC 2447
                    14699.250, #NGC 2420
                    41666.921, #NGC 2236
                    44185.987, #NGC 1907
                    27484.123, #NGC 1245
                    52110.055, #Melotte 71
                    23179.392, #Melotte 66
                    34381.117, #King 5
                    65091.820, #FSR 1252
                    143655.820, #Czernik 37
                    37600.457] #Collinder 110
        
        possumc = possumc + dbin*phi_cpos[i]
        possumf = possumf + dbin*phi_fpos[i]
        pospeakc = peak_cpos[i]
        pospeakf = phi_fpos[i]
        
    Anorm_posC = 1/possumc # to normalize binned function to integral one
    Anorm_posF = 1/possumf
    #print('POS: Ncluster, Nfield, Anorms = ',possumc,possumf,Anorm_posC,Anorm_posF)
      
    # Repeat for parallax
    bin_edges = np.linspace(plx_min,plx_max,num=951,endpoint=True)
    bin_centers = bin_edges[0:-1] + 0.5*( bin_edges[1]-bin_edges[0])
    dbin =  0.5*( bin_edges[1]-bin_edges[0])
    plxsumc = 0.0; plxsumf = 0.0
    for j in range(len(bin_centers)):
        dplx = bin_centers[j]# - selected_cluster_info['Plx'][i]
        phi_cpx = [(75.000*np.exp(-((dplx-.320)**2/(2*0.052**2)))), #Ruprecht 68
                   (43.500*np.exp(-((dplx-0.358)**2/(2*0.040**2)))), #Pismis 18
                   (540.000*np.exp(-((dplx-0.472)**2/(2*0.020**2)))), #NGC 7789
                   (50.000*np.exp(-((dplx-2.293)**2/(2*0.050**2)))), #NGC 752
                   (925.842*np.exp(-((dplx-0.850)**2/(2*0.273**2)))), #NGC 6940
                   (146.938*np.exp(-((dplx-0.518)**2/(2*0.026**2)))), #NGC 6939
                   (247.012*np.exp(-((dplx-0.211)**2/(2*0.073**2)))), #NGC 6791
                   (120.000*np.exp(-((dplx-0.855)**2/(2*0.020**2)))), #NGC 6208
                   (130.000*np.exp(-((dplx-1.144)**2/(2*0.050**2)))), #NGC 2682
                   (57.559*np.exp(-((dplx-0.307)**2/(2*0.070**2)))), #NGC 2660
                   (130.000*np.exp(-((dplx-0.545)**2/(2*0.050**2)))), #NGC 2627
                   (65.731*np.exp(-((dplx-0.360)**2/(2*0.018**2)))), #NGC 2509
                   (186.345*np.exp(-((dplx-0.267)**2/(2*0.051**2)))), #NGC 2506
                   (570.000*np.exp(-((dplx-0.686)**2/(2*0.041**2)))), #NGC 2477
                   (699.896*np.exp(-((dplx-0.800)**2/(2*0.260**2)))), #NGC 2447
                   (65.812*np.exp(-((dplx-0.379)**2/(2*0.033**2)))), #NGC 2420
                   (70.304*np.exp(-((dplx-0.301)**2/(2*0.123**2)))), #NGC 2236
                   (30.000*np.exp(-((dplx-0.623)**2/(2*0.053**2)))), #NGC 1907
                   (112.810*np.exp(-((dplx-0.285)**2/(2*0.070**2)))), #NGC 1245
                   (87.411*np.exp(-((dplx-0.432)**2/(2*0.020**2)))), #Melotte 71
                   (138.563*np.exp(-((dplx-0.188)**2/(2*0.063**2)))), #Melotte 66
                   (29.727*np.exp(-((dplx-0.336)**2/(2*0.100**2)))), #King 5
                   (37.282*np.exp(-((dplx-0.285)**2/(2*0.066**2)))), #FSR 1252
                   (65.000*np.exp(-((dplx-0.380)**2/(2*0.055**2)))), #Czernik 37
                   (144.891*np.exp(-((dplx-0.430)**2/(2*0.060**2))))] #Collinder 110
        
        peak_cpx = [(75.000),
                   (43.500),
                   (540.000),
                   (50.000),
                   (952.959),
                   (146.938),
                   (247.012),
                   (120.000),
                   (130.000),
                   (57.559),
                   (130.000),
                   (65.731),
                   (186.345),
                   (570.000),
                   (699.896),
                   (65.812),
                   (70.304),
                   (30.000),
                   (112.810),
                   (87.411),
                   (138.563),
                   (29.727),
                   (37.282),
                   (65.000),
                   (144.891)]

        phi_fpx = [(200.000*np.exp(-((dplx-0.256)**2/(2*0.232**2))) + 23.672),
                   (140.000*np.exp(-((dplx-0.339)**2/(2*0.200**2))) + 21.057),
                   (663.369*np.exp(-((dplx-0.316)**2/(2*0.234**2))) + 83.308),
                   (493.452*np.exp(-((dplx-0.500)**2/(2*0.900**2))) + 31.763),
                   (4900.910*np.exp(-((dplx-0.263)**2/(2*0.230**2))) + 267.896),
                   (93.331*np.exp(-((dplx-0.350)**2/(2*0.266**2))) + 13.438),
                   (163.732*np.exp(-((dplx-0.273)**2/(2*0.275**2))) + 19.346),
                   (2628.768*np.exp(-((dplx)**2/(2*0.486**2))) + 193.809),
                   (48.367*np.exp(-((dplx-0.563)**2/(2*0.622**2))) + 10.000),
                   (80.000*np.exp(-((dplx-0.267)**2/(2*0.300**2))) + 5.618),
                   (247.682*np.exp(-((dplx-0.240)**2/(2*0.169**2))) + 30.000),
                   (94.471*np.exp(-((dplx-0.297)**2/(2*0.250**2))) + 9.676),
                   (125.332*np.exp(-((dplx-0.313)**2/(2*0.228**2))) + 19.544),
                   (818.705*np.exp(-((dplx-0.245)**2/(2*0.245**2))) + 100.748),
                   (3391.611*np.exp(-((dplx-0.258)**2/(2*0.192**2))) + 145.156),
                   (38.696*np.exp(-((dplx-0.369)**2/(2*0.216**2))) + 5.554),
                   (61.964*np.exp(-((dplx-0.351)**2/(2*0.358**2))) + 9.219),
                   (86.597*np.exp(-((dplx-0.372)**2/(2*0.285**2))) + 14.261),
                   (130.000*np.exp(-((dplx-0.370)**2/(2*0.281**2))) + 20.000),
                   (157.641*np.exp(-((dplx-0.317)**2/(2*0.235**2))) + 19.697),
                   (96.834*np.exp(-((dplx-0.337)**2/(2*0.300**2))) + 9.903),
                   (110.000*np.exp(-((dplx-0.372)**2/(2*0.278**2))) + 17.263),
                   (45.188*np.exp(-((dplx-0.324)**2/(2*0.314**2))) + 6.308),
                   (95.000*np.exp(-((dplx-0.421)**2/(2*0.267**2))) + 15.000),
                   (570.000*np.exp(-((dplx-0.336)**2/(2*0.328**2))) + 100.000) ]
        
        peak_fpx = [(200.000 + 23.672),
                   (140.000 + 21.057),
                   (663.369 + 83.308),
                   (493.452 + 31.763),
                   (4900.910 + 267.896),
                   (93.331 + 13.438),
                   (163.732 + 19.346),
                   (2628.768 + 193.809),
                   (48.367 + 10.000),
                   (80.000 + 5.618),
                   (247.682 + 30.000),
                   (94.471 + 9.676),
                   (125.333 + 19.544),
                   (818.705 + 100.748),
                   (3391.611 + 145.156),
                   (38.698 + 5.554),
                   (61.964 + 9.219),
                   (86.597 + 14.261),
                   (130.000 + 20.000),
                   (157.641 + 19.697),
                   (96.834 + 9.903),
                   (110.000 + 17.263),
                   (45.188 + 6.308),
                   (95.000 + 15.000),
                   (570.000 + 100.000)]
        
        plxsumc = plxsumc + dbin*phi_cpx[i]
        plxsumf = plxsumf + dbin*phi_fpx[i]
        plxpeakc = peak_cpx[i]
        plxpeakf = peak_fpx[i]
        
    Anorm_plxC = 1/plxsumc # to normalize binned function to integral one
    Anorm_plxF = 1/plxsumf
    #print('Parallax: Ncluster, Nfield, Anorms = ',plxsumc,plxsumf,Anorm_plxC,Anorm_plxF)

    # what probability to you get if you use the peak values?
    if selected_cluster_info['Cluster'][i] == 'NGC_6208' or selected_cluster_info['Cluster'][i] == 'NGC_6940' or selected_cluster_info['Cluster'][i] == 'Collinder_110' or selected_cluster_info['Cluster'][i] == 'King_5' or selected_cluster_info['Cluster'][i] == 'NGC_2447':
        zAND = (pmpeakc)/((pmpeakc)+(pmpeakf))
    else:
        zAND = (pmpeakc*pospeakc*plxpeakc)/((pmpeakc*pospeakc*plxpeakc)+(pmpeakf*pospeakf*plxpeakf))
    zOR =1.0 - ( (1.0-pmpeakc/(pmpeakc+pmpeakf))*(1.0-pospeakc/(pospeakc+pospeakf))*(1.0 - plxpeakc/(plxpeakc+plxpeakf)) )
    print('PM  peak prob', pmpeakc/(pmpeakc+pmpeakf) )
    print('POS peak prob', pospeakc/(pospeakc+pospeakf) )
    print('PLX peak prob', plxpeakc/(plxpeakc+plxpeakf) )
    print('AND logic for the center of the distribution', zAND)
    print('OR logic for the center of the distribution',zOR)
    
    #Put each equation for each category and each cluster in a structured array
    #Every equation has a 'c' coefficient: that will act as the constant field star background
    #Prob = (f(x) - c) / (f(x))   [ old "or" logic]
    #x will be dpm, dpos, dplx in each equation  
    ppm = []; ppmc=[];ppmf=[]
    ppos = []; pposc=[];pposf=[]
    pplx = []; pplxc=[];pplxf=[]
    total_probs = []
    total_probs2 = []
    for j in range(len(ra)):
        
        # --------
        # PROPER MOTION
        # --------
      
        dpm = math.sqrt( (pmra[j] - selected_cluster_info['PMra'][i])**2 + (pmdec[j] - selected_cluster_info['PMdec'][i])**2 )
        
        ppm_eqns = [(12429.746 * np.exp(-dpm*11.399)) / (111.261 + 12429.746 * np.exp(-dpm*11.399)), #Ruprecht 68
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
        (28809.904 * np.exp(-dpm*10.447)) / (955.309 + 28809.904 * np.exp(-dpm*10.447))] #Collinder 110
        
        phi_cpm = [(12429.746 * np.exp(-dpm*11.399)) , #Ruprecht 68
        (7889.704 * np.exp(-dpm*10.716) ), #Pismis 18
        (41242.811 * np.exp(-dpm*7.103)) , #NGC 7789
        (820.675* np.exp(-dpm*3.205)) , #NGC 752
        (7707.473 * np.exp(-dpm*6.728)) , #NGC 6940
        (13557.738 * np.exp(-dpm*9.377)) , #NGC 6939
        (47524.902 * np.exp(-dpm*8.087)), #NGC 6791
        (5469.882 * np.exp(-dpm*7.669)) , #NGC 6208
        (6386.179 * np.exp(-dpm*4.596)), #NGC 2682
        (12625.640 * np.exp(-dpm*10.175)), #NGC 2660
        (3115.683 * np.exp(-dpm*1.959)) , #NGC 2627
        (17280.461 * np.exp(-dpm*17.437)) , #NGC 2509
        (33119.650 * np.exp(-dpm*8.828)) , #NGC 2506
        (21370.314 * np.exp(-dpm*5.392)) , #NGC 2477
        (8745.948 * np.exp(-dpm*6.120)) , #NGC 2447
        (12374.103 * np.exp(-dpm*9.113)) , #NGC 2420
        (10819.833 * np.exp(-dpm*9.048)) , #NGC 2236
        (3789.570 * np.exp(-dpm*6.141)) , #NGC 1907
        (13338.112 * np.exp(-dpm*7.465)) , #NGC 1245
        (16787.800 * np.exp(-dpm*9.424)) , #Melotte 71
        (23473.741 * np.exp(-dpm*9.355)) , #Melotte 66
        (7973.347 * np.exp(-dpm*9.869)) , #King 5
        (3992.037 * np.exp(-dpm*9.735)) , #FSR 1252
        (5841.210 * np.exp(-dpm*10.104)) , #Czernik 37
        (28809.904 * np.exp(-dpm*10.447)) ] #Collinder 110 

        phi_fpm = [111.261, 317.237, 1225.596, 15.690, 588.992, 99.831, 
                   130.203, 682.752, 14.800, 155.176, 53.852, 159.432, 
                   117.545, 854.059, 1355.076, 120.148, 428.479, 123.164, 
                   389.433, 136.353, 139.946, 609.280, 143.573, 321.770, 955.309]
        
        ppm.append(ppm_eqns[i]) # add a proper motion probability instance
        #ppmc.append(Anorm_pmC*phi_cpm[i]) # normalized cluster probability instance
        #ppmf.append(Anorm_pmF*phi_fpm[i]) # normalized field probability instance
        ppmc.append(phi_cpm[i])
        ppmf.append(phi_fpm[i])

        # --------
        # POSITION
        # --------
        
        if selected_cluster_info['Cluster'][i] == 'NGC_7789':
            dpos = math.sqrt( ((ra[j] - (selected_cluster_info['Ra'][i]-360.000))*cosdec)**2 + (dec[j] - selected_cluster_info['Dec'][i])**2 )
        else: dpos = math.sqrt( ((ra[j] - selected_cluster_info['Ra'][i])*cosdec)**2 + (dec[j] - selected_cluster_info['Dec'][i])**2 )
        
        ppos_eqns = [(30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))) / (73332.225 + 30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))), #Ruprecht 68
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
                     (30646.548 * np.exp(-dpos*6.359)) / (37600.457 + 30646.548 * np.exp(-dpos*6.359)) ] #Collidner 110


        phi_cpos = [(30000.000 * np.exp(-((dpos-0.006)**2/(2*0.050**2)))), #Ruprecht 68
                     (170521.883 * np.exp(-((dpos-0.002)**2/(2*0.016**2)))), #Pismis 18
                     (55926.096 * np.exp(-((dpos)**2/(2*0.090**2)))), #NGC 7789
                     (2715.388 * np.exp(-dpos*0.048)), #NGC 752
                     (12646.468 * np.exp(-((dpos)**2/(2*0.280**2)))), #NGC 6940
                     (51920.447 * np.exp(-dpos*22.830)), #NGC 6939
                     (479976.011 * np.exp(-dpos*21.117)), #NGC 6791
                     (25208.794 * np.exp(-((dpos+0.091)**2/(2*0.1**2)))), #NGC 6208
                     (20857.247 * np.exp(-dpos*14.331)), #NGC 2682
                     (414912.100 * np.exp(-dpos*59.109)), #NGC 2660
                     (23000.000 * np.exp(-((dpos-0.022)**2/(2*0.039**2)))), #NGC 2627
                     (62734.146 * np.exp(-dpos*39.211)), #NGC 2509
                     (163447.774 * np.exp((-dpos*17.868))), #NGC 2506
                     (84765.416 * np.exp(-dpos*8.043)), #NGC 2477
                     (16000.000 * np.exp(-((dpos-0.030)**2/(2*0.025**2)))), #NGC 2447
                     (139661.988 * np.exp(-dpos*38.210)), #NGC 2420
                     (62730.170 * np.exp(-dpos*20.570)), #NGC 2236
                     (58995.329 * np.exp(-dpos*39.285)), #NGC 1907
                     (77263.108 * np.exp(-dpos*17.225)), #NGC 1245
                     (76703.328 * np.exp(-dpos*18.636)), #Melotte 71
                     (118498.343 * np.exp(-dpos*15.582)), #Melotte 66
                     (70172.124 * np.exp(-dpos*27.637)), #King 5
                     (27350.619*np.exp(-((dpos-0.000)**2/(2*0.030**2)))), #FSR 1252
                     (117912.880*np.exp(-((dpos-0.004)**2/(2*0.019**2)))), #Czernik 37
                     (30646.548 * np.exp(-dpos*6.359))] #Collinder 110 
        
        phi_fpos = [73332.225, #Ruprecht 68
                    261896.000, #Pismis 18
                    66608.971, #NGC 7789
                    2039.455, #NGC 752
                    93792.077, #NGC 6940
                    32503.010, #NGC 6939
                    53557.792, #NGC 6791
                    106234.735, #NGC 6208
                    6456.678, #NGC 2682
                    139757.004, #NGC 2660
                    58828.279, #NGC 2627
                    52748.495, #NGC 2509
                    28867.075, #NGC 2506
                    72807.519, #NGC 2477
                    85553.016, #NGC 2447
                    14699.250, #NGC 2420
                    41666.921, #NGC 2236
                    44185.987, #NGC 1907
                    27484.123, #NGC 1245
                    52110.055, #Melotte 71
                    23179.392, #Melotte 66
                    34381.117, #King 5
                    65091.820, #FSR 1252
                    143655.820, #Czernik 37
                    37600.457] #Collinder 110
        
        ppos.append(ppos_eqns[i]) #distance from cluster center
        #pposc.append(Anorm_posC*phi_cpos[i])
        #pposf.append(Anorm_posF*phi_fpos[i])
        pposc.append(phi_cpos[i])
        pposf.append(phi_fpos[i])
        
        # --------
        # PARALLAX
        # --------
        
        dplx = plx[j] #- selected_cluster_info['Plx'][i] 
        
        pplx_eqns = [(75.000*np.exp(-((dplx-.320)**2/(2*0.052**2)))) / (75.000*np.exp(-((dplx-.320)**2/(2*0.052**2))) + 200.000*np.exp(-((dplx-0.256)**2/(2*0.232**2))) + 23.672), #Ruprecht 68
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

        phi_cpx = [(75.000*np.exp(-((dplx-.320)**2/(2*0.052**2)))), #Ruprecht 68
                   (43.500*np.exp(-((dplx-0.358)**2/(2*0.040**2)))), #Pismis 18
                   (540.000*np.exp(-((dplx-0.472)**2/(2*0.020**2)))), #NGC 7789
                   (50.000*np.exp(-((dplx-2.293)**2/(2*0.050**2)))), #NGC 752
                   (925.842*np.exp(-((dplx-0.850)**2/(2*0.273**2)))), #NGC 6940
                   (146.938*np.exp(-((dplx-0.518)**2/(2*0.026**2)))), #NGC 6939
                   (247.012*np.exp(-((dplx-0.211)**2/(2*0.073**2)))), #NGC 6791
                   (120.000*np.exp(-((dplx-0.855)**2/(2*0.020**2)))), #NGC 6208
                   (130.000*np.exp(-((dplx-1.144)**2/(2*0.050**2)))), #NGC 2682
                   (57.559*np.exp(-((dplx-0.307)**2/(2*0.070**2)))), #NGC 2660
                   (130.000*np.exp(-((dplx-0.545)**2/(2*0.050**2)))), #NGC 2627
                   (65.731*np.exp(-((dplx-0.360)**2/(2*0.018**2)))), #NGC 2509
                   (186.345*np.exp(-((dplx-0.267)**2/(2*0.051**2)))), #NGC 2506
                   (570.000*np.exp(-((dplx-0.686)**2/(2*0.041**2)))), #NGC 2477
                   (699.896*np.exp(-((dplx-0.800)**2/(2*0.260**2)))), #NGC 2447
                   (65.812*np.exp(-((dplx-0.379)**2/(2*0.033**2)))), #NGC 2420
                   (70.304*np.exp(-((dplx-0.301)**2/(2*0.123**2)))), #NGC 2236
                   (30.000*np.exp(-((dplx-0.623)**2/(2*0.053**2)))), #NGC 1907
                   (112.810*np.exp(-((dplx-0.285)**2/(2*0.070**2)))), #NGC 1245
                   (87.411*np.exp(-((dplx-0.432)**2/(2*0.020**2)))), #Melotte 71
                   (138.563*np.exp(-((dplx-0.188)**2/(2*0.063**2)))), #Melotte 66
                   (29.727*np.exp(-((dplx-0.336)**2/(2*0.100**2)))), #King 5
                   (37.282*np.exp(-((dplx-0.285)**2/(2*0.066**2)))), #FSR 1252
                   (65.000*np.exp(-((dplx-0.380)**2/(2*0.055**2)))), #Czernik 37
                   (144.891*np.exp(-((dplx-0.430)**2/(2*0.060**2))))] #Collinder 110


        phi_fpx = [(200.000*np.exp(-((dplx-0.256)**2/(2*0.232**2))) + 23.672), #Ruprecht 68
                   (140.000*np.exp(-((dplx-0.339)**2/(2*0.200**2))) + 21.057), #Pismis 18
                   (663.369*np.exp(-((dplx-0.316)**2/(2*0.234**2))) + 83.308), #NGC 7789
                   (493.452*np.exp(-((dplx-0.500)**2/(2*0.900**2))) + 31.763), #NGC 752
                   (4900.910*np.exp(-((dplx-0.263)**2/(2*0.230**2))) + 267.896), #NGC 6949
                   (93.331*np.exp(-((dplx-0.350)**2/(2*0.266**2))) + 13.438), #NGC 6939
                   (163.732*np.exp(-((dplx-0.273)**2/(2*0.275**2))) + 19.346), #NGC 6791
                   (2628.768*np.exp(-((dplx)**2/(2*0.486**2))) + 193.809), #NGC 6208
                   (48.367*np.exp(-((dplx-0.563)**2/(2*0.622**2))) + 10.000), #NGC 2682
                   (80.000*np.exp(-((dplx-0.267)**2/(2*0.300**2))) + 5.618), #NGC 2660
                   (247.682*np.exp(-((dplx-0.240)**2/(2*0.169**2))) + 30.000), #NGC 2627
                   (94.471*np.exp(-((dplx-0.297)**2/(2*0.250**2))) + 9.676), #NGC 2509
                   (125.332*np.exp(-((dplx-0.313)**2/(2*0.228**2))) + 19.544), #NGC 2506
                   (818.705*np.exp(-((dplx-0.245)**2/(2*0.245**2))) + 100.748), #NGC 2477
                   (3391.611*np.exp(-((dplx-0.258)**2/(2*0.192**2))) + 145.156), #NGC 2447
                   (38.696*np.exp(-((dplx-0.369)**2/(2*0.216**2))) + 5.554), #NGC 2420
                   (61.964*np.exp(-((dplx-0.351)**2/(2*0.358**2))) + 9.219), #NGC 2236
                   (86.597*np.exp(-((dplx-0.372)**2/(2*0.285**2))) + 14.261), #NGC 1907
                   (130.000*np.exp(-((dplx-0.370)**2/(2*0.281**2))) + 20.000), #NGC 1245
                   (157.641*np.exp(-((dplx-0.317)**2/(2*0.235**2))) + 19.697), #Melotte 71
                   (96.834*np.exp(-((dplx-0.337)**2/(2*0.300**2))) + 9.903), #Melotte 66
                   (110.000*np.exp(-((dplx-0.372)**2/(2*0.278**2))) + 17.263), #King 5
                   (45.188*np.exp(-((dplx-0.324)**2/(2*0.314**2))) + 6.308), #FSR 1252
                   (95.000*np.exp(-((dplx-0.421)**2/(2*0.267**2))) + 15.000), #Czernik 37
                   (570.000*np.exp(-((dplx-0.336)**2/(2*0.328**2))) + 100.000) ] #Collinder 110
        
        pplx.append(pplx_eqns[i]) #parallax distribution
        #pplxc.append(Anorm_plxC*phi_cpx[i])
        #pplxf.append(Anorm_plxF*phi_fpx[i])
        pplxc.append(phi_cpx[i])
        pplxf.append(phi_fpx[i])
        
        #calculate total probability from marginals:
        total_probs.append(1 - (1-ppm[j])*(1-pplx[j])*(1-ppos[j]))
        #phic = n_c*ppmc[j]*pposc[j]*pplxc[j]
        #phif = n_f*ppmf[j]*pposf[j]*pplxf[j]
        if selected_cluster_info['Cluster'][i] == 'NGC_6208'  or selected_cluster_info['Cluster'][i] == 'NGC_6940':
            phic = ppmc[j]
            phif = ppmf[j]
        elif  selected_cluster_info['Cluster'][i] == 'NGC_2447' or selected_cluster_info['Cluster'][i] == 'Collinder_110':
            phic = ppmc[j]*pplxc[j]
            phif = ppmf[j]*pplxf[j]
        elif selected_cluster_info['Cluster'][i] == 'King_5':
            phic = ppmc[j]*ppos[j]
            phif = ppmf[j]*ppos[j]
        else:
            phic = ppmc[j]*pposc[j]*pplxc[j]
            phif = ppmf[j]*pposf[j]*pplxf[j]
        total_probs2.append(phic / (phic + phif) )

    counts, edges, plot = hist(total_probs2, range=(min(total_probs2),max(total_probs2)), bins=40, color=orange,rwidth = 1.8)
    title('Membership Probability Distribution: ' + selected_cluster_info['Cluster'][i].replace('_',' '))
    ylim([0.0,300.0])
    yticks(fontsize=15)
    xticks(fontsize=15)
    xlabel('Probability', fontsize=15)
    ylabel('N', fontsize=15)
    #savefig(selected_cluster_info['Cluster'][i] + '_prob_dist.pdf', bbox_inches='tight', format='pdf')
    #close()
    show()

    # plot RA vs Dec, PMRa vs PMDec, CMDs with probabilities from total_probs
    import matplotlib.cm as cm
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6)) 
    
    #Plot CMD using Gaia photometry for all clusters
    bprp_filt = []
    g_filt = []
    total_probs2_filt = []
    for j in range(len(bprp)): #filter out stars of probability less than ___
        if selected_cluster_info['Cluster'][i] == 'Ruprecht_68' or selected_cluster_info['Cluster'][i] == 'King_5' or selected_cluster_info['Cluster'][i] == 'NGC_1907' or selected_cluster_info['Cluster'][i] == 'NGC_2509':
            if total_probs2[j] >= 0.3:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
        if selected_cluster_info['Cluster'][i] == 'NGC_7789' or selected_cluster_info['Cluster'][i] == 'Collinder_110' or selected_cluster_info['Cluster'][i] == 'NGC_2477' or selected_cluster_info['Cluster'][i] == 'NGC_6208' or selected_cluster_info['Cluster'][i] == 'NGC_6940' or selected_cluster_info['Cluster'][i] == 'Melotte_66':
            if total_probs2[j] >= 0.5:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
        if selected_cluster_info['Cluster'][i] == 'NGC_6791':
            if total_probs2[j] >= 0.9:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
        else:
            if total_probs2[j] >= 0.6:
                bprp_filt.append(bprp[j])
                g_filt.append(g[j])
                total_probs2_filt.append(total_probs2[j])
    title('Color Magnitude Diagram (Gaia): '+ selected_cluster_info['Cluster'][i].replace('_',' '))
    xlabel('BP-RP (mag)')
    ylabel('G (mag)')
    im = scatter(bprp_filt,g_filt, s=4, c=total_probs2_filt, cmap=cm.plasma_r)
    xlim([-2.0,6.0])
    ylim([22.0,6.0]) 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar(im, cax=cax, label='Prob')
    show()
    
    #panstarrs clusters
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
    if selected_cluster_info['Cluster'][i] not in ['NGC_2660', 'NGC_2477', 'NGC_6208', 'Ruprecht_68', 'Melotte_66']:
        panstarrs_gi_filt = []
        panstarrs_g_filt = []
        panstarrs_total_probs2_filt = []
        panstarrs_total_probs_filt = []
    for j in range(len(panstarrs_gi)): #filter out stars of probability less than ___
        if selected_cluster_info['Cluster'][i] == 'Ruprecht_68' or selected_cluster_info['Cluster'][i] == 'King_5' or selected_cluster_info['Cluster'][i] == 'NGC_1907' or selected_cluster_info['Cluster'][i] == 'NGC_2509':
            if total_probs2[j] >= 0.3:
                panstarrs_gi_filt.append(panstarrs_gi[j])
                panstarrs_g_filt.append(panstarrs_g[j])
                panstarrs_total_probs2_filt.append(total_probs2[j])
        if selected_cluster_info['Cluster'][i] == 'NGC_7789' or selected_cluster_info['Cluster'][i] == 'NGC_2477' or selected_cluster_info['Cluster'][i] == 'NGC_6208' or selected_cluster_info['Cluster'][i] == 'NGC_6940' or selected_cluster_info['Cluster'][i] == 'Melotte_66':
            if total_probs2[j] >= 0.7:
                panstarrs_gi_filt.append(panstarrs_gi[j])
                panstarrs_g_filt.append(panstarrs_g[j])
                panstarrs_total_probs2_filt.append(total_probs2[j])
        if selected_cluster_info['Cluster'][i] == 'NGC_6791':
            if total_probs2[j] >= 0.9:
                panstarrs_gi_filt.append(panstarrs_gi[j])
                panstarrs_g_filt.append(panstarrs_g[j])
                panstarrs_total_probs2_filt.append(total_probs2[j])
        else:
            if total_probs2[j] >= 0.3:
                panstarrs_gi_filt.append(panstarrs_gi[j])
                panstarrs_g_filt.append(panstarrs_g[j])
                panstarrs_total_probs2_filt.append(total_probs2[j])
    im = scatter(panstarrs_gi_filt,panstarrs_g_filt, s=4, c=panstarrs_total_probs2_filt, cmap=cm.plasma_r)
    xlim([-1.0,4.0])
    ylim([26.0,6.0])   
    title('Color Magnitude Diagram (Pan-STARRS): '+ selected_cluster_info['Cluster'][i].replace('_',' '))
    xlabel('G-I (mag)')
    ylabel('G (mag)')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar(im, cax=cax, label='Prob')
    show()
    
    N_80 = 0
    for j in total_probs2:
        if j >= 0.8:
            N_80 +=1
    N_80_list.append(N_80)
    
    if selected_cluster_info['Cluster'][i] in ['NGC_2660', 'NGC_2477', 'NGC_6208', 'Ruprecht_68', 'Melotte_66']:
        cluster_out2 = Table()
        cluster_out2['BP-RP'] = bprp
        cluster_out2['G'] = g
        cluster_out2['Prob'] = total_probs2
        ascii.write(cluster_out2, selected_cluster_info['Cluster'][i]+'_out2.dat', overwrite=True)
    
    else:
        cluster_out2 = Table()
        cluster_out2['G-I'] = panstarrs_gi
        cluster_out2['G'] = panstarrs_g
        cluster_out2['Prob'] = total_probs2
        ascii.write(cluster_out2, selected_cluster_info['Cluster'][i]+'_out_panstarrs.dat', overwrite=True)

#lists of opening angles and number of stars with proability > 80 for each cluster
print([str(N_80_list[i]) + ' ' + selected_cluster_info['Cluster'][i] for i in [14]])#in range(len(N_80_list))])
print([str(opening_angles[i]) + ' ' + selected_cluster_info['Cluster'][i] for i in [14]]) #range(len(N_80_list))])

print('normal stop')"""