#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 12:38:35 2024

@author: abbychriss
"""

import astropy.units as u
from astropy.io import fits
from astropy.table import Table, vstack
from pylab import *
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import biweight_location
from astropy.stats import biweight_midvariance
import math

num = [26,225,245,368,571]

rcParams['svg.fonttype'] = 'none'
rcParams['text.usetex'] = True

#PARSEC [M/H]=-0.5
isochrone = Table.read('Isochrones/hess_WL11_ev8_955.out',format='ascii',fast_reader=False)
isochrone.rename_column('col8', 'B')
isochrone.rename_column('col9', 'V')
isochrone.rename_column('col11', 'I')
print(isochrone['B','V','I'])
bi = isochrone['B'] - isochrone['I']
imag = isochrone['I']
bi = np.array(bi)
imag = np.array(imag)

iso_shift = [[0.08,0], #smc26
             [0.03,0], #smc225
             [-0.05,0.4], #smc245
             [0.04,0.05], #smc368
             [-0.13,0.05]] #smc571

rc_bi_lim = [[1.40,1.76], #smc26
             [1.467,1.626], #smc225
             [1.45,1.771], #smc245
             [1.40,1.7273], #smc368
             [1.23,1.65]] #smc571

rc_i_lim = [[-0.07,-1.0], #smc26
             [-0.4,-0.78], #smc225
             [-0.233,-0.808], #smc245
             [-0.106,-1.05], #smc368
             [-0.33,-1.1]] #smc571

rc_rgb_bi_lim = [[1.40,3.69], #smc26
             [1.467,3.38], #smc225
             [1.45,3.03], #smc245
             [1.40,2.467], #smc368
             [1.165,2.94]] #smc571

rc_rgb_i_lim = [[-0.07,-4.5], #smc26
             [-0.4,-4.91], #smc225
             [-0.233,-4.0], #smc245
             [-0.106,-3.47], #smc368
             [-0.33,-4.52]] #smc571

#distance in B-I from median RGB location to furthest star on the red side, multiplied by 2, gives an estimate for the total width of the RGB.
#Multiply the half-width by 0.341 to give one sigma
reddest_rgb_bi = [1.927, #smc26 
             2.108, #smc225
             1.966, #smc245
             1.936, #smc368
             1.881] #smc571
             
med_rc_bi_list = []
med_rc_i_list = []
med_rgb_bi_list = []
d_bi_list = []
total_errors = []
rc_numbers = []
rgb_numbers = []
for i in range(len(num)):
    smc_cluster = Table.read('SMC'+str(num[i])+'field_data.csv',format='csv',fast_reader=False)
    B_I = smc_cluster['B-I']
    I = smc_cluster['I']
    smc_cluster = np.sort(smc_cluster, order='I')
    
    #FIND MEDIAN LOCATION OF RED CLUMP
    rc_bi = []
    rc_i = []
    nrgb_rc=0
    #collect all the stars in the red clump by taking the stars within the bp-rp and g limits above
    for j in range(len(B_I)):
        if rc_bi_lim[i][0] <= B_I[j] <= rc_bi_lim[i][1] and rc_i_lim[i][1] <= I[j] <= rc_i_lim[i][0]:
            rc_bi.append(B_I[j])
            rc_i.append(I[j])
        
        if rc_rgb_bi_lim[i][0] <= B_I[j] <= rc_rgb_bi_lim[i][1] and rc_rgb_i_lim[i][1] <= I[j] <= rc_rgb_i_lim[i][0]:
            nrgb_rc+=1
    rc_numbers.append(len(rc_bi))
    n_rgb = nrgb_rc - len(rc_bi)
    rgb_numbers.append(n_rgb)
    
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 10))
    
    #use astropy statistics biweight_location to find the central location of the sample
    #if the median absolute deviation is zero, the biweight location is simply the median
    med_rc_bi = round(biweight_location(rc_bi),6)
    med_rc_i = round(biweight_location(rc_i),6)
    med_rc_bi_list.append(med_rc_bi)
    med_rc_i_list.append(med_rc_i)
    
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx],idx
    
    iso_rgb_bi = np.array(bi[115:150])+iso_shift[i][0] #bp-rp of RGB of isochrone, shifted to fit cluster
    iso_rgb_i = np.array(imag[115:150])+iso_shift[i][1] #gmag of RGB of isochrone, shifted to fit cluster
    
    rgb_g,idx = find_nearest(iso_rgb_i, med_rc_i)
    med_rgb_bi = iso_rgb_bi[idx]
    med_rgb_bi_list.append(round(med_rgb_bi,4))
    
    d_bi = med_rgb_bi - med_rc_bi
    d_bi_list.append(round(d_bi,6))
    
    iso_slope = (iso_rgb_i[idx]-iso_rgb_i[idx-1])/(iso_rgb_bi[idx]-iso_rgb_bi[idx-1])
    #plt.plot(np.linspace(1.4,2.2,80), iso_slope*(np.linspace(1.4,2.2,80)- med_rgb_bi) + med_rc_i,color='orange')
    
    plt.plot(np.linspace(rc_bi_lim[i][0], rc_bi_lim[i][1],80), np.full(80, rc_i_lim[i][0]), color='r',lw=1.5)
    plt.plot(np.linspace(rc_bi_lim[i][0], rc_bi_lim[i][1],80), np.full(80, rc_i_lim[i][1]), color='r',lw=1.5)
    plt.plot(np.full(80, rc_bi_lim[i][0]), np.linspace(rc_i_lim[i][0],rc_i_lim[i][1],80), color='r',lw=1.5)
    plt.plot(np.full(80, rc_bi_lim[i][1]), np.linspace(rc_i_lim[i][1], rc_i_lim[i][0], 80), color='r',lw=1.5)
    
    im = plt.scatter(B_I,I, s=4, color='#020eed',alpha=0.75)
    plt.plot(bi[115:143] + iso_shift[i][0], imag[115:143] + iso_shift[i][1],color='#00de12',lw=1.5,alpha=1,label='0.955 Gyr, [M/H]=-0.5')
    plt.plot(bi[166:234] + iso_shift[i][0], imag[166:234] + iso_shift[i][1],color='#00de12',lw=1.5,alpha=1)
    #plt.legend(frameon=False,fontsize=20,loc='upper right')
    plt.xlim([-1,4])
    plt.ylim([3.5,-5])   
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    #plt.title('Color Magnitude Diagram: Field near SMC'+str(num[i]))
    plt.xlabel(r'$B-I$ (mag)',fontsize=25)
    plt.ylabel(r'$I$ (mag)',fontsize=25)
    
    plt.scatter(med_rc_bi,med_rc_i,s=100,color='white',edgecolor='black',marker='*')
    plt.scatter(med_rgb_bi,med_rc_i,s=100,color='white',edgecolor='black',marker='*')
    #plt.savefig('/Users/abbychriss/Desktop/SMC'+str(num[i])+'_cmd.pdf', bbox_inches='tight', format='pdf')
    #plt.close()
    plt.show()
    
    #use half the width of the RGB*0.341 for 1 standard deviation in RGB color at level of RC
    #photometric error: sigma(B-I) = 0.11 mag, then divide by number of RC and RGB stars
    #total error of averaged values is average of individual errors, divided by sqrt(5)
    d_bi_error = math.sqrt( ((2*((reddest_rgb_bi[i]-med_rgb_bi)*0.2))**2/n_rgb) #uncertainty in RGB color
                           + (biweight_midvariance(rc_bi)/len(rc_bi)) #uncertainty in RC color -> variance = sd^2
                           + (0.11**2/nrgb_rc)  #photometric error in B-I -> sigma(B-I) = sqrt(sigma(B)^2 + sigma(I)^2)
                           + ((0.08**2 / iso_slope)/nrgb_rc )) #photometric error in I divided by the slope of the RGB at the level of the red clump
    
    total_errors.append(round(d_bi_error,6))
    
print('Median B-I red clump: ',med_rc_bi_list)
print('Median I red clump: ',med_rc_i_list)
print('Median B-I red giant branch: ',med_rgb_bi_list)
print('d(B-I): ',d_bi_list)
print('Age of SMC fields: 8.95 +/- 0.07 log years')
print('nRC: ', [('SMC'+str(num[i]), rc_numbers[i]) for i in range(len(rc_numbers))])
print('nRGB: ', [('SMC'+str(num[i]), rgb_numbers[i]) for i in range(len(rgb_numbers))])

average_d_bi = round(np.average(d_bi_list),4)
print(average_d_bi)
average_error_d_bi = round(np.average(total_errors/np.sqrt(4)),4)
print('Average d(B-I): ',average_d_bi,'+/-',average_error_d_bi)
print('Average nRC: ',np.average(rc_numbers))

d_br_list = [0.824*d_bi for d_bi in d_bi_list]
d_br_error_list = [0.824*d_bi_err for d_bi_err in total_errors]
d_br_error_list_new = [np.sqrt(d_br_err**2 + 0.007**2) for d_br_err in d_br_error_list]
d_br_smc = round(np.average(d_br_list),4)
d_br_smc_error = round(np.average(d_br_error_list_new)/np.sqrt(4),4)

d_bprp_list = [d_br/1.77 for d_br in d_br_list]
d_bprp_error_list = [d_br_err/1.77 for d_br_err in d_br_error_list_new]
d_bprp_error_list_new = [np.sqrt(d_bprp_err**2 + 0.004**2) for d_bprp_err in d_bprp_error_list]

print('d(BP-RP): ',[(str(round(d_bprp_list[i],3))+' \pm '+str(round(d_bprp_error_list_new[i],3))) for i in range(len(d_br_list))])
print('d(B-R): ',[(str(round(d_br_list[i],3))+' \pm '+str(round(d_br_error_list_new[i],3))) for i in range(len(d_br_list))])

d_bprp_smc = round(np.average(d_bprp_list),4)
d_bprp_smc_error = round(np.average(d_bprp_error_list_new),4)/np.sqrt(4) #divide by sqrt(4) because there are 5 clusters

print('Average d(B-R):',d_br_smc,'+/-',d_br_smc_error)
print('Average d(BP-RP): ',d_bprp_smc, '+/-',d_bprp_smc_error)
print([round(d_bprp,4) for d_bprp in d_bprp_list])
print([round(err,4) for err in d_bprp_error_list_new])
