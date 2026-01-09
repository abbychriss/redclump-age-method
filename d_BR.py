#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:39:45 2023

@author: abbychriss
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *
from mpl_point_clicker import clicker
from astropy.table import Table
from astropy.io import ascii
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from decimal import Decimal
import sys
from IPython.display import set_matplotlib_formats

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['text.usetex'] = True

#list of B-R values from table
B_R = np.array([-0.014,0.127,0.423,0.836,1.115,1.442,1.819,2.241,2.680],dtype=float)
#list of corresponding BP-RP values
BP_RP = np.array([-0.02543,0.01237,0.3774,0.7023,0.8615,1.0404,1.253,1.6228,2.4106],dtype=float)

#median BP-RP value of the red giant branch at luminosity of red clump for each cluster, sorted by reverse alphabetical order
med_rgb_bprp = [1.297, 1.262, 1.504, 1.24, 1.274, 1.161, 1.092, 1.184, 1.235, 1.226, 1.148, 1.198, 1.254, 1.161]

#median BP-RP value of red clump for each cluster
med_rc_bprp = [1.146086, 1.148982, 1.34344, 1.157591, 1.170835, 1.095183, 0.994505, 1.100035, 1.059469, 1.146977, 1.041511, 1.080467, 1.131439, 1.056953]

med_rgb_br = np.interp(med_rgb_bprp, BP_RP, B_R)
med_rc_br = np.interp(med_rc_bprp, BP_RP, B_R)

d_br = [float(med_rgb_br[i] - med_rc_br[i]) for i in range(len(med_rgb_br))]
d_bprp = [float(med_rgb_bprp[i] - med_rc_bprp[i]) for i in range(len(med_rgb_bprp))]

def func(x, a):#, b):
    return a*x# + b
popt, pcov = curve_fit(func, d_bprp, d_br)
plot(np.linspace(0,0.2,50), func(np.linspace(0,0.2,50), *popt), '#854384', label=r'$y = %5.3fx$'% tuple(popt))#' + %5.3f$' % tuple(popt))
delta_color_slope = tuple(popt)[0]
scatter(d_bprp,d_br)
xlabel(r'$d_{BP-RP}$ ($mag$)')
ylabel(r'$d_{B-R}$ ($mag$)')
xlim([0.04,0.19])
title(r'$d_{B-R}$ vs. $d_{BP-RP}$ around 5000 K')
legend()
show()

cluster_keep = ['King_5','Melotte_66','Melotte_71','NGC_1245','NGC_2420','NGC_2477',
                'NGC_2506','NGC_2509','NGC_2627','NGC_2682','NGC_6208','NGC_6791',
                'NGC_7789','Ruprecht_68']

#list of metallicities from Dias et al 2021MNRAS.504..356D
metallicity_table = Table.read('metallicity_data.txt',format='ascii')
feh_dias = []
err_feh_dias = []
for cluster in np.flip(cluster_keep):
    if cluster in metallicity_table['Cluster']:
        feh_dias.append(metallicity_table['[Fe/H]'][np.where(cluster == metallicity_table['Cluster'])[0][0]])
        err_feh_dias.append(metallicity_table['e_[Fe/H]'][np.where(cluster == metallicity_table['Cluster'])[0][0]])

#from list of clusters we pull the parallax, proper motion, and position from cluster_info.dat table
cluster_info = Table.read('RC_cluster_info.dat',format='ascii',fast_reader=False)

selected_cluster_info = cluster_info[np.where(np.isin(np.array(cluster_info['Cluster']),cluster_keep))]
selected_cluster_info.rename_column('Input file name', 'Gaia file')

#sorted by log age
final_data_d_bprp = Table.read('final_data_d_bprp.dat',format='ascii')

#sorted by alphabetical order
final_data_d_bprp = np.sort(final_data_d_bprp, order='Cluster')

#sorted by reverse alphabetical order
e_bprp = np.flip(final_data_d_bprp['E(BP-RP)'])
d_bprp = np.flip(final_data_d_bprp['d(BP-RP)'])
d_bprp_err = np.flip(final_data_d_bprp['error d(BP-RP)'])
d_br_err = [x*1.77 for x in d_bprp_err]

clusters = np.flip(final_data_d_bprp['Cluster'])
cluster_ages = np.flip(final_data_d_bprp['log age'])
error_cluster_ages = np.flip(final_data_d_bprp['error log age'])
n_rc80 = np.flip(final_data_d_bprp['nRC80'])
fe_h = np.array(np.flip(final_data_d_bprp['[Fe/H]']),dtype='float')
fe_h = [float("%.8f" %x) for x in fe_h]
fe_h[0]=0.13 #Replace value for metallicity for Ruprecht 68 by value reported in (Dias et. al, MNRAS 504, 356–371 (2021))
fe_h[6]=0.082 #Replace value for metallicity for NGC 2509 by value reported in (Dias et. al, MNRAS 504, 356–371 (2021))

print(clusters)
print(fe_h)

smc_feh = -0.6
smc_avg_n_rc = 151.6
smc_log_age = 8.85
smc_log_age_err = 0.07
smc_d_br = 0.1964
smc_d_br_err = 0.007

smc_d_bprp = 0.111 #1.77 = color filter transformation factor
smc_d_bprp_err = 0.00405

#this is sorted by cluster name in reverse alphabetical order
final_data = np.zeros(15, dtype={'names':('cluster', 'd_br', 'error(d_br)','d_bprp', 'error(d_bprp)', 'e_bprp', 'log age', 'error(log age)', 'nRC80', '[Fe/H]','opening angle (deg)'),
                                 'formats':('U20', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8','f8')})
final_data['cluster'] = np.concatenate((clusters,['SMC fields']))
final_data['d_br'] = np.concatenate((d_br,[smc_d_br]))
final_data['error(d_br)'] = np.concatenate((d_br_err,[smc_d_br_err]))
final_data['d_bprp'] = np.concatenate((d_bprp,[smc_d_bprp]))
final_data['error(d_bprp)'] = np.concatenate((d_bprp_err,[smc_d_bprp_err]))
final_data['e_bprp'] = np.concatenate((e_bprp,[99.999])) #99.999 means I don't know
final_data['log age'] = np.concatenate((cluster_ages,[smc_log_age]))
final_data['error(log age)'] = np.concatenate((error_cluster_ages,[smc_log_age_err]))
final_data['nRC80'] = np.concatenate((n_rc80,[smc_avg_n_rc]))
final_data['[Fe/H]'] = np.concatenate((fe_h,[smc_feh]))
final_data['opening angle (deg)'] = np.array([0.22114,0.42180,0.20374,0.48557,0.62159,
                                              0.26039,0.17866,0.26327,0.44970,0.20989,
                                              0.32621,0.22430,0.26235,0.27242,99.999])
age_final_data = np.sort(final_data, order='log age')
d_bprp_data = np.sort(final_data,order='d_bprp')

age_log_err = final_data['error(log age)'] #sorted by reverse alphabetical order
age_log_err_plot = age_final_data['error(log age)'] #sorted by age

#create figure box
size=100
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14,10))
#fit a line to our age calibration ==> not doing because there is no trend
"""def func(x, a, b):
    return (10 ** (a * x)) - b
popt, pcov = curve_fit(func, final_data['log age'],final_data['d_br'])
print(np.linalg.cond(pcov))
plot(range(8,12), func(range(8,12), *popt), '#5B279A', label=r'$d_{B-R} = 10^{%5.3f t} - %5.3f$' % tuple(popt))"""
#scatter(age_final_data['log age'],age_final_data['d_br'],color='#c5001d',s=10,label='Galactic open clusters (This work)') #plot our data in log age

#PLOT SMC POINT
"""scatter(smc_log_age,smc_d_br,label='SMC fields (This work)',s=size,color='white',edgecolor='black',marker='>')
errorbar(smc_log_age,smc_d_br,xerr=smc_log_age_err,yerr=smc_d_br_err,lw=1,capsize=3,fmt='none',color='black') #SMC point"""
#errorbar(age_final_data['log age'],age_final_data['d_br'],xerr=age_log_err_plot,yerr=d_br_err,fmt="o",capsize=4,capthick=1,color='#7D46B7',alpha=0.35) #put in error bars

#four separate structured arrays for Hatzidimitriou data: open, globular, SMC, LMC clusters and their errors
hatz_open_data=np.zeros(6, dtype={'names':('open cluster','log age','err(log age)','d_br','err(d_br)','[Fe/H]'),
                             'formats':('U20','f8','f8','f8','f8','f8')})
hatz_open_data['open cluster'] = ['NGC 2420', 'NGC 2506', 'NGC 2682', 'NGC 6791', 'NGC 7789', 'Mel 66']
hatz_open_data['log age'] = [9.24,9.22,9.63, 9.8, 9.19, 9.63]
hatz_open_data['err(log age)'] = [0.175,0.175,0.15,0.15,0.175,0.15]
hatz_open_data['d_br'] = [0.12,0.13,0.14,0.20,0.13,0.15]
hatz_open_data['err(d_br)'] = [0.02,0.02,0.02,0.03,0.02,0.02]
hatz_open_data['[Fe/H]'] = [-0.6,-0.55,-0.1,0.04,-0.1,-0.36]
hatz_open_data = np.sort(hatz_open_data, order='log age')

hatz_glob_data=np.zeros(4, dtype={'names':('globular cluster', 'd_br', 'err(d_br)', 'log age', 'err(log age)','[Fe/H]'),
                                  'formats':('U20','f8','f8','f8','f8','f8')})
hatz_glob_data['globular cluster'] = ['Pal 4', 'Pal 12', '47 Tuc', 'Eridanus']
hatz_glob_data['d_br'] = [0.34, 0.23, 0.29, 0.30]
hatz_glob_data['err(d_br)'] = [0.04, 0.03, 0.03, 0.02]
hatz_glob_data['log age'] = [np.log10((10.5)*10**9), 10.050, 10.100, np.log10((10.5)*10**9)]
hatz_glob_data['[Fe/H]'] = [-1.41,-0.850,-0.645,-1.4]
hatz_glob_data['err(log age)'] = [0, 0., 0., 0]

hatz_glob_data = np.sort(hatz_glob_data, order='log age')

hatz_smc_data=np.zeros(4, dtype={'names':('SMC cluster', 'd_br', 'err(d_br)', 'log age', 'err(log age)','[Fe/H]'),
                                  'formats':('U20','f8','f8','f8','f8','f8')})
hatz_smc_data['SMC cluster'] = ['L1', 'L113', 'NGC 411', 'K3']
hatz_smc_data['d_br'] = [0.25, 0.20, 0.10, 0.16]
hatz_smc_data['err(d_br)'] = [0.02, 0.02, 0.02, 0.02]
hatz_smc_data['log age'] = [np.log10(9*10**9), np.log10(5.3*10**9), np.log10(1.8*10**9), np.log10(6*10**9)] #new age estimates from https://ui.adsabs.harvard.edu/abs/1998AJ....116.2395M/abstract & https://ui.adsabs.harvard.edu/abs/1998AJ....115.1934D/abstract
hatz_smc_data['err(log age)'] = [0,0,0,0]
#hatz_smc_data['log age'] = [np.log10((0.9*11)*10**9), np.log10((0.9*6)*10**9), np.log10((0.9*1.8)*10**9), np.log10((0.9*8)*10**9)]
#hatz_smc_data['err(log age)'] = [(0.9*1)*10**9, (0.9*1)*10**9, (0.9*0.3)*10**9, (0.9*1)*10**9]
hatz_smc_data['[Fe/H]'] = [-1.3,-1.4,-0.9,-1.3]
hatz_smc_data = np.sort(hatz_smc_data, order='log age')

hatz_lmc_data=np.zeros(5, dtype={'names':('LMC cluster', 'd_br', 'err(d_br)', 'log age', 'err(log age)','[Fe/H]'),
                                  'formats':('U20','f8','f8','f8','f8','f8')})
hatz_lmc_data['LMC cluster'] = ['NGC 1978', 'NGC 2173', 'H4', 'LW47', 'ESO121-SC03']
hatz_lmc_data['d_br'] = [0.12, 0.13, 0.11, 0.09, 0.25]
hatz_lmc_data['err(d_br)'] = [0.01, 0.02, 0.02, 0.01, 0.02]
hatz_lmc_data['log age'] = [np.log10((2)*10**9), np.log10((1.6)*10**9), 9.37, np.log10((2)*10**9), 9.99]
hatz_lmc_data['err(log age)'] = [0, (0.2*10**9)/(np.log(10)*1.6*10**9), 0.02, 0., 0.01] #use error propagation delta(log10(t)) = delta(t)/(t*log10), 
hatz_lmc_data['[Fe/H]'] = [-0.7,-0.4,-0.88,-0.3,-1.4]
hatz_lmc_data = np.sort(hatz_lmc_data, order='log age')

up_err_t_glob = np.log10(np.add(np.power(10,hatz_glob_data['log age']), hatz_glob_data['err(log age)'])) - hatz_glob_data['log age']
low_err_t_glob = hatz_glob_data['log age'] - np.log10(np.subtract(np.power(10,hatz_glob_data['log age']), hatz_glob_data['err(log age)']))
terr_glob = [low_err_t_glob, up_err_t_glob]

up_err_t_smc = np.log10(np.add(np.power(10,hatz_smc_data['log age']), hatz_smc_data['err(log age)'])) - hatz_smc_data['log age']
low_err_t_smc = hatz_smc_data['log age'] - np.log10(np.subtract(np.power(10,hatz_smc_data['log age']), hatz_smc_data['err(log age)']))
terr_smc = [low_err_t_smc, up_err_t_smc]

up_err_t_lmc = np.log10(np.add(np.power(10,hatz_lmc_data['log age']), hatz_lmc_data['err(log age)'])) - hatz_lmc_data['log age']
low_err_t_lmc = hatz_lmc_data['log age'] - np.log10(np.subtract(np.power(10,hatz_lmc_data['log age']), hatz_lmc_data['err(log age)']))
terr_lmc = [low_err_t_lmc, up_err_t_lmc]

hatz_data = np.zeros(19, dtype={'names':('cluster','d_br','err(d_br)','log age','upper err(log age)', 'lower err(log age)', '[Fe/H]'),
                                  'formats':('U20','f8','f8','f8','f8','f8','f8')})
hatz_data['cluster'] = np.concatenate((hatz_open_data['open cluster'], hatz_glob_data['globular cluster'], hatz_smc_data['SMC cluster'], hatz_lmc_data['LMC cluster']))
hatz_data['d_br'] = np.concatenate((hatz_open_data['d_br'], hatz_glob_data['d_br'], hatz_smc_data['d_br'], hatz_lmc_data['d_br']))
hatz_data['err(d_br)'] = np.concatenate((hatz_open_data['err(d_br)'], hatz_glob_data['err(d_br)'], hatz_smc_data['err(d_br)'], hatz_lmc_data['err(d_br)']))
hatz_data['log age'] = np.concatenate((hatz_open_data['log age'], hatz_glob_data['log age'], hatz_smc_data['log age'], hatz_lmc_data['log age']))
hatz_data['upper err(log age)'] = np.concatenate((hatz_open_data['err(log age)'],up_err_t_glob,up_err_t_smc,up_err_t_lmc))
hatz_data['lower err(log age)'] = np.concatenate((hatz_open_data['err(log age)'],low_err_t_glob,low_err_t_smc,low_err_t_lmc))
hatz_data['[Fe/H]'] = np.concatenate((hatz_open_data['[Fe/H]'], hatz_glob_data['[Fe/H]'], hatz_smc_data['[Fe/H]'], hatz_lmc_data['[Fe/H]']))
np.sort(hatz_data, order='log age')

all_data = np.zeros(34, dtype={'names':('cluster','d_br', 'error(d_br)','log age','upper error(log age)','lower error(log age)','[Fe/H]'),
                                  'formats':('U20','f8','f8','f8','f8','f8','f8')})
all_data['cluster'] = np.concatenate((age_final_data['cluster'], hatz_data['cluster']))
all_data['d_br'] = np.concatenate((age_final_data['d_br'], hatz_data['d_br']))
all_data['error(d_br)'] = np.concatenate((age_final_data['error(d_br)'], hatz_data['err(d_br)']))
all_data['log age'] = np.concatenate((age_final_data['log age'],hatz_data['log age']))
all_data['upper error(log age)'] = np.concatenate((age_log_err_plot,hatz_data['upper err(log age)']))
all_data['lower error(log age)'] = np.concatenate((age_log_err_plot,hatz_data['lower err(log age)']))
all_data['[Fe/H]'] = np.concatenate((age_final_data['[Fe/H]'],hatz_data['[Fe/H]']))
all_data = np.sort(all_data, order='log age')
all_data = ascii.write(all_data,'final_data_d_br.dat',overwrite=True)
all_data = Table.read('final_data_d_br.dat',format='ascii')

#plot all data with color bar to get color bar to show up
im = scatter(all_data['log age'],all_data['d_br'],c=all_data['[Fe/H]'],cmap=cm.coolwarm,s=0)

#convert [Fe/H] to a color tuple using the colormap used for scatter
norm = matplotlib.colors.Normalize(vmin=min(all_data['[Fe/H]']), vmax=max(all_data['[Fe/H]']), clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.coolwarm)
feh_color = np.array([(mapper.to_rgba(v)) for v in all_data['[Fe/H]']])

#loop over each data point to plot
j=0
i_thiswork=0
i_hatzopen=0
i_hatzglob=0  
i_hatzsmc=0
i_hatzlmc=0
for x, y, xe1, xe2, ye, color in zip(all_data['log age'], all_data['d_br'], all_data['lower error(log age)'], all_data['upper error(log age)'], all_data['error(d_br)'], feh_color):
    if all_data['cluster'][j] in cluster_keep:
        if i_thiswork==0:
            """if age_final_data['[Fe/H]'][j]>70:
                plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color='black')
                plt.scatter(x, y, marker='o', s=size, color='white', edgecolor='black')
            else:"""
            plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
            plt.scatter(x, y, marker='o', s=size-20, color='white', edgecolor='black', label='Galactic open clusters (This work)')
            plt.scatter(x, y, marker='o', s=size, color=color,edgecolor=color)
        else:
            """if age_final_data['[Fe/H]'][j]>70:
                plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color='black')
                plt.scatter(x, y, marker='o', s=size, color='white', edgecolor='black')
            else:"""
            plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
            plt.scatter(x, y, marker='o', s=size, color=color,edgecolor=color)
        i_thiswork+=1
    
    if all_data['cluster'][j]=='SMC fields':
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
        plt.scatter(x, y, marker='>', s=size-20, color='white', edgecolor='black', label='SMC fields (This work)')
        plt.scatter(x, y, marker='>', s=size, color=color,edgecolor=color)
    
    if all_data['cluster'][j] in hatz_open_data['open cluster']:
        if i_hatzopen==0:
            plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
            plt.scatter(x, y, marker='s', s=size-20, color='white', edgecolor='black',label='Galactic open clusters (H91)')
            plt.scatter(x, y, marker='s', s=size, color=color,edgecolor=color)
        else:
            plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
            plt.scatter(x, y, marker='s',s=size, color=color,edgecolor=color)
        i_hatzopen+=1
      
    if all_data['cluster'][j] in hatz_glob_data['globular cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
        if i_hatzglob==0:
            plt.scatter(x, y, marker='^',s=size-20,color='white', edgecolor='black', label='Galactic globular clusters (H91)')
            plt.scatter(x, y, marker='^',s=size,color=color,edgecolor=color)
        else:
            plt.scatter(x, y, marker='^',s=size, color=color,edgecolor=color)
        i_hatzglob+=1
    
    if all_data['cluster'][j] in hatz_smc_data['SMC cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color)
        if i_hatzsmc==0:
            plt.scatter(x, y, marker='d', s=size-20, color='white', edgecolor='black', label='SMC clusters (H91)')
            plt.scatter(x, y, marker='d', s=size, color=color,edgecolor=color)
        else:
            plt.scatter(x, y, marker='d', s=size, color=color,edgecolor=color)
        i_hatzsmc+=1
        
    if all_data['cluster'][j] in hatz_lmc_data['LMC cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1,fmt='none', capsize=3, color=color)
        if i_hatzlmc==0:
            plt.scatter(x, y, marker='D', s=size-20, color='white', edgecolor='black', label='LMC clusters (H91)')
            plt.scatter(x, y, marker='D', s=size, color=color,edgecolor=color)
        else:
            plt.scatter(x, y, marker='D', s=size, color=color,edgecolor=color)
        i_hatzlmc+=1
    j+=1
    
#fit a curve (exponential with base 10) to Hatzidimitriou's data using scipy.optimize
def func(x, a, b):
    return a*(10**(x)) + b
popt, pcov = curve_fit(func, hatz_data['log age'], hatz_data['d_br'])
plot(np.linspace(8,12,80), func(np.linspace(8,12,80), *popt), color='#2A527C', linestyle='--', alpha = 0.5, label=r'$d_{B-R} = (1.61 \times 10^{-11})10^t + 0.084$')
print(tuple(popt))
print(np.linalg.cond(pcov))

ylim([-0.05,0.65])
xlim([8.1,10.3])
#title('Age Calibration Fit')
xlabel('Age (log years)', fontsize=20)
ylabel(r'$d_{B-R}$ (mag)', fontsize=20)
xticks(fontsize = 20)
yticks(fontsize = 20)
legend(frameon=False,loc=2,fontsize=15)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cax.tick_params(labelsize=20)
cb = colorbar(im,cax=cax)
cb.set_label(label='[Fe/H]', size=20)
savefig('/Users/abbychriss/Desktop/age_calibration.pdf', bbox_inches='tight', format='pdf')
close()
#show()

#------------------------------------------------------------------
#NEW PLOT FOR ALL DATA PLUS THEORETICAL PREDICTIONS
#create figure box
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(14,10))
im = scatter(all_data['log age'],all_data['d_br'],c=all_data['[Fe/H]'],cmap=cm.coolwarm,s=0,alpha=0.35)

#data from isochrones

#marigo 2008
iso_d_br_m7 = [0.24953012, 0.06139919, 0.04410251, 0.0368091, 0.03047921, 0.03183708,
0.0727148, 0.11971314, 0.12480098, 0.09154604, 0.100398, 0.1151616, 0.15479008, 0.25979777]
iso_log_age_m7 = [8.2999429, 8.44994099, 8.59999218, 8.74996808, 8.89998456, 9.04999286,
9.20000187, 9.3499959, 9.50000307, 9.64999651, 9.80000183, 9.94999954,
10.10000158, 10.25000014]

#marigo 2008 (metal poor)
iso_d_br_m8 = [1.37264348, 1.03846112, 0.74298929, 0.4833454 , 0.35738582,
       0.25683363, 0.1937709 , 0.1639044 , 0.15555779, 0.16092502,
       0.17317367, 0.18819244, 0.19196776, 0.19787395, 0.20560033,
       0.22106881, 0.23992555, 0.27705492, 0.348502  , 0.59852846]
iso_log_age_m8 = [ 8.2999429 ,  8.40001964,  8.49996187,  8.59999218,  8.70001106,
        8.80002936,  8.89998456,  9.        ,  9.09999123,  9.20000187,
        9.3000082 ,  9.40000235,  9.50000307,  9.60000309,  9.7000024 ,
        9.80000183,  9.90000097, 10.        , 10.10000158, 10.20000187]

#basti 5
iso_d_br_basti = [0.71450013, 0.14648398, 0.08050823, 0.04782743, 0.13096487,
       0.12808788, 0.10137438, 0.09844984, 0.10341408, 0.10757955,
       0.11214797, 0.12035078, 0.1416016 , 0.18402415]
iso_log_age_basti = [ 7.95424251,  8.30103   ,  8.54406804,  8.69897   ,  8.90308999,
        9.09691001,  9.30103   ,  9.43933269,  9.54406804,  9.65321251,
        9.77815125,  9.87506126,  9.95424251, 10.0211893 ]

#basti 6
iso_d_br_basti6 = [0.41650344, 0.22876971, 0.12991335, 0.06450466, 0.13268641,
       0.13194257, 0.13090103, 0.14111277, 0.12367777, 0.11110069,
       0.1140806 , 0.12302443, 0.13976055, 0.16466348, 0.19873546]
iso_log_age_basti6 = [ 8.20411998,  8.34242268,  8.53147892,  8.66275783,  8.77815125,
        8.87506126,  8.95424251,  9.07918125,  9.25527251,  9.47712125,
        9.65321251,  9.77815125,  9.95424251, 10.07918125, 10.17609126]

#parsec
iso_d_br_parsec = [0.18694578, 0.07047662, 0.03641852, 0.10052414, 0.12967167,
       0.1217676 , 0.13644915, 0.15182921, 0.18251783, 0.20998258]
iso_log_age_parsec = [ 8.43600354,  8.63998425,  8.84397984,  9.04801429,  9.25200302,
        9.45600144,  9.66000172,  9.86400054, 10.07918125, 10.17609126]

#parsec (metal poor)
iso_d_br_parsecmp = [1.88293109, 1.62182664, 0.85892289, 0.39215356, 0.18666677,
       0.13837978, 0.19658041, 0.21376253, 0.2390092 , 0.26536736,
       0.33519734, 0.39500612]
iso_log_age_parsecmp = [ 8.02816442,  8.23197903,  8.43600354,  8.63998425,  8.84397984,
        9.04801429,  9.25200302,  9.45600144,  9.66000172,  9.86400054,
       10.07918125, 10.17609126]


#loop over each data point to plot
j=0
for x, y, xe1, xe2, ye, color in zip(all_data['log age'], all_data['d_br'], all_data['lower error(log age)'], all_data['upper error(log age)'], all_data['error(d_br)'], feh_color):
    if all_data['cluster'][j] in cluster_keep:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color,alpha=0.35)
        plt.scatter(x, y, marker='o', s=size, color=color,edgecolor=color,alpha=0.35)
    
    if all_data['cluster'][j]=='SMC fields':
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color,alpha=0.35)
        plt.scatter(x, y, marker='>', s=size, color=color,edgecolor=color,alpha=0.35)
    
    if all_data['cluster'][j] in hatz_open_data['open cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color,alpha=0.35)
        plt.scatter(x, y, marker='s',s=size, color=color,edgecolor=color,alpha=0.35)
      
    if all_data['cluster'][j] in hatz_glob_data['globular cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color,alpha=0.35)
        plt.scatter(x, y, marker='^',s=size, color=color,edgecolor=color,alpha=0.35)
    
    if all_data['cluster'][j] in hatz_smc_data['SMC cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1, fmt='none', capsize=3, color=color,alpha=0.35)
        plt.scatter(x, y, marker='d', s=size, color=color,edgecolor=color,alpha=0.35)
        
    if all_data['cluster'][j] in hatz_lmc_data['LMC cluster']:
        plt.errorbar(x, y, xerr=[[xe1],[xe2]], yerr=ye, lw=1,fmt='none', capsize=3, color=color,alpha=0.35)
        plt.scatter(x, y, marker='D', s=size, color=color,edgecolor=color,alpha=0.35)
    j+=1

#fit a curve (exponential with base 10) to Hatzidimitriou's data using scipy.optimize
def func(x, a, b):
    return a*(10**(x)) + b
popt, pcov = curve_fit(func, hatz_data['log age'], hatz_data['d_br'])
plot(np.linspace(7,12,80), func(np.linspace(7,12,80), *popt), color='#2A527C', linestyle='--', alpha = 0.25)
print(tuple(popt))
print(np.linalg.cond(pcov))

iso_log_age = np.concatenate((iso_log_age_m7, iso_log_age_m8, iso_log_age_basti, iso_log_age_basti6, iso_log_age_parsec,iso_log_age_parsecmp))
iso_d_br = np.concatenate((iso_d_br_m7,iso_d_br_m8,iso_d_br_basti,iso_d_br_basti6,iso_d_br_parsec,iso_d_br_parsecmp))
iso_colors = np.concatenate(( np.full(14,0), np.full(20,-0.53), np.full(14,0.06), np.full(15,0.06), np.full(10,0), np.full(12,-0.50)))

#im = scatter(iso_log_age,iso_d_br,c=iso_colors,cmap=cm.coolwarm)
scatter(iso_log_age_m7,iso_d_br_m7,color='#e80576',edgecolor='#ffe8f5',marker='v',s=size,label='Padova (2008) [Fe/H]=0.00')
scatter(iso_log_age_m8,iso_d_br_m8,color='#ffe8f5',edgecolor='#e80576',marker='v',s=size,label='Padova (2008) [Fe/H]=–0.53')
scatter(iso_log_age_basti,iso_d_br_basti,color='#55b045',edgecolor='#a7db9e',marker='d',s=size,label='BaSTI (2004) [Fe/H]=+0.06')
scatter(iso_log_age_basti6,iso_d_br_basti6,color='#55b045',edgecolor='#a7db9e',marker='o',s=size,label='BaSTI (2018) [Fe/H]=+0.06')
scatter(iso_log_age_parsec,iso_d_br_parsec,color='#7c36b5',edgecolor='#ddb3ff',marker='s',s=size,label='PARSEC+COLIBRI (2012) [Fe/H]=0.00')
scatter(iso_log_age_parsecmp,iso_d_br_parsecmp,color='#ddb3ff',edgecolor='#7c36b5',marker='s',s=size,label='PARSEC+COLIBRI (2012) [Fe/H]=–0.50') 
ylim([-0.05,0.65])
xlim([8.1,10.3])
#title('Age Calibration')
xlabel('Age (log years)', fontsize=20)
ylabel(r'$d_{B-R}$ (mag)', fontsize=20)
xticks(fontsize = 20)
yticks(fontsize = 20)
legend(frameon=False,loc=9,fontsize=15)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cax.tick_params(labelsize=20)
cb = colorbar(im,cax=cax)
cb.set_label(label='[Fe/H]', size=20)
savefig('/Users/abbychriss/Desktop/theory_age_calibration.pdf', bbox_inches='tight', format='pdf')
close()
#show()

#output table: sorted by reverse alphabetical order
final_output = Table()
final_output['Cluster'] = np.array(final_data['cluster'],dtype=str)
final_output['Log age (years)'] = np.array(final_data['log age'],dtype=float)
final_output['Log age error (years)'] = np.array(final_data['error(log age)'],dtype=float)
final_output['d(B-R) (mag)'] = np.array(final_data['d_br'],dtype=float)
final_output['d(B-R) error (mag)'] = np.array(final_data['error(d_br)'],dtype=float)
final_output['d(BP-RP) (mag)'] = final_data['d_bprp']
final_output['d(BP-RP) error (mag)'] = final_data['error(d_bprp)']
final_output['E(BP-RP) (mag)'] = np.array([round(x,4) for x in final_data['e_bprp']])
final_output['nRC80'] = final_data['nRC80']
#final_output['Opening angle (deg)'] = np.array([round(x,4) for x in final_data['opening angle']])
ascii.write(final_output, 'age_calibrator_table.csv', format='csv', overwrite=True)
final_output = Table.read('age_calibrator_table.csv',format='ascii',fast_reader=False)

print(final_output['Cluster','d(BP-RP) (mag)','d(BP-RP) error (mag)'])
