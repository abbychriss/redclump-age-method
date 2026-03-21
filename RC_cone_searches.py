#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:58:45 2023

@author: abbychriss
"""

import matplotlib.pyplot as plt
import numpy as np
from pylab import *

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

# Uses Astroquery to do a Gaia DR3 cone search
# User should edit "coord", "radius", and "outputfilename" for each cluster.    

# perform a cone search and get the data
radii = np.array([0.157, 0.0417, 0.045, 0.171, 0.23, 0.534, 0.091, 0.088, 0.075, 0.109, 0.075, 0.07, 0.053, 0.202, 0.15, 0.088, 0.06, 0.087, 0.035, 0.208, 0.162, 0.068, 0.123, 0.25, 0.485, 0.211, 0.033, 0.074], dtype=float)
coords = [SkyCoord(ra=099.6770, dec=+02.0690, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=268.3200, dec=-27.3730, unit=(u.degree, u.degree), frame='icrs'),
         SkyCoord(ra=110.6590, dec=-18.7030, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=169.3730, dec=-62.7190, unit=(u.degree, u.degree), frame='icrs'),
         SkyCoord(ra=261.2120, dec=-49.9170, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=279.6490, dec=+05.4350,unit=(u.degree, u.degree), frame='icrs'),
         SkyCoord(ra=048.6820, dec=+52.6950, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=111.5730, dec=-47.6850, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=114.3830, dec=-12.0650, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=048.6910, dec=+47.2350, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=082.0330, dec=+35.3300, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=097.4160, dec=+06.8340, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=114.6020, dec=+21.5750, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=116.1410, dec=-23.8530, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=118.0460, dec=-38.5370, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=120.0100, dec=-10.7730, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=120.2010, dec=-19.0560, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=129.3090, dec=-29.9520, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=130.6670, dec=-47.2010, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=132.8460, dec=+11.8140, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=252.3360, dec=-53.7140, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=290.2210, dec=+37.7780, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=307.9170, dec=+60.6530, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=308.6260, dec=+28.2780, unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=029.2230, dec=+37.7940, unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=359.3340, dec=+56.7260,unit=(u.degree, u.degree), frame='icrs'), SkyCoord(ra=204.2270, dec=-62.0910,unit=(u.degree, u.degree), frame='icrs'), 
         SkyCoord(ra=131.1470, dec=-35.9000,unit=(u.degree, u.degree), frame='icrs') ]
outputfilenames = ['/Users/abbychriss/Desktop/WSU/rup68_cone.csv', '/Users/abbychriss/Desktop/WSU/pismis18_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc7789_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc752_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc6940_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc6939_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc6791_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc6208_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2682_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc2660_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2627_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2509_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc2506_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2477_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2447_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc2420_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc2236_cone.csv', '/Users/abbychriss/Desktop/WSU/ngc1907_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ngc1245_cone.csv', '/Users/abbychriss/Desktop/WSU/mel71_cone.csv', '/Users/abbychriss/Desktop/WSU/mel66_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/king5_cone.csv', '/Users/abbychriss/Desktop/WSU/ic4756_cone.csv', '/Users/abbychriss/Desktop/WSU/ic4651_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/ic2714_cone.csv', '/Users/abbychriss/Desktop/WSU/fsr1252_cone.csv', '/Users/abbychriss/Desktop/WSU/czernik37_cone.csv',
                  '/Users/abbychriss/Desktop/WSU/col110_cone.csv']
outputfilenames.reverse()
clusternames = ['Ruprecht 68', 'Pismis 18', 'NGC 7789', 'NGC 752', 'NGC 6940', 'NGC 6939', 'NGC 6791', 'NGC 6208', 'NGC 2682', 'NGC 2660', 'NGC 2627', 'NGC 2509', 'NGC 2506', 'NGC 2477', 'NGC 2447', 'NGC 2420', 'NGC 2236', 'NGC 1907', 'NGC 1245', 'Melotte 71', 'Melotte 66', 'King 5', 'IC 4756', 'IC 4651', 'IC 2714', 'FSR 1252', 'Czernik 37', 'Collinder 110']
clusternames.reverse()

#add NGC 6791, NGC 2682 (M67), Mel66 290.2210 +37.7780 r=0.0683 plx: 0.192 132.8460 +11.8140 r=0.208 plx: 1.1325 111.5730 -47.6850 r=0.0883 plx 0.183

# Collinder 110 (r=0.157 deg), Czernik 37 (r=0.0417), FSR 1252 (r=0.045), IC 2714 (r=0.171), IC 4651 (r=0.23), IC 4756 (r=0.534), King 5 (r=0.091), Melotte 66 (r=0.0883)
# Melotte 71 (r=0.075), NGC 1245 (r=0.109), NGC 1907 (r=0.075), NGC 2236 (r=0.07), NGC 2420 (r=0.053), NGC 2447 (r=0.202), NGC 2477 (r=0.15),
# NGC 2506 (r=0.088), NGC 2509 (r=0.06), NGC 2627 (r=0.087), NGC 2660 (r=0.035), NGC 2682 (r=0.208), NGC 6208 (r=0.162), NGC 6791 (r=0.068), NGC 6939 (r=0.123), NGC 6940 (r=0.25), NGC 752 (r=0.485),
# NGC 7789 (r=0.211), Pismis 18 (r=0.033), Ruprecht 68 (r=0.074)

i=19
coord = coords[i]
outputfilename = outputfilenames[i]
radius = u.Quantity(3*radii[i], u.deg) #take 2, 3, or 4 times the radius given by simbad

# OK, user can fall asleep, now. Data should be retrieved, then written to output file in csv format.

Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"  # Reselect Data Release 3, 3 is NOT default???
Gaia.ROW_LIMIT = -1  # -1 is "unlimited." <-- be careful with that.
j = Gaia.cone_search_async(coord, radius)
r = j.get_results()
# "r" is a data structure of type "Table" from astropy.table

# Here are some methods for dealing with a "Table"
# More information at https://docs.astropy.org/en/stable/table/

#put in an optional parallax cut for very large clusters

#r.pprint() # print bits of the table
print(r.info) # list the column names (and more, in a pretty format)
#print(r.colnames)  # list the column names crudely

# Extract a column into a 1-d array by name:
name = r['DESIGNATION']
ra = r['ra']
#r['ra_error']
dec = r['dec']
#r['dec_error']
plx = r['parallax']
plxe = r['parallax_error']
pmra = r['pmra']
pmrae = r['pmra_error']
pmdec = r['pmdec']
pmdece = r['pmdec_error']
g = r['phot_g_mean_mag']
#ge = r['phot_g_mean_mag_error']
bp = r['phot_bp_mean_mag']
#bpe = r['phot_bp_mean_mag_error']
rp = r['phot_rp_mean_mag']
#rpe = r['phot_rp_mean_mag_error']
bprp = r['bp_rp']   # BP - RP color, as observed (no dust correction)
# things likely to not exist for every star:
rv = r['radial_velocity']  # km/s
teff= r['teff_gspphot'] # Teff in K
Ag = r['ag_gspphot']  # extinction in G band
Ebr = r['ebpminrp_gspphot']  # color excess E(BP-RP)
vflag = r['phot_variable_flag'] # integer. Probably equals one if Gaia thinks the star is variable

print('Found ',len(pmra),' stars.')

#plot proper motions
plot(r['pmra'],r['pmdec'],'b+')
xlabel(r'$\mu_\alpha$')
ylabel(r'$\mu_\delta$')
title('Proper Motions: ' + str(clusternames[i]))
show()

# filter for NaN (IEEE Not A Number) in parallax (no use in saving these data)
nanfilt = isnan(plx)
distfilt = plx > 0.5

# write cone search results to a file
import csv
with open(outputfilename,'w',newline='') as csvfile:
     cw = csv.writer(csvfile, quoting=csv.QUOTE_MINIMAL)
     cw.writerow(['Gaia ID','RA','DEC','parallax','error in plx','PM(RA)','error in pmra','PM(dec)','error in pmdec','g magnitude','BP','RP','BP-RP','rv','Teff','Ag extinction','E(BP-RP)','variability flag'])
     for i in range(len(ra)):
          if (nanfilt[i] == False):
               cw.writerow([name[i],ra[i],dec[i],plx[i],plxe[i],pmra[i],pmrae[i],pmdec[i],pmdece[i],g[i],bp[i],rp[i],bprp[i],rv[i],teff[i],Ag[i],Ebr[i],vflag[i]])

print( 'normal stop' )
