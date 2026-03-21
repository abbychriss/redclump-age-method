#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 17:04:28 2023

@author: abbychriss
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from mpl_point_clicker import clicker
from astropy.table import Table
from astropy.io import ascii
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

#This program finds the theoretical d_B-R value for the ages in our study

#Isochrone function:

# for convenience, get the list of ages and metallicities (the latter all zero
# in this case.)
#open pops.solar files:
#11
a,feh=np.loadtxt('/Users/abbychriss/Desktop/WSU/pops.solar11',usecols=(0,1),unpack=True)
#basti
abasti,fehbasti=np.loadtxt('pops.solarbasti',usecols=(0,1),unpack=True,skiprows=1)
#print(len(a),len(feh))

def fetchiso(iiso):
    if iiso > 65:
        print('Error! iiso must be 65 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('hess_gaia_ev6res0mixss.out','r') as filehandle:
        filecontents = filehandle.readlines()

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
    rmag = np.array(rmag)
    indx = np.array(indx)
    # nHA is just an integer - the number of point in the isochrone
    return teff,logl,nst,gmag,bmag,rmag,indx,nHA

def fetchbasti(iiso):
    if iiso > 46:
        print('Error! iiso must be 46 or less')
    # seek the "iiso"th isochrone and return the data

    # The hess diagram files are composed of blocks, one block per SSP.
    with open('hess_gaia_ev5res0mixss.out','r') as filehandle:
        filecontents = filehandle.readlines()

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
    rmag = np.array(rmag)
    indx = np.array(indx)
    # nHA is just an integer - the number of point in the isochrone
    return teff,logl,nst,gmag,bmag,rmag,indx,nHA


# program execution: 

# User: change iiso at will
iiso = 50
teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchiso(iiso)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
br = bmag - rmag
im = scatter(br,gmag,label='Age '+str(a[iiso])+' Gyr (Padova)', cmap=cm.plasma_r, c=nst)
legend(frameon=False)
xlim([0.4,5.0])
ylim([17.0,-5.0])
xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
ylabel(r'G (mag)')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(im, cax=cax, label='dN')
show()

iiso = 51
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchiso(iiso)
br = bmag - rmag
im = scatter(br,gmag,label='Age '+str(a[iiso])+' Gyr (Padova)', c=nst, cmap=cm.plasma_r)
legend(frameon=False)
xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
ylabel(r'G (mag)')
divider = make_axes_locatable(ax)
xlim([0.4,5.0])
ylim([15.0,-5.0])
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(im, cax=cax, label='dN')
show()


# plot a BASTI isochrone with fetchbasti
iiso = 25
teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchbasti(iiso)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
br = bmag - rmag
indices = []
for i in range(len(br)):
    if 1 <= br[i] <= 2 and 0 < gmag[i] < 2:
        indices.append(i)
im = scatter([br[i] for i in indices],[gmag[i] for i in indices],label='Age '+str(abasti[iiso])+' Gyr (BASTI)', c=[nst[i] for i in indices], cmap=cm.plasma_r)
xlim([1,2])
ylim([2,-0])
legend(frameon=False)
divider = make_axes_locatable(ax)
xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
ylabel(r'G (mag)')
title('BASTI Isochrone')
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(im, cax=cax, label='dN')
show()

# plot a BASTI isochrone with fetchbasti
iiso = 30  
teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchbasti(iiso)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
br = bmag - rmag
indices = []
for i in range(len(br)):
    if 1 <= br[i] <= 2 and 0 < gmag[i] < 2:
        indices.append(i)
im = scatter([br[i] for i in indices],[gmag[i] for i in indices],label='Age '+str(abasti[iiso])+' Gyr (BASTI)', c=[nst[i] for i in indices], cmap=cm.plasma_r)
xlim([1,2])
ylim([2,-0])
legend(frameon=False)
divider = make_axes_locatable(ax)
xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
ylabel(r'G (mag)')
title('BASTI Isochrone')
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(im, cax=cax, label='dN')
show()


# plot a BASTI isochrone with fetchbasti
iiso = 35
teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchbasti(iiso)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
br = bmag - rmag
indices = []
for i in range(len(br)):
    if 1 <= br[i] <= 2 and 0 < gmag[i] < 2:
        indices.append(i)
im = scatter([br[i] for i in indices],[gmag[i] for i in indices],label='Age '+str(abasti[iiso])+' Gyr (BASTI)', c=[nst[i] for i in indices], cmap=cm.plasma_r)
xlim([1,2])
ylim([2,-0])
legend(frameon=False)
divider = make_axes_locatable(ax)
xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
ylabel(r'G (mag)')
title('BASTI Isochrone')
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(im, cax=cax, label='dN')
show()

# plot a BASTI isochrone with fetchbasti
iiso = 40  
teff,logl,nst,gmag,bmag,rmag,indx,nHA=fetchbasti(iiso)
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(12, 6))
br = bmag - rmag
indices = []
for i in range(len(br)):
    if 1 <= br[i] <= 2 and 0 < gmag[i] < 2:
        indices.append(i)
im = scatter([br[i] for i in indices],[gmag[i] for i in indices],label='Age '+str(abasti[iiso])+' Gyr (BASTI)', c=[nst[i] for i in indices], cmap=cm.plasma_r)
xlim([1,2])
ylim([2,-0])
legend(frameon=False)
divider = make_axes_locatable(ax)
xlabel(r'G$_{BP}$ - G$_{RP}$ (mag)')
ylabel(r'G (mag)')
title('BASTI Isochrone')
cax = divider.append_axes("right", size="5%", pad=0.05)
colorbar(im, cax=cax, label='dN')
show()



