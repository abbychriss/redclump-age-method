#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 19:16:43 2023

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
from decimal import Decimal
import sys
from IPython.display import set_matplotlib_formats

rgb_errors = [0.03, 0.06, 0.01, 0.01, 0.035, 0.015, 0.01, 0.02, 0.01, 0.035, 0.01, 0.04, 0.01, 0.03, 
              0.015, 0.06, 0.11, 0.025, 0.095, 0.025, 0.015, 0.025]
n_rgb = [6,5,73,4,6,18,222,5,32,8,4,9,43,48,0,20,6,13,6,104,9,40]

scatter(n_rgb,rgb_errors)
n=np.linspace(0.7,250)
y = 0.1 * (1 /  sqrt(n) )
plot(n,y)
xlabel('N')
ylabel('RGB Error in BP-RP (mag)')
show()