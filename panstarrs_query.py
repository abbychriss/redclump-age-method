#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 15:56:09 2024

@author: abbychriss
"""

from astropy.io import ascii
from astropy.table import Table

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import json
import requests

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib

cluster_names = np.array(['Collinder_110', 'King_5', 'Melotte_71', 'NGC_1245', 'NGC_1907', 'NGC_2420', 'NGC_2506', 'NGC_2509', 'NGC_2627', 'NGC_2682', 'NGC_6791', 'NGC_6939', 'NGC_6940', 'NGC_752', 'NGC_7789'])
url = f'{https://catalogs.mast.stsci.edu/api/v0.1/panstarrs}'
for name in cluster_names:
    // POST
    /panstarrs/dr2/mean/crossmatch/upload
    file='/Users/abbychriss/Desktop/'+name+'_cross_ref.csv'
    header=True
    resolve=False
    radius=3
    ra_name=ra
    dec_name=dec
    target_name=target
