#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 08:43:46 2021

Loops over 'deptho' files from CMIP6 models and outputs the depth 
near specified lat/lon coordinates

@author: thermans
"""

import xarray as xr
import numpy as np
import os
import fnmatch
from cmip6_find_nearest import cmip6_find_nearest

deptho_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/raw/deptho' #CMIP6 directory with deptho files ordered per model

mindepth = np.array([])

for m,model in enumerate((model for model in os.listdir(deptho_dir) if model!='.DS_Store')): #loop over models

    deptho_model_dir = os.path.join(deptho_dir,model)

    try:
        fn  = fnmatch.filter(os.listdir(deptho_model_dir), "*deptho*gn*nc")[0]
    except:
        print(model + ' has no gn deptho')
        continue
    
    ds = xr.open_dataset(os.path.join(deptho_model_dir,fn)) #open deptho file
    
    #if deptho is 3D variable
    if len(ds.deptho.shape)==3:
        try:
            ds = ds.max(dim='lev')
        except:
            ds = ds.max(dim='deptht')
    
    #want to mask out 0 depth as well, it's land
    ds['deptho'] = ds['deptho'].where(ds.deptho!=0)
    
    min_dist,min_idx,out_lat,out_lon = cmip6_find_nearest(ds.deptho,52,5) #find nearest to query coordinates
    
    mindepth = np.append(mindepth,ds.deptho.values.flatten()[min_idx])
    
    print(model + ', ' + str(ds.deptho.values.flatten()[min_idx][0])+ ' m')
