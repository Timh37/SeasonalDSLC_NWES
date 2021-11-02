#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:11:35 2019

Applies uniform wind-speed change to ERA5 boundary conditions to ROMS model.

@author: thermans
"""

import numpy as np
import os
import xarray as xr
import datetime

in_dir = '/Volumes/Naamloos/PhD_Data/ERA5/ROMS_format/'; #input path bdy conditions
out_dir = '/Volumes/Naamloos/PhD_Data/ERA5/ROMS_format/wind_mod/'; #output path new bdy conditions

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

#prefix to prepend to filename
out_prefix = 'uniform_swWindsq2ms_' #prefix to prepend to output filename

years_to_mod = np.arange(1993,1996) #years to modify in ROMS forcing

for year in years_to_mod: #apply constant offsets to ERA5 wind forcing
    fn = 'ERA5_NorthAtlantic_an_'+str(year)+'_ForROMS.nc'
    ds = xr.open_dataset(os.path.join(in_dir,fn))
    
    with xr.set_options(keep_attrs=True): 
       ds['Uwind'] = ds['Uwind'] + 1 #adapt according to desired magnitude & direction
       ds['Vwind'] = ds['Vwind'] + 1
    
    #overwrite wind variables & save into new file
    ds.attrs['comments'] = 'additional wind component added'
    ds.attrs['mod_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
        
    ds.Uwind.attrs['coordinates'] = 'lon lat'
    ds.Vwind.attrs['coordinates'] = 'lon lat'
    
    with xr.set_options(keep_attrs=True):    
        ds.to_netcdf(os.path.join(out_dir,out_prefix+fn),mode='w',encoding={'Uwind':{'zlib': True,'complevel': 4},'Vwind':{'zlib': True,'complevel': 4}}) #(over)write a new NetCDF
