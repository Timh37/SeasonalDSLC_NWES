#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:56:14 2020

Finds latitudes and longitudes (out_lat, out_lon, min_idx) of unmasked grid 
cells in the grid of the input CMIP6 data array (da) nearest (min_dists) to 
query latitude/longitude points (qlats, qlons).

@author: thermans
"""

import numpy as np
import xarray as xr

def angdist(lats,lons,qlats,qlons): #calculate angular distance
    lat0,lat = np.meshgrid(np.radians(lats),np.radians(qlats))
    lon0,lon = np.meshgrid(np.radians(lons),np.radians(qlons))
    
    temp = np.arctan2(np.sqrt((np.cos(lat)*np.sin(lon-lon0))**2 + (np.cos(lat0)*np.sin(lat) - np.sin(lat0)*np.cos(lat) * np.cos(lon-lon0))**2),
                      (np.sin(lat0)*np.sin(lat) + np.cos(lat0)*np.cos(lat)*np.cos(lon-lon0)))
    return(np.degrees(temp))

def cmip6_find_nearest(da,qlats,qlons):
    #finds nearest unmasked ocean grid cell in xarray dataarray to query latitude
    #and longitude pairs based on minimum angular distance [deg]
    
    #fetch coordinate names
    lonname = np.array(da.coords)[['lon' in x for x in da.coords]][0]
    latname = np.array(da.coords)[['lat' in x for x in da.coords]][0]

    #get lats & lons
    lats = np.array(da[latname])
    lons = np.array(da[lonname])

    if lats.shape!=lons.shape:
        lats,lons = np.meshgrid(lats,lons)
    
    if 'time' in list(da.coords):
        if np.transpose(lats).shape == da.isel(time=0).shape:
            lats = np.transpose(lats)
            lons = np.transpose(lons)
        
    #calculate angular distances between query coordinates and data array grid
    dists=angdist(lats,lons,qlats,qlons)

    #mask land out
    if 'time' in da.dims:
        dists[:,~np.isfinite(da.isel(time=0).values.flatten())] = np.nan
    else:
        dists[:,~np.isfinite(da.values.flatten())] = np.nan
    
    min_dists = np.nanmin(dists,axis=1) #find minimum angular distances in grid to query points
    min_idx = np.nanargmin(dists,axis=1) #indices
    
    out_lat = lats.flatten()[min_idx]
    out_lon = lons.flatten()[min_idx]
    #potentially build in a filter here if unreasonably large distances

    return min_dists, min_idx, out_lat, out_lon

if __name__ == '__main__':
    var_ds= xr.open_dataset('/Volumes/Naamloos/PhD_Data/CMIP6/raw/zos/ACCESS-CM2/zos_Omon_ACCESS-CM2_piControl_r1i1p1f1_gn_095001-144912.nc')
    
    qlats = np.array([33, 39, 40,30])
    qlons = np.array([5, 5, 5,4])
    
    min_dists,min_idx,out_lat,out_lon=cmip6_find_nearest(var_ds.zos,qlats,qlons)
   