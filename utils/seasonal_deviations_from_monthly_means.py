#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 11:09:30 2021

Takes data array (da) with monthly means along time dimension and calculates
DJF, MAM, JJA and SON means relative to annual mean starting in December and
ending in November

@author: thermans
"""
import xarray as xr
import numpy as np

def seasonal_deviations_from_monthly_means(da):
    s = list(da.isel(time=0).shape) #get shape of data array without time
    s.append(int(len(da.time)/12)) #append years and months dimensions
    s.append(12)
    
    my_dims = list(da.isel(time=0).dims) #get dimensions without time
    
    my_coords=[]
    for dim in my_dims:
        my_coords.append(da[dim])
    
    my_dims.append('year') #add year dimension
    
    my_years = np.arange(da.time[0].dt.year.values+1,da.time[-1].dt.year.values+1) #add year coordinates
    my_coords.append(my_years)
    
    ordered_per_year = np.resize(da.transpose(...,'time'),tuple(s)) #order per year & take seasonal means
    
    dec2nov_mean = xr.DataArray(np.mean(ordered_per_year,axis=-1),coords=my_coords,dims=my_dims)
    djf_mean = xr.DataArray(np.mean(ordered_per_year[...,0:3],axis=-1),coords=my_coords,dims=my_dims)
    jja_mean = xr.DataArray(np.mean(ordered_per_year[...,6:9],axis=-1),coords=my_coords,dims=my_dims)
    mam_mean = xr.DataArray(np.mean(ordered_per_year[...,3:6],axis=-1),coords=my_coords,dims=my_dims)
    son_mean = xr.DataArray(np.mean(ordered_per_year[...,9:12],axis=-1),coords=my_coords,dims=my_dims)
    
    #anomalies wrt Dec-Nov annual mean
    djf_mean_anom = djf_mean - dec2nov_mean
    jja_mean_anom = jja_mean - dec2nov_mean
    mam_mean_anom = mam_mean - dec2nov_mean
    son_mean_anom = son_mean - dec2nov_mean
        
    return djf_mean_anom, jja_mean_anom, mam_mean_anom, son_mean_anom, dec2nov_mean