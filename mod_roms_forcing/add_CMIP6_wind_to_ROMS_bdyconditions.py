#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 16:11:35 2019

Applies change in seasonal anomalies of wind-velocity change from CMIP6 model to ERA5 boundary conditions to ROMS model.

@author: thermans
"""

import numpy as np
import cmocean
import os
import fnmatch
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xesmf as xe 
import datetime
from seasonal_deviations_from_monthly_means import seasonal_deviations_from_monthly_means
from cdo import Cdo
cdo = Cdo()

model = 'CNRM-ESM2-1' #CMIP6 model to derive wind from
cmip6_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/' #path to CMIP6 data
ssp = 'ssp585' #SSP

in_dir = '/Volumes/Naamloos/PhD_Data/ERA5/ROMS_format/'; #input path bdy conditions
out_dir = '/Volumes/Naamloos/PhD_Data/ERA5/ROMS_format/wind_mod/'; #output path new bdy conditions
out_prefix = model+'_DJF_wind_' #prefix to prepend to output filename

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

years_to_mod = np.arange(1993,1996) #years to modify in ROMS forcing

#averaging periods
basey = np.arange(1995,2015)
futy = np.arange(2081,2101)

#load variant averaged wind velocity components from cmip6 model
#uas
var_model_dir = os.path.join(cmip6_dir,'uas',model)
    
for f,fn in enumerate(fnmatch.filter(os.listdir(var_model_dir), "*"+ssp+"*nc")): #loop over variants
    #load ssp & historical
    ssp_ds = xr.open_dataset(os.path.join(var_model_dir,fn))
    print(ssp_ds.variant_label)
    hist_ds = xr.open_dataset(os.path.join(var_model_dir,fnmatch.filter(os.listdir(var_model_dir), "*historical*"+ssp_ds.variant_label+"*nc")[0]))
    
    #concatenate historical and ssp and replace model specific time with common time
    hist_ds['time'] = xr.cftime_range(start="1850", periods=len(hist_ds.time), freq="MS", calendar="noleap") #replace model specific time array by common
    ssp_ds['time'] = xr.cftime_range(start="2015", periods=len(ssp_ds.time), freq="MS", calendar="noleap")
    hist_ssp_ds = xr.concat((hist_ds,ssp_ds),dim='time') #concatenate historical & ssp over time dimension
    hist_ssp_ds = xr.decode_cf(hist_ssp_ds)
    if 'time_bounds' in hist_ssp_ds:
        hist_ssp_ds = hist_ssp_ds.rename(time_bounds='time_bnds')
    if f==0: #append variants in 1 ds
        model_ds = []
        model_ds = hist_ssp_ds
    else:
        model_ds = xr.concat((model_ds,hist_ssp_ds),dim='variant') #generating new dimension 'variant'    

if 'variant' in model_ds.dims: #if multiple variants, take variant mean
    uas = model_ds.mean(dim='variant')
else:
    uas = model_ds

#vas, idem
var_model_dir = os.path.join(cmip6_dir,'vas',model)
    
for f,fn in enumerate(fnmatch.filter(os.listdir(var_model_dir), "*"+ssp+"*nc")):
    ssp_ds = xr.open_dataset(os.path.join(var_model_dir,fn))
    print(ssp_ds.variant_label)
    hist_ds = xr.open_dataset(os.path.join(var_model_dir,fnmatch.filter(os.listdir(var_model_dir), "*historical*"+ssp_ds.variant_label+"*nc")[0]))
    
    hist_ds['time'] = xr.cftime_range(start="1850", periods=len(hist_ds.time), freq="MS", calendar="noleap")
    ssp_ds['time'] = xr.cftime_range(start="2015", periods=len(ssp_ds.time), freq="MS", calendar="noleap")
    hist_ssp_ds = xr.concat((hist_ds,ssp_ds),dim='time')
    hist_ssp_ds = xr.decode_cf(hist_ssp_ds)
    if 'time_bounds' in hist_ssp_ds:
        hist_ssp_ds = hist_ssp_ds.rename(time_bounds='time_bnds')
    if f==0:
        model_ds = []
        model_ds = hist_ssp_ds
    else:
        model_ds = xr.concat((model_ds,hist_ssp_ds),dim='variant') #generating new dimension 'model'    

if 'variant' in model_ds.dims:
    vas = model_ds.mean(dim='variant')
else:
    vas = model_ds


#compute change in seasonal anomalies of the wind components
uas = uas.isel(time=np.arange(11+124*12,3011))
vas = vas.isel(time=np.arange(11+124*12,3011))

uas_djf_mean_anom, uas_jja_mean_anom, uas_mam_mean_anom, uas_son_mean_anom, uas_dec2nov_mean = seasonal_deviations_from_monthly_means(uas.uas)
vas_djf_mean_anom, vas_jja_mean_anom, vas_mam_mean_anom, vas_son_mean_anom, vas_dec2nov_mean = seasonal_deviations_from_monthly_means(vas.vas)

uas_djf_mean_anom_d = uas_djf_mean_anom.sel(year=futy).mean(dim='year')-uas_djf_mean_anom.sel(year=basey).mean(dim='year')
vas_djf_mean_anom_d = vas_djf_mean_anom.sel(year=futy).mean(dim='year')-vas_djf_mean_anom.sel(year=basey).mean(dim='year')

#replace wind over land with wind from nearest over-ocean grid cell
lf_path = os.path.join('/Volumes/Naamloos/PhD_Data/CMIP6/raw/sftlf',model) #get land fraction of atmospheric grid cells
lf_file = os.path.join(lf_path,fnmatch.filter(os.listdir(lf_path),'*sftlf*')[0])
lf = xr.open_dataset(lf_file)

if model=='UKESM1-0-LL': #U&V points on different grids, get them on a common grid
    v = np.empty(np.shape(lf.sftlf))
    v[:] = np.nan
    v = (vas_djf_mean_anom_d[0:-1,:].values + vas_djf_mean_anom_d[1:,:].values)/2
    
    uas_djf_mean_anom_d_rep = np.empty((np.shape(lf.sftlf)[0],np.shape(lf.sftlf)[1] + 1))
    uas_djf_mean_anom_d_rep[:,0:-1] = uas_djf_mean_anom_d.values
    uas_djf_mean_anom_d_rep[:,-1] = uas_djf_mean_anom_d.values[:,0]
    
    u = np.empty(np.shape(lf.sftlf))
    u[:] = np.nan
    u = (uas_djf_mean_anom_d_rep[:,0:-1] + uas_djf_mean_anom_d_rep[:,1:])/2
    
    uas_djf_mean_anom_d = xr.DataArray(data=u,dims=["lat", "lon"],coords=[lf.lat, lf.lon])
    vas_djf_mean_anom_d = xr.DataArray(data=v,dims=["lat", "lon"],coords=[lf.lat, lf.lon])

elif model == 'ACCESS-ESM1-5': #U&V points on different grids, get them on a common grid
    v = np.empty((np.shape(vas_djf_mean_anom_d)[0]-1,np.shape(vas_djf_mean_anom_d)[1]))
    v[:] = np.nan
    v = (vas_djf_mean_anom_d[0:-1,:].values + vas_djf_mean_anom_d[1:,:].values)/2
    
    u = np.empty((np.shape(uas_djf_mean_anom_d)[0],np.shape(uas_djf_mean_anom_d)[1]-1))
    u[:] = np.nan
    u = (uas_djf_mean_anom_d[:,0:-1].values + uas_djf_mean_anom_d[:,1:].values)/2
    
    uas_djf_mean_anom_d = xr.DataArray(data=u,dims=["lat", "lon"],coords=[lf.lat, lf.lon[1:]])
    vas_djf_mean_anom_d = xr.DataArray(data=v,dims=["lat", "lon"],coords=[lf.lat[1:-1], lf.lon])
    
elif model == 'MPI-ESM1-2-LR': #same grids u&V and land fraction, but some precision differences
    lf['lon'] = uas_djf_mean_anom_d.lon
    lf['lat'] = uas_djf_mean_anom_d.lat #differences order 10e-14

#apply simple extrapolation correction for land contamination
uas_djf_mean_anom_d_masked = uas_djf_mean_anom_d.where(lf.sftlf==0,drop=False) #mask over-land grid cells
vas_djf_mean_anom_d_masked = vas_djf_mean_anom_d.where(lf.sftlf==0,drop=False)

#extrapolate over-ocean to over-land winds using nearest neighbor using CDO (slightly cumbersome, this could be improved)
if os.path.exists('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_masked.nc'):
    os.remove('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_masked.nc')
if os.path.exists('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_masked.nc'):
    os.remove('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_masked.nc')
    
if os.path.exists('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn.nc'):
    os.remove('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn.nc')
if os.path.exists('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn.nc'):
    os.remove('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn.nc')

if os.path.exists('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn_wrapped.nc'):
    os.remove('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn_wrapped.nc')
if os.path.exists('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn_wrapped.nc'):
    os.remove('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn_wrapped.nc')    
    
uas_djf_mean_anom_d_masked.to_netcdf('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_masked.nc')
vas_djf_mean_anom_d_masked.to_netcdf('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_masked.nc')

cdo.setmisstonn(input='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_masked.nc',
                output='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn.nc') #extrapolate ocean to land
cdo.sellonlatbox(-180,180,-90,90,input='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn.nc',
                output='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn_wrapped.nc') #wrap around lon 

cdo.setmisstonn(input='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_masked.nc',
                output='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn.nc')                                
cdo.sellonlatbox(-180,180,-90,90,input='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn.nc',
                output='/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn_wrapped.nc')

if model == 'MPI-ESM1-2-LR':
    uas_mod = xr.open_dataset('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn.nc') #cdo wrapping somehow doesnt work
    vas_mod = xr.open_dataset('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn.nc')
    
    uas_mod.coords['lon'] = ((uas_mod.coords['lon'] + 180) % 360) - 180 #wrap around 0
    uas_mod = uas_mod.reindex({ 'lon' : np.sort(uas_mod['lon'])})
    vas_mod.coords['lon'] = ((vas_mod.coords['lon'] + 180) % 360) - 180 #wrap around 0
    vas_mod = vas_mod.reindex({ 'lon' : np.sort(vas_mod['lon'])})
else:
    uas_mod = xr.open_dataset('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_uas_d_'+model+'_landfillednn_wrapped.nc') 
    vas_mod = xr.open_dataset('/Users/thermans/Documents/PhD/Phase4_seasonal/djf_anom_vas_d_'+model+'_landfillednn_wrapped.nc') 

#wrap longitude
uas_djf_mean_anom_d.coords['lon'] = ((uas_djf_mean_anom_d.coords['lon'] + 180) % 360) - 180 #wrap around 0
uas_djf_mean_anom_d = uas_djf_mean_anom_d.reindex({ 'lon' : np.sort(uas_djf_mean_anom_d['lon'])})
vas_djf_mean_anom_d.coords['lon'] = ((vas_djf_mean_anom_d.coords['lon'] + 180) % 360) - 180 #wrap around 0
vas_djf_mean_anom_d = vas_djf_mean_anom_d.reindex({ 'lon' : np.sort(vas_djf_mean_anom_d['lon'])})

#extrapolation from over-ocean wind only if leading to higher wind speeds
if 'uas' in uas_mod:
    uas_mod=uas_mod.where(abs(uas_mod.uas)>abs(uas_djf_mean_anom_d),uas_djf_mean_anom_d).uas
    vas_mod=vas_mod.where(abs(vas_mod.vas)>abs(vas_djf_mean_anom_d),vas_djf_mean_anom_d).vas
else:        
    uas_mod=uas_mod.where(abs(uas_mod.__xarray_dataarray_variable__)>abs(uas_djf_mean_anom_d),uas_djf_mean_anom_d).__xarray_dataarray_variable__
    vas_mod=vas_mod.where(abs(vas_mod.__xarray_dataarray_variable__)>abs(vas_djf_mean_anom_d),vas_djf_mean_anom_d).__xarray_dataarray_variable__

#apply the CMIP6 wind velocity change as a constant offset of ERA5 wind forcing
fn = 'ERA5_NorthAtlantic_an_1993_ForROMS.nc' #generate interpolation weights
ds = xr.open_dataset(os.path.join(in_dir,fn))   

#regrid to ERA5 grid
regridder = xe.Regridder(uas_mod,ds,'bilinear')
uas_mod_regrid = regridder(uas_mod)
regridder = xe.Regridder(vas_mod,ds,'bilinear')
vas_mod_regrid = regridder(vas_mod)

for year in years_to_mod: #regrid, add and save
    fn = 'ERA5_NorthAtlantic_an_'+str(year)+'_ForROMS.nc'
    ds = xr.open_dataset(os.path.join(in_dir,fn))
    
    with xr.set_options(keep_attrs=True): 
       ds['Uwind'] = ds['Uwind'] + uas_mod_regrid
       ds['Vwind'] = ds['Vwind'] + vas_mod_regrid
    
    #overwrite wind variables & save into new file
    ds.attrs['comments'] = 'additional wind component added'
    ds.attrs['mod_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
        
    ds.Uwind.attrs['coordinates'] = 'lon lat'
    ds.Vwind.attrs['coordinates'] = 'lon lat'
    
    with xr.set_options(keep_attrs=True):    
        ds.to_netcdf(os.path.join(out_dir,out_prefix+fn),mode='w',encoding={'Uwind':{'zlib': True,'complevel': 4},'Vwind':{'zlib': True,'complevel': 4}}) #(over)write a new NetCDF


