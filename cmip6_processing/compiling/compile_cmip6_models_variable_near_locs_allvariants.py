#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 08:43:46 2021

Compile CMIP6 simulations at native grid cells nearest to query locations into ensemble,
combining  historical and SSP runs per model

@author: thermans
"""

import xarray as xr
import numpy as np
import os
import fnmatch
import datetime
import cftime
from cmip6_find_nearest import cmip6_find_nearest

## user input
target_ssps = ['ssp126']#,'ssp245','ssp370','ssp585'] #desired ssps

variable = 'zos'
var_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/zos_0mean'
out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/'

locations = ['Liverpool','Brest','Newhaven', 'Den Helder', 'Cuxhaven', 'Esbjerg', 'Immingham', 'Aberdeen'] #query location coordinates
qlon = [-3.02,-4.49,0.06,4.8,8.7,8.44,0.187,-2.08]
qlat = [53.45,48.38,50.78,52.95,53.86,55.46,53.63,57.14]
##

for ssp in target_ssps: #loop over desired ssps
    print(ssp) #debugging
    
    ens_models = np.array([]) #model names
    
    first_model=True
    for model in (model for model in os.listdir(var_dir) if model not in ['.DS_Store','AWI-CM-1-1-MR']):
        variants = np.array([])
        print("current model "+model) #debugging
        var_model_dir = os.path.join(var_dir,model) #variable directory

        if not fnmatch.filter(os.listdir(var_model_dir), "*mon*"+ssp+"*nc"): #look for SSP
            print(variable+" for "+ssp+" not available for current model") #if variable not available for current ssp
            continue
        
        for f,fn in enumerate(fnmatch.filter(os.listdir(var_model_dir), "*"+ssp+"*nc")): #loop over variants
            #load ssp & historical
            ssp_ds = xr.open_dataset(os.path.join(var_model_dir,fn))
            try:
                hist_ds = xr.open_dataset(os.path.join(var_model_dir,fnmatch.filter(os.listdir(var_model_dir), "*historical*"+ssp_ds.variant_label+"*nc")[0]))
            except:
                print('no hist available')
                continue
            print(ssp_ds.variant_label)
            variant = ssp_ds.variant_label
            #replace model specific time with common time
            hist_ds['time'] = xr.cftime_range(start="1850", periods=len(hist_ds.time), freq="MS", calendar="noleap") #replace model specific time array by common
            ssp_ds['time'] = xr.cftime_range(start="2015", periods=len(ssp_ds.time), freq="MS", calendar="noleap")
            
            #find nearest coordinates to query locations
            if f==0: #only need to do this once per model
                min_dist,min_idx,out_lat,out_lon = cmip6_find_nearest(ssp_ds[variable],qlat,qlon)
                print(out_lat)
                print(out_lon)
                
            if (len(ssp_ds[variable].shape) > 2): #if lat and lon unraveled in variable data, reshape
                ssp_rs = np.reshape(ssp_ds[variable].values,(np.shape(ssp_ds[variable])[0],-1))
                hist_rs = np.reshape(hist_ds[variable].values,(np.shape(hist_ds[variable])[0],-1))
                
                ssp_ds = xr.DataArray(ssp_rs[:,min_idx],coords=[ssp_ds.time,locations],dims=['time','location'])
                hist_ds = xr.DataArray(hist_rs[:,min_idx],coords=[hist_ds.time,locations],dims=['time','location'])
                
            else: #indix directly
                ssp_rs = ssp_ds[variable].values[:,min_idx]
                hist_rs = hist_ds[variable].values[:,min_idx]
                
                ssp_ds = xr.DataArray(ssp_rs[:,min_idx],coords=[ssp_ds.time,locations],dims=['time','location'])
                hist_ds = xr.DataArray(hist_rs[:,min_idx],coords=[hist_ds.time,locations],dims=['time','location'])
             
            #concatenate historical and ssp
            hist_ssp_ds = xr.concat((hist_ds,ssp_ds),dim='time') #concatenate historical & ssp over time dimension
            
            if model=='CAMS-CSM1-0': #stops in 209912, extrapolate to 210012 by adding average year over 2081-2099
                hist_ssp_ds=hist_ssp_ds.interp(time=xr.cftime_range(start="1850", periods=len(hist_ssp_ds.time)+12, freq="MS", calendar="noleap")
                                   ,kwargs={'fill_value': 'extrapolate'})
                hist_ssp_ds[-12::,...] = hist_ssp_ds[-240::,:].groupby('time.month').mean().values
                
            #hist_ssp_ds = xr.decode_cf(hist_ssp_ds)
            
            if f==0: #if first variant, produce datarray to hold variant data
                model_ds = []
                model_ds = hist_ssp_ds
            else:
                model_ds = xr.concat((model_ds,hist_ssp_ds),dim='variant') #append variant to existing and generating new dimension 'model'    
            
            variants = np.append(variants,variant)
        if f==0:
            model_ds = model_ds.expand_dims('variant')
        model_ds = model_ds.assign_coords(variant=variants)
        
        ### old code for saving mean of variants, now saving all variants instead
        #if ((f>0) and ('variant' in model_ds.dims)):
        #    variant_mean = model_ds.mean(dim='variant') #take mean across variants if more than one variant
        #else:
        #    variant_mean = model_ds #else use single variant
        ###
        
        if first_model: #concatenate variant means along model dimension
            ens_ds=[]
            ens_ds = model_ds #initialise ensemble ds with model ds
            first_model=False
        else: #append next models
            ens_ds = xr.concat((ens_ds,model_ds),dim='model',coords='minimal',compat='override') #generating new dimension 'model'    

        ens_models = np.append(ens_models,model) #store model name in metadata
        
    #add model list as coordinate to ds
    ens_ds = ens_ds.assign_coords(model=ens_models)
    ens_ds = ens_ds.to_dataset(name='zos')

    #save to netcdf
    out_fn = variable+'_CMIP6_'+ssp+'_n'+str(len(ens_ds.model))+'_nslocations_all_variants.nc'
        
    #add attributes
    ens_ds.attrs['comments'] = 'all variants, common calendar used for all models'
    ens_ds.attrs['creation_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    ens_ds.attrs['contact'] = 'tim.hermans@nioz.nl'
    ens_ds.attrs['location'] = 'lat: '+str(qlat)+', lon: '+str(qlon) 
  
    ens_ds.to_netcdf(os.path.join(out_dir,out_fn),mode='w',encoding={variable:{'zlib': True,'complevel': 4}}) #(over)write a new NetCDF file
