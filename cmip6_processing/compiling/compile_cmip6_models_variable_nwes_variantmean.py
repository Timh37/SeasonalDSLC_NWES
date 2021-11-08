#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 08:43:46 2021

Compile CMIP6 simulations on a common grid within the NWES region into ensemble,
combining  historical and SSP runs and averaging across variants per model.

@author: thermans
"""

import xarray as xr
import numpy as np
import os
import fnmatch
import datetime

## user input
target_ssps = ['ssp126']#,'ssp245','ssp370','ssp585'] #desired ssps

variable = 'tauv'
var_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/tauv_1x1/'
out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/'
##

for ssp in target_ssps: #loop over desired ssps
    print(ssp) #debugging
    
    ens_models = np.array([]) #model names
    num_variants = np.array([]) #number of variants averaged over for each model
    
    first_model=True
    
    for model in (model for model in os.listdir(var_dir) if model not in ['.DS_Store']):
        print("current model "+model) #debugging
        var_model_dir = os.path.join(var_dir,model) #variable directory

        if not fnmatch.filter(os.listdir(var_model_dir), "*mon*"+ssp+"*nc"): #look for SSP
            print(variable+" for "+ssp+" not available for current model") #if variable not available for current ssp
            continue
                
        num_variants = np.append(num_variants,len(fnmatch.filter(os.listdir(var_model_dir), "*mon*"+ssp+"*nc"))) #number of variants used for each model
        
        for f,fn in enumerate(fnmatch.filter(os.listdir(var_model_dir), "*"+ssp+"*nc")):
            #load ssp & historical
            ssp_ds = xr.open_dataset(os.path.join(var_model_dir,fn))
            hist_ds = xr.open_dataset(os.path.join(var_model_dir,fnmatch.filter(os.listdir(var_model_dir), "*historical*"+ssp_ds.variant_label+"*nc")[0]))
            
            if variable.startswith('tau'):
                if model=='EC-Earth3':
                    if hist_ds.variant_label=='r11i1p1f1' or hist_ds.variant_label == 'r13i1p1f1': #somehow starts in 184912, generalize this check later
                        hist_ds = hist_ds.isel(time=np.arange(1,1981))
                
            print(ssp_ds.variant_label)
            
            #replace model specific time with common time
            hist_ds['time'] = xr.cftime_range(start="1850", periods=len(hist_ds.time), freq="MS", calendar="noleap") #replace model specific time array by common
            ssp_ds['time'] = xr.cftime_range(start="2015", periods=len(ssp_ds.time), freq="MS", calendar="noleap")
            
            #wrap longitude
            hist_ds.coords['lon'] = ((hist_ds.coords['lon'] + 180) % 360) - 180
            hist_ds = hist_ds.reindex({ 'lon' : np.sort(hist_ds['lon'])})
            
            ssp_ds.coords['lon'] = ((ssp_ds.coords['lon'] + 180) % 360) - 180
            ssp_ds = ssp_ds.reindex({ 'lon' : np.sort(ssp_ds['lon'])})

            #delimit NWES region
            hist_ds = hist_ds.sel(lat=np.arange(40,66),lon=np.arange(-20,15))
            ssp_ds = ssp_ds.sel(lat=np.arange(40,66),lon=np.arange(-20,15))
            
            #concatenate historical & ssp over time dimension
            hist_ssp_ds = xr.concat((hist_ds,ssp_ds),dim='time')
            hist_ssp_ds = xr.decode_cf(hist_ssp_ds)
            
            if f==0: #if first variant, produce datarray to hold variant data
                model_ds = []
                model_ds = hist_ssp_ds
            else:
                model_ds = xr.concat((model_ds,hist_ssp_ds),dim='variant') #append variant to existing and generating new dimension 'model'    
                
        if f>0: #if more than 1 variant, take mean across variants
            variant_mean = model_ds.mean(dim='variant') 
            variant_mean['time_bnds'] = model_ds.time_bnds.isel(variant=0)
        else:
            variant_mean = model_ds #else use single variant
            
        
        if first_model: #concatenate variant means along model dimension
            ens_ds=[]
            ens_ds = variant_mean #initialise ensemble ds with model ds
            first_model=False
        else: #append next models
            ens_ds = xr.concat((ens_ds,variant_mean),dim='model') #generating new dimension 'model'    
            
        ens_models = np.append(ens_models,model) #store model name in metadata
        
    #add model list as coordinate to ds
    ens_ds = ens_ds.assign_coords(model=ens_models)

    #save to netcdf
    out_fn = variable+'_CMIP6_'+ssp+'_n'+str(len(ens_ds.model))+'_nwes_variant_averaged.nc'
        
    #add attributes
    ens_ds['num_variants']=(('model'),num_variants)
    ens_ds = ens_ds.drop_vars(['time_bnds'])
    ens_ds.attrs['comments'] = 'variant mean, common calendar used for all models'
    ens_ds.attrs['creation_date'] = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    ens_ds.attrs['contact'] = 'tim.hermans@nioz.nl'
    
  
    ens_ds.to_netcdf(os.path.join(out_dir,out_fn),mode='w',encoding={variable:{'zlib': True,'complevel': 4}}) #(over)write a new NetCDF file
