'''
Subtracts the area-weighted mean of an ocean variable from each grid cell at each timestep

to-do: implement functionality to calculate area weights for models without 'areacello'

Input parameters:
    variable    ocean variable of interest
    in_dir      path to NetCDF files with ocean variable, ordered by model
    out_dir     path to store the corrected files
    area_dir    path to areacello (ocean grid cell area)

Output:
    awmean   area weighted mean of input variable
    
Created by: Tim Hermans, 29-07-20
'''
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
import fnmatch

#input 
variable = 'zos'
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/' + variable 
out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/zos_0mean/'
area_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/raw/areacello/' #areaweights
skip_existing = True

for i,model in enumerate((model for model in os.listdir(in_dir) if model not in [".DS_Store"])): #loop over models
    model_dir = os.path.join(in_dir,model) #model input directory
    print(model)
    
    if model not in os.listdir(area_dir): #if no directory for model in path to areaweights
        hasAreaweights = False
        print(model+' has no areaweights, moving on')
        continue
    else:
        hasAreaweights = True
        
    for file in (file for file in sorted(os.listdir(model_dir)) if file!='.DS_Store'): #loop over experiments for model
        if skip_existing: #if not processing files already corrected
            if os.path.exists(os.path.join(out_dir,model,file[0:-3]+'_0mean.nc')): #if corrected file already exists
                continue #move on to next
        
        hasAreaweights = True
        
        model_area_dir = os.path.join(area_dir,model) #path to areaweights for model
        areafile=[]
 
        var_ds = xr.open_dataset(os.path.join(model_dir,file),decode_times=False) #load variable file
        coordinates = list(k for k in var_ds[variable].dims if 'time' not in k) #find lon/lat coordinate names
        
        areafile = fnmatch.filter(os.listdir(model_area_dir),"*"+var_ds.grid_label+"*.nc") #look for matching areacello (same grid label)
        
        if not areafile: #if none available
            hasAreaweights=False
            print(model+', '+var_ds.grid_label+' has no areaweights, moving on')
            continue
            
        else:
            areaweights = xr.open_dataset(os.path.join(model_area_dir,areafile[0])) #load first available areaweights with matching grid
            areaweights['areacello'] = areaweights.areacello.where(areaweights.areacello<1e20,np.nan) #some models have 1e35 as missing values, we replace these by nan
         
        if (np.sum(areaweights.areacello)>3.8e14): #sanity check for total ocean surface area [m2]
            #if too large ocean area, land mask may not be contained in areacello, instead extract mask from variable data
            areaweights['areacello'] = areaweights.areacello.values*np.isfinite(var_ds[variable].isel(time=0)) 
        
        if ((np.sum(areaweights.areacello)<3.5e14) or (np.sum(areaweights.areacello)>3.8e14)): #if still invalid total ocean surface, skip file
            print('invalid total ocean area for: ' + file)
            continue
        
        #calculate and subtract areaweighted mean
        awmean = ((var_ds[variable] * areaweights.areacello).sum(dim=coordinates,skipna=True)) / np.sum(areaweights.areacello)
        with xr.set_options(keep_attrs=True):
            var_ds[variable] = var_ds[variable] - awmean #subtract awmean from original variable
        
        #save corrected variable file to netcdf
        if not os.path.exists(os.path.join(out_dir,model)): #make path
            os.mkdir(os.path.join(out_dir,model))
        
        with xr.set_options(keep_attrs=True):    
            var_ds.to_netcdf(os.path.join(out_dir,model,file[0:-3]+'_0mean.nc'),
                                                         mode='w',
                                                         encoding={variable:{'zlib': True,'complevel': 1,'dtype':'float32'}}
                                                         ) #(over)write a new NetCDF file