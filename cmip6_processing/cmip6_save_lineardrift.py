'''
This script stores linearly dedrifted CMIP6 model ocean data (see cmip6_dedrift_linear.py).

Input parameters:
    variable        cmip6 variable to dedrift
    in_dir          input directory with time-merged cmip6 variable data (any frequency), organized by model directory
    out_dir         output directory to store in

Output:
    NetCDF files with linear piControl drift (gridded slopes & intercepts)
    
Created by: Tim Hermans, 05-08-21
'''
import xarray as xr
import numpy as np
import os
import fnmatch

def fit_linear(data, time): #fit polynomial least squares
    pfit = np.polyfit(time, data,1) 
    return np.array([pfit[0],pfit[1]]) #return coefficients, highest power first!

variable = 'zos'
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/' + variable #set raw, time-merged zostoga directory
out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/zos_drift/'

for model in (model for model in os.listdir(in_dir) if model not in ['.DS_Store']):
    print('Current model: '+model)
    model_dir = os.path.join(in_dir,model)
    
    pic_fns = fnmatch.filter(os.listdir(model_dir), "*piControl*nc") #search for piControl runs for current model

    if not pic_fns: #if no piControl available
        print('Model has no piControl')
        continue
    
    #loop over piControl runs for current model
    for pic_fn in pic_fns: 
        
        #open piControl run (decode false, xarray can't handle piC calendar years)
        pic = xr.open_dataset(os.path.join(model_dir,pic_fn),decode_times=False) 
        
        print('Current piControl variant: '+pic.variant_label+', grid: '+pic.grid_label)

        #fit polynomial to piControl
        control_pfit = xr.apply_ufunc(
                fit_linear, pic[variable] , pic.time,
                input_core_dims=[["time"], ["time"]], #core dimension: time, loop over the others
                output_core_dims=[["coefs"]], #outputs 1st degree and intercept
                vectorize=True, 
                dask='allowed', #allow calculating in chunks (dask='parallelized' doesn't work)
                output_dtypes=[float],
                output_sizes={"coefs": 2}, #output must be numpy array
                )
        
        #save pfit
        if not os.path.exists(os.path.join(out_dir,model)): #make path
            os.mkdir(os.path.join(out_dir,model))
        
        control_pfit = control_pfit.to_dataset(name='pfit')
        control_pfit.to_netcdf(os.path.join(out_dir,model,'pfit_linear_'+pic_fn),
                                         mode='w') #(over)write a new NetCDF file