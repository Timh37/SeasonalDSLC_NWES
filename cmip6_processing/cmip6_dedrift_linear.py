'''
This script linearly dedrifts CMIP6 model ocean data using piControl. The drift
is estimated by fitting a linear polynomial to the full piControl run and subtracting
the polynomial from the historical and scenario runs. The output is saved to NetCDFs,
separately for experiment, variant and grid.

Input parameters:
    variable        cmip6 variable to dedrift
    in_dir          input directory with time-merged cmip6 variable data (any frequency), organized by model directory
    out_dir         output directory to store dedrifted cmip6 variable data 

Output:
    NetCDF files with dedrifted timeseries
    
Created by: Tim Hermans, 29-07-20
'''
import xarray as xr
import numpy as np
import os
import fnmatch

def fit_linear(data, time): #fit polynomial least squares
    pfit = np.polyfit(time, data,1) 
    return np.array([pfit[0],pfit[1]]) #return coefficients, highest power first!

def evaluate_drift(pfit, time): #evaluate polynomial
    drift = np.polyval(pfit,time)
    return drift

def make_dedrifted_ds(ds,dedrifted_data,variable): #generate dedrifted dataset from raw dataset
    dedrifted_data.attrs = ds[variable].attrs #copy attributes from original data
    dedrifted_ds = ds.copy(deep=True) #copy original dataset
    
    dedrifted_ds[variable] = dedrifted_data #replace raw with dedrifted data
    dedrifted_ds[variable].attrs['dedrifted'] = 'Dedrifted using linear polynomial fitted to full piControl'
    return dedrifted_ds

def cmip6_dedrift_linear(in_dir,out_dir,variable,skip_existing):
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
             
            #search for historical & ssp filenames with this grid
            hist_fns = fnmatch.filter(os.listdir(model_dir), "*historical*"+pic.grid_label+"*nc")
            ssp_fns = fnmatch.filter(os.listdir(model_dir), "*ssp*"+pic.grid_label+"*nc")
            
            if not hist_fns: #if no historical file, we cannot dedrift, so move on to next piControl run
                print('Current piControl has no historical')
                continue #go to next piControl
                
            if skip_existing:
                ssp_fns_mirror = ssp_fns
                ssp_fns = []
                
                for ssp_fn in ssp_fns_mirror: #generate list of ssp filenames not yet dedrifted
                    if not os.path.exists(os.path.join(out_dir,model,ssp_fn[0:-3]+'_dedrifted.nc')):
                        ssp_fns = np.append(ssp_fns,ssp_fn)
                ssp_fns = list(ssp_fns)

                if not ssp_fns: #if no ssps to dedrift
                    all_hist_exist = 1
                    for hist_fn in hist_fns:
                        if not os.path.exists(os.path.join(out_dir,model,hist_fn[0:-3]+'_dedrifted.nc')):
                            all_hist_exist = 0
                            break
                        
                    if all_hist_exist:
                        print('no new files to dedrift')
                        continue #no files to dedrift, move on to next pic
                        
            #to-do: implement smth to check here if ssp_fns is empty and also historical already done, then we may move on to next pic version

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
            
            print('Fitted polynomial to piControl')
    
            #dedrift historical        
            for hist_fn in hist_fns: #loop over historical files corresponding to current piControl
                dedrift_hist = True
                
                #open historical file
                hist = xr.open_dataset(os.path.join(model_dir,hist_fn),decode_times=False)
                
                if skip_existing: #"beta-version 22-02-21" assumes ssp variant label always same as hist variant label, this is probably always true?
                    if os.path.exists(os.path.join(out_dir,model,hist_fn[0:-3]+'_dedrifted.nc')): #if dedrifted historical file exists
                        dedrift_hist = False #don't have to dedrift again
                        
                        if not fnmatch.filter(ssp_fns,"*"+hist.variant_label+"*"+hist.grid_label+"*"): #if no hist & no new matching ssps to be dedrifted
                            continue #skip file altogether
                    
                #if current historical is  branched from current piControl
                if (((hist.parent_experiment_id=='piControl') or (hist.parent_experiment_id=='p i C o n t r o l')) and (hist.parent_variant_label == pic.variant_label)): 
                    
                    #evaluate polynomial of piControl for historical
                    hist_drift = xr.apply_ufunc(
                        evaluate_drift,control_pfit,hist.time, 
                        input_core_dims=[["coefs"], ["time"]], 
                        output_core_dims=[["time"]], 
                        vectorize=True,
                        dask='allowed',
                        output_dtypes=[hist[variable].dtype],
                        )
                    hist_drift = hist_drift - hist_drift.isel(time=0) #we're interested in subtracting the slope, not absolute values
                    
                    if dedrift_hist:
                        
                        hist_dedrifted = hist[variable] - hist_drift #subtract drift
                        hist_dedrifted_ds = make_dedrifted_ds(hist,hist_dedrifted,variable)

                        if not os.path.exists(os.path.join(out_dir,model)): #make path
                                os.mkdir(os.path.join(out_dir,model))
                
                        #save to netcdf, use overwrite and compress variable
                        hist_dedrifted_ds.to_netcdf(os.path.join(out_dir,model,hist_fn[0:-3]+'_dedrifted.nc'),
                                                     mode='w',
                                                     encoding={variable:{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
                
                        print('Dedrifted '+hist.experiment_id+' ('+hist.variant_label+')')  
                    
                        hist_dedrifted_ds.close()
        
                else: #current historical not branched from current piControl    
                    print('Current ' +hist.experiment_id+' ('+hist.variant_label+') not branched from current piControl')
                    continue #move on to next historical (we cannot dedrift SSPs either in this case)
    
                #dedrift SSPs
                for ssp_fn in ssp_fns:
                    ssp = xr.open_dataset(os.path.join(model_dir,ssp_fn),decode_times=False)
                
                    #error in EC Earth metadata SSP parent variant label (14/7/2021: this occurs for more models, can use ncatted to edit)
                    if (model.startswith('EC-Earth3') and (ssp.variant_label==hist.variant_label)):
                        ssp.attrs['parent_variant_label'] = hist.variant_label
                    #if current ssp is branched from current historical
                    if (((ssp.parent_experiment_id=='historical') or (ssp.parent_experiment_id=='h i s t o r i c a l')) and (ssp.parent_variant_label == hist.variant_label)):
                            
                        #evaluate quadratic polynomial of piControl for ssp time
                        ssp_drift = xr.apply_ufunc(  
                            evaluate_drift,control_pfit,ssp.time,
                            input_core_dims=[["coefs"], ["time"]], 
                            output_core_dims=[["time"]], 
                            vectorize=True,
                            dask='allowed',
                            output_dtypes=[ssp[variable].dtype],
                            )
                        #replace first value of drift with last value of historical drift & add 1 month of drift (this is necessary because ssp time might have a different reference date from hist time)
                        ssp_drift = ssp_drift - ssp_drift.isel(time=0) + hist_drift.isel(time=-1) + ( ssp_drift.isel(time=1)-ssp_drift.isel(time=0))
                        ssp_dedrifted = ssp[variable] - ssp_drift
        
                        ssp_dedrifted_ds = make_dedrifted_ds(ssp,ssp_dedrifted,variable)
        
                        if not os.path.exists(os.path.join(out_dir,model)): #make path
                                os.mkdir(os.path.join(out_dir,model))
                        
                        #save to netcdf, use overwrite and compress variable
                        ssp_dedrifted_ds.to_netcdf(os.path.join(out_dir,model,ssp_fn[0:-3]+'_dedrifted.nc'),
                                                     mode='w',
                                                     encoding={variable:{'zlib': True,'complevel': 1}}
                                                     ) #(over)write a new NetCDF file
                       
                        print('Dedrifted '+ssp.experiment_id+' ('+ssp.variant_label+')')  
                        
                        ssp_dedrifted_ds.close()
                    
                    else:
                        print('Current ' +ssp.experiment_id+' ('+ssp.variant_label+') not branched from current historical, but from '+ssp.parent_experiment_id+' ('+ssp.parent_variant_label+')')
                        continue
                        
                    ssp.close()
                
            hist.close()
        pic.close()
    
#example
if __name__ == '__main__':
    variable = 'zos'
    in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/' + variable #set raw, time-merged zostoga directory
    out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/' + variable #set output directory
    
    skip_existing  = True
    
    cmip6_dedrift_linear(in_dir,out_dir,variable,skip_existing)