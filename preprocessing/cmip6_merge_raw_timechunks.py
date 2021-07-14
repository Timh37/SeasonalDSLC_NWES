import xarray as xr
import numpy as np
import os
import fnmatch
from cdo import Cdo
'''
This script merges time chunks of raw ESGF data to complete time series and 
saves output to NetCDFs, separately for experiment, variant and grid.

Note: this script assumes input is in raw monthly means with ESGF file naming
conventions.

Note: could consider filling data gaps using interpolation rather than skipping the data

Input parameters:
    in_dir          input directory with monthly mean raw cmip6 variable data, organized by model directory
    out_dir         output directory to store time-merged monthly mean cmip6 variable data 
    variable        cmip6 variable for which to concatenate time chunks

Output:
    NetCDF files with concatenated timeseries
    
Created by: Tim Hermans, 10-07-20
'''
def cmip6_merge_raw_timechunks(in_dir,out_dir,variable,skip_existing):
    #loop over models
    cdo = Cdo()
    for model in (model for model in os.listdir(in_dir) if model!=".DS_Store"):
        print('Current model: '+model)
    
        model_dir = os.path.join(in_dir,model) #model input directory
    
        for experiment in ['piControl','historical','ssp']: #loop over experiment types
            runs=list(set([i.split('_')[3] for i in fnmatch.filter(os.listdir(model_dir), "*"+experiment+"*.nc")])) #different runs for this experiment (e.g., ssp126 for ssp)
        
            for run in runs:
                variants=list(set([i.split('_')[4] for i in fnmatch.filter(os.listdir(model_dir), "*"+run+"*.nc")])) #different variants for this run (e.g., r[n]i[n]p[n]f[n])
                
                for variant in variants:    
                    grids=list(set([i.split('_')[5] for i in fnmatch.filter(os.listdir(model_dir), "*"+run+"*"+variant+"*.nc")])) #different grids for this variant (e.g., 'gm')
                    
                    for grid in grids:
                        try: #open multifile dataset for this grid
                            concat_ds = xr.open_mfdataset(model_dir+"/*mon*"+run+"*"+variant+"*"+grid+"*.nc",combine='by_coords',decode_times=False)
                        except:
                            print('Cannot open files or find common time dimension for '+model+', '+run+', '+variant+', '+grid)
                            continue
                        
                        if concat_ds.frequency != 'mon': #check if raw data is monthly means
                            print('no monthly means input for '+model+', '+run+', '+variant+', '+grid)
                            continue
                        elif max(np.diff(concat_ds.time))>31: #check if no gaps larger than 31 days (1 month)
                            print('data gaps for '+model+', '+run+', '+variant+', '+grid)
                            #possibly add interpolaton methods for filling missing data
                            continue
                        else: #if data can be concatenated
                            filenames = sorted(fnmatch.filter(os.listdir(model_dir), "*mon*"+run+"*"+variant+'*'+grid+"*.nc")) #get filenames
                            
                            fileyears = [fn[-16:-3] for fn in filenames] #extract years for new filename
                            
                            #save concatenated time series to netcdf
                            if not os.path.exists(out_dir): #generate output directory
                                os.mkdir(out_dir)
                                
                            if not os.path.exists(os.path.join(out_dir,model)): #generate output model directory
                                os.mkdir(os.path.join(out_dir,model))
        
                            concat_ds.close() #close multifile dataset
                            
                            if skip_existing: #if files already concatenated before, skip
                                if os.path.exists(os.path.join(out_dir,model,filenames[0][0:-16]+fileyears[0][0:6]+'-'+fileyears[-1][7::]+'.nc')):
                                    continue
                                
                            try: #concatenate using CDO (fast & easy) & compress
                                cdo.cat(input=model_dir+'/*mon*'+run+'*'+variant+'*'+grid+'*.nc',
                                        output=os.path.join(out_dir,model,filenames[0][0:-16]+fileyears[0][0:6]+'-'+fileyears[-1][7::]+'.nc'),
                                        options='-z zip_1')
                                print('saved: '+filenames[0][0:-16]+fileyears[0][0:6]+'-'+fileyears[-1][7::]+'.nc')
                            except:
                                print('could not save: '+filenames[0][0:-16]+fileyears[0][0:6]+'-'+fileyears[-1][7::])
                            
                            
#example
if __name__ == '__main__':
    variable = 'zos'
    in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/raw/' + variable #+'/'+model+'/'
    out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/'+variable
    skip_existing = True
    
    cmip6_merge_raw_timechunks(in_dir,out_dir,variable,skip_existing)