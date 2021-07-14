import numpy as np
import os
import iris
from RegridCMIP6 import regrid

'''
Regrids CMIP6 model variable fields to common 1x1 grid using bilinear interpolation
with ESMValTool (iris internal regridding for regular grids and ESMF for irregular)
https://doi.org/10.5194/gmd-13-3383-2020

Note: no extrapolation where source is land, land mask is larger than on native grid and varies per model on the common grid,
'extrap_method' argument to ESMF currently not implemented in ESMValTool 
some weird wrapping error around poles, differs per model where exactly-maybe solve later

Input parameters:
    variable    ocean variable of interest
    in_dir      path to NetCDF files with ocean variable, ordered by model
    out_dir     path to store the regridded files

Output:
    regridded files
    
GG;    
Note: In order to run this code, activate the "esmpy" conda environment
conda activate esmpy

Created by: Tim Hermans, 02-11-20
'''

#input
variable = 'zos'
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/zos_0mean/' #+'/'+
out_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/zos_0mean_1x1/' #+'/'+
skip_existing = True

for i,model in enumerate((model for model in os.listdir(in_dir) if model not in [".DS_Store"])): #loop over models
    model_dir = os.path.join(in_dir,model) #model input directory
    print(model)
    
    for file in sorted((file for file in os.listdir(model_dir) if file!='.DS_Store')): #loop over experiments for model
    
        if skip_existing: #if not re-regridding already regridded files
            if os.path.exists(os.path.join(out_dir,model,file[0:-3]+'_1x1.nc')):
                continue #if regridded file already exists, move on to next
        print(file)        
        
       	# Load in a CF-compliant netCDF file (i.e. from CMIP6 archive)
       	# Note: iris.load returns a list of cubes. There's should be only one
        # cube in this file, so extract that cube from the list.
       	cubes = iris.load(os.path.join(model_dir,file))
        if len(cubes)>1: #file contains multiple cubes, we want the zos cube
            var_names = [cube.var_name for cube in cubes]   
            cube=cubes[var_names.index(variable)]
        else:
            cube=cubes[0]
           
        if len(cube.coords(axis='x'))==len(cube.coords(axis='y'))==2: #if dealing with duplicate coordinate entries (BCC models for example)
            #we need coordinate with either bounds or circularity
            for k in np.arange(0,2):
                if cube.coords(axis='x')[k].bounds is None:
                    
                    cube.remove_coord(cube.coords(axis='x')[k].standard_name) #remove that coordinate
                    cube.remove_coord(cube.coords(axis='y')[k].standard_name)
                    cube.coord(axis='x').rename('longitude') #rename remaining coordinates to standard
                    cube.coord(axis='y').rename('latitude')
                    break
        
        if ((cube.coord('longitude').bounds is None) & (isinstance(cube.coord('longitude'),iris.coords.DimCoord)==False)): #if no bounds, cannot determine if circular
            print('Cannot determine circularity of longitudes without bounds for: '+file+', moving on')
            continue
       	
       	# Run the regridding method on the loaded data
       	# Usage: regrid(cube, resolution (degree x degree), method, lat_offset, lon_offset)
       	new_cube = regrid(cube, "1x1", "linear", False, False)
   	
        if not os.path.exists(os.path.join(out_dir,model)): #make path
            os.mkdir(os.path.join(out_dir,model))
           
        # Write the regridded data to a netCDF file
        iris.save(new_cube,os.path.join(out_dir,model,file[0:-3]+'_1x1.nc'), zlib=True, complevel=4, fill_value=np.nan)