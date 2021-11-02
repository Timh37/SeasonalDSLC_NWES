'''
compute grid resolution using areacell(a/o) with dedicated CMIP6 script:

uses largest great-circle distance between vertices of each grid cell,
then computes the area-weighted mean of these

Nominal resolutions CMIP6: https://github.com/WCRP-CMIP/CMIP6_CVs/blob/master/CMIP6_source_id.json
Documentation routines: https://docs.google.com/document/d/1h0r8RZr_f3-8egBMMh7aqLwy3snpD6_MrDz1q8n5XUk/edit#
PCMDI code: https://github.com/PCMDI/nominal_resolution/blob/master/lib/api.py

Loops over CMIP6 models incorporated in the ensemble, assuming folder-structure
organized by variable by model

Created by: Tim Hermans, 13-10-21
'''
import xarray as xr
import numpy
import numpy as np
import os
import fnmatch
import warnings 

### user input
realm = 'atmosphere'
target_area='nwes' #'NWES' or 'global'
model_list = list(np.genfromtxt('/Users/thermans/Documents/PhD/Phase4_seasonal/Analysis/compiling_ensembles/ens_model_list_'+'ssp585'+'.txt',dtype='str'))


#routines for computing max distances and nominal resolution from CMIP6
try:
    import cdms2  # noqa
    has_cdms = True
except BaseException:
    has_cdms = False
    
def degrees2radians(data, radianLimit):
    """Converts data from deg to radians
    If no data is outside the range -radianLimit to + radianLimit, assume data
    is already in radians and do nothing
    :Example:
        .. doctest::
            >>> data = nominal_resolution.degrees2radians(data)
    :param data: array to potentially convert from degrees to radians
    :type data: cdms2.tvariable.TransientVariable
    :param radianLimit: If the magnitude of any data value exceeds radianLimit,
                        convert data from deg to rad.
    :type radianLimit: `float`_
    :return: The same data converted from degrees to radians
    :rtype: `cdms2.tvariable.TransientVariable`_
    """
    if numpy.ma.absolute(data).max() > radianLimit:
        data = (data / 180.) * numpy.pi
    else:
        warnings.warn("It appears your data is already in radians," +
                      " nothing done by function degrees2radians")
    return data


def mean_resolution(cellarea, latitude_bounds, longitude_bounds,
                    convertdeg2rad, returnMaxDistance=False):
    """Computes mean nominal resolution
    formula from: https://en.wikipedia.org/wiki/Great-circle_distance
    :Example:
        .. doctest::
            >>> mean = nominal_resolution.mean_resolution(cellarea, latitude_bounds, longitude_bounds,
                                                          convertdeg2rad, returnMaxDistance=False)
    :param cellarea: simple array, python masked array, or cdms2 variable containing area of each cell
    :type cellarea: `cdms2.tvariable.TransientVariable`_
    :param latitude_bounds: 2D numpy-like array containing latitudes vertices (ncell, nvertices). If None is passed
                            and cdms2 is available, will try to obtain from cellarea grid
    :type latitude_bounds: `numpy.ndarray`_
    :param longitude_bounds: 2D numpy-like array containing longitudes vertices (ncell, nvertices). If None is passed
                            and cdms2 is available, will try to obtain from cellarea grid
    :type longitude_bounds: `numpy.ndarray`_
    :param convertdeg2rad: set to True if lat/lon in degrees; set to false if in radians
    :type convertdeg2rad: `bool`_
    :param returnMaxDistance: Returns an array representing the maximum distance (in km) between vertices for each cell
    :type returnMaxDistance: `bool`_
    :return: the mean nominal resolution in km and optionally the maximum distance array
    :rtype: `float [,cdms2.tvariable.TransientVariable]`_
    """

    if longitude_bounds is None and latitude_bounds is None and has_cdms:
        try:  # Ok maybe we can get this info from cellarea data
            mesh = cellarea.getGrid().getMesh()
            latitude_bounds = mesh[:, 0]
            longitude_bounds = mesh[:, 1]
            warnings.warn(
                "You did not pass lat/lon bounds but we inferred them from cellarea")
        except BaseException:
            pass

    if longitude_bounds is None or latitude_bounds is None:
        raise RuntimeError(
            "You did not pass lat/lon bounds and couldn't infer them from cellarea")

    if convertdeg2rad:
        radianLimit = 0.
        latitude_bounds = degrees2radians(latitude_bounds, radianLimit)
        longitude_bounds = degrees2radians(longitude_bounds, radianLimit)

    # distance between successive corners
    nverts = latitude_bounds.shape[-1]
    sh = list(latitude_bounds.shape[:-1])
    del_lats = numpy.ma.zeros(sh, dtype=numpy.float)
    del_lons = numpy.ma.zeros(sh, dtype=numpy.float)
    max_distance = numpy.ma.zeros(sh, dtype=numpy.float)

    for i in range(nverts - 1):
        for j in range(i + 1, nverts):
            del_lats = numpy.ma.absolute(
                latitude_bounds[:, i] - latitude_bounds[:, j])
            del_lons = numpy.ma.absolute(
                longitude_bounds[:, i] - longitude_bounds[:, j])
            del_lons = numpy.ma.where(
                numpy.ma.greater(
                    del_lons, numpy.pi), numpy.ma.absolute(
                    del_lons - 2. * numpy.pi), del_lons)

                        # formula from: https://en.wikipedia.org/wiki/Great-circle_distance
            distance = 2. * numpy.ma.arcsin(numpy.ma.sqrt(numpy.ma.sin(del_lats / 2.)**2 + numpy.ma.cos(
                latitude_bounds[:, i]) * numpy.ma.cos(latitude_bounds[:, j]) * numpy.ma.sin(del_lons / 2.)**2))

            max_distance = numpy.ma.maximum(max_distance, distance.filled(0.0))

    radius = 6371.  # in km
    accumulation = numpy.ma.sum(
        cellarea * max_distance) * radius / numpy.ma.sum(cellarea)
    #accumulation = numpy.nansum(
    #    cellarea * max_distance) * radius / numpy.nansum(cellarea)
    if returnMaxDistance:
        return accumulation, max_distance * radius
    else:
        return accumulation    


if realm == 'ocean': #set paths and variable names depending on realm
    areavar = 'areacello'
    area_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/raw/areacello/'
    
    datavar = 'zos' #used to determine grid label & land mask
    data_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/zos_0mean/'
elif realm == 'atmosphere':
    areavar = 'areacella'
    area_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/raw/areacella/'
    
    datavar = 'tauu' #used to determine grid label
    data_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/tauu/'



mres_native=np.array([])
for i,model in enumerate(model_list): #loop over models

    model_data_dir = os.path.join(data_dir,model) #model input directory
    data_file = fnmatch.filter(os.listdir(model_data_dir), "*ssp585*nc")[0] #grab arbitrary file
    data_ds = xr.open_dataset(os.path.join(model_data_dir,data_file),decode_times=False) #load variable file
    
    print(model+', '+data_ds.grid_label) #print current model & grid label
    
    model_area_dir = os.path.join(area_dir,model) #model input directory
    try:
        area_file = fnmatch.filter(os.listdir(model_area_dir), "*"+data_ds.grid_label+"*nc")[0] #first area file available with the right grid label
    except:
        print('no area available, skipping model')
        continue
    
    area_ds = xr.open_dataset(os.path.join(model_area_dir,area_file),decode_times=False) #load areavar file
    
    lat_name = list(c for c in area_ds.coords if 'lat' in c)[0] #fetch lon and lat variable names
    lon_name = list(c for c in area_ds.coords if 'lon' in c)[0]
    
    try: #find lon/lat bounds
        lat_bnds_name = list(c for c in area_ds.variables if ((('bnds' in c) or ('bounds' in c) or ('vert' in c)) and ('lat' in c)))[0]
        lon_bnds_name = list(c for c in area_ds.variables if ((('bnds' in c) or ('bounds' in c) or ('vert' in c)) and ('lon' in c)))[0]
        bounds_available = True
    except:
        bounds_available = False

    if ((target_area == 'NWES') or (target_area == 'nwes')): #mask out cell areas outside the region of interest
        if np.any(area_ds[lon_name]>180): #lon defined from 0 to 360 instead of -180 to 180, so area to constrain between 0 and 12 and 360-15 to 360
            cellarea = area_ds[areavar].where((area_ds[lat_name]>=45) & (area_ds[lat_name]<=62) & ~((area_ds[lon_name]>10) & (area_ds[lon_name]<345)) ).values
        else:
            cellarea = area_ds[areavar].where((area_ds[lat_name]>=45) & (area_ds[lat_name]<=62) & (area_ds[lon_name]>=-15) & (area_ds[lon_name]<=10) ).values
            
    elif target_area == 'global':    
        cellarea=area_ds[areavar].values

    if realm == 'ocean':
        cellarea[cellarea>1e20] = np.nan #some models use large values as mask
        
        if (np.nansum(area_ds[areavar])>3.8e14): #sanity check for total ocean surface area [m2]
            print('replacing land mask')
            #if too large ocean area, land mask may not be contained in areacello, instead extract mask from variable data:
            cellarea[~np.isfinite(data_ds[datavar].isel(time=0))] = np.nan
    
    #generate 1D masked array of cell areas
    cellarea=numpy.ma.array(cellarea.flatten(),mask=~np.isfinite(cellarea.flatten()))
        
    if not bounds_available: #if no bounds, grid may be regular, and we can create our own bounds
        if ((len(area_ds[lat_name].shape)==1) & (len(area_ds[lon_name].shape)==1)): #regular grid
            dlat = np.mean(np.diff(area_ds[lat_name]))    
            dlon = np.mean(np.diff(area_ds[lon_name]))
            
            blats = np.column_stack([area_ds[lat_name].values,area_ds[lat_name].values])
            blats[:,0] = blats[:,0] - dlat/2
            blats[:,1] = blats[:,1] + dlat/2
            
            blons = np.column_stack([area_ds[lon_name].values,area_ds[lon_name].values])
            blons[:,0] = blons[:,0] - dlon/2
            blons[:,1] = blons[:,1] + dlon/2
            
        else:
            print('bounds not available and no regular grid, cannot compute grid resolution')
            continue
    else:
        #if bounds available, use CMIP6 routine to compute the area-weighted mean resolution (from maximum great-circle distance between vertices)   
        if len(area_ds[lon_bnds_name].shape)>2: #reduce to 1 dimension
            blats = area_ds[lat_bnds_name].values.reshape(area_ds[lat_bnds_name].shape[0]*area_ds[lat_bnds_name].shape[1],area_ds[lat_bnds_name].shape[-1])
            blons = area_ds[lon_bnds_name].values.reshape(area_ds[lon_bnds_name].shape[0]*area_ds[lon_bnds_name].shape[1],area_ds[lon_bnds_name].shape[-1])
        elif len(area_ds[lon_bnds_name].shape)==2:
            blats = area_ds[lat_bnds_name].values
            blons = area_ds[lon_bnds_name].values
            
    if blats.shape[-1]==2: #regular grid, bounds only given at sides of each coordinate? (BCC) -> necessary to repeat for mesh
        blats = np.repeat(blats,len(area_ds[lon_name]),axis=0)
        blons = np.repeat(blons,len(area_ds[lat_name]),axis=0)
    
    mresolution = mean_resolution(cellarea,blats,blons,convertdeg2rad=True,returnMaxDistance=False ) #run the CMIP6 routine
    
    if area_ds.grid_label=='gn':
        mres_native = np.append(mres_native,mresolution)
    print(mresolution)
print('mean native grid resolution:')
print(np.mean(mres_native))
    
