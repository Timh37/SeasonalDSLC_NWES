import numpy as np
import os
import re
import iris
from iris.analysis import AreaWeighted, Linear, Nearest, UnstructuredNearest
from RegridESMpy import ESMF_REGRID_METHODS
from RegridESMpy import regrid as esmpy_regrid


#TH: see https://github.com/ESMValGroup/ESMValCore/blob/master/esmvalcore/preprocessor/_regrid.py

# Regular expression to parse a "MxN" cell-specification.
_CELL_SPEC = re.compile(
    r'''\A
        \s*(?P<dlon>\d+(\.\d+)?)\s*
        x
        \s*(?P<dlat>\d+(\.\d+)?)\s*
        \Z
     ''', re.IGNORECASE | re.VERBOSE)

# Default fill-value.
_MDI = 1e+20

# Stock cube - global grid extents (degrees).
_LAT_MIN = -90.0
_LAT_MAX = 90.0
_LAT_RANGE = _LAT_MAX - _LAT_MIN
_LON_MIN = 0.0
_LON_MAX = 360.0
_LON_RANGE = _LON_MAX - _LON_MIN


# A cached stock of standard horizontal target grids.
_CACHE = dict()


# Supported horizontal regridding schemes.
HORIZONTAL_SCHEMES = {
    'linear': Linear(extrapolation_mode='mask'),
    'linear_extrapolate': Linear(extrapolation_mode='extrapolate'),
    'nearest': Nearest(extrapolation_mode='mask'),
    'area_weighted': AreaWeighted(),
    'unstructured_nearest': UnstructuredNearest(),
}


def _attempt_irregular_regridding(cube, scheme):
    """Check if irregular regridding with ESMF should be used."""
    if scheme in ESMF_REGRID_METHODS:
        try:
            lat_dim = cube.coord('latitude').ndim
            lon_dim = cube.coord('longitude').ndim
            if lat_dim == lon_dim == 2:
                return True
        except iris.exceptions.CoordinateNotFoundError:
            pass
    return False



def regrid(cube, target_grid, scheme, lat_offset=True, lon_offset=True):
    """
    Perform horizontal regridding.
    Parameters
    ----------
    cube : cube
        The source cube to be regridded.
    target_grid : cube or str
        The cube that specifies the target or reference grid for the regridding
        operation. Alternatively, a string cell specification may be provided,
        of the form 'MxN', which specifies the extent of the cell, longitude by
        latitude (degrees) for a global, regular target grid.
    scheme : str
        The regridding scheme to perform, choose from
        'linear',
        'linear_extrapolate',
        'nearest',
        'area_weighted',
        'unstructured_nearest'.
    lat_offset : bool
        Offset the grid centers of the latitude coordinate w.r.t. the
        pole by half a grid step. This argument is ignored if `target_grid`
        is a cube or file.
    lon_offset : bool
        Offset the grid centers of the longitude coordinate w.r.t. Greenwich
        meridian by half a grid step.
        This argument is ignored if `target_grid` is a cube or file.
    Returns
    -------
    cube
    See Also
    --------
    extract_levels : Perform vertical regridding.
    """
    if HORIZONTAL_SCHEMES.get(scheme.lower()) is None:
        emsg = 'Unknown regridding scheme, got {!r}.'
        raise ValueError(emsg.format(scheme))

    if isinstance(target_grid, str):
        if os.path.isfile(target_grid):
            target_grid = iris.load_cube(target_grid)
        else:
            # Generate a target grid from the provided cell-specification,
            # and cache the resulting stock cube for later use.
            target_grid = _CACHE.setdefault(
                target_grid,
                _stock_cube(target_grid, lat_offset, lon_offset),
            )
            # Align the target grid coordinate system to the source
            # coordinate system.
            src_cs = cube.coord_system()
            xcoord = target_grid.coord(axis='x', dim_coords=True)
            ycoord = target_grid.coord(axis='y', dim_coords=True)
            xcoord.coord_system = src_cs
            ycoord.coord_system = src_cs

    if not isinstance(target_grid, iris.cube.Cube):
        raise ValueError('Expecting a cube, got {}.'.format(target_grid))

    # Unstructured regridding requires x2 2d spatial coordinates,
    # so ensure to purge any 1d native spatial dimension coordinates
    # for the regridder.
    if scheme == 'unstructured_nearest':
        for axis in ['x', 'y']:
            coords = cube.coords(axis=axis, dim_coords=True)
            if coords:
                [coord] = coords
                cube.remove_coord(coord)

    # Perform the horizontal regridding.
    
    if _attempt_irregular_regridding(cube, scheme): #fails for AWI model
        cube = esmpy_regrid(cube, target_grid, scheme) #somehow fails for atmospheric variables
        #print("Irregular grid found")
    else:
        cube = cube.regrid(target_grid, HORIZONTAL_SCHEMES[scheme])
    
    ##we'll just use esmpy for all #TH: 02-11-2020
    #cube = esmpy_regrid(cube, target_grid, scheme)
    return cube
    



def parse_cell_spec(spec):
    """
    Parse an MxN cell specification string.
    Parameters
    ----------
    spec: str
    Returns
    -------
    tuple
        tuple of (float, float) of parsed (lon, lat)
    Raises
    ------
    ValueError
        if the MxN cell specification is malformed.
    ValueError
        invalid longitude and latitude delta in cell specification.
    """
    cell_match = _CELL_SPEC.match(spec)
    if cell_match is None:
        emsg = 'Invalid MxN cell specification for grid, got {!r}.'
        raise ValueError(emsg.format(spec))

    cell_group = cell_match.groupdict()
    dlon = float(cell_group['dlon'])
    dlat = float(cell_group['dlat'])

    if (np.trunc(_LON_RANGE / dlon) * dlon) != _LON_RANGE:
        emsg = ('Invalid longitude delta in MxN cell specification '
                'for grid, got {!r}.')
        raise ValueError(emsg.format(dlon))

    if (np.trunc(_LAT_RANGE / dlat) * dlat) != _LAT_RANGE:
        emsg = ('Invalid latitude delta in MxN cell specification '
                'for grid, got {!r}.')
        raise ValueError(emsg.format(dlat))

    return dlon, dlat



def _stock_cube(spec, lat_offset=True, lon_offset=True):
	"""
	Create a stock cube.
	Create a global cube with M degree-east by N degree-north regular grid
	cells.
	The longitude range is from 0 to 360 degrees. The latitude range is from
	-90 to 90 degrees. Each cell grid point is calculated as the mid-point of
	the associated MxN cell.
	Parameters
	----------
	spec : str
		Specifies the 'MxN' degree cell-specification for the global grid.
	lat_offset : bool
		Offset the grid centers of the latitude coordinate w.r.t. the
		pole by half a grid step. This argument is ignored if `target_grid`
		is a cube or file.
	lon_offset : bool
		Offset the grid centers of the longitude coordinate w.r.t. Greenwich
		meridian by half a grid step.
		This argument is ignored if `target_grid` is a cube or file.
	Returns
	-------
		A :class:`~iris.cube.Cube`.
	"""
	dlon, dlat = parse_cell_spec(spec)
	mid_dlon, mid_dlat = dlon / 2, dlat / 2

	# Construct the latitude coordinate, with bounds.
	if lat_offset:
		latdata = np.linspace(_LAT_MIN + mid_dlat, _LAT_MAX - mid_dlat,
							  int(_LAT_RANGE / dlat))
	else:
		latdata = np.linspace(_LAT_MIN, _LAT_MAX, int(_LAT_RANGE / dlat) + 1)

	# Construct the longitude coordinat, with bounds.
	if lon_offset:
		londata = np.linspace(_LON_MIN + mid_dlon, _LON_MAX - mid_dlon,
							  int(_LON_RANGE / dlon))
	else:
		londata = np.linspace(_LON_MIN, _LON_MAX - dlon,
							  int(_LON_RANGE / dlon))

	lats = iris.coords.DimCoord(
		latdata,
		standard_name='latitude',
		units='degrees_north',
		var_name='lat')
	lats.guess_bounds()

	lons = iris.coords.DimCoord(
		londata,
		standard_name='longitude',
		units='degrees_east',
		var_name='lon')
	lons.guess_bounds()

	# Construct the resultant stock cube, with dummy data.
	shape = (latdata.size, londata.size)
	dummy = np.empty(shape, dtype=np.dtype('int8'))
	coords_spec = [(lats, 0), (lons, 1)]
	cube = iris.cube.Cube(dummy, dim_coords_and_dims=coords_spec)

	return cube

#example
if __name__ == "__main__":
	cube = iris.load("~/Downloads/zos_Omon_CanESM5_piControl_r1i1p1f1_gn_520101-540012.nc")[0]
	
	new_cube = regrid(cube, "1x1", "area_weighted", False, False)
	
	iris.save(new_cube, "test_iris_cube.nc")
	
	exit()
	