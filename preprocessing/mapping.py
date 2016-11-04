# -*- coding: utf-8 -*-
import netCDF4
import numpy as np
import pyproj
import rasterio
import scipy.spatial

def calc_grid_mapping(lons, lats, search_lons, search_lats):
    """
    Calculate mapping (nearest neighbors) between two grids.

    Args:
        lons: 2-D array of reference longitudes.
        lats: 2-D array of reference latitudes.
        search_lons: 2-D array of longitudes for which their nearest neighbor (within the lons/lats arrays) should be calculated.
        search_lats: 2-D array of latitudes for which their nearest neighbor (within the lons/lats arrays) should be calculated.

    Returns:
        Two arrays containing the x and y indices of the nearest neighbors for all search coordinates.
    """
    lonlats = np.array([lons.ravel(), lats.ravel()]).T
    x_index, y_index = np.meshgrid(np.arange(lons.shape[1]), np.arange(lons.shape[0]))
    indices = np.array([x_index.ravel(), y_index.ravel()]).T
    kd = scipy.spatial.KDTree(lonlats)
    query = kd.query(np.array([search_lons.ravel(), search_lats.ravel()]).T)
    nearest_points = indices[query[1]]
    map_x = nearest_points[:, 0].reshape(search_lons.shape)
    map_y = nearest_points[:, 1].reshape(search_lons.shape)

    return map_x, map_y


basin_file = '../examples/data/basin_kirchbichl_1k.tif'
reference_model_file = '../examples/data/HISTALP_temperature_1780-2008_inn.nc'
forecast_model_file = '/Users/flo/tmp/CFS/flxf.01.1998050100.nc'
reference_model_mapping_file = '/tmp/histalp_mapping.npz'
forecast_model_mapping_file = '/tmp/cfs_mapping.npz'

with rasterio.open(basin_file) as ds:
    basin = ds.read(1, masked=True)
    basin = basin.astype(int).filled(-1) # nodata values in uint8 files are often set to 255 -> set them to -1

proj_basin = pyproj.Proj(ds.crs)
proj_wgs84 = pyproj.Proj('+init=EPSG:4326')

x = np.linspace(ds.bounds.left + ds.res[0] / 2, ds.bounds.right - ds.res[0] / 2, ds.width) # (cell centers)
y = np.linspace(ds.bounds.top - ds.res[0] / 2, ds.bounds.bottom + ds.res[0] / 2, ds.height)
X, Y = np.meshgrid(x, y)
basin_lon, basin_lat = pyproj.transform(proj_basin, proj_wgs84, X, Y)

nc_reference = netCDF4.Dataset(reference_model_file, 'r')
nc_forecast = netCDF4.Dataset(forecast_model_file, 'r')

reference_lons = nc_reference.variables['lon'][:]
reference_lats = nc_reference.variables['lat'][:]
forecast_lons = nc_forecast.variables['longitude'][:]
forecast_lats = nc_forecast.variables['latitude'][:]

# check for longitudes greater 180 degrees
reference_lons[reference_lons > 180] -= 360
forecast_lons[forecast_lons > 180] -= 360

reference_lon_grid, reference_lat_grid = np.meshgrid(reference_lons, reference_lats)
forecast_lon_grid, forecast_lat_grid = np.meshgrid(forecast_lons, forecast_lats)

# forecast model -> reference model
hc_map_x, hc_map_y = calc_grid_mapping(forecast_lon_grid, forecast_lat_grid, reference_lon_grid, reference_lat_grid)

# reference model -> model grid
model_map_x, model_map_y = calc_grid_mapping(reference_lon_grid, reference_lat_grid, basin_lon, basin_lat)

np.savez(reference_model_mapping_file, mapping_x=model_map_x, mapping_y=model_map_y)
np.savez(forecast_model_mapping_file, mapping_x=hc_map_x, mapping_y=hc_map_y)
