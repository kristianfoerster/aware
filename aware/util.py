from __future__ import print_function, division, absolute_import
import dateutil
import netCDF4
import numpy as np
import pandas as pd
import rasterio

def read_gdal_file(filename, fill_value=None, return_prj_settings=False):
    with rasterio.open(filename) as ds:
        data = ds.read(1, masked=(fill_value is not None))

        if fill_value is not None:
            data = data.filled(fill_value)

        if return_prj_settings:
            meta = ds.meta
            meta = {k: v for k, v in meta.items() if k != 'transform'} # avoid "GDAL-style transforms are deprecated and will not be supported in Rasterio 1.0" warning
            return data, meta
        else:
            return data

def write_gdal_file(filename, data, prj_settings):
    with rasterio.open(filename, 'w', **prj_settings) as ds:
        ds.write_band(1, data)

def num2date(ncvar):
    time_vals = ncvar[:]

    # special case for handling "months since" which is not supported by netCDF4-python:
    if 'months since' in ncvar.units:
        time_vals_i = np.array(time_vals, dtype = int)
        basedate = pd.Timestamp(ncvar.units.split('since')[1])

        # dr = pd.date_range(start=basedate, periods=np.ceil(time_vals_i.max()) + 1, freq='MS') # caution: "MS" sets the date to the beginning of the month (whereas "M" sets it to the end)!
        # s = pd.Series(dr)
        s = pd.Series(index=np.floor(time_vals_i).astype(int)) # take floor because relativedelta can only handle integers
        for timediff in s.index:
            if type(timediff) == int:
                td = timediff
            else:
                td = np.asscalar(timediff)
            s[timediff] = basedate + dateutil.relativedelta.relativedelta(months=td)

        s = s.apply(lambda t: (t - pd.datetime(1970, 1, 1)).total_seconds()) # convert to seconds
        s = s.reindex(sorted(list(set(time_vals_i) | set(s.index)))).interpolate(method='linear')
        s = pd.to_datetime(s, unit='s') 
        s = s[time_vals_i]

        dates = pd.DatetimeIndex(s).to_period('M').to_timestamp('M') # take last day of month as date
        return dates
    else:
        return netCDF4.num2date(time_vals, ncvar.units)

def center_lon_lat(filename):
    with rasterio.open(filename) as ds:
        return ds.lnglat()
