#!/usr/bin/env python

import xarray as xr
import pandas as pd
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def open_qfed(fname):
    from glob import glob
    from numpy import sort
    vrs = ['BC','CH4','CO','CO2','NH3','NOx','OC','PM2.5','SO2']
    qfed_vars = ['bc','ch4','co','co2','nh3','no','oc','pm25','so2']
    das = []
    
    print('')
    print('Opening QFED Files...')
    if len(fname) > 1:
        files = sort(fname)
    else:
        files = sort(glob(fname))

    for v in files:
        good = [ i for i in qfed_vars if i in v]
        if good:
            print('  opening:',v,good)
            dset = xr.open_dataset(v,decode_cf=False)
            var_index = vrs[qfed_vars.index(good[0])] # variable index
            da = dset['biomass']
            das.append(da)
    dset_dict = {}
    for index,v in enumerate(vrs):
        dset_dict[v] = das[index]
    dset = xr.Dataset(dset_dict)
    return dset

def open_climatology(fname):
    from glob import glob
    from numpy import sort
    
    # array to house datasets
    das = [] 
    print('')
    print('Opening Climatology Files...')
    
    if len(fname) > 1:
        files = sort(fname)
    else:
        files = sort(glob(fname))
    #print(files)
    xr.open_dataset(files[0])
    for i,f in enumerate(files):
        print('  opening:',f)
        das.append(xr.open_dataset(f, engine='netcdf4'))

    return xr.concat(das, dim='time')


def write_ncf(dset, outfile):
    """
    Write the given dataset to a NetCDF file with specified encoding.

    Parameters:
    - dset (xarray.Dataset): The dataset to be written to the NetCDF file.
    - outfile (str): The path and filename of the output NetCDF file.

    Returns:
    None
    """
    print('Output File:', outfile)
    encoding = {}
    for v in dset.data_vars:
        encoding[v] = dict(zlib=True, complevel=4)
    if 'latitude' in dset:
        encoding['latitude'] = dict(zlib=True, complevel=4)
        encoding['longitude'] = dict(zlib=True, complevel=4)
    if 'lat_b' in dset:
        encoding["lat_b"] = dict(zlib=True, complevel=4)
        encoding["lon_b"] = dict(zlib=True, complevel=4)
    if 'time' in dset:
        encoding['time'] = dict(dtype='i4')
    dset.load().to_netcdf(outfile, encoding=encoding)

def create_climatology(emissions, climatology, lat_coarse=50, lon_coarse=50):
    """
    Create scaled climatology data based on emission data.

    Parameters:
    emissions (xarray.DataArray): Emission data.
    climatology (xarray.Dataset): Input climatology data.
    lat_coarse (int, optional): Coarsening factor for latitude. Defaults to 50.
    lon_coarse (int, optional): Coarsening factor for longitude. Defaults to 50.

    Returns:
    xarray.Dataset: Scaled climatology data.

    """
    # Create a copy of the climatology
    clim = climatology.copy()

    # Coarsen the climatology
    clim_coarse = climatology.coarsen(lat=lat_coarse, lon=lon_coarse, boundary='trim').sum()

    # Calculate the ratio of emissions to climatology and handle NaN values
    ratio = (emissions.squeeze().data / clim_coarse.where(clim_coarse > 0)).fillna(0)

    # Interpolate the ratio to match the coordinates of the climatology
    ratio_interp = ratio.sel(lat=clim.lat, lon=clim.lon, method='nearest')

    # Loop through each time slice and scale the climatology
    for index, time_slice in enumerate(clim.time):
        # Get the current time slice of the climatology
        clim_slice = clim.data[index, :, :]

        # Calculate the weighted alpha ratio parameter
        alpha = 1. - 1. / (index + 1)

        # Scale the current time slice
        scaled_slice = clim_slice * ratio_interp[index, :, :]

        # Update the climatology with the scaled time slice
        clim.data[index, :, :] = scaled_slice.squeeze().data

    return clim.compute()


def make_fire_emission(d=None, climos=None, ratio=0.9, scale_climo=True, n_forecast_days=35, obsfile='GBBEPx_all01GRID.emissions_v004_20190601.nc',climo_directory='climMean'):   
    """
    Generate fire emissions data for a given date and forecast period.

    Parameters:
    - d (str or pd.Timestamp): The date for which fire emissions are generated.
    - climos (dict): Dictionary containing pre-calculated climatology data for scaling.
    - ratio (float): The ratio of original data to climatology data for blending.
    - scale_climo (bool): Flag indicating whether to scale the climatology data.
    - n_forecast_days (int): Number of forecast days.
    - obsfile (str): Path to the file containing observed fire emissions data.
    - climo_directory (str): Directory containing climatology files.

    Returns:
    - list: A list of xarray.Dataset objects representing fire emissions data for each forecast day.
    """
    import pandas as pd
    import numpy as np
    from glob import glob
    # get the timestamp
    dd = pd.Timestamp(d) 
    
    #open fire emission
    if len(obsfile) > 1:
        if "QFED" in obsfile[0]:
            g = open_qfed(obsfile)
        else:
            g = xr.open_mfdataset(obsfile, decode_cf=False)
    else:
        if 'QFED' in obsfile:
            g = open_qfed(obsfile)
        else:
            g = xr.open_mfdataset(obsfile, decode_cf=False)
    
    # climo files to open
    final = dd + pd.Timedelta('{} days'.format(n_forecast_days))
    if dd.is_leap_year:
        dates = pd.DatetimeIndex([ d + pd.DateOffset(years=1) + pd.DateOffset(days=n) for n in range(n_forecast_days)])
    else:
        dates = pd.date_range(start=dd,end=final,freq='D')
    files = [t.strftime('{}/GBBEPx-all01GRID_v4r0_climMean_%m%d.nc'.format(climo_directory)) for t in dates]

    # open climo file 
    #climo = xr.open_mfdataset(files)
    climo = open_climatology(files)
    climo = climo.sel(lat=g['lat'],lon=g['lon'],method='nearest')

    # make weighted climo 
    gc = g.coarsen(lat=150,lon=150,boundary='trim').sum()
    
    dsets = []
    climos = {}
    for tslice in np.arange(n_forecast_days):
        #print(tslice)
        #make copy of original data 
        if tslice==0:
            dset = g.copy()
        else:
            dset = dsets[tslice-1].copy()
        
        # fix time in copy 
        
        import pandas as pd

        dset.update({'time': [float(tslice*24)]})
        dset.time.attrs = g.time.attrs

        for v in g.data_vars:
            if scale_climo == False:
                if tslice > 5:
                    dset[v].data = (ratio * dset[v] + (1 - ratio) * climo[v].data[tslice, :, :])
            else:
                if tslice == 0:
                    print('creating climatology scaling for', v)
                    climos[v] = create_climatology(gc[v], climo[v], lon_coarse=150, lat_coarse=150)
                else:
                    # cn = create_climdata(gc[v],climo[v])
                    # print(cn)
                    if tslice > 5:
                        dset[v].data = (ratio * dset[v] + (1 - ratio) * climos[v].data[tslice, :, :])
                    else:
                        dset[v] = dset[v]
        dsets.append(dset)
    return dsets

if __name__ == '__main__':
    parser = ArgumentParser(description='Regrid MERRA2 to FV3 native netcdf input files', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-s',
        '--start_date',
        type=str,
        default=None,
        help='start day of processing. *Note: leave blank for Real-time format YYYY-mm-dd')
    parser.add_argument(
        '-n',
        '--n_forecast_days',
        type=int,
        default=None,
        help='number of days to process')
    parser.add_argument('-c',
                        '--climo_directory',
                        default='/scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/nexus/GBBEPx/v4/climMean',
                        help='Directory of the climate data')
    parser.add_argument('-f',
                        '--fire_file',
                        default='/scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/nexus/QFED/2021/',
                        help='input fire emission file/files', 
                        nargs='+')
    parser.add_argument('-o',
                        '--output_filename',
                        default='/scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/nexus/QFED/2021/',
                        help='output file name')
    parser.add_argument('-r',
                        '--ratio',
                        default=0.95,
                        type=float,
                        help='weighting ratio')
    args = parser.parse_args()

    start_date = args.start_date
    n_forecast_days = args.n_forecast_days
    climate_directory = args.climo_directory
    observation_file = args.fire_file
    outfname = args.output_filename
    ratio = float(args.ratio)

    start = pd.Timestamp(start_date)

    dsets = make_fire_emission(d=start,
                     n_forecast_days=n_forecast_days,
                     climo_directory=climate_directory,
                               obsfile=observation_file,
                               ratio=ratio)
    dset = xr.concat(dsets, dim='time').fillna(0.)
    dset['time'] = dset.time.astype('int')
    dset['time'].attrs['units'] = start.strftime('hours since %Y-%m-%d 12:00:00')
    
    # write netcdf file out 
    write_ncf(dset,outfname)


    

