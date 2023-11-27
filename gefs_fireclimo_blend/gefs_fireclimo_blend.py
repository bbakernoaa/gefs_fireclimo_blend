#!/bin/bash/env python

import xarray as xr
import pandas as pd
import numpy as np

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
    ratio_interp = ratio.interp(lat=clim.lat, lon=clim.lon, method='nearest')

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
def make_fire_emission(d=None,climos=None, ratio=0.9,scale_climo=True, n_forecast_days=35, obsfile='GBBEPx_all01GRID.emissions_v004_20190601.nc',climo_directory='climMean'):   
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
    
    # get the timestamp
    dd = pd.Timestamp(d) 
    
    #open fire emission
    g = xr.open_mfdataset(obsfile, decode_cf=False)
    
    # climo files to open
    final = dd + pd.Timedelta('{} days'.format(n_forecast_days))
    dates = pd.date_range(start=dd,end=final,freq='D')
    files = [t.strftime('{}}/GBBEPx-all01GRID_v4r0_climMean_%m%d.nc'.format(climo_directory)) for t in dates]

    # open climo file 
    climo = xr.open_mfdataset(files)

    # make weighted climo 
    gc = g.coarsen(lat=150,lon=150,boundary='trim').sum()
    
    dsets = []
    if scale_climo==True:
        climos = {}
    else:
        print('Not scaling climotology')
    for tslice in np.arange(n_forecast_days):
        print(tslice)
        #make copy of original data 
        if tslice==0:
            dset = g.copy()
        else:
            dset = dsets[tslice-1].copy()
        
        # fix time in copy 
        
        import pandas as pd

        dset.update({'time': [tslice]})
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
            dsets.append(dset)
    return dsets

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='blending GEFs with fire emissions',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-s',
        '--start_date',
        type=str,
        default=None,
        help='start day of processing. *Note: leave blank for Real-time format YYYY-mm-dd')
    parser.add_argument(
        '-n',
        '--n_forecast_days',
        type=str,
        default=None,
        help='number of days to process')
    parser.add_argument('-c',
                        '--climo_directory',
                        default='/scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/nexus/GBBEPx/v4/climMean',
                        help='Directory of the climate data')
    parser.add_argument('-o',
                        '--observation_file',
                        default='/scratch1/RDARCH/rda-arl-gpu/Barry.Baker/emissions/nexus/QFED/2021/',
                        help='NASA PASSWORD')
    args = parser.parse_args()

    start_date = args.start_doy
    n_forecast_days = args.n_forecast_days
    climate_directory = args.climo_directory
    observation_file = args.observation_file

    start = pd.Timestamp(start_date)

    dsets = make_fire_emission(d=start,
                     n_forecast_days=n_forecast_days,
                     climate_directory=climate_directory,
                     obsfile=observation_file)
    dset = xr.concat(dsets,dim='time')
    dset['time'] = dset.time.astype('int')
    #set.OC.where(dset.OC > 0).isel(time=-1).monet.quick_imshow(norm=LogNorm(),cmap='turbo',vmin=1e-11,vmax=1e-8)
    #for i in arange(35):
    #    dset.OC.where(dset.OC > 0).isel(time=i).monet.quick_imshow(norm=LogNorm(),cmap='turbo',vmin=1e-11,vmax=1e-8)
    #    outname = "/scratch2/data_untrusted/Barry.Baker/blended_output_{}.jpg".format(i)
    #    savefig(outname,dpi=200)
    #close('all')
    
    def write_ncf(dset,outfile):
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
        dset.load().to_netcdf(outfile, encoding=encoding)

