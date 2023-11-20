# coding: utf-8
import xarray as xr
import monet
import warnings 
from matplotlib.colors import LogNorm
def create_climdata(emis_coarse,clim,lat_coarse=50,lon_coarse=50):
    # create copy of climatology
    cn = clim.copy()
    
    # coarsen climatology
    clim_coarse = clim.coarsen(lat=lat_coarse,lon=lon_coarse,boundary='trim').sum()
    
    ratio = (emis_coarse.squeeze().data / clim_coarse.where(clim_coarse > 0)).fillna(0)

    ratio_i = ratio.interp(lat=cn.lat,lon=cn.lon,method='nearest')
    
    # loop through and create 
    for n,t in enumerate(cn.time):
        # current time slice of climo
        c = cn.data[n,:,:]
        
        # weighted alpha ratio param 
        alpha = 1. -  1. / (n + 1) 
        
        # ratio on coarsen grid 
    #        ratio = emis_coarse / clim_coarse[n,:,:]
        
        # ratio interpolated back to fine grid
        #ratio_i = ratio.interp(lat=cn.lat,lon=cn.lon,method='nearest')    
        
        # current time slice scaling 
        #ans = c * ( alpha * ( ratio_i) )
        ans = c * ratio_i[n,:,:] # * ( 1 - ( 1 / (36 -n ) ))
        # climo scaling for period
        cn.data[n,:,:] = ans.squeeze().data
    
    return cn.compute() #, ratio_i
    
def make_fire_emission(d='20190601',climos=None, ratio=0.9,scale_climo=True):
    import pandas as pd
    
    # get the timestamp
    dd = pd.Timestamp(d) 
    
    #open fire emission
    g = xr.open_dataset(dd.strftime('GBBEPx_all01GRID.emissions_v004_%Y%m%d.nc'), decode_cf=False)
    
    # climo files to open
    final = dd + pd.Timedelta('35 days')
    dates = pd.date_range(start=dd,end=final,freq='D')
    files = [t.strftime('climMean/GBBEPx-all01GRID_v4r0_climMean_%m%d.nc') for t in dates]
    # open climo file 
    climo = xr.open_mfdataset(files)
    # make weighted climo 
    gc = g.coarsen(lat=150,lon=150,boundary='trim').sum()
    
    dsets = []
    if scale_climo==True:
        climos = {}
    else:
        print('Not scaling climotology')
    for tslice in np.arange(35):
        print(tslice)
        #make copy of original data 
        if tslice==0:
            dset = g.copy()
        else:
            dset = dsets[tslice-1].copy()
        
        # fix time in copy 
        
        dset.update({'time':[tslice]})
        dset.time.attrs = g.time.attrs
        
        for v in g.data_vars:
            if scale_climo == False:
                if tslice > 5:
                    dset[v].data = (ratio * dset[v] + (1 - ratio) * climo[v].data[tslice,:,:])
            else:
                if tslice == 0:
                    print('creating climotology scaling for ',v)
                    climos[v] = create_climdata(gc[v],climo[v],lon_coarse=150,lat_coarse=150)
                else:
                    #cn = create_climdata(gc[v],climo[v])
                    #print(cn)
                    if tslice > 5:
                        dset[v].data = (ratio * dset[v] + (1 - ratio) * climos[v].data[tslice,:,:])
        dsets.append(dset)
    return dsets
    
dsets = make_fire_emission()
dset = xr.concat(dsets,dim='time')
dset.OC.where(dset.OC > 0).isel(time=-1).monet.quick_imshow(norm=LogNorm(),cmap='turbo',vmin=1e-11,vmax=1e-8)
for i in arange(35):
    dset.OC.where(dset.OC > 0).isel(time=i).monet.quick_imshow(norm=LogNorm(),cmap='turbo',vmin=1e-11,vmax=1e-8)
    outname = "/scratch2/data_untrusted/Barry.Baker/blended_output_{}.jpg".format(i)
    savefig(outname,dpi=200)
    
close('all')
