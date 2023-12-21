import numpy as np
import xarray as xr
from windspharm.standard import VectorWind
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import os
from plot_tl_sf import interpolate_to_TL_ndim, transform_velocities_to_TL_interp

parser = argparse.ArgumentParser()
parser.add_argument("output_dir", type=str)
args = parser.parse_args()

Rp = 1.66e7
g = 12.4

ds = xr.open_dataset('atmos_daily_interp.nc', decode_times=False)

lon_TL = ds.grid_xt.data
lat_TL = ds.grid_yt.data

u = ds.ucomp[-10:].mean('time')
v = ds.vcomp[-10:].mean('time')

utl, vtl = transform_velocities_to_TL_interp(u.data,v.data,
                                             np.radians(ds.grid_yt.data), np.radians(ds.grid_xt.data),
                                             np.radians(lat_TL), np.radians(lon_TL),
                                             lon_ss = 180.)

temp = ds.temp[-10:].mean('time')

ttl = interpolate_to_TL_ndim(np.radians(lat_TL), np.radians(lon_TL),
                             np.radians(lat_TL[:,None]), np.radians(lon_TL[None,:]),
                             temp.data, lon_ss=180.)


fig,ax = plt.subplots(3,1,figsize = [5,10], sharex=True, sharey=True)
for n,p in enumerate([100,150,200]):
    i = np.searchsorted(ds.level, p)

    ax[n].contourf(lon_TL, lat_TL, ttl[i], levels=15, cmap=mpl.cm.coolwarm)
    ax[n].quiver(lon_TL[::4], lat_TL[::4], utl[i,::4,::4], vtl[i,::4,::4])

plt.savefig('tlwind.pdf')
