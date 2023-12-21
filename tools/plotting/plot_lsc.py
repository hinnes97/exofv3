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

u = ds.ucomp[-1]
v = ds.vcomp[-1]


data = ds.t_dt_lsc[-1]

fig,ax = plt.subplots(3,1,figsize = [5,10], sharex=True, sharey=True)
for n,p in enumerate([50,75,200]):
    i = np.searchsorted(ds.level, p)

    im = ax[n].contourf(lon_TL, lat_TL, data[i], levels=15, cmap=mpl.cm.coolwarm, norm=mpl.colors.CenteredNorm())
    ax[n].quiver(lon_TL[::4], lat_TL[::4], u.data[i,::4,::4], v.data[i,::4,::4])
    plt.colorbar(im, ax=ax[n])
plt.savefig('lsc.pdf')
