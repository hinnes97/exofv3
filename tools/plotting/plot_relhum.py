import xarray as xr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os


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


data = ds.relhum[-10:].mean('time')


fig,ax = plt.subplots(3,1,figsize = [5,10], sharex=True, sharey=True)
for n,p in enumerate([5, 50, 500]):
    i = np.searchsorted(ds.level, p)

    im = ax[n].contourf(lon_TL, lat_TL, data[i], levels=15, cmap=mpl.cm.coolwarm, vmin=0, vmax=1.0)
    ax[n].quiver(lon_TL[::4], lat_TL[::4], u.data[i,::4,::4], v.data[i,::4,::4])

fig.colorbar(mpl.cm.ScalarMappable(cmap=mpl.cm.coolwarm, norm=mpl.colors.Normalize(vmin=0,vmax=1.0)), ax=ax[:])
plt.savefig('relhum.pdf')
