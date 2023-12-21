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

w = VectorWind(np.moveaxis(np.flip(ds.ucomp[-1].data,axis=1),0,-1), np.moveaxis(np.flip(ds.vcomp[-1].data,axis=1),0,-1),rsphere=1.66e7)

udiv,vdiv,urot,vrot = w.helmholtz(truncation=20)

udiv = np.flip(urot, axis=0)
vdiv = np.flip(vrot, axis=0)

fig,ax = plt.subplots(3,1,figsize = [5,10], sharex=True, sharey=True)
for n,p in enumerate([50,75,100]):
    i = np.searchsorted(ds.level, p)

    quiv = ax[n].quiver(ds.grid_xt[::4], ds.grid_yt[::4], udiv[::4,::4,i], vdiv[::4,::4,i])
    ax[n].quiverkey(quiv, 0.9,1.05, np.amax(udiv[:,:,i]), label=f'{np.amax(udiv[:,:,i]):.1f} [m/s]',coordinates='axes')
    plt.savefig('rotu.pdf')
