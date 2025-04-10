import xarray as xr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
from scipy.integrate import cumulative_trapezoid
parser = argparse.ArgumentParser()
parser.add_argument("output_dir", type=str)
args = parser.parse_args()
f = 'atmos_daily.nc'

cmap = mpl.cm.coolwarm

with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:

    dp = ds.phalf.diff('phalf').swap_dims(dict(phalf='pfull'))    
    rp = 1.66e7
    g = 12.4
    
    strf = 2*np.pi*rp/g * np.cumsum((ds.vcomp[0,:-2].mean(['grid_xt'])*np.cos(np.deg2rad(ds.grid_yt))).data * dp.data[:-2,np.newaxis],axis=0)
    
    #data = ds.ucomp[-10:].mean(dim=['time', 'grid_xt'])

    plt.figure()

    plt.contourf(ds.grid_yt, ds.pfull[:-2]*100., strf, cmap=cmap, norm =mpl.colors.CenteredNorm(), levels=15)
    plt.colorbar()

    plt.xlabel('Latitude')
    plt.ylabel('Pressure [Pa]')
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')

    plt.tight_layout()

    plt.savefig(os.path.join(args.output_dir,'sf_cutoff.pdf'))
