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
f = 'atmos_daily_interp.nc'

cmap = mpl.cm.coolwarm

with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:

    dp = np.r_[ds.level[0]*100, np.diff(ds.level.data*100)]
    rp = 1.66e7
    g = 12.4
    
    strf = 2*np.pi*rp/g * np.cumsum((ds.vcomp[0].mean(['grid_xt'])*np.cos(np.deg2rad(ds.grid_yt))).data * dp[:,np.newaxis],axis=0)
    
    #data = ds.ucomp[-10:].mean(dim=['time', 'grid_xt'])

    plt.figure()

    plt.contourf(ds.grid_yt, ds.level[:]*100., strf, cmap=cmap, norm =mpl.colors.CenteredNorm(), levels=15)
    plt.colorbar()

    plt.xlabel('Latitude')
    plt.ylabel('Pressure [Pa]')
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')

    plt.tight_layout()

    plt.savefig(os.path.join(args.output_dir,'sf_interp.pdf'))
