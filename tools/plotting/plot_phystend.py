import xarray as xr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os
import seaborn as sns

mpl.rcParams['font.size'] = 10
sns.set_palette('colorblind')
lw2 = 7.126625 # Exact double-column width in inches
lw = 3.365 # single column in inches
phi = 0.5*(5**0.5 + 1)
hgt = 9.129


# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("output_dir", type=str)
parser.add_argument("--plevs", nargs='+',type=float)
args = parser.parse_args()

cmap = mpl.cm.magma



with xr.open_dataset(os.path.join(args.output_dir, 'atmos_daily_interp.nc'), decode_times=False) as ds:
    for p in args.plevs:
        fig, ax = plt.subplots(3,2, figsize=[lw2, lw2*1.3])


        tmin = 1.e20
        tmax = -1.e20
        
        data = ds.t_dt_ding[-10:].mean('time').sel(level=p,method='nearest')
        tmax = max(tmax, np.amax(data))
        tmin = min(tmin, np.amin(data))

        data = ds.t_dt_lsc[-10:].mean('time').sel(level=p,method='nearest')
        tmax = max(tmax, np.amax(data))
        tmin = min(tmin, np.amin(data))

        data = ds.t_dt_rainout[-10:].mean('time').sel(level=p,method='nearest')
        tmax = max(tmax, np.amax(data))
        tmin = min(tmin, np.amin(data))

        norm = mpl.colors.Normalize(vmin=tmin, vmax=tmax)
        data = ds.t_dt_ding[-10:].mean('time').sel(level=p,method='nearest')
        ax[0,0].contourf(data.grid_xt, data.grid_yt, data,
                         norm = norm, cmap=cmap)
        

        data = ds.t_dt_lsc[-10:].mean('time').sel(level=p,method='nearest')
        ax[1,0].contourf(data.grid_xt, data.grid_yt, data,
                         norm=norm, cmap=cmap)

        data = ds.t_dt_rainout[-10:].mean('time').sel(level=p,method='nearest')
        imt = ax[2,0].contourf(data.grid_xt, data.grid_yt, data,
                               norm=norm, cmap=cmap)

        qmin = 1.e20
        qmax = -1.e20
        
        data = ds.q_dt_ding[-10:].mean('time').sel(level=p,method='nearest')
        qmax = max(tmax, np.amax(data))
        qmin = min(tmin, np.amin(data))

        data = ds.q_dt_lsc[-10:].mean('time').sel(level=p,method='nearest')
        qmax = max(tmax, np.amax(data))
        qmin = min(tmin, np.amin(data))

        data = ds.q_dt_rainout[-10:].mean('time').sel(level=p,method='nearest')
        qmax = max(tmax, np.amax(data))
        qmin = min(tmin, np.amin(data))

        norm = mpl.colors.Normalize(vmin=tmin, vmax=tmax)
        data = ds.q_dt_ding[-10:].mean('time').sel(level=p,method='nearest')
        ax[0,0].contourf(data.grid_xt, data.grid_yt, data,
                         norm=norm, cmap=cmap)
        
        data = ds.q_dt_lsc[-10:].mean('time').sel(level=p,method='nearest')
        ax[1,0].contourf(data.grid_xt, data.grid_yt, data,
                         norm=norm, cmap=cmap)

        data = ds.q_dt_rainout[-10:].mean('time').sel(level=p,method='nearest')
        imq = ax[2,0].contourf(data.grid_xt, data.grid_yt, data,
                               norm=norm, cmap=cmap)

        cbar1 = fig.colorbar(imt, ax=ax[:,0], location='left')
        cbar2 = fig.colorbar(imq, ax=ax[:,1], location='right')

        plt.tight_layout()
        
        plt.savefig(os.path.join(args.output_dir, f"phystend_p{p:.0f}.pdf"))
        
        
