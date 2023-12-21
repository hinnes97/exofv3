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
f = 'atmos_daily_interp.nc'

cmap = mpl.cm.coolwarm

with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:
    data = ds.t_dt_lsc[-1:,:,28:36,:].weighted(np.cos(np.deg2rad(ds.grid_yt[28:36]))).mean(('time', 'grid_yt'))

    plt.figure()

    plt.contourf(data.grid_xt, data.level*100., data, cmap=cmap)
    plt.colorbar()

    plt.xlabel('Longitude')
    plt.ylabel('Pressure [Pa]')
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')

    plt.tight_layout()

    plt.savefig(os.path.join(args.output_dir,'eqtdtlsc.pdf'))
