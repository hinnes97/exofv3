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
    data = ds.ucomp[-1].mean(dim=['grid_xt'])

    plt.figure()

    plt.contourf(data.grid_yt, data.level*100., data, cmap=cmap, levels=15, norm=mpl.colors.CenteredNorm())
    plt.colorbar(label='Zonal wind [m/s]')

    plt.xlabel('Latitude')
    plt.ylabel('Pressure [Pa]')
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')

    plt.tight_layout()

    plt.savefig(os.path.join(args.output_dir,'zwind.pdf'))

