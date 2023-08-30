import xarray as xr
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import os

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("output_dir", type=str)
parser.add_argument("--lons", nargs='+',type=float)
parser.add_argument("--lats", nargs='+', type=float)
args = parser.parse_args()

files = ['atmos_daily_interp.nc', 'atmos_average_interp.nc']

cmap = mpl.cm.coolwarm
for f in files:
    with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:
        plt.figure()
        for x in args.lons:
            for y in args.lats:

                data = ds.temp.sel(grid_xt=x, method='nearest').sel(grid_yt=y,method='nearest')



                plt.plot(data[-1], data.level, label=f"Lon = {x:.0f} deg, Lat = {y:.0f} deg")



    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()

    plt.plot(np.maximum(500*(data.level/data.level[-1])**(2./7.), 200.), data.level)
    plt.legend()
    plt.tight_layout()
    file_root='_'.join(f.split('_')[:2])
    plt.savefig(os.path.join(args.output_dir, f"{file_root}_pt.pdf"))
