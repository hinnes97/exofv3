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
parser.add_argument("--plevs", nargs='+',type=float)
args = parser.parse_args()

files = ['atmos_daily_interp.nc', 'atmos_average_interp.nc']

cmap = mpl.cm.coolwarm
for f in files:
    with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:
        for p in args.plevs:
            data = ds.temp.sel(level=p, method='nearest')

            plt.figure()
            
            plt.contourf(data.grid_xt, data.grid_yt, data[-1], cmap=cmap)
            plt.colorbar()

            plt.xlabel("Longitude")
            plt.ylabel("Latitude")
            plt.title(f"p = {p:.0f} mbar")

            plt.tight_layout()
            file_root='_'.join(f.split('_')[:2])
            plt.savefig(os.path.join(args.output_dir,f"{file_root}_temp_p{p:.0f}.pdf"))
            
