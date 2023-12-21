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
args = parser.parse_args()

lats=[0,45,70]
lons = [0,90,180,270]
files = ['atmos_daily_interp.nc', 'atmos_average_interp.nc']

cmap = mpl.cm.coolwarm
for f in files:
    with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:
        plt.figure()
        for x in lons:
            for y in lats:

                data = ds.relhum.sel(grid_xt=x, method='nearest').sel(grid_yt=y,method='nearest')



                plt.plot(data[-1], data.level, label=f"Lon = {x:.0f} deg, Lat = {y:.0f} deg")



    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()

    #plt.plot(np.maximum(500*(data.level/data.level[-1])**(2./7.), 200.), data.level)
    plt.legend()
    plt.tight_layout()
    file_root='_'.join(f.split('_')[:2])
    plt.savefig(os.path.join(args.output_dir, f"{file_root}_rhprof.pdf"))
