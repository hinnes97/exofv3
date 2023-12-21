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
parser.add_argument('q', type=str)
args = parser.parse_args()
lons = [0,180]
lats= [0,45,70]
files = ['atmos_daily_interp.nc', 'atmos_average_interp.nc']

cmap = mpl.cm.coolwarm
for f in files:
    with xr.open_dataset(os.path.join(args.output_dir, f), decode_times=False) as ds:
        plt.figure()
        for x in lons:
            for y in lats:

                data = ds.temp.sel(grid_xt=x, method='nearest').sel(grid_yt=y,method='nearest')



                plt.plot(data[-1], data.level*100, label=f"Lon = {x:.0f} deg, Lat = {y:.0f} deg")
                


    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    #plt.plot(np.maximum(500*(data.level/data.level[-1])**(2./7.), 200.), data.level)
    ds_2 = xr.open_dataset(f'/network/group/aopp/planetary/RTP026_INNES_MSN/moist_subnep/soc-rad-conv-fort/fortran/prelim_runs/ic_q0.1.nc')
    plt.plot(ds_2.Tf, ds_2.pfull, label='Initial condition')
    plt.legend()
    plt.tight_layout()
    file_root='_'.join(f.split('_')[:2])
    plt.savefig(os.path.join(args.output_dir, f"{file_root}_pt_new.pdf"))
