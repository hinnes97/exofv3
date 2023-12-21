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

parser = argparse.ArgumentParser()
parser.add_argument("output_dir", type=str)
args = parser.parse_args()

Rp = 1.66e7
g = 12.4

def transform_latlon_to_TL(lat_tmp,lon_tmp,lon_ss=0.):
    if lat_tmp.ndim==1 or lon_tmp.ndim==1:
        lon,lat = np.meshgrid(lon_tmp,lat_tmp)
    else:
        lat,lon = lat_tmp,lon_tmp
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)

    lon_ss = -lon_ss
    xprime = np.cos(lon_ss)*x - np.sin(lon_ss)*y
    yprime = np.sin(lon_ss)*x + np.cos(lon_ss)*y

    lat_TL = np.arcsin(xprime)
    lon_TL = np.arctan2(yprime,z)
    # change from lon in {-pi,pi} to lon in {0,2pi}:
    lon_TL[lon_TL<0.] = lon_TL[lon_TL<0.] + 2.*np.pi

    return lat_TL, lon_TL

def transform_velocities_to_TL(u,v,lat_tmp,lon_tmp,lon_ss=0.):
    # This returns TL velocities but on irregular grid.
    if lat_tmp.ndim==1 or lon_tmp.ndim==1:
        lon,lat = np.meshgrid(lon_tmp,lat_tmp)
    else:
        lat,lon = lat_tmp,lon_tmp
    lat_TL,lon_TL = transform_latlon_to_TL(lat,lon,lon_ss)

    lon_ss = -lon_ss

    # (tedious algebra - got these via Mathematica)
    Dlon_TL_Dlon = 1./( 1./np.cos(lon+lon_ss)*np.tan(lat)+np.cos(lat)/np.sin(lat)*np.sin(lon+lon_ss)*np.tan(lon+lon_ss) )
    Dlon_TL_Dlat = 8.*np.sin(lon+lon_ss) / ( -6.+2.*np.cos(2.*lat)+np.cos(2.*(lat-lon-lon_ss))+\
                                                  2.*np.cos(2*(lon+lon_ss))+np.cos(2.*(lat+lon+lon_ss)) )
    Dlat_TL_Dlon = -np.cos(lat)*np.sin(lon+lon_ss)/( np.sqrt(1.-np.cos(lat)**2*np.cos(lon+lon_ss)**2) )
    Dlat_TL_Dlat = -np.cos(lon+lon_ss)*np.sin(lat)/( np.sqrt(1.-np.cos(lat)**2*np.cos(lon+lon_ss)**2) )

    u_TL = Dlon_TL_Dlon * np.cos(lat_TL)/np.cos(lat)*u + Dlon_TL_Dlat*np.cos(lat_TL)*v
    v_TL = Dlat_TL_Dlon * u/np.cos(lat) + Dlat_TL_Dlat*v
    return u_TL,v_TL

def transform_velocities_to_TL_interp( u,v,lat,lon,lat_TL,
                                      lon_TL,lon_ss=0.,method="nearest"):
    # This returns TL velocities, interpolated onto regular grid.
    u_TL,v_TL = transform_velocities_to_TL(u,v,lat,lon,lon_ss)
    u_TL_i = interpolate_to_TL_ndim(lat,lon,lat_TL,lon_TL,u_TL,lon_ss,method=method)
    v_TL_i = interpolate_to_TL_ndim(lat,lon,lat_TL,lon_TL,v_TL,lon_ss,method=method)
    return u_TL_i,v_TL_i

def interpolate_to_TL(lat,lon,lat_TL,lon_TL,data,
                          lon_ss=0.,method="nearest"):
    # Assume data is 2D!
    # Interpolate data, given on an earth-like lat-lon grid, to a TL lat-lon grid.
    # First, transform the given lat-lon points into TL coords.
    # Second, use the given points as basis for interpolation.
    lat_TL_given,lon_TL_given = transform_latlon_to_TL(lat,lon,lon_ss)
    data_interp = griddata( (lat_TL_given.ravel(),lon_TL_given.ravel()),\
                                data.ravel(),(lat_TL,lon_TL),method=method )
    return data_interp

def interpolate_to_TL_ndim(lat,lon,lat_TL,lon_TL,data,lon_ss=0.,method="nearest"):
    # Uses above, but for N-dim data.
    # Also assume lat/lon are last two dims.
    if data.ndim==2:
        data_interp = np.squeeze(interpolate_to_TL(lat,lon,\
                                            lat_TL[None,:],\
                                            lon_TL[:,None],data,lon_ss,method))
    elif data.ndim==3:
        data_interp = np.zeros((data.shape[0],lat_TL.size,lon_TL.size))
        for i in range(data.shape[0]):
            data_interp[i] = \
            np.squeeze(interpolate_to_TL(lat,lon,lat_TL[:,None],lon_TL[None,:],\
                                      data[i],lon_ss,method))
    elif data.ndim==4:
        data_interp = np.zeros((data.shape[0],data.shape[1], \
                                lat_TL.size,lon_TL.size))
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                data_interp[i,j,...] = \
                                     interpolate_to_TL(lat,lon,\
                                        lat_TL[None,:],lon_TL[:,None],\
                                        data[i,j,...],lon_ss,method)
    else:
        # fix this later...x
        print("(interpolate_to_TL_ndim) Error! expected 4 or fewer dimensions")
        pass

    return data_interp 

with xr.open_dataset('atmos_daily.nc', decode_times=False) as ds:
    lon_TL = ds.grid_xt.data*np.pi/180
    lat_TL = ds.grid_yt.data*np.pi/180

    u = ds.ucomp[0]
    v = ds.vcomp[0]
    
    u_TL, v_TL = transform_velocities_to_TL_interp(u.data, v.data,
                                                   np.radians(ds.grid_yt.data), np.radians(ds.grid_xt.data),
                                                   lat_TL, lon_TL, lon_ss = np.pi)

    dp = ds.phalf.diff('phalf').swap_dims(dict(phalf='pfull')).data

    integrand = dp[:,np.newaxis]/g * np.pi * 2* Rp * np.average(v_TL, axis=-1) * np.cos(lat_TL)[np.newaxis,:]
    strf = np.cumsum(integrand, axis=0)
    
    plt.figure()
    cmap = mpl.cm.coolwarm
    plt.contourf(np.degrees(lat_TL), ds.pfull*100., strf, cmap=cmap, norm =mpl.colors.CenteredNorm())
    plt.colorbar()

    plt.xlabel('Latitude')
    plt.ylabel('Pressure [Pa]')
    plt.gca().invert_yaxis()
    plt.gca().set_yscale('log')

    plt.tight_layout()

    plt.savefig(os.path.join(args.output_dir,'tlsf.pdf'))
