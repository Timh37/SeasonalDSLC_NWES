#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:40:08 2020

#plot the winter dSSLA of ACCESS-ESM1-5 and composite plots of the response of sea level & barotropic currents in ROMS
to winter dSWSpA from ACCESS-ESM1-5 with a closed and open English Channel

@author: thermans
"""
import xarray as xr
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import os
import numpy as np
import fnmatch
import cmocean
import cartopy.crs as ccrs
plt.close('all')

out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#initialize figure
my_cmap = cmocean.cm.balance
my_cmap.set_bad('grey')

fig=plt.figure(figsize=(8, 4.25))
gs = fig.add_gridspec(1,3)
gs.update(hspace = .2,wspace=.2,bottom=0.12)

#load & plot ACCESS-ESM1-5 change data (similar to cmip6_plot_dSWSA_dSSLA_native_im.py)
ssp = 'ssp585'
cmip6_zos_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/'
cmip6_wind_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/'
dzos_model_path = os.path.join(cmip6_zos_dir,'dzos_0mean','DJF_anom','ACCESS-ESM1-5')

zos_djf = xr.open_mfdataset(dzos_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',coords='minimal',data_vars=['zos'])
zos_djf_vm = zos_djf.zos.mean(dim='variant') #take mean across variants

unmasked = np.isfinite(zos_djf_vm.values[:])
masked = ~np.isfinite(zos_djf_vm.values[:])

ax = plt.subplot(gs[0,0], projection=ccrs.Orthographic(0, 50))
#scatterplots
p=plt.scatter(zos_djf_vm.longitude.values[:][unmasked],zos_djf_vm.latitude.values[:][unmasked],c=100*zos_djf_vm.values[:][unmasked],vmin=-6, vmax=6,cmap=my_cmap,s=80,transform=ccrs.PlateCarree(),marker='s')
ma=plt.scatter(zos_djf_vm.longitude.values[:][masked],zos_djf_vm.latitude.values[:][masked],c='grey',s=85,transform=ccrs.PlateCarree(),marker='s')

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax,levels=[-200],colors=['white'],add_colorbar=False) #shelf break

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}
gl.xlocator = mticker.FixedLocator([-5, 0, 5])

ax.coastlines(resolution='50m',color='black')
ax.set_extent([-8.5,7.5,47,59.5])
ax.set_title('(a) ACCESS-ESM1-5')

#load & plot ROMS experiments (similar to plot_cmip6_roms_experiments.py)
plot_titles = ['(b) Exp_ACC_cc','(c) Exp_ACC']

exps_dir = '/Users/thermans/Documents/Modeling/ROMS/northsea8/'
exps = ['ACCESS-ESM1-5_DJF_wind_chclosed','ACCESS-ESM1-5_DJF_wind'] #experiments to plot

#open reference experiment (with channel closed or open)
org_dir = '/Users/thermans/Documents/Modeling/ROMS/northsea8/'
org_exps = ['Exp70_1993_1995_era5_gl12v1_v3e10_rx013_chclosed','Exp39_1993_2018_era5_gl12v1_v3e10_rx013']

for e,exp in enumerate(exps):
    dns = fnmatch.filter(os.listdir(exps_dir),"*"+exp)
    ds = xr.open_dataset(os.path.join(exps_dir,dns[0],'NorthSea8_avg_timeseries_monthly.nc'))
    org_ds = xr.open_dataset(os.path.join(org_dir,org_exps[e],'NorthSea8_avg_timeseries_monthly.nc'))
    
    season='DJF'
   
    if season == 'DJF':
        time_i = [1,11,12,13,23,24,25,35]
    elif season == 'MAM':
        time_i = [2,3,4,14,15,16,26,27,28]    
    elif season == 'JJA':
        time_i = [5,6,7,17,18,19,29,30,31]
    elif season == 'SON':
        time_i = [8,9,10,20,21,22,32,33,34]

    
    diff_zeta = (ds-org_ds).zeta.isel(ocean_time=time_i).mean(dim='ocean_time')
    diff_vbar = (ds-org_ds).vbar.isel(ocean_time=time_i).mean(dim='ocean_time')
    diff_ubar = (ds-org_ds).ubar.isel(ocean_time=time_i).mean(dim='ocean_time')
    
    #interpolate u,v points to rho points
    diff_vbar_rho = np.empty((218,242))
    diff_vbar_rho[:] = np.nan
    diff_vbar_rho[1:-1,:] = (diff_vbar[0:-1,:] + diff_vbar[1:,:])/2
    diff_ubar_rho = np.empty((218,242))
    diff_ubar_rho[:] = np.nan
    diff_ubar_rho[:,1:-1] = (diff_ubar[:,0:-1] + diff_ubar[:,1:])/2

    ax = fig.add_subplot(gs[0,e+1], projection=ccrs.Orthographic(0, 50))
    
    
    im=(100*diff_zeta).plot.pcolormesh("lon_rho", "lat_rho", ax=ax,vmin=-6,vmax=6,cmap=my_cmap,label='dSSH [m]',transform=ccrs.PlateCarree(),add_colorbar=False,rasterized=True)
    bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax,levels=[-200],colors=['white'],add_colorbar=False)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
    gl.xlocator = mticker.FixedLocator([-5, 0, 5])
    
    gl.top_labels = gl.right_labels = gl.left_labels = False #don't label top and right axes
    
    gl.xlabel_style = {'color': 'black','rotation':0}
    gl.ylabel_style = {'color': 'black','rotation':0}
                                                               
    ax.coastlines(resolution='10m',color='black',zorder=4)
    
    if e==0:
        plt.scatter(1.55,51.05,edgecolor='lightgreen',facecolor='none',transform=ccrs.PlateCarree(),marker='o',s=250,linewidth=1.5,zorder=5)

    
    ax.set_extent([-8.5,7.5,47,59.5])
    q=ax.quiver(ds.lon_rho.values[::3,::3],ds.lat_rho.values[::3,::3],diff_ubar_rho[::3,::3],diff_vbar_rho[::3,::3],
              scale=.5,color='black',width=.005,edgecolors='k',transform=ccrs.PlateCarree())

    ax.set_title(plot_titles[e])
    
    if e==1:
        cax = fig.add_axes([0.25, .135, .5, .03])
        fig.colorbar(im, cax=cax,orientation='horizontal',label='Sea-level response [cm]',extend='neither')
        qk = ax.quiverkey(q, 0.81, 0.14, 0.05, label='5 cm/s', labelpos='E',coordinates='figure')
#fig.savefig(os.path.join(out_dir,'Figure8_access.pdf'),dpi=300)