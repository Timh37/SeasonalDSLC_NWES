#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot composite plots sea level & barotropic currents for ROMS sensitivity experiments with uniform wind speed change and open english channel

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

#initialize figure
my_cmap = cmocean.cm.balance
my_cmap.set_bad('grey')
fig=plt.figure(figsize=(6, 4.25))
gs = fig.add_gridspec(1,2)
gs.update(hspace = .2,wspace=.2,bottom=0.15)

plot_titles = ['(a) Exp_SW, DJF','(b) Exp_NE, JJA']

exps_dir = '/Users/thermans/Documents/Modeling/ROMS/northsea8/'
exps = ['swWindUniform_sq2ms','neWindUniform_sq2ms'] #experiments to plot

#open reference experiment (original)
org_dir = '/Users/thermans/Documents/Modeling/ROMS/northsea8/Exp39_1993_2018_era5_gl12v1_v3e10_rx013'
org_ds = xr.open_dataset(os.path.join(org_dir,'NorthSea8_avg_timeseries_monthly.nc'))#.isel(ocean_time=np.arange(0,36))

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

for e,exp in enumerate(exps):
    dns = fnmatch.filter(os.listdir(exps_dir),"*"+exp)
    ds = xr.open_dataset(os.path.join(exps_dir,dns[0],'NorthSea8_avg_timeseries_monthly.nc'))
    
    if e<1: #select months from desired season
        season='DJF'
    else:
        season='JJA'
        
    if season == 'DJF':
        time_i = [1,11,12,13,23,24,25,35]
    elif season == 'MAM':
        time_i = [2,3,4,14,15,16,26,27,28]    
    elif season == 'JJA':
        time_i = [5,6,7,17,18,19,29,30,31]
    elif season == 'SON':
        time_i = [8,9,10,20,21,22,32,33,34]

    
    diff_zeta = (ds-org_ds).zeta.isel(ocean_time=time_i).mean(dim='ocean_time') #calculate diff w.r.t. reference
    diff_vbar = (ds-org_ds).vbar.isel(ocean_time=time_i).mean(dim='ocean_time')
    diff_ubar = (ds-org_ds).ubar.isel(ocean_time=time_i).mean(dim='ocean_time')
    
    #interpolate u,v points to rho points
    diff_vbar_rho = np.empty((218,242))
    diff_vbar_rho[:] = np.nan
    diff_vbar_rho[1:-1,:] = (diff_vbar[0:-1,:] + diff_vbar[1:,:])/2
    diff_ubar_rho = np.empty((218,242))
    diff_ubar_rho[:] = np.nan
    diff_ubar_rho[:,1:-1] = (diff_ubar[:,0:-1] + diff_ubar[:,1:])/2

    ax = fig.add_subplot(gs[0,e], projection=ccrs.Orthographic(0, 50))
    
    im=(100*diff_zeta).plot.pcolormesh("lon_rho", "lat_rho", ax=ax,vmin=-12,vmax=12,cmap=my_cmap,label='dSSH [m]',transform=ccrs.PlateCarree(),add_colorbar=False,rasterized=True)
    bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
    gl.xlocator = mticker.FixedLocator([-5, 0, 5])
    
    if e>0:
        gl.top_labels = gl.right_labels = gl.left_labels = False #don't label top and right axes
    else:
        gl.top_labels = gl.right_labels = False #don't label top and right axes
    
    
    gl.xlabel_style = {'color': 'black','rotation':0}
    gl.ylabel_style = {'color': 'black','rotation':0}
                                                               
    ax.coastlines(resolution='10m',color='black')
    ax.set_extent([-8.5,7.5,47,59.5])
    q=ax.quiver(ds.lon_rho.values[::3,::3],ds.lat_rho.values[::3,::3],diff_ubar_rho[::3,::3],diff_vbar_rho[::3,::3],
              scale=.5,color='black',width=.005,edgecolors='k',transform=ccrs.PlateCarree())

    ax.set_title(plot_titles[e])
    
    if e==1:
        cax = fig.add_axes([0.25, .125, .5, .03])
        fig.colorbar(im, cax=cax,orientation='horizontal',label='Sea-level response [cm]',extend='neither')
        qk = ax.quiverkey(q, 0.81, 0.14, 0.05, label='5 cm/s', labelpos='E',coordinates='figure')
fig.savefig(os.path.join(out_dir,'Figure6_uniform_wind_roms.pdf'),dpi=300)