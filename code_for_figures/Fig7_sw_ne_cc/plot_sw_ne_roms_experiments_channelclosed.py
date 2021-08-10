#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:40:08 2020

plot composite plots sea level & barotropic currents for ROMS sensitivity experiments with uniform wind speed change and closed english channel,
and difference in responses w.r.t. same experiments with an open english channel

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
fig=plt.figure(figsize=(6, 7.5))
gs = fig.add_gridspec(2,2)
gs.update(hspace = .1,wspace=.2,bottom=0.1,top=.98)

plot_titles = ['(a) Exp_SW_cc, DJF','(b) Exp_NE_cc, JJA','(c) Closed minus open','(d) Closed minus open']

exps_dir = '/Users/thermans/Documents/Modeling/ROMS/northsea8/'
exps = ['swWindUniform_sq2ms_chclosed','neWindUniform_sq2ms_chclosed','swWindUniform_sq2ms','neWindUniform_sq2ms'] #experiments to plot

#reference experiments (with open and closed channels)
org_dir = '/Users/thermans/Documents/Modeling/ROMS/northsea8/'
org_exps = ['Exp70_1993_1995_era5_gl12v1_v3e10_rx013_chclosed','Exp70_1993_1995_era5_gl12v1_v3e10_rx013_chclosed',
            'Exp39_1993_2018_era5_gl12v1_v3e10_rx013','Exp39_1993_2018_era5_gl12v1_v3e10_rx013']

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

for e,exp in enumerate(exps):
    dns = fnmatch.filter(os.listdir(exps_dir),"*"+exp)
    ds = xr.open_dataset(os.path.join(exps_dir,dns[0],'NorthSea8_avg_timeseries_monthly.nc'))
    org_ds = xr.open_dataset(os.path.join(org_dir,org_exps[e],'NorthSea8_avg_timeseries_monthly.nc'))
    
    if np.mod(e,2)==0: #select months from desired season
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

    if e==0: #store outcomes for closed channel to subtract outcomes from open channel
        diff_zeta_sw_cc = diff_zeta
        diff_ubar_rho_sw_cc = diff_ubar_rho
        diff_vbar_rho_sw_cc = diff_vbar_rho
    elif e==1:
        diff_zeta_ne_cc = diff_zeta
        diff_ubar_rho_ne_cc = diff_ubar_rho
        diff_vbar_rho_ne_cc = diff_vbar_rho        
        
    if e<2:
        ax = fig.add_subplot(gs[0,e], projection=ccrs.Orthographic(0, 50))
    else:
        ax = fig.add_subplot(gs[1,e-2], projection=ccrs.Orthographic(0, 50))
    

    sl = 100*diff_zeta
    u = diff_ubar_rho[::3,::3]
    v = diff_vbar_rho[::3,::3]
    
    if e==2:
        sl = 100*diff_zeta_sw_cc - sl
        u = diff_ubar_rho_sw_cc[::3,::3]-u
        v = diff_vbar_rho_sw_cc[::3,::3]-v
        
    elif e==3:
        sl = 100*diff_zeta_ne_cc - sl
        u = diff_ubar_rho_ne_cc[::3,::3]-u
        v = diff_vbar_rho_ne_cc[::3,::3]-v
    
    if e>1:
        im=sl.plot.pcolormesh("lon_rho", "lat_rho", ax=ax,vmin=-4,vmax=4,cmap=my_cmap,label='dSSH [m]',transform=ccrs.PlateCarree(),add_colorbar=False,rasterized=True)
    else:
        im=sl.plot.pcolormesh("lon_rho", "lat_rho", ax=ax,vmin=-12,vmax=12,cmap=my_cmap,label='dSSH [m]',transform=ccrs.PlateCarree(),add_colorbar=False,rasterized=True)
    bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax,levels=[-200],colors=['white'],add_colorbar=False)

    plt.scatter(1.55,51.05,edgecolor='lightgreen',facecolor='none',transform=ccrs.PlateCarree(),marker='o',s=250,linewidth=1.5,zorder=5) #mark closed channel with green circle
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
    gl.xlocator = mticker.FixedLocator([-5, 0, 5])
    
    gl.top_labels = gl.right_labels = gl.left_labels = gl.bottom_labels = False
    if np.mod(e,2)==0:
        gl.left_labels = True #don't label top and right axes
    if e>1:
        gl.bottom_labels = True
    
    gl.xlabel_style = {'color': 'black','rotation':0}
    gl.ylabel_style = {'color': 'black','rotation':0}
                                                               
    ax.coastlines(resolution='10m',color='black',zorder=4)
    ax.set_extent([-8.5,7.5,47,59.5])
    q=ax.quiver(ds.lon_rho.values[::3,::3],ds.lat_rho.values[::3,::3],u,v,
              scale=.5,color='black',width=.005,edgecolors='k',transform=ccrs.PlateCarree())

    ax.set_title(plot_titles[e])

    if e==1:
        cax = fig.add_axes([0.26, .57, .5, .02])
        fig.colorbar(im, cax=cax,orientation='horizontal',label='Sea-level response [cm]',extend='neither')
        qk = ax.quiverkey(q, 0.82, 0.59, 0.05, label='5 cm/s', labelpos='E',coordinates='figure')
        
    if e==3:
        cax = fig.add_axes([0.26, .09, .5, .02])
        fig.colorbar(im, cax=cax,orientation='horizontal',label='Sea-level difference [cm]',extend='neither')
        qk = ax.quiverkey(q, 0.82, 0.11, 0.05, label='5 cm/s', labelpos='E',coordinates='figure')
fig.savefig(os.path.join(out_dir,'Figure7_uniform_wind_cc_roms.pdf'),dpi=300)