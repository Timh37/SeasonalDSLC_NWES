#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:53:08 2021

Plot composites of SSLA and SWSA, from NCEP (wind) and AVISO (sea level) in winter 2007 and summer 1995
In these years SWSA has a direction similar to dSWSA in CMIP6 models

@author: thermans
"""
import os
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean
import cartopy.crs as ccrs
from seasonal_deviations_from_monthly_means import seasonal_deviations_from_monthly_means

plt.close('all') 

#load data
aviso_dir = '/Users/thermans/Documents/Data/AVISO/'
ncep_dir = '/Users/thermans/Documents/Data/NCEP/'
out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

aviso = xr.open_dataset(os.path.join(aviso_dir,'dt_global_allsat_msla_h_y1993m01-2019m04.nc'))
aviso = aviso.isel(time=np.arange(11,311)) #monthly means 1993-2018

aviso.coords['lon'] = ((aviso.coords['lon'] + 180) % 360) - 180 #wrap around 0
aviso = aviso.reindex({ 'lon' : np.sort(aviso['lon'])})

aviso_djf_mean_anom, aviso_jja_mean_anom, aviso_mam_mean_anom, aviso_son_mean_anom, aviso_dec2nov_mean = seasonal_deviations_from_monthly_means(aviso.sla) #get SSLA
#
ncep = xr.open_mfdataset(ncep_dir+'*')
ncep = ncep.isel(time=np.arange(11,875))

ncep.coords['lon'] = ((ncep.coords['lon'] + 180) % 360) - 180 #wrap around 0
ncep = ncep.reindex({ 'lon' : np.sort(ncep['lon'])})

taux = 1.2 * 0.0015 * ncep.uwnd*np.sqrt(ncep.uwnd**2+ncep.vwnd**2) #convert wind speed to wind stress using a simple approximation
tauy = 1.2 * 0.0015 * ncep.vwnd*np.sqrt(ncep.uwnd**2+ncep.vwnd**2)

taux_djf_mean_anom, taux_jja_mean_anom, taux_mam_mean_anom, taux_son_mean_anom, taux_dec2nov_mean = seasonal_deviations_from_monthly_means(taux) #get SWSA
tauy_djf_mean_anom, tauy_jja_mean_anom, tauy_mam_mean_anom, tauy_son_mean_anom, tauy_dec2nov_mean = seasonal_deviations_from_monthly_means(tauy)

#plotting
my_cmap = cmocean.cm.balance
my_cmap.set_bad('grey')
fig = plt.figure(figsize=(6,4.25))
gs = fig.add_gridspec(1,2)

ax = fig.add_subplot(gs[0,0],projection=ccrs.Orthographic(0, 50)) #winter 2007

p=(100*aviso_djf_mean_anom.sel(year=2007)).plot(transform=ccrs.PlateCarree(),vmin=-12,vmax=12,cmap=my_cmap,ax=ax,add_colorbar=False,rasterized=True) #SSLA
p.set_edgecolor('face')
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax,levels=[-200],colors=['white'],add_colorbar=False)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}
                                                               
ax.coastlines(resolution='50m',color='black')
ax.set_extent([-9.5,9,45,60])
x=taux_djf_mean_anom.lon.values
y=tauy_djf_mean_anom.lat.values
u = taux_djf_mean_anom.sel(year=2007).values #SWSA
v = tauy_djf_mean_anom.sel(year=2007).values
ax.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.5,color='black',width=.007,edgecolors='k',zorder=5) 
ax.set_title('(a)')


ax = fig.add_subplot(gs[0,1],projection=ccrs.Orthographic(0, 50)) #summer 1995
p=(100*aviso_jja_mean_anom).sel(year=1995).plot(transform=ccrs.PlateCarree(),vmin=-12,vmax=12,cmap=my_cmap,ax=ax,add_colorbar=False,rasterized=True)
p.set_edgecolor('face')
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax,levels=[-200],colors=['white'],add_colorbar=False)
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}
                                                               
ax.coastlines(resolution='50m',color='black')
ax.set_extent([-9.5,9,45,60])
x=taux_jja_mean_anom.lon.values
y=tauy_jja_mean_anom.lat.values
u = taux_jja_mean_anom.sel(year=1995).values
v = tauy_jja_mean_anom.sel(year=1995).values
q = ax.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.5,color='black',width=.007,edgecolors='k',zorder=5) 
qk = ax.quiverkey(q, 0.8, 0.11, 0.05, label='0.05 N/m'r'$^{2}$', labelpos='E',
                   coordinates='figure')
ax.set_title('(b)')

cax = fig.add_axes([0.29, .12, .45, .03])

fig.colorbar(p, cax=cax,orientation='horizontal',label='SSLA [cm]',extend='neither')
#fig.savefig(os.path.join(out_dir,'suppFigure6_obs.pdf'),dpi=300)