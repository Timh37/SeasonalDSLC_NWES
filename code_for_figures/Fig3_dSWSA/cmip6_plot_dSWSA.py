#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#plot maps of dSWSA

@author: thermans
"""
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import numpy as np
import fnmatch
import cmocean
from seasonal_deviations_from_monthly_means import seasonal_deviations_from_monthly_means
import cartopy.crs as ccrs
plt.close('all')

#settings
ssp = 'ssp585' #desired ssp
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/' #input directory ensemble data
out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure
model_list = list(np.genfromtxt('/Users/thermans/Documents/PhD/Phase4_seasonal/Analysis/compiling_ensembles/ens_model_list_'+ssp+'.txt',dtype='str'))

baseyrs = np.arange(1995,2015)
futyrs = np.arange(2081,2101)

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#load zos to apply its land mask to wind stress
zos_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc")[0]))
numfin = np.sum(np.isfinite(zos_ds.zos.isel(time=0)),axis=0) #get number of models with ocean value for each grid point
min_num = 5

#open ensemble dataset wind stress zonal
taux_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "tauu_CMIP6_ssp585_n33_nwes_variant_averaged.nc")[0]))
taux_ds = taux_ds.isel(time=np.arange(11+124*12,3011))
taux_ds = taux_ds.sel(model=model_list)

#mask land values
taux_ds = taux_ds.where(numfin>min_num-1)

taux_ds.coords['lon'] = ((taux_ds.coords['lon'] + 180) % 360) - 180 #wrap around 0
taux_ds = taux_ds.reindex({ 'lon' : np.sort(taux_ds['lon'])})

modi = np.where(np.isin(taux_ds.model.values,['CIESM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CAS-ESM2-0', 'FIO-ESM-2-0']))[0]
taux_ds['tauu'][modi,...] = -1*taux_ds['tauu'][modi,...] #these models appear positive westward, so reverse
taux_djf_mean_anom, taux_jja_mean_anom, taux_mam_mean_anom, taux_son_mean_anom, taux_dec2nov_mean = seasonal_deviations_from_monthly_means(taux_ds.tauu)

#open ensemble dataset wind stress meridional
tauy_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "tauv_CMIP6_ssp585_n33_nwes_variant_averaged.nc")[0]))
tauy_ds = tauy_ds.isel(time=np.arange(11+124*12,3011))
tauy_ds = tauy_ds.sel(model=model_list)
#mask land values
tauy_ds = tauy_ds.where(numfin>min_num-1)

tauy_ds.coords['lon'] = ((tauy_ds.coords['lon'] + 180) % 360) - 180 #wrap around 0
tauy_ds = tauy_ds.reindex({ 'lon' : np.sort(tauy_ds['lon'])})

modi = np.where(np.isin(tauy_ds.model.values,['CIESM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CAS-ESM2-0', 'FIO-ESM-2-0']))[0]
tauy_ds['tauv'][modi,...] = -1*tauy_ds['tauv'][modi,...]
tauy_djf_mean_anom, tauy_jja_mean_anom, tauy_mam_mean_anom, tauy_son_mean_anom, tauy_dec2nov_mean = seasonal_deviations_from_monthly_means(tauy_ds.tauv)

#calculate change future relative to baseline
taux_djf_mean_anom_d = taux_djf_mean_anom.sel(year=futyrs).mean(dim='year') - taux_djf_mean_anom.sel(year=baseyrs).mean(dim='year')
tauy_djf_mean_anom_d = tauy_djf_mean_anom.sel(year=futyrs).mean(dim='year') - tauy_djf_mean_anom.sel(year=baseyrs).mean(dim='year')

taux_jja_mean_anom_d = taux_jja_mean_anom.sel(year=futyrs).mean(dim='year') - taux_jja_mean_anom.sel(year=baseyrs).mean(dim='year')
tauy_jja_mean_anom_d = tauy_jja_mean_anom.sel(year=futyrs).mean(dim='year') - tauy_jja_mean_anom.sel(year=baseyrs).mean(dim='year')

taux_mam_mean_anom_d = taux_mam_mean_anom.sel(year=futyrs).mean(dim='year') - taux_mam_mean_anom.sel(year=baseyrs).mean(dim='year')
tauy_mam_mean_anom_d = tauy_mam_mean_anom.sel(year=futyrs).mean(dim='year') - tauy_mam_mean_anom.sel(year=baseyrs).mean(dim='year')

taux_son_mean_anom_d = taux_son_mean_anom.sel(year=futyrs).mean(dim='year') - taux_son_mean_anom.sel(year=baseyrs).mean(dim='year')
tauy_son_mean_anom_d = tauy_son_mean_anom.sel(year=futyrs).mean(dim='year') - tauy_son_mean_anom.sel(year=baseyrs).mean(dim='year')

taux_dec2nov_mean_d = taux_dec2nov_mean.sel(year=futyrs).mean(dim='year') - taux_dec2nov_mean.sel(year=baseyrs).mean(dim='year')
tauy_dec2nov_mean_d = tauy_dec2nov_mean.sel(year=futyrs).mean(dim='year') - tauy_dec2nov_mean.sel(year=baseyrs).mean(dim='year')
 

### plotting
my_cmap = cmocean.cm.tempo
my_cmap.set_bad('grey')

fig=plt.figure(figsize=(6, 7.5))
gs = fig.add_gridspec(2,2)
gs.update(hspace = 0.1)
gs.update(wspace = 0.1)

#a) winter
abs_d = np.sqrt(taux_djf_mean_anom_d.mean(dim='model')**2+tauy_djf_mean_anom_d.mean(dim='model')**2) #absolute change
ax1 = plt.subplot(gs[0,0], projection=ccrs.Orthographic(0, 50))

im = abs_d.plot(transform=ccrs.PlateCarree(),ax=ax1,vmin=0,vmax=.03,cmap=my_cmap,add_colorbar=False)
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
im.set_edgecolor('face') #avoid whitespace between grid cells when rendering
x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_djf_mean_anom_d.mean(dim='model').values
v = tauy_djf_mean_anom_d.mean(dim='model').values
ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.2,color='black',width=.007,edgecolors='k',zorder=5)   
    
ax1.coastlines(resolution='50m',color='black')
ax1.set_extent([-9.5,9,45,60])
ax1.set_title('(a) Winter (DJF)')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

#b) spring
abs_d = np.sqrt(taux_mam_mean_anom_d.mean(dim='model')**2+tauy_mam_mean_anom_d.mean(dim='model')**2)
ax1 = plt.subplot(gs[0,1], projection=ccrs.Orthographic(0, 50))

im = abs_d.plot(transform=ccrs.PlateCarree(),ax=ax1,vmin=0,vmax=.03,cmap=my_cmap,add_colorbar=False)
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False)
im.set_edgecolor('face')
x=taux_djf_mean_anom.lon.values
y=taux_djf_mean_anom.lat.values
u = taux_mam_mean_anom_d.mean(dim='model').values
v = tauy_mam_mean_anom_d.mean(dim='model').values
ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.2,color='black',width=.007,edgecolors='k',zorder=5)                                                
ax1.coastlines(resolution='50m',color='black')
ax1.set_extent([-9.5,9,45,60])
ax1.set_title('(b) Spring (MAM)')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels= gl.bottom_labels= False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

#c) summer
abs_d = np.sqrt(taux_jja_mean_anom_d.mean(dim='model')**2+tauy_jja_mean_anom_d.mean(dim='model')**2)
ax1 = plt.subplot(gs[1,0], projection=ccrs.Orthographic(0, 50))

im = abs_d.plot(transform=ccrs.PlateCarree(),ax=ax1,vmin=0,vmax=.03,cmap=my_cmap,add_colorbar=False)
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False)
im.set_edgecolor('face')
x=taux_djf_mean_anom.lon.values
y=taux_djf_mean_anom.lat.values
u = taux_jja_mean_anom_d.mean(dim='model').values
v = tauy_jja_mean_anom_d.mean(dim='model').values
ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.2,color='black',width=.007,edgecolors='k',zorder=5)          
ax1.coastlines(resolution='50m',color='black')
ax1.set_extent([-9.5,9,45,60])
ax1.set_title('(c) Summer (JJA)')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

#d) autumn
abs_d = np.sqrt(taux_son_mean_anom_d.mean(dim='model')**2+tauy_son_mean_anom_d.mean(dim='model')**2)
ax1 = plt.subplot(gs[1,1], projection=ccrs.Orthographic(0, 50))

im = abs_d.plot(transform=ccrs.PlateCarree(),ax=ax1,vmin=0,vmax=.03,cmap=my_cmap,add_colorbar=False)
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False)
im.set_edgecolor('face')
x=taux_djf_mean_anom.lon.values
y=taux_djf_mean_anom.lat.values
u = taux_son_mean_anom_d.mean(dim='model').values
v = tauy_son_mean_anom_d.mean(dim='model').values
q = ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.2,color='black',width=.007,edgecolors='k',zorder=5)          
ax1.coastlines(resolution='50m',color='black')
ax1.set_extent([-9.5,9,45,60])
ax1.set_title('(d) Autumn (SON)')

qk = ax1.quiverkey(q, 0.8, 0.08, 0.02, label='0.02 N/m'r'$^{2}$', labelpos='E',
                   coordinates='figure') #quiver key on scale

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels= False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax1.legend(bbox_to_anchor=(1.1, -1),frameon=False)
cax = fig.add_axes([0.31, .065, .4, .02])

fig.colorbar(im, cax=cax,orientation='horizontal',label=r'$\Delta$' 'SWSA [N/m'r'$^{2}$' ']')

#fig.savefig(os.path.join(out_dir,'Figure3_ens_wind_v2.pdf'),dpi=300) #store figure