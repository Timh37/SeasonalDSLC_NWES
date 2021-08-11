#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:40:08 2020

#plot maps of historical and future mean SWSA

@author: thermans
"""
import xarray as xr
import matplotlib.pyplot as plt
import os
import numpy as np
import fnmatch
import cmocean
import cartopy.crs as ccrs
import matplotlib
import matplotlib.gridspec as gridspec
from seasonal_deviations_from_monthly_means import seasonal_deviations_from_monthly_means
plt.close('all')

#settings
ssp = 'ssp585' #desired ssp
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/' #input directory ensemble data
out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure
model_list = list(np.genfromtxt('/Users/thermans/Documents/PhD/Phase4_seasonal/Analysis/compiling_ensembles/ens_model_list_'+ssp+'.txt',dtype='str'))

histmeanyrs = np.arange(1995,2015)
futmeanyrs = np.arange(2081,2101)

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#fetch land mask from zos data
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

modi = np.where(np.isin(taux_ds.model.values,['CIESM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CAS-ESM2-0', 'FIO-ESM-2-0']))[0] #defined positive westward
taux_ds['tauu'][modi,...] = -1*taux_ds['tauu'][modi,...]
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

####plotting
my_cmap = cmocean.cm.tempo #continuous colormap
my_cmap.set_bad('grey')

fig=plt.figure(figsize=(8.5,11.5))   
gs = fig.add_gridspec(8,3)
gs.update(top=.97,bottom=.1,hspace=.32,wspace=.1)

#historical
#winter
ax3 = plt.subplot(gs[0:2,0], projection=ccrs.Orthographic(0, 50)) #maps
absw = np.sqrt(taux_djf_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_djf_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False) #absolute magnitude
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_djf_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_djf_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   
    
gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(a) Winter (DJF)') #subplot title


#future
#winter
ax3 = plt.subplot(gs[0:2,1], projection=ccrs.Orthographic(0, 50)) #maps
absw = np.sqrt(taux_djf_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_djf_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_djf_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_djf_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = gl.left_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(e) Winter (DJF)') #subplot title

#historical
#spring
ax3 = plt.subplot(gs[2:4,0], projection=ccrs.Orthographic(0, 50))
absw = np.sqrt(taux_mam_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_mam_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_mam_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_mam_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.left_labels = True
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(b) Spring (MAM)') #subplot title

#future
#spring
ax3 = plt.subplot(gs[2:4,1], projection=ccrs.Orthographic(0, 50))

absw = np.sqrt(taux_mam_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_mam_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_mam_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_mam_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(f) Spring (MAM)') #subplot title


#historical
#summer
ax3 = plt.subplot(gs[4:6,0], projection=ccrs.Orthographic(0, 50))

absw = np.sqrt(taux_jja_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_jja_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_jja_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_jja_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(c) Summer (JJA)') #subplot title

#future
#summer
ax3 = plt.subplot(gs[4:6,1], projection=ccrs.Orthographic(0, 50))
absw = np.sqrt(taux_jja_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_jja_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_jja_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_jja_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(g) Summer (JJA)') #subplot title

#historical
#autumn
ax3 = plt.subplot(gs[6:8,0], projection=ccrs.Orthographic(0, 50))
absw = np.sqrt(taux_son_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_son_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_son_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_son_mean_anom.sel(year=histmeanyrs).mean(dim='year').mean(dim='model').values
ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(d) Autumn (SON)') #subplot title

#future
#autumn
ax3 = plt.subplot(gs[6:8,1], projection=ccrs.Orthographic(0, 50))
absw = np.sqrt(taux_son_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2+tauy_son_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model')**2) #absolute change

p = absw.plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=0,vmax=.1,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid whitespace between grid cells when rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_son_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
v = tauy_son_mean_anom.sel(year=futmeanyrs).mean(dim='year').mean(dim='model').values
q=ax3.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.65,color='black',width=.007,edgecolors='k',zorder=5)   

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels= False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(h) Autumn (SON)') #subplot title

#add colorbar
cax = fig.add_axes([0.23, .05, .3, .02])
fig.colorbar(p, cax=cax,orientation='horizontal',label='SWSA [N/m'r'$^{2}$'']')

#add quiverkey
qk = ax3.quiverkey(q, 0.575, 0.06, 0.05, label='0.05 N/m'r'$^{2}$', labelpos='E',
                   coordinates='figure') #quiver key on scale

#fig.savefig(os.path.join(out_dir,'suppFigure5_seasonal_windstress_ssp585.pdf'),dpi=300)
