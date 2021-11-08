#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

#plot maps of ensemble mean annual mean dslc and wind stress change

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

ssp = 'ssp585'
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/' #input directory ensemble data
out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure
model_list = list(np.genfromtxt('/Users/thermans/Documents/PhD/Phase4_seasonal/Analysis/compiling_ensembles/ens_model_list_'+ssp+'.txt',dtype='str'))

baseyears = np.arange(1995,2015) #reference periods
futyears = np.arange(2081,2101)

#load variant averaged ensemble for NWES (1x1)
nwes_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc")[0]))  
nwes_ds = nwes_ds.sel(model=model_list)
nwes_ds = nwes_ds.isel(time=np.arange(11,3011)) #start in December end in November

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#get seasonal anomalies
ens_djf_mean_anom, ens_jja_mean_anom, ens_mam_mean_anom, ens_son_mean_anom, ens_dec2nov_mean = seasonal_deviations_from_monthly_means(nwes_ds.zos)
ens_dec2nov_mean_d = ens_dec2nov_mean - ens_dec2nov_mean.sel(year=baseyears).mean(dim='year')#change relative to base period

#prepare data for plotting
numfin = np.sum(np.isfinite(ens_dec2nov_mean_d.isel(year=0)),axis=0) #get number of models with ocean value for each grid point
min_num = 5

a = ens_dec2nov_mean_d #winter
a = a.where(np.isfinite(ens_dec2nov_mean_d.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
a = a.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

#load wind stress
taux_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "tauu_CMIP6_ssp585_n33_nwes_variant_averaged.nc")[0]))
taux_ds = taux_ds.isel(time=np.arange(11+124*12,3011))
taux_ds = taux_ds.sel(model=model_list)
taux_ds = taux_ds.where(numfin>min_num-1) #apply mask from zos

taux_ds.coords['lon'] = ((taux_ds.coords['lon'] + 180) % 360) - 180 #wrap around 0
taux_ds = taux_ds.reindex({ 'lon' : np.sort(taux_ds['lon'])})

modi = np.where(np.isin(taux_ds.model.values,['CIESM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CAS-ESM2-0', 'FIO-ESM-2-0']))[0] #defined positive westward
taux_ds['tauu'][modi,...] = -1*taux_ds['tauu'][modi,...] #change sign
taux_djf_mean_anom, taux_jja_mean_anom, taux_mam_mean_anom, taux_son_mean_anom, taux_dec2nov_mean = seasonal_deviations_from_monthly_means(taux_ds.tauu)

#open ensemble dataset wind stress y
tauy_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "tauv_CMIP6_ssp585_n33_nwes_variant_averaged.nc")[0]))
tauy_ds = tauy_ds.isel(time=np.arange(11+124*12,3011))
tauy_ds = tauy_ds.sel(model=model_list)
tauy_ds = tauy_ds.where(numfin>min_num-1) #apply mask from zos

tauy_ds.coords['lon'] = ((tauy_ds.coords['lon'] + 180) % 360) - 180 #wrap around 0
tauy_ds = tauy_ds.reindex({ 'lon' : np.sort(tauy_ds['lon'])})

modi = np.where(np.isin(tauy_ds.model.values,['CIESM', 'CMCC-CM2-SR5', 'CMCC-ESM2', 'CAS-ESM2-0', 'FIO-ESM-2-0']))[0]
tauy_ds['tauv'][modi,...] = -1*tauy_ds['tauv'][modi,...]
tauy_djf_mean_anom, tauy_jja_mean_anom, tauy_mam_mean_anom, tauy_son_mean_anom, tauy_dec2nov_mean = seasonal_deviations_from_monthly_means(tauy_ds.tauv)

#calculate change future relative to baseline
taux_dec2nov_mean_d = taux_dec2nov_mean.sel(year=futyears).mean(dim='year') - taux_dec2nov_mean.sel(year=baseyears).mean(dim='year')
tauy_dec2nov_mean_d = tauy_dec2nov_mean.sel(year=futyears).mean(dim='year') - tauy_dec2nov_mean.sel(year=baseyears).mean(dim='year')

####plotting
my_cmap = cmocean.cm.balance #continuous colormap
my_cmap.set_bad('grey')

fig=plt.figure(figsize=(8,5))   
gs = fig.add_gridspec(1,2)
gs.update(hspace = 0.75,wspace=0.3,right=.84,top=.95,bottom=.05)

#annual mean slc
ax3 = plt.subplot(gs[0,0], projection=ccrs.Orthographic(0, 50)) #maps
p=(100*a).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-25,vmax=25,cmap=my_cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(a)') #subplot title

cax = fig.add_axes([0.125, .15, .312, .02])
fig.colorbar(p, cax=cax,orientation='horizontal',label='SLC [cm]')

#annual mean wind-stress change
my_cmap = cmocean.cm.tempo
my_cmap.set_bad('grey')
abs_d = np.sqrt(taux_dec2nov_mean_d.mean(dim='model')**2+tauy_dec2nov_mean_d.mean(dim='model')**2) #absolute change
ax1 = plt.subplot(gs[0,1], projection=ccrs.Orthographic(0, 50)) #maps

im = abs_d.plot(transform=ccrs.PlateCarree(),ax=ax1,vmin=0,vmax=.03,cmap=my_cmap,add_colorbar=False)
bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False)
im.set_edgecolor('face') #avoid whitespace between grid cells when rendering
x=taux_djf_mean_anom.lon.values #quivers
y=taux_djf_mean_anom.lat.values
u = taux_dec2nov_mean_d.mean(dim='model').values
v = tauy_dec2nov_mean_d.mean(dim='model').values
q=ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.2,color='black',width=.007,edgecolors='k',zorder=5)   
qk = ax1.quiverkey(q, 0.73, 0.239, 0.02, label='0.02 N/m'r'$^{2}$', labelpos='E',
                   coordinates='figure') #quiver key on scale
    
ax1.coastlines(resolution='50m',color='black')
ax1.set_extent([-9.5,9,45,60])
ax1.set_title('(b)')

gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}
cax = fig.add_axes([0.529, .15, .312, .02])
fig.colorbar(im, cax=cax,orientation='horizontal',label='Wind-stress change [N/m'r'$^{2}$'']')

#fig.savefig(os.path.join(out_dir,'suppFigure1_annualmean_change_ssp585.pdf'),dpi=300)
