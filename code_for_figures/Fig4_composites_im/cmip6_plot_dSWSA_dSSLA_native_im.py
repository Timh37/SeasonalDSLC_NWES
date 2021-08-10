#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:40:08 2020

#plot composite plots of winter and summer dSSLA and dSWSA on raw/native model grids, 
and plot timeseries of winter and summerdSSLA and dSWSA near Esbjerg
@author: thermans
"""
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import numpy as np
import fnmatch
import cmocean
from seasonal_deviations_from_monthly_means import seasonal_deviations_from_monthly_means
import cartopy.crs as ccrs
import xesmf as xe
plt.close('all')

#settings
ssp = 'ssp585' #desired ssp
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/' #input directory ensemble data
out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure
cmip6_wind_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/time_merged/' #windstress directory
cmip6_zos_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/dedrifted_linear/' #sea level directory

models = ['CanESM5','UKESM1-0-LL','IPSL-CM6A-LR']

baseyrs = np.arange(1995,2015)

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#initialize plot
my_cmap = cmocean.cm.balance
my_cmap.set_bad('grey')
fig=plt.figure(figsize=(8, 7.5))
gs = fig.add_gridspec(9,3) #for 3 models
gs.update(hspace = .8,bottom=.05,top=.95)

abc = 'abcdefghij'
for m,model in enumerate(models): #loop over desired models and load dSSLA and dSWSA from model-ordered directories
    #DJF winter
    #load native grid change data
    dzos_model_path = os.path.join(cmip6_zos_dir,'dzos_0mean','DJF_anom',model)
    dtauu_model_path = os.path.join(cmip6_wind_dir,'dtauu','DJF_anom',model)
    dtauv_model_path = os.path.join(cmip6_wind_dir,'dtauv','DJF_anom',model)
    
    #open all variants at once and take the average
    tauu_djf = xr.open_mfdataset(dtauu_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',coords='minimal')
    tauv_djf = xr.open_mfdataset(dtauv_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',coords='minimal')
    zos_djf = xr.open_mfdataset(dzos_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',data_vars=['zos'],coords='minimal',compat='override')

    taux_djf_vm = tauu_djf.tauu.mean(dim='variant')
    tauy_djf_vm = tauv_djf.tauv.mean(dim='variant')
    zos_djf_vm = zos_djf.zos.mean(dim='variant')
    
    if taux_djf_vm.shape is not tauy_djf_vm.shape: #if u & v on different grids, then regrid y to x
        
        lon = [s for s in list(taux_djf_vm.coords) if "lon" in s][0]
        lat = [s for s in list(taux_djf_vm.coords) if "lat" in s][0]
        
        taux_djf_vm.coords[lon] = ((taux_djf_vm.coords[lon] + 180) % 360) - 180 #wrap around 0
        taux_djf_vm = taux_djf_vm.reindex({ lon : np.sort(taux_djf_vm[lon])})
        tauy_djf_vm.coords[lon] = ((tauy_djf_vm.coords[lon] + 180) % 360) - 180 #wrap around 0
        tauy_djf_vm = tauy_djf_vm.reindex({ lon : np.sort(tauy_djf_vm[lon])})
            
        regridder = xe.Regridder(tauy_djf_vm,taux_djf_vm,'bilinear') #regrid v to u points bilinearly (approximation, but only used for Figure 5)
        
        tauy_djf_vm = regridder(tauy_djf_vm)
    
    zos_lon_name = [s for s in list(zos_djf.coords) if "lon" in s][0] #get zos coordinate names for plotting
    zos_lat_name = [s for s in list(zos_djf.coords) if "lat" in s][0]
    
    #determine masked and unmasked ocean grid cells
    unmasked = np.isfinite(zos_djf_vm.values[:])
    masked = ~np.isfinite(zos_djf_vm.values[:])

    ax1 = plt.subplot(gs[0:3,np.mod(m,3)], projection=ccrs.Orthographic(0, 50))
    ax1.set_facecolor('grey')
    #scatterplots
    p=plt.scatter(zos_djf_vm[zos_lon_name].values[:][unmasked],zos_djf_vm[zos_lat_name].values[:][unmasked],c=100*zos_djf_vm.values[:][unmasked],vmin=-12, vmax=12,cmap=my_cmap,s=50,transform=ccrs.PlateCarree(),marker='s',rasterized=True)
    ma=plt.scatter(zos_djf_vm[zos_lon_name].values[:][masked],zos_djf_vm[zos_lat_name].values[:][masked],c='grey',s=50,transform=ccrs.PlateCarree(),marker='s',rasterized=True)
    plt.scatter(8.44,55.46,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=7) #location of Esbjerg
    
    bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False)
    
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
    gl.top_labels = gl.right_labels = gl.bottom_labels = False #don't label top and right axes
    if m>0:
        gl.left_labels = False
    
    gl.xlabel_style = {'color': 'black','rotation':0}
    gl.ylabel_style = {'color': 'black','rotation':0}

    #quivers
    x=taux_djf_vm.lon.values
    y=taux_djf_vm.lat.values
    u = taux_djf_vm.values
    v = tauy_djf_vm.values
    q=ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.4,color='black',width=.007,edgecolors='k',zorder=6)       
    
    ax1.coastlines(resolution='50m',color='black',zorder=5)
    ax1.set_extent([-9.5,9,45,60])
    ax1.set_title('('+str(abc[m])+') '+model+' (n='+str(len(zos_djf.variant))+')')
    
    #####
    #JJA summer similar to above
    dzos_model_path = os.path.join(cmip6_zos_dir,'dzos_0mean','JJA_anom',model)
    dtauu_model_path = os.path.join(cmip6_wind_dir,'dtauu','JJA_anom',model)
    dtauv_model_path = os.path.join(cmip6_wind_dir,'dtauv','JJA_anom',model)

    tauu_jja = xr.open_mfdataset(dtauu_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',coords='minimal')
    tauv_jja = xr.open_mfdataset(dtauv_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',coords='minimal')
    zos_jja = xr.open_mfdataset(dzos_model_path+'/*'+ssp+'*.nc',combine='nested',concat_dim='variant',data_vars=['zos'],coords='minimal',compat='override')

    taux_jja_vm = tauu_jja.tauu.mean(dim='variant')
    tauy_jja_vm = tauv_jja.tauv.mean(dim='variant')
    zos_jja_vm = zos_jja.zos.mean(dim='variant')
    
    if taux_jja_vm.shape is not tauy_jja_vm.shape: #u & v on different grids
        
        lon = [s for s in list(taux_jja_vm.coords) if "lon" in s][0]
        lat = [s for s in list(taux_jja_vm.coords) if "lat" in s][0]
        
        taux_jja_vm.coords[lon] = ((taux_jja_vm.coords[lon] + 180) % 360) - 180 #wrap around 0
        taux_jja_vm = taux_jja_vm.reindex({ lon : np.sort(taux_jja_vm[lon])})
        tauy_jja_vm.coords[lon] = ((tauy_jja_vm.coords[lon] + 180) % 360) - 180 #wrap around 0
        tauy_jja_vm = tauy_jja_vm.reindex({ lon : np.sort(tauy_jja_vm[lon])})
            
        regridder = xe.Regridder(tauy_jja_vm,taux_jja_vm,'bilinear') #regrid v to u points bilinearly (approximation, but only used for Figure 5)
        
        tauy_jja_vm = regridder(tauy_jja_vm)

    unmasked = np.isfinite(zos_jja_vm.values[:])
    masked = ~np.isfinite(zos_jja_vm.values[:])

    ax1 = plt.subplot(gs[3:6,np.mod(m,3)], projection=ccrs.Orthographic(0, 50))
    ax1.set_facecolor('grey')
    #scatterplots
    p=plt.scatter(zos_jja_vm[zos_lon_name].values[:][unmasked],zos_jja_vm[zos_lat_name].values[:][unmasked],c=100*zos_jja_vm.values[:][unmasked],vmin=-12, vmax=12,cmap=my_cmap,s=50,transform=ccrs.PlateCarree(),marker='s',rasterized=True)
    ma=plt.scatter(zos_jja_vm[zos_lon_name].values[:][masked],zos_jja_vm[zos_lat_name].values[:][masked],c='grey',s=50,transform=ccrs.PlateCarree(),marker='s',rasterized=True)
        
    plt.scatter(8.44,55.46,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=7)
    
    bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax1,levels=[-200],colors=['white'],add_colorbar=False)
    
    gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
    gl.top_labels = gl.right_labels = False #don't label top and right axes
    if m>0:
        gl.left_labels = False
    gl.xlabel_style = {'color': 'black','rotation':0}
    gl.ylabel_style = {'color': 'black','rotation':0}
    
    x=taux_jja_vm.lon.values
    y=taux_jja_vm.lat.values
    u = taux_jja_vm.values
    v = tauy_jja_vm.values
    q = ax1.quiver(x,y,u,v,transform=ccrs.PlateCarree(),scale=.4,color='black',width=.007,edgecolors='k',zorder=6)       

    if m == len(models)-1:
        qk = ax1.quiverkey(q, 0.815, 0.3, 0.05, label='0.05 N/m'r'$^{2}$', labelpos='N',coordinates='figure')
        
    ax1.coastlines(resolution='50m',color='black',zorder=5)
    ax1.set_extent([-9.5,9,45,60])
    ax1.set_title('('+str(abc[m+3])+') '+model+' (n='+str(len(zos_djf.variant))+')')

    if m==2: #add colorbar
        cax = fig.add_axes([0.31, .31, .4, .02])
        fig.colorbar(p, cax=cax,orientation='horizontal',label=r'$\Delta$' 'SSLA [cm]')
    
    
    ########### timeseries at Esbjerg in winter and summer
    loc_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n39_nslocations_all_variants.nc")[0])).sel(location='Esbjerg')
    loc_ds = loc_ds.isel(time=np.arange(11+124*12,3011))
    
    loc_djf_mean_anom, loc_jja_mean_anom, loc_mam_mean_anom, loc_son_mean_anom, loc_dec2nov_mean = seasonal_deviations_from_monthly_means(loc_ds.zos)
    seas_cmap = matplotlib.colors.ListedColormap(cmocean.cm.balance((20,95,170,235)), name='my_colormap_name') #generate 4 level discrete colormap
    
    #timeseries
    ax3 = plt.subplot(gs[7:9,np.mod(m,3)])
 
    sl_djf = (loc_djf_mean_anom.sel(model=model).mean(dim='variant')-loc_djf_mean_anom.sel(model=model).mean(dim='variant').sel(year=baseyrs).mean(dim='year'))
    sl_jja = (loc_jja_mean_anom.sel(model=model).mean(dim='variant')-loc_jja_mean_anom.sel(model=model).mean(dim='variant').sel(year=baseyrs).mean(dim='year'))
    
    p1=(100*sl_djf).plot(ax=ax3,label='DJF',color=seas_cmap(0))    
    p2=(100*sl_jja).plot(ax=ax3,label='JJA',color=seas_cmap(3))     
    
    if m==0:
        ax3.legend(ncol=2,loc='upper left',frameon=False)
         
    ax3.grid()
    ax3.set_ylim([-20,20])
    ax3.set_xlim([1990,2100])
    ax3.set_yticks([-20,-10,0,10,20])
    ax3.set_title('('+str(abc[m+6])+') '+model+' (n='+str(len(zos_djf.variant))+')')
    ax3.set_xlabel('')
    
    if m==0 or m==3:
        ax3.set_ylabel('')
    else:
        ax3.set_ylabel('')
        ax3.set_yticklabels('')
    if m==0:
        ax3.set_ylabel(r'$\Delta$' 'SSLA [cm]')

#fig.savefig(os.path.join(out_dir,'Figure4_composite_models.pdf'),dpi=300)
