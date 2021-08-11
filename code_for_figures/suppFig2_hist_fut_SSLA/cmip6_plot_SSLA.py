#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 15:40:08 2020

#plot maps of historical and future mean SSLA, and at 8 example locations

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

baseyears = np.arange(1995,2015)
futyears = np.arange(2081,2101)

locations = ['Aberdeen','Immingham','Liverpool','Newhaven','Brest','Den Helder','Cuxhaven','Esbjerg']
no_locations = ['1','2','3','4','5','6','7','8']

ns_coords_lon = [-2.08,0.187,-3.02,0.06,-4.49,4.75,8.717,8.44]
ns_coords_lat = [57.14,53.63,53.45,50.78,48.38,52.964,53.87,55.46]

#load variant averaged ensemble for NWES (1x1) and for 8 example locations
if ssp=='ssp126':
    nwes_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp126_n35_nwes_variant_averaged.nc")[0]))
    locs_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp126_n36_nslocations_all_variants.nc")[0])) 
elif ssp=='ssp585':
    nwes_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc")[0]))
    locs_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n39_nslocations_all_variants.nc")[0]))     

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#select models
nwes_ds = nwes_ds.sel(model=model_list)
locs_ds = locs_ds.sel(model=model_list)
locs_ds = locs_ds.mean(dim='variant')

#extract seasonal & shifted annual means
nwes_ds = nwes_ds.isel(time=np.arange(11,3011)) #start in December end in November
locs_ds = locs_ds.isel(time=np.arange(11,3011))

ens_djf_mean_anom, ens_jja_mean_anom, ens_mam_mean_anom, ens_son_mean_anom, ens_dec2nov_mean = seasonal_deviations_from_monthly_means(nwes_ds.zos)

#prepare data for plotting
numfin = np.sum(np.isfinite(ens_djf_mean_anom.isel(year=0)),axis=0) #get number of models with ocean value for each grid point
min_num = 5

a = ens_djf_mean_anom #winter
a = a.where(np.isfinite(ens_djf_mean_anom.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
a = a.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

b = ens_mam_mean_anom #spring
b = b.where(np.isfinite(ens_djf_mean_anom.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
b = b.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

c = ens_jja_mean_anom #summer
c = c.where(np.isfinite(ens_djf_mean_anom.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
c = c.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded


d = ens_son_mean_anom #autumn
d = d.where(np.isfinite(ens_djf_mean_anom.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
d = d.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

####plotting
my_cmap = matplotlib.colors.ListedColormap(cmocean.cm.balance((20,95,170,235)), name='my_colormap_name') #generate 4 level discrete colormap
my_cmap.set_bad('grey') #set land mask color

cma = cmocean.cm.balance #continuous colormap
cma.set_bad('grey')

fig=plt.figure(figsize=(8.5,11.5))   
gs = fig.add_gridspec(8,3)
gs.update(top=.97,bottom=.1,hspace=.32,wspace=.1)

#historical
#winter
ax3 = plt.subplot(gs[0:2,0], projection=ccrs.Orthographic(0, 50)) #maps
p=(100*a).sel(year=baseyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False) #ssla
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5) #example locations
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(a) Winter (DJF)') #subplot title

#future
#winter
ax3 = plt.subplot(gs[0:2,1], projection=ccrs.Orthographic(0, 50)) #maps
p=(100*a).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = gl.left_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(e) Winter (DJF)') #subplot title

#historical
#spring
ax3 = plt.subplot(gs[2:4,0], projection=ccrs.Orthographic(0, 50))
p=(100*b).sel(year=baseyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)
    

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.left_labels = True
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(b) Spring (MAM)') #subplot title

#future
#spring
ax3 = plt.subplot(gs[2:4,1], projection=ccrs.Orthographic(0, 50))
p=(100*b).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)
    

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(f) Spring (MAM)') #subplot title


#historical
#summer
ax3 = plt.subplot(gs[4:6,0], projection=ccrs.Orthographic(0, 50))
p=(100*c).sel(year=baseyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)#cbar_kwargs={'label': 'Seasonal SLC [cm]'}
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(c) Summer (JJA)') #subplot title

#future
#summer
ax3 = plt.subplot(gs[4:6,1], projection=ccrs.Orthographic(0, 50))
p=(100*c).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)#cbar_kwargs={'label': 'Seasonal SLC [cm]'}
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels = gl.bottom_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(g) Summer (JJA)') #subplot title

#historical
#autumn
ax3 = plt.subplot(gs[6:8,0], projection=ccrs.Orthographic(0, 50))
p=(100*d).sel(year=baseyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)#cbar_kwargs={'label': 'Seasonal SLC [cm]'}
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.5,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(d) Autumn (SON)') #subplot title

#future
#autumn
ax3 = plt.subplot(gs[6:8,1], projection=ccrs.Orthographic(0, 50))
p=(100*d).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-15,vmax=15,cmap=cma,add_colorbar=False)#cbar_kwargs={'label': 'Seasonal SLC [cm]'}
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.5,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels= False #don't label top and right axes
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(h) Autumn (SON)') #subplot title

#add colorbar
cax = fig.add_axes([0.23, .05, .3, .02])
fig.colorbar(p, cax=cax,orientation='horizontal',label='SSLA [cm]')

plt.annotate('(i)',(.665,.979),xycoords='figure fraction',fontsize=12) 

########## plot SSLA as seasonal cycle at example locations
#cycles
#North Sea locations
locs_djf_mean_anom, locs_jja_mean_anom, locs_mam_mean_anom, locs_son_mean_anom, locs_dec2nov_mean = seasonal_deviations_from_monthly_means(locs_ds.zos)

row = [0,0,1,1,2,2,3,3] #subplot arrangement
for l,loc in enumerate(locations):
    ax = plt.subplot(gs[l,2])
    ax.axhline(y=0,color='black',linestyle='dashed')
    ax.plot([0,1,2,3],[100*locs_djf_mean_anom.sel(year=baseyears,location=loc).mean(dim='year').mean(dim='model'),
                    100*locs_mam_mean_anom.sel(year=baseyears,location=loc).mean(dim='year').mean(dim='model'),
                    100*locs_jja_mean_anom.sel(year=baseyears,location=loc).mean(dim='year').mean(dim='model'),
                    100*locs_son_mean_anom.sel(year=baseyears,location=loc).mean(dim='year').mean(dim='model')],marker='o',linestyle='solid',label='1995-2014',markersize=4)
    
    ax.plot([0,1,2,3],[100*locs_djf_mean_anom.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),
                    100*locs_mam_mean_anom.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),
                    100*locs_jja_mean_anom.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),
                    100*locs_son_mean_anom.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model')],marker='s',linestyle='solid',label='2081-2100',markersize=4)
    ax.set_xticks([0,1,2,3])
    ax.set_title(no_locations[l]+' '+loc,x=0.5,y=1.05,ha='center',va='center',fontsize=12)  
    ax.set_ylim(-15,15)
    ax.set_xlim(-.1,3.1)
    ax.set_xticklabels([])
    ax.grid()
    ax.set_yticks([-15,-10,-5,0,5,10,15])
    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    
    ax.set_ylabel('SSLA [cm]')
    ax.set_yticklabels(['','-10','','0','','10',''])
        
    if l==7:
        ax.legend(bbox_to_anchor=(.85, -.3))
        ax.set_xticklabels(['DJF','MAM','JJA','SON'])
    
#fig.savefig(os.path.join(out_dir,'suppFigure2_seasonal_zos_ssp585.pdf'),dpi=300)