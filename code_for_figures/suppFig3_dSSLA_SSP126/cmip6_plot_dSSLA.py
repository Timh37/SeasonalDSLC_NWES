#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#plot maps of dSSLA and multi-model distributions nearest to 8 example locations

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

#load variant averaged ensembles for NWES (1x1) and closest point to coordinates on native grids
if ssp=='ssp126':
    nwes_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp126_n35_nwes_variant_averaged.nc")[0]))
    locs_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp126_n36_nslocations_all_variants.nc")[0])) 
elif ssp=='ssp585':
    nwes_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc")[0]))
    locs_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n39_nslocations_all_variants.nc")[0]))     

nwes_ds = nwes_ds.sel(model=model_list)
locs_ds = locs_ds.sel(model=model_list)
locs_ds = locs_ds.mean(dim='variant')

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#extract seasonal & shifted annual means
nwes_ds = nwes_ds.isel(time=np.arange(11,3011)) #start in December end in November
locs_ds = locs_ds.isel(time=np.arange(11,3011))

#NWES:
ens_djf_mean_anom, ens_jja_mean_anom, ens_mam_mean_anom, ens_son_mean_anom, ens_dec2nov_mean = seasonal_deviations_from_monthly_means(nwes_ds.zos)

#change in seasonal anomalies relative to baseyears
ens_djf_mean_anom_d = ens_djf_mean_anom - ens_djf_mean_anom.sel(year=baseyears).mean(dim='year')
ens_jja_mean_anom_d = ens_jja_mean_anom - ens_jja_mean_anom.sel(year=baseyears).mean(dim='year')
ens_mam_mean_anom_d = ens_mam_mean_anom - ens_mam_mean_anom.sel(year=baseyears).mean(dim='year')
ens_son_mean_anom_d = ens_son_mean_anom - ens_son_mean_anom.sel(year=baseyears).mean(dim='year')

#North Sea locations
locs_djf_mean_anom, locs_jja_mean_anom, locs_mam_mean_anom, locs_son_mean_anom, locs_dec2nov_mean = seasonal_deviations_from_monthly_means(locs_ds.zos)

#change in seasonal anomalies relative to baseyears
locs_djf_mean_anom_d = locs_djf_mean_anom - locs_djf_mean_anom.sel(year=baseyears).mean(dim='year')
locs_jja_mean_anom_d = locs_jja_mean_anom - locs_jja_mean_anom.sel(year=baseyears).mean(dim='year')
locs_mam_mean_anom_d = locs_mam_mean_anom - locs_mam_mean_anom.sel(year=baseyears).mean(dim='year')
locs_son_mean_anom_d = locs_son_mean_anom - locs_son_mean_anom.sel(year=baseyears).mean(dim='year')


#prepare data for plotting
numfin = np.sum(np.isfinite(ens_djf_mean_anom_d.isel(year=0)),axis=0) #get number of models with ocean value for each grid point
min_num = 5

a = ens_djf_mean_anom_d #winter
a = a.where(np.isfinite(ens_djf_mean_anom_d.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
a = a.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

b = ens_mam_mean_anom_d #spring
b = b.where(np.isfinite(ens_djf_mean_anom_d.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
b = b.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

c = ens_jja_mean_anom_d #summer
c = c.where(np.isfinite(ens_djf_mean_anom_d.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
c = c.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

d = ens_son_mean_anom_d #autumn
d = d.where(np.isfinite(ens_djf_mean_anom_d.sel(year=futyears).mean(dim='year').mean(dim='model',skipna=True) )) #mask where slc masked
d = d.where(numfin>min_num-1) #mask where minimum number of models with ocean value not exceeded

####plotting
discrete_cmap = matplotlib.colors.ListedColormap(cmocean.cm.balance((20,95,170,235)), name='my_colormap_name') #generate 4 level discrete colormap
discrete_cmap.set_bad('grey') #set land mask color

cmap = cmocean.cm.balance #continuous colormap
cmap.set_bad('grey')

fig=plt.figure(figsize=(6.5,11.5)) #generate figure  
gs = fig.add_gridspec(10,4)
gs.update(hspace = 0.75,wspace=0.3,right=.84,top=.95,bottom=.05)

#a) winter
ax3 = plt.subplot(gs[0:3,0:2], projection=ccrs.Orthographic(0, 50)) #maps
p=(100*a).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-5,vmax=5,cmap=cmap,add_colorbar=False)
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False) #shelf break
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5) #example locations
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.bottom_labels = False 
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(a) Winter (DJF)') #subplot title

#b) spring
ax3 = plt.subplot(gs[0:3,2:4], projection=ccrs.Orthographic(0, 50))
p=(100*b).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-5,vmax=5,cmap=cmap,add_colorbar=False)
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
ax3.set_title('(b) Spring (MAM)') #subplot title

#c) summer
ax3 = plt.subplot(gs[3:6,0:2], projection=ccrs.Orthographic(0, 50))
p=(100*c).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-5,vmax=5,cmap=cmap,add_colorbar=False)#cbar_kwargs={'label': 'Seasonal SLC [cm]'}
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.7,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = False 
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(c) Summer (JJA)') #subplot title

#d) autumn
ax3 = plt.subplot(gs[3:6,2:4], projection=ccrs.Orthographic(0, 50))
p=(100*d).sel(year=futyears).mean(dim='year').mean(dim='model').plot(transform=ccrs.PlateCarree(),ax=ax3,vmin=-5,vmax=5,cmap=cmap,add_colorbar=False)#cbar_kwargs={'label': 'Seasonal SLC [cm]'}
p.set_edgecolor('face') #avoid white lines in pdf rendering

bathy.Band1.plot.contour(transform=ccrs.PlateCarree(),ax=ax3,levels=[-200],colors=['white'],add_colorbar=False)
ax3.coastlines(resolution='50m',color='black',linewidth=1)

plt.scatter(ns_coords_lon,ns_coords_lat,facecolor='white',transform=ccrs.PlateCarree(),s=26,marker='o',edgecolor='black',zorder=5)
for i, no in enumerate(no_locations): #annotate coastal points
    plt.text(ns_coords_lon[i]+.5,ns_coords_lat[i]-.5,no_locations[i],color='white',transform=ccrs.PlateCarree(),size=14,fontweight='bold',zorder=7)

gl = ax3.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1, color='lightgrey', alpha=0, linestyle='-')
gl.top_labels = gl.right_labels = gl.left_labels= False 
gl.xlabel_style = {'color': 'black','rotation':0}
gl.ylabel_style = {'color': 'black','rotation':0}

ax3.set_extent([-9.5,9,45,60]) #set longitude & latitude extent
ax3.set_title('(d) Autumn (SON)') #subplot title

#add colorbar
cax = fig.add_axes([0.85, .515, .02, .35])
fig.colorbar(p, cax=cax,orientation='vertical',label=r'$\Delta$' 'SSLA [cm]')

#plot multi-model distributions in panel e
row = [6,6,7,7,8,8,9,9] #subplot arrangement

for l,loc in enumerate(locations): #loop over example locations
    if np.mod(l,2)==1:
       plt.subplot(gs[row[l],2:4])
    else:
       plt.subplot(gs[row[l],0:2])
    plt.grid()
    
    plt.scatter(100*locs_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'),2*np.ones(len(100*locs_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'))),facecolors='none', edgecolors=discrete_cmap(3))
    plt.scatter(100*locs_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),2,facecolor='white',edgecolor='black')
            
    plt.scatter(100*locs_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'),np.ones(len(100*locs_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'))),facecolors='none', edgecolors=discrete_cmap(1))
    plt.scatter(100*locs_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),1,facecolor='white',edgecolor='black')

    plt.scatter(100*locs_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'),4*np.ones(len(100*locs_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'))),facecolors='none', edgecolors=discrete_cmap(0))
    plt.scatter(100*locs_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),4,facecolor='white',edgecolor='black')
    
    plt.scatter(100*locs_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'),3*np.ones(len(100*locs_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year'))),facecolors='none', edgecolors=discrete_cmap(2))
    plt.scatter(100*locs_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),3,facecolor='white',edgecolor='black')
    
    plt.xticks([-15,-10,-5,0,5,10,15],[])
        
    if (l % 2) == 0:
        plt.yticks([1,2,3,4],['SON','JJA','MAM','DJF'])
    else:
        plt.yticks([1,2,3,4],[])
        
    if l>5:
        plt.xlabel(r'$\Delta$' 'SSLA [cm]')
        plt.xticks([-15,-10,-5,0,5,10,15],['-15','-10','-5','0','5','10','15'])
    else:
        plt.xticks([-15,-10,-5,0,5,10,15],[])
        
    plt.ylim([0.25,4.75])
    plt.xlim([-15,15])
    plt.title(no_locations[l]+' '+loc,x=0.5,y=1.05,ha='center',va='center',fontsize=12)  
    
plt.annotate('(e)',(.465,.405),xycoords='figure fraction',size=12) 

#save figure
'''
if ssp=='ssp126':
    fig.savefig(os.path.join(out_dir,'suppFigure3_seasonal_dzos_ssp126.pdf'),dpi=300)
elif ssp=='ssp585':
    fig.savefig(os.path.join(out_dir,'Figure1_seasonal_dzos_ssp585.pdf'),dpi=300)
'''