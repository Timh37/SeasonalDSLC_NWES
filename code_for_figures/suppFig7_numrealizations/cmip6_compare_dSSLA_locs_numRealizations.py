#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare dSSLA at coastal locations using a single v.s. using all available realizations per CMIP6 model

@author: thermans
"""
import xarray as xr
import matplotlib.pyplot as plt
import os
import numpy as np
import fnmatch
import cmocean
import matplotlib
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

#load multi-variant multi-model ensembles for closest point to coordinates on native grids

locs_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n39_nslocations_all_variants.nc")[0]))     
locs_ds = locs_ds.sel(model=model_list)
locs_ds = locs_ds.isel(time=np.arange(11,3011))

#generate dataset with only first variants in there
locs_single_ds = xr.DataArray(data=locs_ds.zos.values[np.arange(len(locs_ds.zos.model)),np.argmax(np.isfinite(locs_ds.zos[:,:,0,0]).values,axis=1),:,:],dims=['model','time','location'],coords=dict(model=locs_ds.model,time=locs_ds.time,location=locs_ds.location))
locs_single_ds = locs_single_ds.to_dataset(name='zos')       

locs_allvar_ds=locs_ds

bathy = xr.open_dataset('/Users/thermans/Documents/Modeling/ROMS/preprocessing/bathymetry/etopo1_bedrock_bigNorthSea1.nc')

#get seasonal anomaly changes
locs_allvar_djf_mean_anom, locs_allvar_jja_mean_anom, locs_allvar_mam_mean_anom, locs_allvar_son_mean_anom, locs_allvar_dec2nov_mean = seasonal_deviations_from_monthly_means(locs_allvar_ds.zos)

#change in seasonal anomalies relative to baseyears
locs_allvar_djf_mean_anom_d = locs_allvar_djf_mean_anom - locs_allvar_djf_mean_anom.sel(year=baseyears).mean(dim='year')
locs_allvar_jja_mean_anom_d = locs_allvar_jja_mean_anom - locs_allvar_jja_mean_anom.sel(year=baseyears).mean(dim='year')
locs_allvar_mam_mean_anom_d = locs_allvar_mam_mean_anom - locs_allvar_mam_mean_anom.sel(year=baseyears).mean(dim='year')
locs_allvar_son_mean_anom_d = locs_allvar_son_mean_anom - locs_allvar_son_mean_anom.sel(year=baseyears).mean(dim='year')

locs_single_djf_mean_anom, locs_single_jja_mean_anom, locs_single_mam_mean_anom, locs_single_son_mean_anom, locs_single_dec2nov_mean = seasonal_deviations_from_monthly_means(locs_single_ds.zos)

#change in seasonal anomalies relative to baseyears
locs_single_djf_mean_anom_d = locs_single_djf_mean_anom - locs_single_djf_mean_anom.sel(year=baseyears).mean(dim='year')
locs_single_jja_mean_anom_d = locs_single_jja_mean_anom - locs_single_jja_mean_anom.sel(year=baseyears).mean(dim='year')
locs_single_mam_mean_anom_d = locs_single_mam_mean_anom - locs_single_mam_mean_anom.sel(year=baseyears).mean(dim='year')
locs_single_son_mean_anom_d = locs_single_son_mean_anom - locs_single_son_mean_anom.sel(year=baseyears).mean(dim='year')

####plotting
discrete_cmap = matplotlib.colors.ListedColormap(cmocean.cm.balance((20,95,170,235)), name='my_colormap_name') #generate 4 level discrete colormap
discrete_cmap.set_bad('grey') #set land mask color

cmap = cmocean.cm.amp #continuous colormap
cmap.set_bad('grey')

fig=plt.figure(figsize=(6.5,11.5)) #generate figure  
gs = fig.add_gridspec(10,4)
gs.update(hspace = 0.5,wspace=0.3,right=.84,top=.95,bottom=.05)

#plot multi-model distributions in panel e
row = [1,1,2,2,3,3,4,4] #subplot arrangement

for l,loc in enumerate(locations): #loop over example locations
    if np.mod(l,2)==1:
       plt.subplot(gs[row[l],2:4])
    else:
       plt.subplot(gs[row[l],0:2])
    plt.grid()

    plt.errorbar(100*locs_allvar_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').mean(dim='model'),2+.25,xerr=100*locs_allvar_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').std(dim='model'),fmt='o',capsize=2,ecolor='black',mfc='black',mec='black',label='All realizations')
    plt.errorbar(100*locs_allvar_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').mean(dim='model'),1+.25,xerr=100*locs_allvar_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').std(dim='model'),fmt='o',capsize=2,ecolor='black',mfc='black',mec='black')
    plt.errorbar(100*locs_allvar_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').mean(dim='model'),4+.25,xerr=100*locs_allvar_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').std(dim='model'),fmt='o',capsize=2,ecolor='black',mfc='black',mec='black')
    plt.errorbar(100*locs_allvar_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').mean(dim='model'),3+.25,xerr=100*locs_allvar_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='variant').std(dim='model'),fmt='o',capsize=2,ecolor='black',mfc='black',mec='black')

    plt.errorbar(100*locs_single_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),2-.25,xerr=100*locs_single_jja_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').std(dim='model'),fmt='o',capsize=2,ecolor='red',mfc='red',mec='red',label='Single realization')
    plt.errorbar(100*locs_single_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),1-.25,xerr=100*locs_single_son_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').std(dim='model'),fmt='o',capsize=2,ecolor='red',mfc='red',mec='red')
    plt.errorbar(100*locs_single_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),4-.25,xerr=100*locs_single_djf_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').std(dim='model'),fmt='o',capsize=2,ecolor='red',mfc='red',mec='red')
    plt.errorbar(100*locs_single_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').mean(dim='model'),3-.25,xerr=100*locs_single_mam_mean_anom_d.sel(year=futyears,location=loc).mean(dim='year').std(dim='model'),fmt='o',capsize=2,ecolor='red',mfc='red',mec='red')

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
    if l==7:
        plt.legend(bbox_to_anchor=(.35, -.75))


#save figure
#fig.savefig(os.path.join(out_dir,'suppFigure7_single_multivar.pdf'),dpi=300)
