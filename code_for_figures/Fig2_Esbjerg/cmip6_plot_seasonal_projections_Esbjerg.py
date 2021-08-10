#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 08:43:46 2021

Plot dynamic sea-level projections at Esbjerg: a) seasonal & annual means,
b) change in SSLA

@author: thermans
"""

import xarray as xr
import numpy as np
import os
import fnmatch
import matplotlib
import matplotlib.pyplot as plt
import cmocean
plt.close('all')

baseyears = np.arange(1995,2015) #1995-2014
futyears = np.arange(2081,2101) #2081-2100

ssp = 'ssp585'
in_dir = '/Volumes/Naamloos/PhD_Data/CMIP6/ensemble_netcdfs/' #input directory ensemble data
out_dir = '/Users/thermans/Documents/PhD/Phase4_seasonal/Figures' #where to store the figure
model_list = list(np.genfromtxt('/Users/thermans/Documents/PhD/Phase4_seasonal/Analysis/compiling_ensembles/ens_model_list_'+ssp+'.txt',dtype='str'))

if ssp=='ssp126':
    zos_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp126_n36_nslocations_all_variants.nc")[0])) 
elif ssp=='ssp585':
    zos_ds = xr.open_dataset(os.path.join(in_dir,fnmatch.filter(os.listdir(in_dir), "zos_CMIP6_ssp585_n39_nslocations_all_variants.nc")[0]))  

zos_ds = zos_ds.sel(model=model_list,location='Esbjerg') #select Esbjerg
zos_ds = zos_ds.isel(time=np.arange(11,3011)) #select Dec->Nov
zos_ds['zos'] = zos_ds.zos*100 #convert from [m] to [cm]

#calculate seasonal means 
ordered_per_year = np.resize(zos_ds.zos,(len(zos_ds.model),len(zos_ds.variant),int(len(zos_ds.time)/12),12))
am_dec2nov = np.mean(ordered_per_year,axis=3)

ens_dec2nov_mean = xr.DataArray(am_dec2nov,coords=[zos_ds.model,zos_ds.variant,np.arange(1851,2101)],dims=['model','variant','year'])#zos_ds_am.copy(deep=True)
ens_djf_mean = xr.DataArray(np.mean(ordered_per_year[:,:,:,0:3],axis=3),coords=[zos_ds.model,zos_ds.variant,np.arange(1851,2101)],dims=['model','variant','year'])#zos_ds_am.copy(deep=True)
ens_mam_mean = xr.DataArray(np.mean(ordered_per_year[:,:,:,3:6],axis=3),coords=[zos_ds.model,zos_ds.variant,np.arange(1851,2101)],dims=['model','variant','year'])#zos_ds_am.copy(deep=True)
ens_jja_mean = xr.DataArray(np.mean(ordered_per_year[:,:,:,6:9],axis=3),coords=[zos_ds.model,zos_ds.variant,np.arange(1851,2101)],dims=['model','variant','year'])#zos_ds_am.copy(deep=True)
ens_son_mean = xr.DataArray(np.mean(ordered_per_year[:,:,:,9:12],axis=3),coords=[zos_ds.model,zos_ds.variant,np.arange(1851,2101)],dims=['model','variant','year'])#zos_ds_am.copy(deep=True)

#seasonal anomalies relative to annual mean
ens_djf_mean_anom = ens_djf_mean-ens_dec2nov_mean
ens_jja_mean_anom = ens_jja_mean-ens_dec2nov_mean
ens_mam_mean_anom = ens_mam_mean-ens_dec2nov_mean
ens_son_mean_anom = ens_son_mean-ens_dec2nov_mean

#relate to baseyears 1995-2014
ens_djf_mean_d = (ens_djf_mean - ens_djf_mean.sel(year=baseyears).mean(dim='year'))
ens_jja_mean_d = (ens_jja_mean - ens_jja_mean.sel(year=baseyears).mean(dim='year'))
ens_mam_mean_d = (ens_mam_mean - ens_mam_mean.sel(year=baseyears).mean(dim='year'))
ens_son_mean_d = (ens_son_mean - ens_son_mean.sel(year=baseyears).mean(dim='year'))

ens_djf_mean_anom_d = ens_djf_mean_anom - ens_djf_mean_anom.sel(year=baseyears).mean(dim='year')
ens_jja_mean_anom_d = ens_jja_mean_anom - ens_jja_mean_anom.sel(year=baseyears).mean(dim='year')
ens_mam_mean_anom_d = ens_mam_mean_anom - ens_mam_mean_anom.sel(year=baseyears).mean(dim='year')
ens_son_mean_anom_d = ens_son_mean_anom - ens_son_mean_anom.sel(year=baseyears).mean(dim='year')
ens_dec2nov_mean_d = ens_dec2nov_mean - ens_dec2nov_mean.sel(year=baseyears).mean(dim='year')
 

####plotting
my_cmap = matplotlib.colors.ListedColormap(cmocean.cm.balance((20,95,170,235)), name='my_colormap_name') #generate 4 level discrete colormap

fig=plt.figure(figsize=(6.5,4.5))
gs = fig.add_gridspec(2,3)
ax = fig.add_subplot(gs[0,0:3]) #0:2 to leave space for bars

#plot timeseries of absolute change
ax.axhline(0,color='grey',linestyle='dashed',lw=1)
ens_dec2nov_mean_d.mean(dim='variant',skipna=True).median(dim='model').plot(color='black',label='Annual mean',zorder=5)
ens_djf_mean_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(0),label='Winter (DJF)',zorder=4)
ens_mam_mean_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(2),label='Spring (MAM)',zorder=1)
ens_jja_mean_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(3),label='Summer (JJA)',zorder=2)
ens_son_mean_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(1),label='Autumn (SON)',zorder=3)

ax.axvline(2100,color='black',linestyle='solid',lw=.75)
ax.legend(ncol=2,loc='upper left')

#plot distribution during 2081-2100 on the right hand side of the plot:
ax.fill_between([2103,2104], ens_dec2nov_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_dec2nov_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'),
                facecolor='black',linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2102,2105], ens_dec2nov_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), 
                ens_dec2nov_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), 
                facecolor='black',linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2102,2105],[ens_dec2nov_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),
                     ens_dec2nov_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],
        color='white',linewidth=2,clip_on=False,zorder=5) 


ax.fill_between([2108,2109], ens_djf_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_djf_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(0),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2107,2110], ens_djf_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), ens_djf_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(0),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2107,2110],[ens_djf_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),ens_djf_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5) 

ax.fill_between([2112,2113], ens_mam_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_mam_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(2),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2111,2114], ens_mam_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), ens_mam_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(2),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2111,2114],[ens_mam_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),ens_mam_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5)

ax.fill_between([2116,2117], ens_jja_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_jja_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(3),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2115,2118], ens_jja_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), ens_jja_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(3),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2115,2118],[ens_jja_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),ens_jja_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5)

ax.fill_between([2120,2121], ens_son_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_son_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(1),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2119,2122], ens_son_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), ens_son_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(1),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2119,2122],[ens_son_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),ens_son_mean_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5)

ax.text(2103,36,'2081-2100')

#axes customization
ax.grid()
ax.set_xlim([1980,2125])
ax.set_ylim([-5,35])
ax.set_ylabel('SLC [cm]')
ax.get_yaxis().set_label_coords(-0.08,0.5)

ax.set_title('(a)')
ax.set_xticks(np.linspace(1980,2100,7))
ax.set_xticklabels('')
ax.set_xlabel('')

#repeat the above for changes in SSLA
gs = fig.add_gridspec(2,3)
ax = fig.add_subplot(gs[1,0:3])

ax.axhline(0,color='grey',linestyle='dashed',lw=1)


ens_djf_mean_anom_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(0),zorder=4)
ens_mam_mean_anom_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(2),zorder=1)
ens_jja_mean_anom_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(3),zorder=2)
ens_son_mean_anom_d.mean(dim='variant',skipna=True).median(dim='model').plot(color=my_cmap(1),zorder=3)

ax.axvline(2100,color='black',linestyle='solid',lw=.75)

ax.fill_between([2107,2110], ens_djf_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), 
                ens_djf_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(0),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2108,2109], ens_djf_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_djf_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(0),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2107,2110],[ens_djf_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),
                    ens_djf_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5) 

ax.fill_between([2111,2114], ens_mam_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), 
                ens_mam_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(2),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2112,2113], ens_mam_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_mam_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(2),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2111,2114],[ens_mam_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),
                    ens_mam_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5)

ax.fill_between([2115,2118], ens_jja_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'), 
                ens_jja_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(3),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2116,2117], ens_jja_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'), 
                ens_jja_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(3),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2115,2118],[ens_jja_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),
                ens_jja_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5)


ax.fill_between([2119,2122], ens_son_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.17,dim='model'),
                ens_son_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.83,dim='model'), facecolor=my_cmap(1),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.fill_between([2120,2121], ens_son_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.05,dim='model'),
                ens_son_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.95,dim='model'), facecolor=my_cmap(1),linewidth=1.0,clip_on=False,zorder=5) #plot bar
ax.plot([2119,2122],[ens_son_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model'),
                     ens_son_mean_anom_d.sel(year=futyears).mean(dim='variant',skipna=True).mean(dim='year').quantile(.5,dim='model')],color='white',linewidth=2,clip_on=False,zorder=5)

ax.plot(2103.5,0,color='black',marker='x',clip_on=False,zorder=5) 

ax.text(2103,10.5,'2081-2100')

ax.grid()
ax.set_xlim([1980,2125])
ax.set_ylim([-10,10])
ax.set_ylabel(r'$\Delta$' 'SSLA [cm]')
ax.get_yaxis().set_label_coords(-0.08,0.5)
ax.set_title('(b)')
ax.set_xlabel('Year')
ax.set_xticks(np.linspace(1980,2100,7))
ax.set_xticklabels(['1980','2000','2020','2040','2060','2080','2100'])
fig.tight_layout()

'''
if ssp=='ssp126':
    fig.savefig(os.path.join('suppFigure4_projections_Esbjerg_ssp126.pdf'),dpi=300)
elif ssp=='ssp585':
    fig.savefig(os.path.join('Figure2_projections_Esbjerg_ssp585.pdf'),dpi=300)
'''