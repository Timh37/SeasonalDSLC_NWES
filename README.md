# SeasonalSLC_NWES
Python code to analyze future seasonal SLC on the NWES using CMIP6 data and high-resolution regional ocean model experiments.

Please cite
```
Publication
```
when using this repository for your own studies.

## Downloading CMIP6 data
Using [**wget scripts**](https://github.com/Timh37/SeasonalSLC_NWES/tree/main/cmip6_downloading), e.g.:

```
bash wget-20210723114729_ssp585_tauu.sh
```

Wget scripts can be generated using search URLs.

## Processing CMIP6 data
The typical process chain starting with raw monthly mean data, organized by variable by model:
1. Merge separate time-chunks ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/cmip6_merge_raw_timechunks.py))
2. (*For 'zos'*) Dedrift using a linear fit to piControl ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/cmip6_dedrift_linear.py))
3. (*For 'zos'*) Subtract area-weighted mean at each timestep ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/cmip6_subtract_areawmean_ocean.py))
4. (*Optional*) Regrid to a common grid, for example 1 by 1 degrees ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/regridding/cmip6_regrid_to_common.py))

## ROMS model experiments
Code to prepare and run the ROMS model are available from:
'''
Hermans, Tim; Le Bars, D. (Dewi); Katsman, C.A. (Caroline); Carolina M.L. Camargo; Gerkema, Theo; Calafat, F. M. (Francisco); et al. (2020): Model input and output accompanying Drivers of interannual sea-level variability on the Northwestern European Shelf. 4TU.ResearchData. Dataset. https://doi.org/10.4121/uuid:d9656541-ff40-45d0-8859-ac644b155dfb 
'''
Scripts to add dSWSpA from CMIP6 models to the ERA5-based wind-speed forcing used by [**Hermans et al. (2020) (JGRo)**](https://doi.org/10.1029/2020JC016325) are available here *to-do*.

## List of figures in the manuscript and supplementary information
*The input data is available at the 4TU Research.Data Repository, DOI:* [To-do](http://github.com)

| Figure | Code | Input data | Brief description |
| ------------- |:-------------:| -----| -----|
| Fig. 1 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig1_dSSLA/cmip6_plot_dSSLA.py) | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, zos_CMIP6_ssp585_n39_nslocations_all_variants.nc, ens_model_list_ssp585.txt | Ensemble mean dSSLA for SSP5-8.5, 2081-2100 relative to 1995-2014, maps; and multi-model distributions at 8 example coastal locations |
| Fig. 2 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig2_Esbjerg/cmip6_plot_seasonal_projections_Esbjerg.py) | zos_CMIP6_ssp585_n39_nslocations_all_variants.nc, ens_model_list_ssp585.txt | Probabilistic projections of seasonal and annual mean SLC, and of SSLA, near Esbjerg, for SSP5-8.5 |
| Fig. 3 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig3_dSWSA/cmip6_plot_dSWSA.py) | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, tauu_CMIP6_ssp585_n33_nwes_variant_averaged.nc, tauv_CMIP6_ssp585_n33_nwes_variant_averaged.nc, ens_model_list_ssp585.txt | Ensemble mean dSWSA for SSP5-8.5, 2081-2100 relative to 1995-2014, maps |
| Fig. 4 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig4_composites_im/cmip6_plot_dSWSA_dSSLA_native_im.py) | dzos_Omon_CanESM5_ssp585_r(1-25)i1p1f1_gn_DJF(JJA)_anom.nc, dtauu(v)_Amon_CanESM5_ssp585_r(1-25)i1p1f1_gn_DJF(JJA).nc, dzos_Omon_UKESM1-0-LL_ssp585_r(1-4,8)i1p1f2_gn_DJF(JJA)_anom.nc, dtauu(v)_Amon_UKESM1-0-LL_ssp585_r(1-4,8)i1p1f2_gn_DJF(JJA).nc, dzos_Omon_IPSL-CM6A-LR_ssp585_r(1-4,6,14)i1p1f1_gn_DJF(JJA)_anom.nc, dtauu(v)_Amon_IPSL-CM6A-LR_ssp585_r(1-4,6,14)i1p1f1_gr_DJF(JJA).nc, zos_CMIP6_ssp585_n39_nslocations_all_variants.nc | Composite plots of dSSLA and dSWSA in 3 individual CMIP6 models and timeseries near Esbjerg in these models |
| Fig. 5 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig5_composites_roms/plot_cmip6_roms_experiments.py) | NorthSea8_avg_timeseries_monthly.nc *to-do* | Sea level and barotropic currents response to winter and summer dSWSpA from 3 individual CMIP6 models in ROMS |
| Fig. 6 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig6_sw_ne/plot_sw_ne_roms_experiments.py) | NorthSea8_avg_timeseries_monthly.nc *to-do* | Sea level and barotropic currents response to uniform southwesterly and northeasterly wind-stress change in winter and summer in ROMS |
| Fig. 7 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig7_sw_ne_cc/plot_sw_ne_roms_experiments_channelclosed.py) | NorthSea8_avg_timeseries_monthly.nc *to-do* | Sea level and barotropic currents response to uniform southwesterly and northeasterly wind-stress change in winter and summer in ROMS, with a closed English Channel |
| Fig. 8 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig8_access/plot_ACCESS-ESM1-5_roms_experiments.py) | dzos_Amon_ACCESS-ESM1-5_ssp585_r(1-10)i1p1f1_gn_DJF_anom.nc, NorthSea8_avg_timeseries_monthly.nc *to-do* | Winter dSSLA in ACCESS-ESM1-5, and sea level and barotropic currents response to winter dSWSpA from ACCESS-ESM1-5 in ROMS, with a closed and an open English Channel |
| Suppl. Fig. 1 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/suppFig1_annualmean/cmip6_plot_annualmean_change.py) | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, tauu_CMIP6_ssp585_n33_nwes_variant_averaged.nc, tauv_CMIP6_ssp585_n33_nwes_variant_averaged.nc, ens_model_list_ssp585.txt | Ensemble mean dynamic SLC and wind-stress change for SSP5-8.5, 2081-2100 relative to 1995-2014, maps |
| Suppl. Fig. 2 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/suppFig2_hist_fut_SSLA/cmip6_plot_SSLA.py) | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, zos_CMIP6_ssp585_n39_nslocations_all_variants.nc, ens_model_list_ssp585.txt | Maps of historical and future SSLA for SSP5-8.5, and at 8 example locations |
| Suppl. Fig. 3 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/suppFig3_dSSLA_SSP126/cmip6_plot_dSSLA.py) | zos_CMIP6_ssp126_n35_nwes_variant_averaged.nc, zos_CMIP6_ssp126_n36_nslocations_all_variants.nc, ens_model_list_ssp126.txt | Fig. 1 for SSP1-2.6 |
| Suppl. Fig. 4 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/suppFig4_Esbjerg_SSP126/cmip6_plot_seasonal_projections_Esbjerg.py) | zos_CMIP6_ssp126_n36_nslocations_all_variants.nc, ens_model_list_ssp126.txt | Fig. 2 for SSP1-2.6 |
| Suppl. Fig. 5 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/suppFig5_hist_fut_SWSA/cmip6_plot_SWSA.pyy) | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, tauu_CMIP6_ssp585_n33_nwes_variant_averaged.nc, tauv_CMIP6_ssp585_n33_nwes_variant_averaged.nc, ens_model_list_ssp585.txt | Maps of historical and future SWSA for SSP5-8.5 |
| Suppl. Fig. 6 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/suppFig6_obs/plot_ncep_aviso_SWSA_SSLA.py) | dt_global_allsat_msla_h_y1993m01-2019m04.nc, uwnd.mon.mean.nc, vwnd.mon.mean.nc | Composite plots of SSLA (Aviso) and SWSA (derived from NCEP wind speed) in winter 2007 and summer 1995 |
