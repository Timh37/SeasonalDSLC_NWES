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
The typical process chain starting with raw monthly mean data:
1. Merge separate time-chunks ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/cmip6_merge_raw_timechunks.py))
2. (*For 'zos'*) Dedrift using a linear fit to piControl ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/cmip6_dedrift_linear.py))
3. (*For 'zos'*) Subtract area-weighted mean at each timestep ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/cmip6_subtract_areawmean_ocean.py))
4. (*Optional*) Regrid to a common grid, for example 1 by 1 degrees ([**Code**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/cmip6_processing/regridding/cmip6_regrid_to_common.py))

## List of figures in the manuscript and supplementary information
*The input data is available at the 4TU Research.Data Repository, DOI:* [To-do](http://github.com)

| Figure | Code | Input data | Brief description |
| ------------- |:-------------:| -----| -----|
| Fig. 1 | [**Link**](https://github.com/Timh37/SeasonalSLC_NWES/blob/main/code_for_figures/Fig1_dSSLA/cmip6_plot_dSSLA.py) | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, zos_CMIP6_ssp585_n39_nslocations_all_variants.nc, ens_model_list_ssp585.txt | ensemble mean dSSLA for SSP5-8.5, 2081-2100 relative to 1995-2014, maps; and multi-model distributions at 8 example coatal locations |
| Fig. 2 | link | zos_CMIP6_ssp585_n38_nwes_variant_averaged.nc, zos_CMIP6_ssp585_n39_nslocations_all_variants.nc, ens_model_list_ssp585.txt | ensemble mean dSSLA for SSP5-8.5, 2081-2100 relative to 1995-2014, maps; and multi-model distributions at 8 example coatal locations |
