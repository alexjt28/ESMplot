# ESMplot

Welcome to ESMplot: the Earth System Model plotting package on Python! <br/>
Author: Alex Thompson (ajthompson@wustl.edu) <br/>
Date: 2023-01-24 <br/>
<br/>
This is the prototype of ESMplot, a Python package designed for flexible visualization of Earth system model netCDF output. This ReadMe file displays the types of modules and scripts currently available in this package, with a brief description of each included.

# Installation
For current prototype version, download entire directory 'ESMplot' (size: ~92MB) from GitHub and place in your own working directory.

# Conda environment
Use /conda_envs/environment_ESMplot_py310.yml to create a conda environment with python 3.10 to run this package.<br/>

With 'environment_ESMplot_py310.yml' in your working directory, type "conda env create -f environment_ESMplot_py310.yml"<br/>

Then type "conda activate ESMplot_py3.10" to activate the conda environment.<br/>

NOTE: ESMplot currently requires Python version 3.10 to run properly.<br/>

# Table of Contents

====================================== <br/>
Example script: **calculate_seasavg.py** <br/>
====================================== <br/>

Example script for calculating a seasonal average from any number of netCDF files using the user-defined features of this package's functions.

====================================== <br/>
Example script: **calculate_seascyc.py** <br/>
====================================== <br/>

Example script for calculating the seasonal cycle (line and map plots) from any number of netCDF files using the user-defined features of this package's functions. <br/>

======================================= <br/>
Example script: **calculate_watertags.py** <br/>
======================================= <br/>

Example script for calculating seasonally-averaged precipitation and isotopic composition from water tagged climate model simulations. <br/>

============================== <br/>
directory **climate_analysis** <br/>
============================== <br/>

#---**climatology.py**--- : *list of functions for calculating climatologies from time series* <br/>

*function* : **clmMonTLL()** - calculate 12 month climatology from time series with dimensions (time x lat x lon) <br/>
*function* : **clmMonTLLL()** - calculate 12 month climatology from time series with dimensions (time x atm level x lat x lon) <br/>
*function* : **clmMonTSLL()** - calculate 12 month climatology from time series with dimensions (time x soil level x lat x lon) <br/>

#---**coordinate_functions.py**--- : *list of functions for processing and indexing coordinates* <br/>

*function* : **lat_lon_index_array()** - calculate index arrays for lat and Lon based on given value ranges <br/>

#---**mon_wgt_avg.py**--- : *list of functions for weighting 12 month climatology variables by the fractional length of each month* <br/>
  
*function* : **mon_wgt_avg()** - calculate seasonal average that is weighted by the proportional length of each included month <br/>

#---**seas_avg_LL.py**--- : *list of functions for calculating seasonally averaged values on a map with final dimensions (lat x lon; hence LL)* <br/>
  
*function* : **seasavg_var_LL()** - calculates seasonally averaged global map of a specified variable that may or may not include an atmospheric level dimension <br/>
*function* : **seasavg_prect_LL()** - calculates seasonally averaged global map of precipitation <br/>
*function* : **seasavg_soilvar_LL()** - calculates seasonally averaged global map of variable that includes a soil level dimension <br/>
*function* : **seasavg_rainiso_LL()** - calculates seasonally averaged global map of precipitation isotopes <br/>
*function* : **seasavg_soiliso_LL()** - calculates seasonally averaged global map of soil water isotopes <br/>
*function* : **seasavg_vaporiso_LL()** - calculates seasonally averaged global map of water vapor isotopes at a specified atmospheric level <br/>
*function* : **seasavg_isoroot_LL()** - calculates seasonally averaged global map of soil water isotopes weighted by rooting depth fraction of vegetation <br/>
*function* : **seasavg_wind_vec_LL()** - calculates seasonally averaged global map of U and V wind vectors <br/>
*function* : **seasavg_IVT_vec_LL()** - calculates seasonally averaged global map of U and V integrated vapor transport (IVT) vectors <br/>

#---**seas_cycle_TLL.py**--- : *list of functions for calculating seasonal cycle variable on a map with final dimensions (time = 12 months x lat x lon; hence TLL)* <br/>

*function* : **seascyc_var_TLL()** - calculates seasonal cycle global map of a specified variable that may or may not include an atmospheric level dimension <br/>
*function* : **seascyc_prect_TLL()** - calculates seasonal cycle global map of precipitation <br/>
*function* : **seascyc_soilvar_TLL()** - calculates seasonal cycle global map of variable that includes a soil level dimension <br/>
*function* : **seascyc_rainiso_TLL()** - calculates seasonal cycle global map of precipitation isotopes <br/>
*function* : **seascyc_soiliso_TLL()** - calculates seasonal cycle global map of soil water isotopes <br/>
*function* : **seascyc_vaporiso_TLL()** - calculates seasonal cycle global map of water vapor isotopes at a specified atmospheric level <br/>
*function* : **seascyc_isoroot_TLL()** - calculates seasonal cycle global map of soil water isotopes weighted by rooting depth fraction of vegetation <br/>
*function* : **seascyc_wind_vec_TLL()** - calculates seasonal cycle global map of U and V wind vectors <br/>
*function* : **seascyc_IVT_vec_TLL()** - calculates seasonal cycle global map of U and V integrated vapor transport (IVT) vectors <br/>

#---**seas_cycle_TSLL.py**--- : *list of functions for calculating seasonal cycle variable on a map with final dimensions (time = 12 months x soil level x lat x lon; hence TSLL)* <br/>

*function* : **seascyc_soilvar_TSLL()** - calculates seasonal cycle global map of a specified variable, final variable retains the soil level dimension <br/>
*function* : **seascyc_soiliso_TSLL()** - calculates seasonal cycle global map of soil water isotopes, final variable retains the soil level dimension <br/>

============================== <br/>
directory **plotting** <br/>
============================== <br/>

#---**plot_functions.py**--- : *list of functions that are called for making plots elsewhere in this package* <br/>

*function* : **save_multi_image()** - saves multiple images within the same PDF file as separate pages <br/>
*function* : **map_ticks_and_labels()** - plot map ticks and labels for latitude and longitude coordinates in a variety of flexible, user-defined ways <br/>
*function* : **draw_region_box()** - draw a region box on a map given latitude and longitude limits <br/>

#---**plot_map_avg_functions.py**--- : *functions for plotting contour maps of netCDF output* <br/>

*function* : **plot_contour_map_avg()** - makes panel plot with option for individual plots of seasonally averaged global variable, many user-defined options are available for custom plotting <br/>
*function* : **plot_diff_contour_map_avg()** - makes panel plot with option for individual plots of seasonally averaged global variable differences, many user-defined options are available for custom plotting <br/>

#---**plot_seascycle_functions.py**--- : *functions for plotting seasonal cycles of netCDF output in many possible ways* <br/>

*function* : **plot_seasonal_cycle()** - plot seasonal cycle with options for 1) line plot, 2) map plots of each month, and 3) animated GIF of map plot for each month <br/>

============================== <br/>
directory **print_values** <br/>
============================== <br/>

#---**print_spatial_average.py**--- : *functions for printing a spatial average from a global variable* <br/>

*function* : **print_global_average()** - prints spatially weighted global average of a global variable <br/>
*function* : **print_region_average()** - prints spatially weighted averaged over a specified region from a global variable <br/>
*function* : **print_point_average()** - prints average at a specified point (grid cell) from a global variable <br/>

============================== <br/>
directory **watertagging** <br/>
============================== <br/>

#---**print_watertag_values.py**--- : *functions for printing values for each water tagged region to the command line or to an Excel file* <br/>

*function* : **print_watertag_values()** - prints values for each tag region for a specified case to the main screen <br/>
*function* : **monthly_watertag_values_to_excel()** - prints monthly values for precipitation and d18Op for each tag region for all cases to an Excel file <br/>

#---**seas_avg_LL_watertags.py**--- : *functions for taking water tagged isotopic inputs and returning seasonally averaged variables* <br/>

*function* : **seasavg_watertagging_vars()** - creates seasonally averaged global map of precipitation/d18Op for each specified tag region <br/>

#---**tagged_regions.py**--- : *functions for defining/drawing tagged regions on a map* <br/>

*function* : **draw_land_tags()** - draw land tagged regions on a map <br/>
*function* : **draw_ocean_tags()** - draw ocean tagged regions on a map <br/>

#---**watertag_plots.py**--- : *functions for creating map plot visualizations of water tagged climate model simulations* <br/>

*function* : **watertagging_values_on_map()** - creates three map plots showing values of 1) precipitation, 2) precipitation percentage, and 3) d18Op for land and ocean water tag regions <br/>
*function* : **plot_tagged_precip_and_d18Op()** - creates a multipage PDF file with a precipitation and d18Op map for each tagged region <br/>

============================== <br/>
directory **conda_envs** <br/>
============================== <br/>

Contains 'environment_ESMplot_py310.yml' file for creating a conda environment from which to run ESMplot

============================== <br/>
directory **pdfs** <br/>
============================== <br/>

Directory to which output files are sent
