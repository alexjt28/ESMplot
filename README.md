# ESMplot

Welcome to ESMplot: the Earth System Model plotting package on Python! <br/>
Author: Alex Thompson (ajthompson@wustl.edu) <br/>
Date: 2022-12-06 <br/>
<br/>
This is the prototype of ESMplot, a Python package designed for flexible vizualization of Earth system model netCDF output. This ReadMe file displays the types of modules and scripts currently available in this package, with a brief description of each included.

# Installation
For current prototype version, download entire directory 'ESMplot' (size: ~42MB) from GitHub and place in your own working directory.

# Table of Contents

====================================== <br/>
Example script: **calculate_seasavg.py** <br/>
====================================== <br/>

Example script for calculating a seasonal average from any number of netCDF files using the user-defined features of this package's functions.

====================================== <br/>
Example script: **calculate_seascyc.py** <br/>
====================================== <br/>

Example script for calculating the seasonal cycle (line and map plots) from any number ofnetCDF files using the user-defined features of this package's functions. <br/>

============================== <br/>
directory **climate_analysis** <br/>
============================== <br/>

#---**climatology.py**--- : *list of functions for calculating climatologies from time series* <br/>

*function* : **clmMonTLL()** - calculate 12 month climatology from time series with dimensions (time x lat x lon) <br/>
*function* : **clmMonTLLL()** - calculate 12 month climatology from time series with dimensions (time x atm level x lat x lon) <br/>
*function* : **clmMonTSLL()** - calculate 12 month climatology from time series with dimensions (time x soil level x lat x lon) <br/>

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

#---**seas_cycle_TLL.py**--- : *list of functions for calculating seasonal cycle variable on a map with final dimensions (time = 12 months x lat x lon; hence TLL)* <br/>

*function* : **seascyc_var_TLL()** - calculates seasonal cycle global map of a specified variable that may or may not include an atmospheric level dimension <br/>
*function* : **seascyc_prect_TLL()** - calculates seasonal cycle global map of precipitation <br/>
*function* : **seascyc_soilvar_TLL()** - calculates seasonal cycle global map of variable that includes a soil level dimension <br/>
*function* : **seascyc_rainiso_TLL()** - calculates seasonal cycle global map of precipitation isotopes <br/>
*function* : **seascyc_soiliso_TLL()** - calculates seasonal cycle global map of soil water isotopes <br/>
*function* : **seascyc_vaporiso_TLL()** - calculates seasonal cycle global map of water vapor isotopes at a specified atmospheric level <br/>
*function* : **seascyc_isoroot_TLL()** - calculates seasonal cycle global map of soil water isotopes weighted by rooting depth fraction of vegetation <br/>
*function* : **seascyc_wind_vec_TLL()** - calculates seasonal cycle global map of U and V wind vectors <br/>

#---**seas_cycle_TSLL.py**--- : *list of functions for calculating seasonal cycle variable on a map with final dimensions (time = 12 months x soil level x lat x lon; hence TSLL)* <br/>

*function* : **seascyc_soilvar_TSLL()** - calculates seasonal cycle global map of a specified variable, final variable retains the soil level dimension <br/>
*function* : **seascyc_soiliso_TSLL()** - calculates seasonal cycle global map of soil water isotopes, final variable retains the soil level dimension <br/>

============================== <br/>
directory **plotting** <br/>
============================== <br/>

#---**plot_functions.py**--- : *list of functions that are called for making plots elsewhere in this package* <br/>

*function* : **save_multi_image()** - saves multiple images within the same PDF file as separate pages <br/>
*function* : **map_ticks_and_labels()** - plot map ticks and labels for latitude and longitude coordinates in a variety of flexible, user-defined ways <br/>

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
