# ESMplot

Welcome to ESMplot: the Earth System Model plotting package on Python! \
Author: Alex Thompson (ajthompson@wustl.edu) \
Date: 2022-12-06

# Installation
For current prototype version, download entire directory 'ESMplot' (size: ~42MB) from GitHub and place in your own working directory.

# Table of Contents

======================================<br/>
Example script: calculate_seasavg.py:<br/>
====================================== <br/>

 Example script for calculating a seasonal average from any number of netCDF files using the
  user-defined features of this package's functions.


======================================<br/>
Example script: calculate_seascyc.py:<br/>
====================================== <br/>

 Example script for calculating the seasonal cycle (line and map plots) from any number of
  netCDF files using the user-defined features of this package's functions.<br/>


==============================<br/>
directory 'climate_analysis':<br/>
============================== <br/>
 --------------<br/>
 climatology.py : list of functions for calculating climatologies from time series<br/>
 --------------<br/>
  function : clmMonTLL() - calculate 12 month climatology from time series with dimensions
                           (time x lat x lon)<br/>
  function : clmMonTLLL() - calculate 12 month climatology from time series with dimensions
                            (time x atm level x lat x lon)<br/>
  function : clmMonTSLL() - calculate 12 month climatology from time series with dimensions
                            (time x soil level x lat x lon)<br/>
