#######################################################################################################
#
# Functions for calculating climatologies from global time series variables 
#
#######################################################################################################

import numpy as np
import xarray as xr

#######################################################################################################
# Calculate 12 month climatology variables 
#######################################################################################################

def clmMonTLL(var_time: xr.DataArray,
              time_dim: str = 'time'):
   return xr.DataArray(np.array(var_time).reshape(int(len(var_time[time_dim])/12),
                                                  12,
                                                  len(var_time.lat),
                                                  len(var_time.lon)          
                                                  ),
                       dims=['year',time_dim,'lat','lon'],
                       coords={time_dim: np.arange(12), 'lat': var_time.lat, 'lon': var_time.lon}
                       ).mean('year')

   ''' Designed to emulate https://www.ncl.ucar.edu/Document/Functions/Contributed/clmMonTLL.shtml
   Calculates a 12 month climatology variable from a time series variable with dimensions 
   (time x lat x lon), hence TLL. Order of dimensions in xarray input variable must be (time,lat,lon). 
   Returned variable will be of size (12 x lat x lon), with 12 corresponding to the mean value for
   each month (first element = Jan avg., etc.) from the original time series. Dimension for months
   is named 'time' for ease of future weighting the variable by time, but can be specified as a
   different name if desired (ex. 'month') '''


def clmMonTLLL(var_time: xr.DataArray,
              time_dim: str = 'time',
              lev_dim: str = 'plev'):
   return xr.DataArray(np.array(var_time).reshape(int(len(var_time[time_dim])/12),
                                                  12,
                                                  len(var_time[lev_dim]),
                                                  len(var_time.lat),
                                                  len(var_time.lon)
                                                  ),
                       dims=['year',time_dim,lev_dim,'lat','lon'],
                       coords={time_dim:np.arange(12),lev_dim:var_time.lev,'lat':var_time.lat,'lon':var_time.lon}
                       ).mean('year')

   ''' Designed to emulate https://www.ncl.ucar.edu/Document/Functions/Contributed/clmMonTLLL.shtml
   Identical to 'clmMonTLL' above but requires an antmospheric level dimension (here 'plev' which 
   indicates pressure levels) such that the input has dimensions (time,lev or plev,lat,lon). 
   Returned variable will be of size (12 x # of atmospheric levels x lat x lon), hence TLLL. '''

def clmMonTSLL(var_time: xr.DataArray,
               time_dim: str = 'time',
               soil_dim: str = 'levgrnd'):
   return xr.DataArray(np.array(var_time).reshape(int(len(var_time[time_dim])/12),
                                                  12,
                                                  len(var_time[soil_dim]),
                                                  len(var_time.lat),
                                                  len(var_time.lon)
                                                  ),
                       dims=['year',time_dim,soil_dim,'lat','lon'],
                       coords={time_dim:np.arange(12),soil_dim:var_time[soil_dim],'lat':var_time.lat,'lon':var_time.lon}
                       ).mean('year')

   ''' Designed to emulate https://www.ncl.ucar.edu/Document/Functions/Contributed/clmMonTLLL.shtml
   Identical to 'clmMonTLL' above but requires a 'lev' dimension such that the input has dimensions
   (time,lev,lat,lon). Returned variable will be of size (12 x lev x lat x lon), hence TLLL. '''

