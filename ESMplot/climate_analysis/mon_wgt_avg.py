#######################################################################################################
# 
# Functions for calculating a seasonal average that includes weighting by month 
#
#######################################################################################################

import numpy as np
import xarray as xr

# Default arrays
default_wgt_mon = xr.DataArray(np.ones(12),dims=['time']).astype(float)

#######################################################################################################

def mon_wgt_avg(var: xr.DataArray,
                months: list,
                wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in a climatology variable with dimensions (time: 12 months x lat x lon) and calculates
  a seasonal average that is weighted by proportional length of each included month.

  Required parameters:
  ----------------------------
  var,months

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  var: :class:'xarray.DataArray'
       A variable with dimensions (time: 12 months x lat x lon)             

  months: :class:'list'
          A list of indices that defines the months to be included in the seasonal average.
          Ex. Annual average = [0,1,2,3,4,5,6,7,8,9,10,11] or JJA average = [5,6,7]

  wgt_mon: :class:'xarray.DataArray', Optional
           Array of monthly weights (size: 12). This feature is available for paleoclimate analysis when
           the relative weight of each month differed from the values of today. 
           Default = all months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  # Weight by each month and calculate seasonal average
  var_wgt = var[months,:,:].weighted(wgt_mon)
  var_avg = var_wgt.mean('time')

  return var_avg

