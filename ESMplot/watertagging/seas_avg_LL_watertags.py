#######################################################################################################
# 
# These functions take isotopic inputs from netCDF output and return a seasonally averaged variable 
# (weighted by a corresponding climate variable if allowable) of size (lat x lon). The functions
# contained here are specifically used for water tagged CESM output. 
#
#######################################################################################################

import numpy as np
import xarray as xr
from ESMplot.climate_analysis.climatology import clmMonTLL
from ESMplot.climate_analysis.mon_wgt_avg import mon_wgt_avg

# Default arrays
default_wgt_mon = xr.DataArray(np.ones(12),dims=['time']).astype(float)

#######################################################################################################
# Precipitation isotopes: d18Op, dDp, or dexcess of precip 
#######################################################################################################

def seasavg_watertagging_vars(tagcode:str,
                              months: list,
                              path: str,
                              begi: str = 'beg',
                              endi: str = 'end',
                              ptiny: float = 1.E-18,
                              mult: float = 86400000.,
                              wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in precipitation isotopes and precpitation (here defined as the sum of 16O variables) from a 
  netCDF file (that includes a dimension for time) and creates a seasonally averaged global map for a
  specified tag region (lat x lon). Weighting by precipitation amount and fractional length of each 
  month occurs if parameters 'months' includes two or more values.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
  tagcode
   class: 'string', Code name of tagged region. This parameter is used for reading in water tag variables.  

  months
   class: 'list', A list of indices that defines the months to be included in the seasonal average.
                  Ex. Annual average = [0,1,2,3,4,5,6,7,8,9,10,11] or JJA average = [5,6,7]

  path
   class: 'string', The file path to the netCDF file containing the output variables to be read in. In 
                    addition to the variable being read in by var1, file at path must include hyam, hybm, 
                    and ps if var1 contains dimension 'lev'.

  Optional Parameters                       
  --------------------------------------------------------------------------------------------------------
  begi, endi: class 'string' or 'int', Determines the first and last time element to be included in 
                                       analysis. Defaults 'beg' and 'end' read in the entire length of
                                       the file. Specifying an integer will read in time from begi as int 
                                       to endi as int. If integers, (endi-begi)/12 must be multiple of 12.  
  ptiny: class 'float', Minimum value for light isotopic variable. This is required in the calculation of 
                        isotopes so as not to have artificially high delta values due to extremely low 
                        light isotopic values. Default = 1.E-18
  mult: class 'float', A float value that becomes a scalar to the variable's values. Assumes original 
                       variable is in m/s. Default = 86400000., which converts m/s to mm/day
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged variables for a specific tag: Pi [lat,lon], d18Opsink [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  tagcode = str(tagcode)
  months = list(months)
  path = str(path)
  ptiny = float(ptiny)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #-----------------------------------
  # Read in dataset 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Assign values for begi and endi
  begi = 0 if begi == 'beg' else begi
  endi = len(data.time) if endi == 'end' else endi

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of precipitation 
  #----------------------------------------------------------------------------------

  # Read in light isotopic values
  liso_d18O_init = xr.DataArray(data['PRECRC_'+tagcode+'r'][begi:endi,:,:] + \
                                data['PRECSC_'+tagcode+'s'][begi:endi,:,:] + \
                                data['PRECRL_'+tagcode+'R'][begi:endi,:,:] + \
                                data['PRECSL_'+tagcode+'S'][begi:endi,:,:])

  # Seasonal cycle of liso_d18O
  prect_clim = clmMonTLL(liso_d18O_init) 

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of d18Op
  #----------------------------------------------------------------------------------

  # Remove extremely small values from light isotopic variable
  liso_d18O = xr.where(liso_d18O_init < ptiny,ptiny,liso_d18O_init)

  # Read in heavy isotopic values
  hiso_d18O = xr.DataArray(data['PRECRC_'+tagcode+'18Or'][begi:endi,:,:] + \
                           data['PRECSC_'+tagcode+'18Os'][begi:endi,:,:] + \
                           data['PRECRL_'+tagcode+'18OR'][begi:endi,:,:] + \
                           data['PRECSL_'+tagcode+'18OS'][begi:endi,:,:])

  # Calculate delta values
  d18O = ( (hiso_d18O / liso_d18O) - 1. ) * 1000.

  # Perform seasonal averaging
  d18O_clim = clmMonTLL(d18O)

  #----------------------------------------------------------------------------------
  # Process variables before performing weighting
  #----------------------------------------------------------------------------------

  # Remove any negative or fill values from wtvar and set to NaN
  wtvar_proc = xr.where((prect_clim <= 0.) | (prect_clim > 1E29),float('NaN'),prect_clim)

  # Set weight variable (prect) as standardized xr.DataArray
  wtvar = xr.DataArray(wtvar_proc,dims=['time','lat','lon'])

  # Set iso variable (d18O) as standardized xr.DataArray
  iso = xr.DataArray(d18O_clim,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------------------------
  # Weight precipitation monthly values by multiplying them by fractional length of each month
  # NOTE: 'wgt_mon' variable will already account for how many months are present in the calculation,
  #       but 'wtvar' will not, so need to index 'months' for 'wtvar' here
  #----------------------------------------------------------------------------------------------------

  # result -> [months x lat x lon], precip weighted by month
  wtvar_mon_wgt_calc = wgt_mon*wtvar[months,:,:]

  # Set as standardized xr.DataArray
  wtvar_mon_wgt = xr.DataArray(wtvar_mon_wgt_calc,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------------------------
  # Weight isotopic values by seasonal rainfall
  # NOTE: 'wtvar_mon_wgt' variable already accounts for how many months are present, but 'iso' does not 
  #----------------------------------------------------------------------------------------------------

  d18Opwt_tag = xr.DataArray( xr.DataArray(wtvar_mon_wgt*iso[months,:,:]).sum(dim='time') / \
                              xr.DataArray(wtvar_mon_wgt).sum(dim='time')                   )

  #----------------------------------------------------------------------------------
  # Process final variables
  #----------------------------------------------------------------------------------

  # PRECT 
  Pi_monwgt = mon_wgt_avg(var=prect_clim,months=months,wgt_mon=wgt_mon)   
  Pi = xr.DataArray(Pi_monwgt,dims=['lat','lon'],coords=dict(lat=data.lat,lon=data.lon))*mult

  # d18Op
  d18Opsink = xr.DataArray(d18Opwt_tag,dims=['lat','lon'],coords=dict(lat=data.lat,lon=data.lon)) 

  return Pi, d18Opsink

