#######################################################################################################
# 
# These functions take isotopic inputs from netCDF output and return a seasonally averaged variable 
# (weighted by a corresponding climate variable if allowable) of size (lat x lon). 
#
#######################################################################################################

import numpy as np
import xarray as xr
from geocat.comp.interpolation import interp_hybrid_to_pressure
from ESMplot.climate_analysis.climatology import clmMonTSLL
from ESMplot.climate_analysis.mon_wgt_avg import mon_wgt_avg
from ESMplot.climate_analysis import seas_cycle_TLL as seascyc
from ESMplot.climate_analysis import seas_cycle_TSLL as seassoil

# Default arrays
default_levels  = np.arange(0,1050,50) # 0 hPa to 1000 hPa by 50 
default_wgt_mon = xr.DataArray(np.ones(12),dims=['time']).astype(float)

#======================================================================================================
# File contains the following functions:
#
# Contours (output: lat x lon)
# seasavg_var_LL, seasavg_prect_LL, seasavg_soilvar_LL, seasavg_rainiso_LL, seasavg_soiliso_LL,
# seasavg_vaporiso_LL, seasavg_isoroot_LL    
#
# Vectors (output: lat x lon)
# seasavg_wind_vec_LL, seasavg_IVT_vec_LL
#
#======================================================================================================

#######################################################################################################
# Generic variable at surface or pressure level of the atmosphere 
#######################################################################################################

def seasavg_var_LL(var1: str,
                   months: list,
                   path: str,
                   begi: str = 'beg',
                   endi: str = 'end',
                   var2: str = None,
                   math: str = None,
                   level: int = None,
                   plev: np.ndarray = default_levels,
                   mult: float = 1.,
                   wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in up to two variables from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) output variable.

  Required Parameters:
  --------------------------------------------------------------------------------------------------------
  var1 
   class: 'string', A string specifying the variable to be read in from the netCDF file.
                    Ex. 'FSNT' (CESM output variable for net shortwave flux at top of model, as found at 
                    https://www.cesm.ucar.edu/models/cesm2/atmosphere/docs/ug6/hist_flds_f2000.html)

  months
   class: 'list', A list of indices that defines the months to be included in the seasonal average.
                  Ex. Annual average = [0,1,2,3,4,5,6,7,8,9,10,11] or JJA average = [5,6,7]

  path
   class: 'string', The file path to the netCDF file containing the output variables to be read in. In 
                    addition to the variable being read in by var1, file at path must include hyam, hybm, 
                    and ps if var1 contains dimension 'lev'.

  Optional Parameters:
  --------------------------------------------------------------------------------------------------------
  begi, endi: class 'string' or 'int', Determines the first and last time element to be included in 
                                       analysis. Defaults 'beg' and 'end' read in the entire length of
                                       the file. Specifying an integer will read in time from begi as int 
                                       to endi as int. If integers, (endi-begi)/12 must be multiple of 12.  
  var2: class 'string', A string for a second variable to read in from the netCDF file. If 'var2' is 
                        called, 'math' must also be specified. Ex. 'FLNT' (CESM output variable for net 
                        longwave flux at top of model)
  math: class 'string', Specified mathematical operation to be applied to 'var1' and 'var2'. Choices are: 
                        'add' (var1+var2), 'sub' (var1-var2), 'mul' (var1*var2), 'div' (var1/var2), and 
                        'exp' (var1**var2).
  level: class 'int', The integer value with units hPa/mb of the atmospheric level chosen from which to 
                      calculate the variable's values. Requires plev to be specified for the range of 
                      pressure levels. Ex. 850 -> calculating variable's values at 850 hPa 
  plev: class 'numpy.ndarray', The new range of pressure levels with units hPa/mb, of which 'level' must 
                               be an available value. Default = np.arange(0,1050,50) 
                               This script will interpolate climate model output from hybrid to pressure 
                               levels using geocat.comp's interp_hybrid_to_pressure function. NOTE: If 
                               NaN's present in variable while interpolating from hybrid to pressure, may 
                               see this error, which comes from metpy and can be ignored: 'UserWarning: 
                               Interpolation point out of data bounds encountered 
                               warnings.warn('Interpolation point out of data bounds encountered')'
                               Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with 
                               stride of 50 hPa. Default = np.arange(0,1050,50) 
  mult: class 'float', A float value that becomes a scalar to the variable's values. Default = 1.0. Will 
                       revert to default if not specified. 
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  var1 = str(var1)
  months = list(months)
  path = str(path)
  plev = np.array(plev)
  mult = float(mult)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time']) 
  if var2 != None:
   var2 = str(var2)
  if math != None:
   math = str(math)
  if level != None:
   level = int(level)
  
  #-----------------------------
  # Calculate seasonal cycle
  #-----------------------------

  var_clim = seascyc.seascyc_var_TLL(var1=var1,path=path,begi=begi,endi=endi,level=level,plev=plev,
                                     mult=mult,var2=var2,math=math)

  #--------------------------------------------
  # Weight by fractional length of each month
  #--------------------------------------------

  var_avg = mon_wgt_avg(var=var_clim,months=months,wgt_mon=wgt_mon)

  return var_avg

#######################################################################################################
# Precipitation variable 
#######################################################################################################

def seasavg_prect_LL(months: list,
                     path: str,
                     begi: str = 'beg',
                     endi: str = 'end',
                     mult: float = 86400000.,
                     wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in precipitation variables from a netCDF file (that includes a dimension for time) and 
  creates a seasonally averaged global map (lat x lon) output variable.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
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
  mult: class 'float', A float value that becomes a scalar to the variable's values. Assumes original 
                       variable is in m/s. Default = 86400000., which converts m/s to mm/day
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  ________________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  months = list(months)
  path = str(path)
  mult = float(mult)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #-----------------------------
  # Calculate seasonal cycle
  #-----------------------------

  var_clim = seascyc.seascyc_prect_TLL(path=path,begi=begi,endi=endi,mult=mult)

  #--------------------------------------------
  # Weight by fractional length of each month
  #--------------------------------------------

  var_avg = mon_wgt_avg(var=var_clim,months=months,wgt_mon=wgt_mon)

  return var_avg

#######################################################################################################
# Precipitation-Evaporation variable (Effective moisture)
#######################################################################################################

def seasavg_PminE_LL(months: list,
                     path: str,
                     begi: str = 'beg',
                     endi: str = 'end',
                     Pmult: float = 86400000.,
                     Emult: float = 86400.,
                     wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in precipitation and evaporation variables from a netCDF file (that includes a dimension 
  for time) and creates a seasonally averaged global map (lat x lon) output variable of
  Precipitation minus (-) Evaporation, otherwise known as effective moisture.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
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
  Pmult: class 'float', A float value that becomes a scalar to precipitation values. Assumes original 
                        variable is in m/s. Default = 86400000., which converts m/s to mm/day. 
  Emult: class 'float', A float value that becomes a scalar to evaporation values. Assumes original 
                        variable is in m/s. Default = 86400, which converts m/s to mm/day.
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  ________________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  months = list(months)
  path = str(path)
  Pmult = float(Pmult)
  Emult = float(Emult)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #-----------------------------
  # Calculate seasonal cycle
  #-----------------------------

  var_clim = seascyc.seascyc_PminE_TLL(path=path,begi=begi,endi=endi,Pmult=Pmult,Emult=Emult)

  #--------------------------------------------
  # Weight by fractional length of each month
  #--------------------------------------------

  var_avg = mon_wgt_avg(var=var_clim,months=months,wgt_mon=wgt_mon)

  return var_avg

#######################################################################################################
# Land variable with soil levels 
#######################################################################################################

def seasavg_soilvar_LL(var1: str,
                       months: list,
                       path: str,
                       begi: str = 'beg',
                       endi: str = 'end',
                       var2: str = None,
                       math: str = None,
                       soillev: list = [0],
                       soildimname: str = 'levgrnd',
                       mult: float = 1.,
                       wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in up to two land variables from a netCDF file (that includes a dimension for time and soil
  level) and creates a seasonally averaged global map (lat x lon) output variable.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
  var1
   class: 'string', A string of the variable to read in from the netCDF file. Ex. 'FSNT' (CESM output 
                    variable for net shortwave flux at top of model, as found at 
                    https://www.cesm.ucar.edu/models/cesm2/atmosphere/docs/ug6/hist_flds_f2000.html)

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
  var2: class 'string', A string for a second variable to read in from the netCDF file. If 'var2' is 
                        called, 'math' must also be specified. Ex. 'FLNT' (CESM output variable for net 
                        longwave flux at top of model)
  math: class 'string', Specified mathematical operation to be applied to 'var1' and 'var2'. Choices are: 
                        'add' (var1+var2), 'sub' (var1-var2), 'mul' (var1*var2), 'div' (var1/var2), and 
                        'exp' (var1**var2).
  soillev: class 'list', List of integers corresponding to the indices of the soil level dimension (ex. 
                         levgrnd) over which an average will be calculated. If not specified, first soil 
                         level will be indexed. Ex. [0,1,2] calculates the average of the first three soil 
                         level indices.
  soildimname: class 'str', Name of soil level dimension. Default is 'levgrnd', which is standard for CESM 
                            variables. This must correspond to the correct dimension name in var for proper 
                            averaging to occur.
  mult: class 'float', A float value that becomes a scalar to the variable's values. Default = 1.0. 
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  ________________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  var1 = str(var1)
  months = list(months)
  path = str(path)
  soillev = list(soillev)
  soildimname = str(soildimname)
  mult = float(mult)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])
  if var2 != None:
   var2 = str(var2)
  if math != None:
   math = str(math)

  #-----------------------------
  # Calculate seasonal cycle
  #-----------------------------

  var_clim = seascyc.seascyc_soilvar_TLL(var1=var1,path=path,begi=begi,endi=endi,soillev=soillev,
                                          mult=mult,var2=var2,math=math)

  #--------------------------------------------
  # Weight by fractional length of each month
  #--------------------------------------------

  var_avg = mon_wgt_avg(var=var_clim,months=months,wgt_mon=wgt_mon)

  return var_avg

#######################################################################################################
# Precipitation isotopes: d18Op, dDp, or dexcess of precip 
#######################################################################################################

def seasavg_rainiso_LL(iso_type:str,
                       months: list,
                       path: str,
                       begi: str = 'beg',
                       endi: str = 'end',
                       ptiny: float = 1.E-18,
                       wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in precipitation isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) output variable. Weighting by precipitation amount and 
  fractional length of each month occurs if parameters 'months' includes two or more values.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
  iso_type 
   class: 'string', Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

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
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  ________________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  iso_type = str(iso_type)
  months = list(months)
  path = str(path)
  ptiny = float(ptiny)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time']) 

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of precipitation 
  #----------------------------------------------------------------------------------

  # Seasonal cycle of precipitation (months x lat x lon)
  wtvar_init = seascyc.seascyc_prect_TLL(path=path,begi=begi,endi=endi) # no need to specify 'mult' here

  # Remove any negative or fill values from wtvar and set to NaN
  wtvar_proc = xr.where((wtvar_init <= 0.) | (wtvar_init > 1E29),float('NaN'),wtvar_init)

  # Set as standardized xr.DataArray
  wtvar = xr.DataArray(wtvar_proc,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of unweighted precipitation isotopes
  #----------------------------------------------------------------------------------

  # Calculation
  iso_init = seascyc.seascyc_rainiso_TLL(iso_type=iso_type,path=path,begi=begi,endi=endi,ptiny=ptiny)

  # Set as standardized xr.DataArray
  iso = xr.DataArray(iso_init,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------------------------
  # Weight precipitation monthly values by multiplying them by fractional length of each month
  # NOTE: 'wgt_mon' variable will already account for how many months are present in the calculation,
  #       but 'wtvar' will not, so need to index 'months' for 'wtvar' here
  #----------------------------------------------------------------------------------------------------
  
  # result -> [months x lat x lon], precip weighted by month
  wtvar_mon_wgt_calc = wgt_mon*wtvar[months,:,:]  
 
  # Set as standardized xr.DataArray
  wtvar_mon_wgt = xr.DataArray(wtvar_mon_wgt_calc,dims=['time','lat','lon'])
 
  #---------------------------------------------------------------------------------------------------------
  # Weight isotopic values by seasonal rainfall
  # NOTE: 'wtvar_mon_wgt' variable will already account for how many months are present, but 'iso' will not 
  #---------------------------------------------------------------------------------------------------------

  iso_wt = xr.DataArray( xr.DataArray(wtvar_mon_wgt*iso[months,:,:]).sum(dim='time') / \
                         xr.DataArray(wtvar_mon_wgt).sum(dim='time')                   ) 

  return iso_wt

#######################################################################################################
# Soil water isotopes: d18Op, dDp, or dexcess of soil water 
#######################################################################################################

def seasavg_soiliso_LL(iso_type:str,
                       months: list,
                       path: str,
                       begi: str = 'beg',
                       endi: str = 'end',
                       ptiny: float = 1.E-18,
                       soillev: list = [0],
                       soildimname: str = 'levgrnd',
                       wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in soil water isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) output variable. Weighting by soil water amount and 
  fractional length of each month occurs if parameter 'months' includes two or more values.

  Required Parameters
  --------------------------------------------------------------------------------------------------------
  iso_type
   class 'string', Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

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
  soillev: class 'list', List of integers corresponding to the indices of the soil level dimension (ex. 
                         levgrnd) over which an average will be calculated. If not specified, first soil 
                         level will be indexed. Ex. [0,1,2] calculates the average of the first three soil 
                         level indices.
  soildimname: class 'str', Name of soil level dimension. Default is 'levgrnd', which is standard for CESM 
                            variables. This must correspond to the correct dimension name in var for 
                            proper averaging to occur.
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  ________________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-----------------------------
  # Ensure type of parameters
  #-----------------------------

  iso_type = str(iso_type)
  months = list(months)
  path = str(path)
  ptiny = float(ptiny)
  soillev = list(soillev)
  soildimname = str(soildimname)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of soil water 
  #----------------------------------------------------------------------------------

  # Seasonal cycle of soil water (months x lat x lon)
  wtvar_init = seascyc.seascyc_soilvar_TLL(var1='H2OSOI',path=path,begi=begi,endi=endi,
                                           soillev=soillev)

  # Remove any negative or fill values from wtvar and set to NaN
  wtvar_proc = xr.where((wtvar_init <= 0.) | (wtvar_init > 1E29),float('NaN'),wtvar_init)

  # Set as standardized xr.DataArray
  wtvar = xr.DataArray(wtvar_proc,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of unweighted soil water isotopes
  #----------------------------------------------------------------------------------

  # Calculation
  iso_init = seascyc.seascyc_soiliso_TLL(iso_type=iso_type,path=path,begi=begi,endi=endi,ptiny=ptiny,
                                         soillev=soillev)

  # Set as standardized xr.DataArray
  iso = xr.DataArray(iso_init,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------------------------
  # Weight soil water monthly values by multiplying them by fractional length of each month
  # NOTE: 'wgt_mon' variable will already account for how many months are present in the calculation,
  #       but 'wtvar' will not, so need to index 'months' for 'wtvar' here
  #----------------------------------------------------------------------------------------------------

  # result -> [months x lat x lon], soil water weighted by month
  wtvar_mon_wgt_calc = wgt_mon*wtvar[months,:,:]

  # Set as standardized xr.DataArray
  wtvar_mon_wgt = xr.DataArray(wtvar_mon_wgt_calc,dims=['time','lat','lon'])

  #---------------------------------------------------------------------------------------------------------
  # Weight isotopic values by seasonal soil water amount
  # NOTE: 'wtvar_mon_wgt' variable will already account for how many months are present, but 'iso' will not 
  #---------------------------------------------------------------------------------------------------------

  iso_wt = xr.DataArray( xr.DataArray(wtvar_mon_wgt*iso[months,:,:]).sum(dim='time') / \
                         xr.DataArray(wtvar_mon_wgt).sum(dim='time')                   )

  return iso_wt

#######################################################################################################
# Vapor isotopes at specified pressure level: d18Op, dDp, or dexcess of water vapor
#######################################################################################################

def seasavg_vaporiso_LL(iso_type: str,
                        months: list,
                        path: str,
                        level: int,
                        begi: str = 'beg',
                        endi: str = 'end',
                        plev: np.ndarray = default_levels,
                        wtvar_name: str = 'Q',
                        ptiny: float = 1.E-18,
                        wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in water vapor isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) output variable. Weighting by specific humidity and 
  fractional length of each month occurs if parameter 'months' includes two or more values.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
  iso_type
   class 'string', Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

  months
   class: 'list', A list of indices that defines the months to be included in the seasonal average.
                  Ex. Annual average = [0,1,2,3,4,5,6,7,8,9,10,11] or JJA average = [5,6,7]

  path
   class: 'string', The file path to the netCDF file containing the output variables to be read in. In 
                    addition to the variable being read in by var1, file at path must include hyam, hybm, 
                    and ps if var1 contains dimension 'lev'.

  level: class 'int', The integer value with units hPa/mb of the atmospheric level chosen from which to 
                      calculate the variable's values. Requires plev to be specified for the range of 
                      pressure levels. Ex. 850 -> calculating variable's values at 850 hPa 

  Optional Parameters                       
  --------------------------------------------------------------------------------------------------------
  begi, endi: class 'string' or 'int', Determines the first and last time element to be included in 
                                       analysis. Defaults 'beg' and 'end' read in the entire length of
                                       the file. Specifying an integer will read in time from begi as int 
                                       to endi as int. If integers, (endi-begi)/12 must be multiple of 12.  
  plev: class 'numpy.ndarray', The new range of pressure levels with units hPa/mb, of which 'level' must 
                               be an available value. Default = np.arange(0,1050,50) 
                               This script will interpolate climate model output from hybrid to pressure 
                               levels using geocat.comp's interp_hybrid_to_pressure function. NOTE: If 
                               NaN's present in variable while interpolating from hybrid to pressure, may 
                               see this error, which comes from metpy and can be ignored: 'UserWarning: 
                               Interpolation point out of data bounds encountered 
                               warnings.warn('Interpolation point out of data bounds encountered')'
                               Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with 
                               stride of 50 hPa. Default = np.arange(0,1050,50) 
  wtvar_name: class 'string', Variable to weight the isotopic values by. Default is 'Q', which for CESM is 
                              specific humidity. 'Q' is the typical variable by which to weight vapor 
                              isotopes. 
  ptiny: class 'float', Minimum value for light isotopic variable. This is required in the calculation of 
                        isotopes so as not to have artificially high delta values due to extremely low 
                        light isotopic values. Default = 1.E-18
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________
  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  iso_type = str(iso_type)
  months = list(months)
  path = str(path)
  level = int(level)
  plev = np.array(plev)
  wtvar_name = str(wtvar_name)
  ptiny = float(ptiny)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of variable to weight by (ex. Q: specific humidty) 
  #----------------------------------------------------------------------------------

  # Seasonal cycle of wtvar (months x lat x lon)
  wtvar_init = seascyc.seascyc_var_TLL(var1=wtvar_name,path=path,begi=begi,endi=endi,
                                                     level=level,plev=plev)

  # Remove any negative or fill values from wtvar and set to NaN
  wtvar_proc = xr.where((wtvar_init <= 0.) | (wtvar_init > 1E29),float('NaN'),wtvar_init)

  # Set as standardized xr.DataArray
  wtvar = xr.DataArray(wtvar_proc,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of unweighted water vapor isotopes   
  #----------------------------------------------------------------------------------

  # Calculation
  iso_init = seascyc.seascyc_vaporiso_TLL(iso_type=iso_type,path=path,begi=begi,endi=endi,ptiny=ptiny,
                                            level=level,plev=plev)

  # Set as standardized xr.DataArray
  iso = xr.DataArray(iso_init,dims=['time','lat','lon'])

  #----------------------------------------------------------------------------------------------------
  # Weight wtvar monthly values by multiplying them by fractional length of each month
  # NOTE: 'wgt_mon' variable will already account for how many months are present in the calculation,
  #       but 'wtvar' will not, so need to index 'months' for 'wtvar' here
  #----------------------------------------------------------------------------------------------------

  # result -> [months x lat x lon], wtvar weighted by month
  wtvar_mon_wgt_calc = wgt_mon*wtvar[months,:,:]

  # Set as standardized xr.DataArray
  wtvar_mon_wgt = xr.DataArray(wtvar_mon_wgt_calc,dims=['time','lat','lon'])

  #---------------------------------------------------------------------------------------------------------
  # Weight isotopic values by seasonal wtvar (ex. Q)
  # NOTE: 'wtvar_mon_wgt' variable will already account for how many months are present, but 'iso' will not 
  #---------------------------------------------------------------------------------------------------------

  iso_wt = xr.DataArray( xr.DataArray(wtvar_mon_wgt*iso[months,:,:]).sum(dim='time') / \
                         xr.DataArray(wtvar_mon_wgt).sum(dim='time')                   )

  return iso_wt

#######################################################################################################
# Soil water isotopes weighted by root depth fraction: d18Op, dDp, or dexcess of soil water 
#######################################################################################################

def seasavg_isoroot_LL(iso_type:str,
                       months: list,
                       path: str,
                       begi: str = 'beg',
                       endi: str = 'end', 
                       ptiny: float = 1.E-18,
                       rootwgt: str = 'ROOTR_COLUMN',
                       soildimname: str = 'levgrnd',
                       wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in soil water isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) output variable weighted by root depth fraction. 
  Weighting by soil water amount and fractional length of each month occurs if parameter 'months' 
  includes two or more values. This variable is used for comparison to leaf wax proxy records.

  Required Parameters
  --------------------------------------------------------------------------------------------------------
  iso_type
   class 'string', Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

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
  rootwgt: class 'str', Name for variable used for weighting isotopic values by root depth fraction. In 
                        CESM, possible options are 'ROOTR_COLUMN', 'ROOTFR', or 'ROOTR'. Default = 
                        'ROOTR_COLUMN'.
  soildimname: class 'str', Name of soil level dimension. Default is 'levgrnd', which is standard for CESM 
                            variables. This must correspond to the correct dimension name in var for 
                            proper averaging to occur.
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable [lat,lon]

  '''

  #------------------------------
  # Ensure type of parameters
  #------------------------------

  #------------------------------
  # Ensure type of parameters
  #------------------------------
  
  iso_type = str(iso_type)
  months = list(months)
  path = str(path)
  ptiny = float(ptiny)
  rootwgt = str(rootwgt)
  soildimname = str(soildimname)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of soil water 
  #----------------------------------------------------------------------------------

  # Seasonal cycle of soil water (months x lat x lon)
  wtvar_init = seassoil.seascyc_soilvar_TSLL(var1='H2OSOI',path=path,begi=begi,endi=endi)

  # Remove any negative or fill values from wtvar and set to NaN
  wtvar_proc = xr.where((wtvar_init <= 0.) | (wtvar_init > 1E29),float('NaN'),wtvar_init)

  # Set as standardized xr.DataArray
  wtvar = xr.DataArray(wtvar_proc,dims=['time',soildimname,'lat','lon'])

  #----------------------------------------------------------------------------------
  # Calculate seasonal cycle of unweighted soil water isotopes
  #----------------------------------------------------------------------------------

  # Calculation
  iso_init = seassoil.seascyc_soiliso_TSLL(iso_type=iso_type,path=path,begi=begi,endi=endi,ptiny=ptiny)

  # Set as standardized xr.DataArray
  iso = xr.DataArray(iso_init,dims=['time',soildimname,'lat','lon'])

  #----------------------------------------------------------------------------------------------------
  # Weight soil water monthly values by multiplying them by fractional length of each month
  # NOTE: 'wgt_mon' variable will already account for how many months are present in the calculation,
  #       but 'wtvar' will not, so need to index 'months' for 'wtvar' here
  #----------------------------------------------------------------------------------------------------

  # result -> [months x lat x lon], soil water weighted by month
  wtvar_mon_wgt_calc = wgt_mon*wtvar[months,:,:,:]

  # Set as standardized xr.DataArray
  wtvar_mon_wgt = xr.DataArray(wtvar_mon_wgt_calc,dims=['time',soildimname,'lat','lon'])

  #---------------------------------------------------------------------------------------------------------
  # Weight isotopic values by seasonal soil water amount
  # NOTE: 'wtvar_mon_wgt' variable will already account for how many months are present, but 'iso' will not 
  #---------------------------------------------------------------------------------------------------------

  iso_wtvar = xr.DataArray( xr.DataArray(wtvar_mon_wgt*iso[months,:,:,:]).sum(dim='time') / \
                            xr.DataArray(wtvar_mon_wgt).sum(dim='time')                   )
 
  #---------------------------------------------------------
  # Weight soil water isotopes by fractional rooting depth
  #---------------------------------------------------------

  # Read in climatology of rooting depth variable
  root_clim = seassoil.seascyc_soilvar_TSLL(var1=rootwgt,path=path,begi=begi,endi=endi)

  # Weight root_clim by month
  root_prewgt = root_clim[months,:,:,:].weighted(wgt_mon)
  root_mon_wgt = root_prewgt.mean(dim='time')

  # Remove values < 0 in average root depth fraction variable
  root_mon_wgt_proc = xr.where(root_mon_wgt < 0, float('NaN'), root_mon_wgt) # remove values < 0 and set them to NaN

  # Weight month- and soil water amount-weighted isotopic variable by fractional rooting depth
  iso_wt = xr.DataArray( xr.DataArray(iso_wtvar*root_mon_wgt_proc).sum(dim=soildimname) /                      \
                         xr.DataArray(root_mon_wgt_proc,dims=[soildimname,'lat','lon']).sum(dim=soildimname) )
  
  return iso_wt

#######################################################################################################
# Wind U and V vectors
#######################################################################################################

def seasavg_wind_vec_LL(months: list,
                        path: str,
                        begi: str = 'beg',
                        endi: str = 'end',
                        level: int = None,
                        plev: np.ndarray = default_levels,
                        wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in wind (U,V) variables from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) of vectors output variable.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
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
  level: class 'int', The integer value with units hPa/mb of the atmospheric level chosen from which to 
                      calculate the variable's values. Requires plev to be specified for the range of 
                      pressure levels. Ex. 850 -> calculating variable's values at 850 hPa 
  plev: class 'numpy.ndarray', The new range of pressure levels with units hPa/mb, of which 'level' must 
                               be an available value. Default = np.arange(0,1050,50) 
                               This script will interpolate climate model output from hybrid to pressure 
                               levels using geocat.comp's interp_hybrid_to_pressure function. NOTE: If 
                               NaN's present in variable while interpolating from hybrid to pressure, may 
                               see this error, which comes from metpy and can be ignored: 'UserWarning: 
                               Interpolation point out of data bounds encountered 
                               warnings.warn('Interpolation point out of data bounds encountered')'
                               Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with 
                               stride of 50 hPa. Default = np.arange(0,1050,50) 
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable of vectors [lat,lon]
  
  '''

  #----------------------------
  # Ensure type of parameters
  #----------------------------

  months = list(months)
  path = str(path)
  level = int(level)
  plev = np.array(plev)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #-----------------------------
  # Calculate seasonal cycle
  #-----------------------------

  u_clim, v_clim = seascyc.seascyc_wind_vec_TLL(path=path,begi=begi,endi=endi,level=level,plev=plev)

  #--------------------------------------------
  # Weight by fractional length of each month
  #--------------------------------------------

  u_avg = mon_wgt_avg(var=u_clim,months=months,wgt_mon=wgt_mon)
  v_avg = mon_wgt_avg(var=v_clim,months=months,wgt_mon=wgt_mon)

  return u_avg, v_avg

#######################################################################################################
# Integrated vapor transport (IVT) vectors
#######################################################################################################

def seasavg_IVT_vec_LL(months: list,
                       path: str,
                       begi: str = 'beg',
                       endi: str = 'end',
                       ptop: float = 0.,
                       pbot: float = 1000., 
                       plev: np.ndarray = default_levels,
                       wgt_mon: xr.DataArray = default_wgt_mon) -> xr.DataArray:

  '''Reads in U, V, and Q variables from a netCDF file (that includes a dimension for time) and creates
  a seasonally averaged global map (lat x lon) of integrated vapor transport vectors. 

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
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
  ptop: class 'float', Top boundary value of averaging region for vertical integration, in pressure (hPa).
                       Ex. 50. -> vertical integration will sum from pbot to 50 hPa.
  pbot: class 'float', Bottom boundary value of averaging region for vertical integration, in pressure 
                       (hPa). Ex. 1000. -> vertical integration will sum from 1000 hPa to ptop. 
  plev: class 'numpy.ndarray', The new range of pressure levels with units hPa/mb, of which 'level' must 
                               be an available value. Default = np.arange(0,1050,50) 
                               This script will interpolate climate model output from hybrid to pressure 
                               levels using geocat.comp's interp_hybrid_to_pressure function. NOTE: If 
                               NaN's present in variable while interpolating from hybrid to pressure, may 
                               see this error, which comes from metpy and can be ignored: 'UserWarning: 
                               Interpolation point out of data bounds encountered 
                               warnings.warn('Interpolation point out of data bounds encountered')'
                               Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with 
                               stride of 50 hPa. Default = np.arange(0,1050,50) 
  wgt_mon: class 'xarray.DataArray', Array of monthly weights (size: 12). This feature is available for 
                                     paleoclimate analysis when the relative weight of each month differed 
                                     from the values of today. Must be same size as 'months'. Default = all 
                                     months have the same weight. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonally averaged global variable of IVT vectors [lat,lon]
  
  '''

  #----------------------------
  # Ensure type of parameters
  #----------------------------

  #----------------------------
  # Ensure type of parameters
  #----------------------------

  months = list(months)
  path = str(path)
  ptop = float(ptop)
  pbot = float(pbot) 
  plev = np.array(plev)
  wgt_mon = xr.DataArray(wgt_mon,dims=['time'])

  #-----------------------------
  # Calculate seasonal cycle
  #-----------------------------

  u_clim, v_clim = seascyc.seascyc_IVT_vec_TLL(path=path,begi=begi,endi=endi,ptop=ptop,pbot=pbot,plev=plev)

  #--------------------------------------------
  # Weight by fractional length of each month
  #--------------------------------------------

  u_avg = mon_wgt_avg(var=u_clim,months=months,wgt_mon=wgt_mon)
  v_avg = mon_wgt_avg(var=v_clim,months=months,wgt_mon=wgt_mon)

  return u_avg, v_avg

