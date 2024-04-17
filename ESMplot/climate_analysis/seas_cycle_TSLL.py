#######################################################################################################
#
# These functions take an input from netCDF output and return a seasonal cycle variable that retains
# the soil level dimension. Final variable will be of size (12 months x soil level x lat x lon). 
#
#######################################################################################################

import numpy as np
import xarray as xr
from geocat.comp.interpolation import interp_hybrid_to_pressure 
from ESMplot.climate_analysis.climatology import clmMonTLL,clmMonTLLL,clmMonTSLL

# Default arrays
default_levels  = np.arange(0,1050,50) # 0 hPa to 1000 hPa by 50   

#======================================================================================================
# File contains the following functions:
#
# Contours (output: 12 months x soil level x lat x lon)
# seascyc_soilvar_TSLL,seascyc_soiliso_TSLL                                                                      
#
#======================================================================================================

#######################################################################################################
# Land variable with soil levels 
#######################################################################################################

def seascyc_soilvar_TSLL(var1: str,
                         path: str,
                         begi: str = 'beg',
                         endi: str = 'end',
                         var2: str = None,
                         math: str = None,
                         soildimname: str = 'levgrnd',
                         mult: float = 1.) -> xr.DataArray:

  '''Reads in up to two land variables from a netCDF file (that includes a dimension for time and soil 
  level) and creates a seasonal cycle global map (12 months x soil level x lat x lon) output variable.

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
  var1 
   class: 'string', A string specifying the variable to be read in from the netCDF file.
                    Ex. 'FSNT' (CESM output variable for net shortwave flux at top of model, as found at 
                    https://www.cesm.ucar.edu/models/cesm2/atmosphere/docs/ug6/hist_flds_f2000.html)

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
  soildimname: class 'str', Name of soil level dimension. Default is 'levgrnd', which is standard for CESM 
                            variables. This must correspond to the correct dimension name in var for proper
                            averaging to occur.
  mult: class 'float', A float value that becomes a scalar to the variable's values. Default = 1.0. Will 
                       revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,soillevel,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  var1 = str(var1)
  path = str(path)
  soildimname = str(soildimname)
  mult = float(mult)
  if var2 != None:
   var2 = str(var2)
  if math != None:
   math = str(math)

  #-----------------------------------
  # Read in variables and scale them 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Assign values for begi and endi
  begi = 0 if begi == 'beg' else begi
  endi = len(data.time) if endi == 'end' else endi

  # Read in variable(s)
  var1_time = data[var1][begi:endi,:,:,:]*mult
  if var2 != None:
   var2_time = data[var2][begi:endi,:,:,:]*mult

  #--------------------------------------
  # Perform math operator if necessary
  #--------------------------------------

  if var2 != None:
   if math == 'add':
    var_time = var1_time + var2_time
   elif math == 'sub':
    var_time = var1_time - var2_time
   elif math == 'mul':
    var_time = var1_time * var2_time
   elif math == 'div':
    var_time = var1_time / var2_time
   elif math == 'exp':
    var_time = var1_time ** var2_time
  else:
   var_time = var1_time

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = clmMonTSLL(var_time,soil_dim=soildimname)

  return var_clim

#######################################################################################################
# Soil water isotopes: d18Op, dDp, or dexcess of soil water 
#######################################################################################################

def seascyc_soiliso_TSLL(iso_type:str,
                         path: str,
                         begi: str = 'beg',
                         endi: str = 'end',
                         ptiny: float = 1.E-18,
                         soildimname: str = 'levgrnd') -> xr.DataArray:

  '''Reads in soil water isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x soil level x lat x lon) output variable.

  Required Parameters
  --------------------------------------------------------------------------------------------------------
  iso_type
   class: 'string', Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

  path
   class: 'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  Optional Parameters
  --------------------------------------------------------------------------------------------------------
  begi, endi: class 'string' or 'int', Determines the first and last time element to be included in 
                                       analysis. Defaults 'beg' and 'end' read in the entire length of
                                       the file. Specifying an integer will read in time from begi as int 
                                       to endi as int. If integers, (endi-begi)/12 must be multiple of 12. 
  ptiny: class 'float', Minimum value for light isotopic variable. This is required in the calculation of 
                        isotopes so as not to have artificially high delta values due to extremely low 
                        light isotopic values. Default = 1.E-18
  soildimname: class 'str', Name of soil level dimension. Default is 'levgrnd', which is standard for CESM 
                            variables. This must correspond to the correct dimension name in var for proper 
                            averaging to occur.
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,soillevel,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  iso_type = str(iso_type)
  path = str(path)
  soildimname = str(soildimname)

  #-----------------------------------
  # Read in dataset
  #-----------------------------------

  # Read in dataset
  data = xr.open_dataset(path)

  # Assign values for begi and endi
  begi = 0 if begi == 'beg' else begi
  endi = len(data.time) if endi == 'end' else endi

  #------------------------------------
  # Soil water isotopic variable   
  #------------------------------------

  # Read in light isotopic values
  liso_d18O_init = xr.DataArray(data.H2OSOI_H2OTR[begi:endi,:,:,:])
  liso_dHDO_init = xr.DataArray(data.H2OSOI_H2OTR[begi:endi,:,:,:])

  # Remove extremely small values from light isotopic variable
  liso_d18O = xr.where(liso_d18O_init < ptiny,ptiny,liso_d18O_init)
  liso_dHDO = xr.where(liso_dHDO_init < ptiny,ptiny,liso_dHDO_init)

  # Read in heavy isotopic values
  hiso_d18O = xr.DataArray(data.H2OSOI_H218O[begi:endi,:,:,:])
  hiso_dHDO = xr.DataArray(data.H2OSOI_HDO[begi:endi,:,:,:])

  # Calculate delta values
  d18O = ( (hiso_d18O / liso_d18O) - 1. ) * 1000.
  dHDO = ( (hiso_dHDO / liso_dHDO) - 1. ) * 1000.

  #---------------------
  # Apply isotope type
  #---------------------

  if iso_type == 'd18O':
   var_time = xr.DataArray(d18O,dims=['time',soildimname,'lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dHDO':
   var_time = xr.DataArray(dHDO,dims=['time',soildimname,'lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dexcess':
   var_time = xr.DataArray(dHDO - 8 * d18O,dims=['time',soildimname,'lat','lon'],coords=dict(lat=data.lat,lon=data.lon))

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = xr.DataArray(clmMonTSLL(var_time,soil_dim=soildimname),dims=['time',soildimname,'lat','lon'],
                          coords={'time':np.arange(12),soildimname:data[soildimname],'lat':data.lat,'lon':data.lon})

  return var_clim

