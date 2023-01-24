#######################################################################################################
#
# These functions take an input from netCDF output and return a seasonal cycle variable of size 
# (12 months x lat x lon). 
#
#######################################################################################################

import numpy as np
import xarray as xr
from geocat.comp.interpolation import interp_hybrid_to_pressure 
from climate_analysis.climatology import clmMonTLL,clmMonTLLL,clmMonTSLL
from geocat.f2py.dpres_plevel_wrapper import dpres_plevel
from warnings import simplefilter

# Default arrays
default_levels  = np.arange(0.,1050.,50.) # 0 hPa to 1000 hPa by 50   

#======================================================================================================
# File contains the following functions:
#
# Contours (output: 12 months x lat x lon)
# seascyc_var_TLL, seascyc_prect_TLL, seascyc_soilvar_TLL, seascyc_rainiso_TLL, seascyc_soiliso_TLL,
# seascyc_vaporiso_TLL, seascyc_isoroot_TLL    
#
# Vectors (output: 12 months x lat x lon)
# seascyc_wind_vec_TLL, seascyc_IVT_vec_TLL
#
#======================================================================================================

#######################################################################################################
# Generic variable at surface or pressure level of the atmosphere 
#######################################################################################################

def seascyc_var_TLL(var1: str,
                    path: str,
                    begi: int,
                    endi: int,
                    var2: str = None, 
                    math: str = None, 
                    level: int = None, 
                    plev: np.ndarray = default_levels,
                    mult: float = 1.) -> xr.DataArray: 

  '''Reads in up to two variables from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  var1,path,begi,endi

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  var1: :class:'string'
        A string of the variable to read in from the netCDF file.
        Ex. 'FSNT' (CESM output variable for net shortwave flux at top of model, as found at 
        https://www.cesm.ucar.edu/models/cesm2/atmosphere/docs/ug6/hist_flds_f2000.html)

  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  var2: :class:'string', Optional 
        Optional: A string for a second variable to read in from the netCDF file. If 'var2' is called,
        'math' must also be specified.
        Ex. 'FLNT' (CESM output variable for net longwave flux at top of model)
 
  math: :class:'string', Optional 
        Specified mathematical operation to be applied to 'var1' and 'var2'.
        Choices are: 'add' (var1+var2), 'sub' (var1-var2), 'mul' (var1*var2), 'div' (var1/var2),     
        and 'exp' (var1**var2).

  level: :class:'int', Optional 
         The integer value with units hPa/mb of the atmospheric level chosen from which to calculate 
         the variable's values. Requires plev to be specified for the range of pressure levels.
         Ex. 850 -> calculating variable's values at 850 hPa 

  plev: :class:'numpy.ndarray', Optional - default = np.arange(0,1050,50) 
         The new range of pressure levels with units hPa/mb, of which 'level' must be an available value. 
         This script will interpolate climate model output from hybrid to pressure levels using 
         geocat.comp's interp_hybrid_to_pressure function.
         NOTE: If NaN's present in variable while interpolating from hybrid to pressure, may see this
         error, which comes from metpy and can be ignored: 
          'UserWarning: Interpolation point out of data bounds encountered 
           warnings.warn('Interpolation point out of data bounds encountered')'
         Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with stride of 50 hPa  

  mult: :class:'float', Optional - default = 1.
        A float value that becomes a scalar to the variable's values. 
        Default = 1.0. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  var1 = str(var1)
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  plev = np.array(plev)
  mult = float(mult)
  if var2 != None:
   var2 = str(var2)
  if math != None:
   math = str(math)

  # Assign index for level in plev
  if level != None:
   lev_ind = int(np.where(np.array(plev) == level)[0])

  #-----------------------------------
  # Read in variables and scale them 
  #-----------------------------------
 
  # Read in dataset 
  data = xr.open_dataset(path)
 
  # Read in variable(s)
  var1_data = data[var1][begi:endi,:,:]*mult
  if var2 != None:
   var2_data = data[var2][begi:endi,:,:]*mult
 
  # Convert sigma coordinates to pressure coordinates
  if level != None:
   simplefilter(action='ignore', category=UserWarning) # ignore unnecessary error from metpy
   hyam = data.hyam 
   hybm = data.hybm
   psrf = data.PS
   p0Pa  = 100000.
   var1_lev = interp_hybrid_to_pressure(data=var1_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.) # converts plev to units=Pa
   var1_time = xr.DataArray(var1_lev[:,lev_ind,:,:], dims=['time','lat','lon'],coords=dict(time=var1_data.time,lat=var1_data.lat,lon=var1_data.lon),
                                                     attrs=dict(description=str(var1)+' at '+str(level)+' hPa'))
   if var2 != None:
    var2_lev = interp_hybrid_to_pressure(data=var2_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.) # converts plev to units=Pa
    var2_time = xr.DataArray(var2_lev[:,lev_ind,:,:], dims=['time','lat','lon'],coords=dict(time=var1_data.time,lat=var1_data.lat,lon=var1_data.lon),
                                                      attrs=dict(description=str(var2)+' at '+str(level)+' hPa'))
  else:
   var1_time = xr.DataArray(var1_data)
   if var2 != None:
    var2_time = xr.DataArray(var2_data)

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
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Precipitation variable 
#######################################################################################################

def seascyc_prect_TLL(path: str,
                      begi: int,
                      endi: int,
                      mult: float = 86400000.) -> xr.DataArray:

  '''Reads in precipitation variables from a netCDF file (that includes a dimension for time) 
  and creates a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  path,begi,endi

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  mult: :class:'float', Optional - default = 86400000. which converts m/s to mm/day
        A float value that becomes a scalar to the variable's values. Assumes original variable is in m/s.
        Default = 86400000. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  mult = float(mult)

  #-----------------------------------
  # Read in variables and scale them 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Read in variable(s)
  var_time = (data['PRECC'][begi:endi,:,:] + data['PRECL'][begi:endi,:,:])*mult

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Precipitation-Evaporation variable (Effective moisture)
#######################################################################################################

def seascyc_PminE_TLL(path: str,
                      begi: int,
                      endi: int,
                      Pmult: float = 86400000.,
                      Emult: float = 86400.) -> xr.DataArray:

  '''Reads in precipitation and evaporation variables from a netCDF file (that includes a dimension 
  for time) and creates a seasonal cycle global map (12 months x lat x lon) output variable of 
  Precipitation minus (-) Evaporation, otherwise known as effective moisture.

  Required parameters:
  ----------------------------
  path,begi,endi

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  Pmult: :class:'float', Optional - default = 86400000. which converts m/s to mm/day
         A float value that becomes a scalar to precipitation values. Assumes original variable is in m/s.
         Default = 86400000. Will revert to default if not specified. 

  Emult: :class:'float', Optional - default = 86400. which converts m/s to mm/day
         A float value that becomes a scalar to evaporation values. Assumes original variable is in mm/s.
         Default = 86400. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  Pmult = float(Pmult)
  Emult = float(Emult)

  #-----------------------------------
  # Read in variables and scale them 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Read in variables
  p_time = (data['PRECC'][begi:endi,:,:] + data['PRECL'][begi:endi,:,:])*Pmult
  e_time = (data['QFLX'][begi:endi,:,:])*Emult

  # Calculate P-E
  PminE_time = p_time - e_time

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = clmMonTLL(PminE_time)

  return var_clim

#######################################################################################################
# Land variable with soil levels 
#######################################################################################################

def seascyc_soilvar_TLL(var1: str,
                        path: str,
                        begi: int,
                        endi: int,
                        var2: str = None,
                        math: str = None,
                        soillev: list = [0],
                        soildimname: str = 'levgrnd',
                        mult: float = 1.) -> xr.DataArray:

  '''Reads in up to two land variables from a netCDF file (that includes a dimension for time and soil 
  level) and creates a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  var1,path,begi,endi

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  var1: :class:'string'
        A string of the variable to read in from the netCDF file.
        Ex. 'FSNT' (CESM output variable for net shortwave flux at top of model, as found at 
        https://www.cesm.ucar.edu/models/cesm2/atmosphere/docs/ug6/hist_flds_f2000.html)

  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  var2: :class:'string', Optional 
        Optional: A string for a second variable to read in from the netCDF file. If 'var2' is called,
        'math' must also be specified.
        Ex. 'FLNT' (CESM output variable for net longwave flux at top of model)
 
  math: :class:'string', Optional 
        Specified mathematical operation to be applied to 'var1' and 'var2'.
        Choices are: 'add' (var1+var2), 'sub' (var1-var2), 'mul' (var1*var2), 'div' (var1/var2),     
        and 'exp' (var1**var2).

  soillev: :class:'list', Optional 
           List of integers corresponding to the indices of the soil level dimension (ex. levgrnd) over
           which an average will be calculated. If not specified, first soil level will be indexed.
           Ex. [0,1,2] calculates the average of the first three soil level indices.

  soildimname: :class:'str', Optional
               Name of soil level dimension. Default is 'levgrnd', which is standard for CESM variables.
               This must correspond to the correct dimension name in var for proper averaging to occur.

  mult: :class:'float', Optional - default = 1.
        A float value that becomes a scalar to the variable's values. 
        Default = 1.0. Will revert to default if not specified. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  var1 = str(var1)
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  soillev = list(soillev)
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

  # Read in variable(s)
  var1_time = data[var1][begi:endi,soillev,:,:].mean(dim=soildimname)*mult
  if var2 != None:
   var2_time = data[var2][begi:endi,soillev,:,:].mean(dim=soildimname)*mult

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
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Precipitation isotopes: d18Op, dDp, or dexcess of precip 
#######################################################################################################

def seascyc_rainiso_TLL(iso_type:str,
                        path: str,
                        begi: int,
                        endi: int,
                        ptiny: float = 1.E-18) -> xr.DataArray:

  '''Reads in precipitation isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  iso_type,path,begi,endi

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  iso_type: :class:'string'
            Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  ptiny: :class:'float' - Optional, default = 1.E-18
         Minimum value for light isotopic variable. This is required in the calculation of isotopes so as
         not to have artificially high delta values due to extremely low light isotopic values. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  iso_type = str(iso_type)
  path = str(path)
  begi = int(begi)
  endi = int(endi)

  #-----------------------------------
  # Read in dataset 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  #------------------------------------
  # Precipitation isotopic variable   
  #------------------------------------

  # Read in light isotopic values
  liso_d18O_init = xr.DataArray(data.PRECRC_H216Or[begi:endi,:,:] + data.PRECSC_H216Os[begi:endi,:,:] + \
                                data.PRECRL_H216OR[begi:endi,:,:] + data.PRECSL_H216OS[begi:endi,:,:])
  liso_dHDO_init = xr.DataArray(data.PRECRC_H2Or[begi:endi,:,:] + data.PRECSC_H2Os[begi:endi,:,:] + \
                                data.PRECRL_H2OR[begi:endi,:,:] + data.PRECSL_H2OS[begi:endi,:,:])

  # Remove extremely small values from light isotopic variable
  liso_d18O = xr.where(liso_d18O_init < ptiny,ptiny,liso_d18O_init)
  liso_dHDO = xr.where(liso_dHDO_init < ptiny,ptiny,liso_dHDO_init)

  # Read in heavy isotopic values
  hiso_d18O = xr.DataArray(data.PRECRC_H218Or[begi:endi,:,:] + data.PRECSC_H218Os[begi:endi,:,:] + \
                           data.PRECRL_H218OR[begi:endi,:,:] + data.PRECSL_H218OS[begi:endi,:,:])
  hiso_dHDO = xr.DataArray(data.PRECRC_HDOr[begi:endi,:,:] + data.PRECSC_HDOs[begi:endi,:,:] + \
                           data.PRECRL_HDOR[begi:endi,:,:] + data.PRECSL_HDOS[begi:endi,:,:])

  # Calculate delta values
  d18O = ( (hiso_d18O / liso_d18O) - 1. ) * 1000.
  dHDO = ( (hiso_dHDO / liso_dHDO) - 1. ) * 1000.

  #---------------------
  # Apply isotope type
  #---------------------

  if iso_type == 'd18O':
   var_time = xr.DataArray(d18O,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dHDO':
   var_time = xr.DataArray(dHDO,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dexcess':
   var_time = xr.DataArray(dHDO - 8 * d18O,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Soil water isotopes: d18Op, dDp, or dexcess of soil water 
#######################################################################################################

def seascyc_soiliso_TLL(iso_type:str,
                        path: str,
                        begi: int,
                        endi: int,
                        ptiny: float = 1.E-18,
                        soillev: list = [0],
                        soildimname: str = 'levgrnd') -> xr.DataArray:

  '''Reads in soil water isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  iso_type,path,begi,endi

  Parameters
  --------------------------------------------------------------------------------------------------------
  iso_type: :class:'string'
            Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January.

  endi: :class:'int'
        The index of the last time to be included in analysis.
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  ptiny: :class:'float' - Optional, default = 1.E-18
         Minimum value for light isotopic variable. This is required in the calculation of isotopes so as
         not to have artificially high delta values due to extremely low light isotopic values.

  soillev: :class:'list', Optional
           List of integers corresponding to the indices of the soil level dimension (ex. levgrnd) over
           which an average will be calculated. If not specified, first soil level will be indexed.
           Ex. [0,1,2] calculates the average of the first three soil level indices.

  soildimname: :class:'str', Optional
               Name of soil level dimension. Default is 'levgrnd', which is standard for CESM variables.
               This must correspond to the correct dimension name in var for proper averaging to occur.
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  iso_type = str(iso_type)
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  soillev = list(soillev)
  soildimname = str(soildimname)

  #-----------------------------------
  # Read in dataset
  #-----------------------------------

  # Read in dataset
  data = xr.open_dataset(path)

  #------------------------------------
  # Soil water isotopic variable   
  #------------------------------------

  # Read in light isotopic values
  liso_d18O_init = xr.DataArray(data.H2OSOI_H2OTR[begi:endi,soillev,:,:]).mean(dim=soildimname) 
  liso_dHDO_init = xr.DataArray(data.H2OSOI_H2OTR[begi:endi,soillev,:,:]).mean(dim=soildimname)

  # Remove extremely small values from light isotopic variable
  liso_d18O = xr.where(liso_d18O_init < ptiny,ptiny,liso_d18O_init)
  liso_dHDO = xr.where(liso_dHDO_init < ptiny,ptiny,liso_dHDO_init)

  # Read in heavy isotopic values
  hiso_d18O = xr.DataArray(data.H2OSOI_H218O[begi:endi,soillev,:,:]).mean(dim=soildimname)
  hiso_dHDO = xr.DataArray(data.H2OSOI_HDO[begi:endi,soillev,:,:]).mean(dim=soildimname)

  # Calculate delta values
  d18O = ( (hiso_d18O / liso_d18O) - 1. ) * 1000.
  dHDO = ( (hiso_dHDO / liso_dHDO) - 1. ) * 1000.

  #---------------------
  # Apply isotope type
  #---------------------

  if iso_type == 'd18O':
   var_time = xr.DataArray(d18O,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dHDO':
   var_time = xr.DataArray(dHDO,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dexcess':
   var_time = xr.DataArray(dHDO - 8 * d18O,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Vapor isotopes at specified pressure level: d18Op, dDp, or dexcess of water vapor
#######################################################################################################

def seascyc_vaporiso_TLL(iso_type: str,
                         path: str,
                         begi: int,
                         endi: int,
                         level: int,
                         plev: np.ndarray = default_levels,
                         ptiny: float = 1.E-18) -> xr.DataArray:

  '''Reads in water vapor isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  iso_type,path,begi,endi,level

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  iso_type: :class:'string'
            Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  level: :class:'int'
         The integer value with units hPa/mb of the atmospheric level chosen from which to calculate 
         the variable's values. Requires plev to be specified for the range of pressure levels.
         Ex. 850 -> calculating variable's values at 850 hPa 

  plev: :class:'numpy.ndarray', Optional - default = np.arange(0,1050,50) 
         The new range of pressure levels with units hPa/mb, of which 'level' must be an available value. 
         This script will interpolate climate model output from hybrid to pressure levels using 
         geocat.comp's interp_hybrid_to_pressure function.
         NOTE: If NaN's present in variable while interpolating from hybrid to pressure, may see this
         error, which comes from metpy and can be ignored: 
          'UserWarning: Interpolation point out of data bounds encountered 
           warnings.warn('Interpolation point out of data bounds encountered')'
         Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with stride of 50 hPa 

  ptiny: :class:'float' - Optional, default = 1.E-18
         Minimum value for light isotopic variable. This is required in the calculation of isotopes so as
         not to have artificially high delta values due to extremely low light isotopic values. 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  iso_type = str(iso_type)
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  level = int(level)

  # Time variables for reading in values 
  numyrs = int(((endi-begi)+1)/12) # number of years calculated from begi,endi

  # Assign index for level in plev
  lev_ind = int(np.where(np.array(plev) == level)[0])

  #-----------------------------------
  # Read in dataset 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Read in variables for converting sigma to pressure coordinates
  hyam = data.hyam
  hybm = data.hybm
  psrf = data.PS[begi:endi,:,:]
  p0Pa  = 100000.

  #------------------------------------
  # Water vapor isotopic variable   
  #------------------------------------

  # Read in light isotopic values
  liso_d18O_init = data.H216OV[begi:endi,:,:,:]
  liso_dHDO_init = data.H2OV[begi:endi,:,:,:]

  # Ignore unnecessary error from metpy
  simplefilter(action='ignore', category=UserWarning) 

  # Convert light isotopic values from sigma to pressure coordinates
  liso_d18O_lev = interp_hybrid_to_pressure(data=liso_d18O_init,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,
                                            new_levels=plev*100.)
  liso_dHDO_lev = interp_hybrid_to_pressure(data=liso_dHDO_init,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,
                                            new_levels=plev*100.)
  liso_d18O_time = xr.DataArray(liso_d18O_lev[:,lev_ind,:,:],dims=['time','lat','lon'],
                            coords=dict(time=liso_d18O_lev.time,lat=liso_d18O_lev.lat,lon=liso_d18O_lev.lon))
  liso_dHDO_time = xr.DataArray(liso_dHDO_lev[:,lev_ind,:,:],dims=['time','lat','lon'],
                            coords=dict(time=liso_dHDO_lev.time,lat=liso_dHDO_lev.lat,lon=liso_dHDO_lev.lon))

  # Remove extremely small values from light isotopic variable
  liso_d18O = xr.where(liso_d18O_time < ptiny,ptiny,liso_d18O_time)
  liso_dHDO = xr.where(liso_dHDO_time < ptiny,ptiny,liso_dHDO_time)

  # Read in heavy isotopic values
  hiso_d18O_init = data.H218OV[begi:endi,:,:,:]
  hiso_dHDO_init = data.HDOV[begi:endi,:,:,:]

  # Convert light isotopic values from sigma to pressure coordinates
  hiso_d18O_lev = interp_hybrid_to_pressure(data=hiso_d18O_init,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,
                                            new_levels=plev*100.)
  hiso_dHDO_lev = interp_hybrid_to_pressure(data=hiso_dHDO_init,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,
                                            new_levels=plev*100.)
  hiso_d18O = xr.DataArray(hiso_d18O_lev[:,lev_ind,:,:],dims=['time','lat','lon'],
                            coords=dict(time=hiso_d18O_lev.time,lat=hiso_d18O_lev.lat,lon=hiso_d18O_lev.lon))
  hiso_dHDO = xr.DataArray(hiso_dHDO_lev[:,lev_ind,:,:],dims=['time','lat','lon'],
                            coords=dict(time=hiso_dHDO_lev.time,lat=hiso_dHDO_lev.lat,lon=hiso_dHDO_lev.lon))

  # Calculate delta values
  d18O = ( (hiso_d18O / liso_d18O) - 1. ) * 1000.
  dHDO = ( (hiso_dHDO / liso_dHDO) - 1. ) * 1000.

  #---------------------
  # Apply isotope type
  #---------------------

  if iso_type == 'd18O':
   var_time = xr.DataArray(d18O,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dHDO':
   var_time = xr.DataArray(dHDO,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))
  if iso_type == 'dexcess':
   var_time = xr.DataArray(dHDO - 8 * d18O,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Soil water isotopes weighted by root depth fraction: d18Op, dDp, or dexcess of soil water 
#######################################################################################################

def seascyc_isoroot_TLL(iso_type:str,
                        path: str,
                        begi: int,
                        endi: int,
                        ptiny: float = 1.E-18,
                        rootwgt: str = 'ROOTR_COLUMN',
                        soildimname: str = 'levgrnd') -> xr.DataArray:

  '''Reads in soil water isotopes from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable weighted by root depth fraction. 
  This variable is used for comparison to leaf wax proxy records.

  Required parameters:
  ----------------------------
  iso_type,path,begi,endi

  Parameters
  --------------------------------------------------------------------------------------------------------
  iso_type: :class:'string'
            Type of isotopic variable to plot. Options are: 'd18O', 'dD', or 'dexcess'

  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January.

  endi: :class:'int'
        The index of the last time to be included in analysis.
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  ptiny: :class:'float' - Optional, default = 1.E-18
         Minimum value for light isotopic variable. This is required in the calculation of isotopes so as
         not to have artificially high delta values due to extremely low light isotopic values.

  rootwgt: :class:'str' - Optional, default = 'ROOTR_COLUMN'
           Name for variable used for weighting isotopic values by root depth fraction. In CESM, possible
           options are 'ROOTR_COLUMN', 'ROOTFR', or 'ROOTR' 

  soildimname: :class:'str', Optional
               Name of soil level dimension. Default is 'levgrnd', which is standard for CESM variables.
               This must correspond to the correct dimension name in var for proper averaging to occur.
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  iso_type = str(iso_type)
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  ptiny = float(ptiny)
  rootwgt = str(rootwgt)
  soildimname = str(soildimname)

  # Time variables for reading in values
  numyrs = int(((endi-begi)+1)/12) # number of years calculated from begi,endi

  #-----------------------------------
  # Read in dataset
  #-----------------------------------

  # Read in dataset
  data = xr.open_dataset(path)

  #------------------------------------
  # Soil water isotopic variable   
  #------------------------------------

  # Read in light isotopic values
  liso_d18O_init = data.H2OSOI_H2OTR[begi:endi,:,:,:]
  liso_dHDO_init = data.H2OSOI_H2OTR[begi:endi,:,:,:]

  # Remove extremely small values from light isotopic variable
  liso_d18O = xr.where(liso_d18O_init < ptiny,ptiny,liso_d18O_init)
  liso_dHDO = xr.where(liso_dHDO_init < ptiny,ptiny,liso_dHDO_init)

  # Read in heavy isotopic values
  hiso_d18O = data.H2OSOI_H218O[begi:endi,:,:,:]
  hiso_dHDO = data.H2OSOI_HDO[begi:endi,:,:,:]

  # Calculate delta values
  d18Os = ( (hiso_d18O / liso_d18O) - 1. ) * 1000.
  dHDOs = ( (hiso_dHDO / liso_dHDO) - 1. ) * 1000.

  #-----------------------
  # Apply isotope type
  #-----------------------

  # Variable now has shape (levgrnd x lat x lon)

  if iso_type == 'd18O':
   isoslev = d18Os
  if iso_type == 'dHDO':
   isoslev = dHDOs
  if iso_type == 'dexcess':
   isoslev = xr.DataArray(dHDOs - 8 * d18Os,dims=['time',soildimname,'lat','lon'],
                          coords={'time':data.time,soildimname:data[soildimname],
                                  'lat':data.lat,'lon':data.lon})

  #---------------------------------------------------------
  # Weight soil water isotopes by fractional rooting depth
  #---------------------------------------------------------

  # Read in rooting depth and average years together
  root_init = data[rootwgt][begi:endi,:,:,:]

  # Remove values < 0 in average root depth fraction variable
  root_proc = xr.where(root_init < 0, float('NaN'), root_init) # remove values < 0 and set them to NaN

  # Weight the isotopes by the fractional rooting depth
  iso_time = xr.DataArray( xr.DataArray(isoslev*root_proc).sum(dim=soildimname) /                      \
                           xr.DataArray(root_proc,dims=['time',soildimname,'lat','lon']).sum(dim=soildimname) )

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Set root_seas as xarray DataArray
  var_time = xr.DataArray(iso_time,dims=['time','lat','lon'],coords=dict(lat=data.lat,lon=data.lon))

  # Calculate averages for each month
  var_clim = clmMonTLL(var_time)

  return var_clim

#######################################################################################################
# Wind U and V vectors
#######################################################################################################

def seascyc_wind_vec_TLL(path: str,
                         begi: int,
                         endi: int,
                         level: int = None,
                         plev: np.ndarray = default_levels) -> xr.DataArray:

  '''Reads in wind (U,V) variables from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable.

  Required parameters:
  ----------------------------
  path,begi,endi,level,plev

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  level: :class:'int'
         The integer value with units hPa/mb of the atmospheric level chosen from which to calculate 
         the variable's values. Requires plev to be specified for the range of pressure levels.
         Ex. 850 -> calculating variable's values at 850 hPa 

  plev: :class:'numpy.ndarray'
         The new range of pressure levels with units hPa/mb, of which 'level' must be an available value. 
         This script will interpolate climate model output from hybrid to pressure levels using 
         geocat.comp's interp_hybrid_to_pressure function.
         NOTE: If NaN's present in variable while interpolating from hybrid to pressure, may see this
         error, which comes from metpy and can be ignored: 
          'UserWarning: Interpolation point out of data bounds encountered 
           warnings.warn('Interpolation point out of data bounds encountered')'
         Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with stride of 50 hPa  
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]
  
  '''
  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  level = int(level)
  plev = np.array(plev)

  # Assign index for level in plev
  lev_ind = int(np.where(np.array(plev) == level)[0])

  #-----------------------------------
  # Read in variables and scale them 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Read in variable(s)
  u_data = data['U'][begi:endi,:,:,:]
  v_data = data['V'][begi:endi,:,:,:]

  # Convert sigma coordinates to pressure coordinates
  simplefilter(action='ignore', category=UserWarning) # ignore unnecessary error from metpy
  hyam = data.hyam
  hybm = data.hybm
  psrf = data.PS
  p0Pa  = 100000.
  u_lev = interp_hybrid_to_pressure(data=u_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.) # converts plev to units=Pa
  v_lev = interp_hybrid_to_pressure(data=v_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.) # converts plev to units=Pa
  u_time = xr.DataArray(u_lev[:,lev_ind,:,:], dims=['time','lat','lon'],coords=dict(time=u_data.time,lat=u_data.lat,lon=u_data.lon),
                                                     attrs=dict(description='U winds at '+str(level)+' hPa'))
  v_time = xr.DataArray(v_lev[:,lev_ind,:,:], dims=['time','lat','lon'],coords=dict(time=v_data.time,lat=v_data.lat,lon=v_data.lon),
                                                     attrs=dict(description='V winds at '+str(level)+' hPa'))

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  u_clim = clmMonTLL(u_time)
  v_clim = clmMonTLL(v_time)

  # Set NaN values
  u_clim = u_clim.where(u_clim < 1E10)
  v_clim = v_clim.where(v_clim < 1E10)

  return u_clim, v_clim

#######################################################################################################
# Integrated vapor transport (IVT) vectors 
#######################################################################################################

def seascyc_IVT_vec_TLL(path: str,
                        begi: int,
                        endi: int,
                        ptop: float = 0.,
                        pbot: float = 1000.,
                        plev: np.ndarray = default_levels) -> xr.DataArray:

  '''Reads in U,V, and Q variables from a netCDF file (that includes a dimension for time) and creates
  a seasonal cycle global map (12 months x lat x lon) output variable of integrated vapor transport vectors.

  NOTE: To use this function, geocat.f2py must be version 2022.4.0. Later versions do not work currently.

  Required parameters:
  ----------------------------
  path,begi,endi,ptop,pbot,plev

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  path: :class:'string'
        The file path to the netCDF file containing the output variables to be read in. In addition to
        the variable being read in by var1, file at path must include hyam, hybm, and ps if var1 contains
        dimension 'lev'.

  begi: :class:'int'
        The index of the first time to be included in analysis
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=0 starts the variable that is read
            at January. 

  endi: :class:'int'
        The index of the last time to be included in analysis. 
        NOTE: (endi-begi) / 12 must be a multiple of 12
        Ex. If netCDF file contains a 12-month climatology, begi=12 ends the variable that is read
            at December. NOTE: This follows python standards, so 12 here actually means index 11.

  ptop: :class:'float'
        Top boundary value of averaging region for vertical integration, in pressure (hPa).
        Ex. 50. -> vertical integration will sum from pbot to 50 hPa.

  pbot: :class:'float'
        Bottom boundary value of averaging region for vertical integration, in pressure (hPa).
        Ex. 1000. -> vertical integration will sum from 1000 hPa to ptop. 

  plev: :class:'numpy.ndarray'
         The new range of pressure levels with units hPa/mb, of which 'level' must be an available value. 
         This script will interpolate climate model output from hybrid to pressure levels using 
         geocat.comp's interp_hybrid_to_pressure function.
         NOTE: If NaN's present in variable while interpolating from hybrid to pressure, may see this
         error, which comes from metpy and can be ignored: 
          'UserWarning: Interpolation point out of data bounds encountered 
           warnings.warn('Interpolation point out of data bounds encountered')'
         Ex. np.arange(0,1050,50) -> pressure levels between 0 hPa and 1000 hPa with stride of 50 hPa  
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: :class:'xarray.DataArray'
          Seasonal cycle global variable [12 months,lat,lon]

  '''

  #-------------------------------------------------------------------------
  # Process parameters to ensure they will run within this function 
  #-------------------------------------------------------------------------

  # Ensure type of parameters
  path = str(path)
  begi = int(begi)
  endi = int(endi)
  ptop = float(ptop)
  pbot = float(pbot)
  plev = np.array(plev)

  #-----------------------------------
  # Read in variables and scale them 
  #-----------------------------------

  # Read in dataset 
  data = xr.open_dataset(path)

  # Read in variables
  u_data = data.U[begi:endi,:,:,:]
  v_data = data.V[begi:endi,:,:,:]
  q_data = data.Q[begi:endi,:,:,:]

  # Convert sigma coordinates to pressure coordinates
  simplefilter(action='ignore', category=UserWarning) # ignore unnecessary error from metpy
  hyam = data.hyam
  hybm = data.hybm
  psrf = data.PS
  p0Pa  = 100000.
  u_lev = xr.DataArray(interp_hybrid_to_pressure(data=u_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.),
                       dims=['time','plev','lat','lon'],coords=dict(time=data.time,plev=plev,lat=data.lat,lon=data.lon))
  v_lev = xr.DataArray(interp_hybrid_to_pressure(data=v_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.),
                       dims=['time','plev','lat','lon'],coords=dict(time=data.time,plev=plev,lat=data.lat,lon=data.lon))
  q_lev = xr.DataArray(interp_hybrid_to_pressure(data=q_data,ps=psrf,hyam=hyam,hybm=hybm,p0=p0Pa,new_levels=plev*100.),
                       dims=['time','plev','lat','lon'],coords=dict(time=data.time,plev=plev,lat=data.lat,lon=data.lon))

  # Calculate climatologies
  u_time = clmMonTLLL(u_lev)
  v_time = clmMonTLLL(v_lev)
  q_time = clmMonTLLL(q_lev)

  # Calculate UQ and VQ
  uq = u_time * q_time
  vq = v_time * q_time

  #------------------------------
  # Perform calculation for IVT 
  #------------------------------

  # Calculate pressure layer thickness (dp) with function from geocat.f2py, this requires two specific operations: 
  # 1) Finding index of plev that corresponds to ptop, this avoids error where ptop is not <= min(plev)
  ptop_ind = np.where(plev == ptop)[0][0] 
  # 2) Converting pressure to Pa and then converting it right back to hPa, this is a workaround for dpres_plevel to work
  dp = (dpres_plevel(pressure_levels=plev[ptop_ind:]*100.,pressure_surface=pbot*100.,pressure_top=ptop*100.)/100.) 

  # Weight each atmospheric layer by dp/g
  g = 9.8 # m/s^2
  UQ_wgt = np.zeros((12,len(plev),len(data.lat),len(data.lon)),dtype=float)
  VQ_wgt = np.zeros((12,len(plev),len(data.lat),len(data.lon)),dtype=float)
  for v in range(len(dp)):
   UQ_wgt[:,v,:,:] = (100.*dp[v]/g)*uq[:,v,:,:] # multiply by 100 for hPa to Pa conversion
   VQ_wgt[:,v,:,:] = (100.*dp[v]/g)*vq[:,v,:,:] # and this makes the final units = kg/(m*s) 

  # Sum across pressure levels
  UQ_time = xr.DataArray(UQ_wgt,dims=['time','plev','lat','lon'],coords=dict(time=np.arange(12),lat=data.lat,lon=data.lon),
                         ).sum(dim='plev')
  VQ_time = xr.DataArray(VQ_wgt,dims=['time','plev','lat','lon'],coords=dict(time=np.arange(12),lat=data.lat,lon=data.lon),
                         ).sum(dim='plev')

  #-------------------------------
  # Perform seasonal averaging
  #-------------------------------

  # Calculate averages for each month
  u_clim = clmMonTLL(UQ_time)
  v_clim = clmMonTLL(VQ_time)

  # Set NaN values
  u_clim = u_clim.where(u_clim < 1E10)
  v_clim = v_clim.where(v_clim < 1E10)

  return u_clim, v_clim
