#######################################################################################################
# 
# These functions print values for each water tagged region to the command line or to an Excel file 
#
#######################################################################################################

import numpy as np
import xarray as xr
import pandas as pd
import climate_analysis.seas_avg_LL as seasavg
from watertagging.seas_avg_LL_watertags import seasavg_watertagging_vars
from climate_analysis.coordinate_functions import lat_lon_index_array
from warnings import simplefilter

# Default arrays
default_wgt_mon = xr.DataArray(np.ones(12),dims=['time']).astype(float)

#######################################################################################################

def print_watertag_values(precip: xr.DataArray, d18Op: xr.DataArray,
                          precip_sum: float, d18Op_sum: float,
                          precip_reg_gbl: float, d18Op_reg_gbl: float,
                          lat: xr.DataArray, lon: xr.DataArray,
                          case: str, tagnames: list, season: str, 
                          slat: float = None, nlat: float = None, 
                          wlon: float = None, elon: float = None):

  '''Print values for each tag region for a specified case to the main screen 

  Parameters
  --------------------------------------------------------------------------------------------------
  precip
   class: 'xarray.DataArray', variable containing averaged precipitation values for each tag region

  d18Op
   class: 'xarray.DataArray', variable containing averaged d18Op values for each tag region

  precip_sum
   class: 'float', sum of all values in 'precip' as a single float value 

  d18Op_sum
   class: 'float', sum of all values in 'd18Op' as a single float value 

  precip_reg_gbl
   class: 'float', average precipitation value at defined region/point from global variable
                   (rather than tag region variable)

  d18Op_reg_gbl
   class: 'float', average 18Op value at defined region/point from global variable (rather than 
                   tag region variable)

  lat
   class: 'xarray.DataArray', latitude coordinate array

  lon
   class: 'xarray.DataArray', longitude coordinate array

  case
   class: 'string', name of the case corresponding to the printed values

  tagnames
   class: 'list', list of strings indicating each tag region's name

  slat, nlat, wlon, elon
   class: 'float', southern/northern latitude and western/eastern longitude defining region/point 
  
  '''

  #----------------------------
  # Define region/point 
  #----------------------------

  # Process coordinates
  latadj = abs(np.float64((lat[2]-lat[1])/2))
  lonadj = abs(np.float64((lon[2]-lon[1])/2))

  # Discover direction of latitude array
  if lat[0] > 0:
   lat_direction = 'N->S'
  elif lat[0] < 0:
   lat_direction = 'S->N'

  # Discover treatment of longitude array
  if any(i > 180. for i in lon):
   lon_treatment = 'degW = 180-360'
  elif any(i < 0. for i in lon):
   lon_treatment = 'degW = -180-0'

  # If longitude array treats degW as 180-360, modify lonw,lone
  if lon_treatment == 'degW = 180-360':
   if wlon < 0.:
    wlon = wlon + 360
   if elon < 0.:
    elon = elon + 360

  # Select coordinate indices
  lats = list(lat.values).index(lat.sel(lat=slat, method='nearest').lat)
  latn = list(lat.values).index(lat.sel(lat=nlat, method='nearest').lat)
  lonw = list(lon.values).index(lon.sel(lon=wlon, method='nearest').lon)
  lone = list(lon.values).index(lon.sel(lon=elon, method='nearest').lon)

  # Define perimeter
  rlats = np.float64(lat[lats]-latadj)
  rlatn = np.float64(lat[latn]+latadj)
  rlonw = np.float64(lon[lonw]-lonadj)
  rlone = np.float64(lon[lone]+lonadj)
  if rlonw > 180:
    rlonw = rlonw - 360
  if rlone > 180:
    rlone = rlone - 360

  #-----------------------------------------------------------------------------------
  # Convert to pandas df and then print values in neat columns 
  # This needs to happen after defining region/point 
  #-----------------------------------------------------------------------------------

  # First, defining precipitation and d18Op values for each tag
  df1 = pd.DataFrame({'Tag': tagnames, 'Precip': precip, 'd18Op': d18Op})

  # Next, adding a row for the sum of all tagged regions together
  df2 = pd.concat([df1,pd.DataFrame({'Tag': 'Region summed', 'Precip': np.float64([precip_sum]), 
                                     'd18Op': np.float64([d18Op_sum])})]).reset_index(drop=True)

  # Lastly, adding a row for the average calculated from the global variable for comparison with the sum of all tagged regions
  df  = pd.concat([df2,pd.DataFrame({'Tag': 'Region avg with global variable', 'Precip': [float(precip_reg_gbl)], 
                                     'd18Op': [float(d18Op_reg_gbl)]})]).reset_index(drop=True)

  # Print case and region bounds
  print(season+' '+str(case)+' values for ('+str(rlats)+' to '+str(rlatn)+'°N, '+str(rlonw)+' to '+str(rlone)+'°E) are ...')

  # Print dataframe (all printed values are set to 6 decimal  places)
  pd.options.display.float_format = '{:.6f}'.format  # setting all printed values to 6 decimal places
  print(df.to_string())

#######################################################################################################################

def monthly_watertag_values_to_excel(CASES: str, cases: str, begi: int, endi: int, tagnames: list,tagcodes: list, 
                                     slat: float = None, nlat: float = None,   
                                     wlon: float = None, elon: float = None,
                                     folderpath: str = '', reg_name: str = '', filesuf: str = '.xlsx'):

  '''Print monthly values for precipitation and d18Op for each tag region for all cases to an Excel file. 

  Parameters
  ---------------------------------------------------------------------------------------------------------------
  CASES 
   class: 'list', List of strings with each element containing the full file path to each water tagged 
                  simulations output file. Example: CASES = ['/usr/local/watertag1.nc','/usr/local/watertag2.nc']
 
  cases  
   class: 'list', List of strings containing the name of case(s). Example: cases = ['case1','case2']

  begi
   class: 'int', The index of the first time to be included in analysis. NOTE: (endi-begi) / 12 must be a multiple 
                 of 12. Example: If netCDF file contains a 12-month climatology, begi=0 starts the variable that is 
                 read at January. 

  endi
   class: 'int', The index of the last time to be included in analysis. 
 
  tagnames
   class: 'list', List of strings containing the long-form name of each tagged region. Example: tagnames = 
                  ['Antarctica','North America','South America']

  tagcodes
   class: 'list', List of strings containing code name of each tagged region. Example: ['ANTA','NAM','SAM']

  slat, nlat, wlon, elon
   class: 'float', Bounding latitude and longitude coordinates of region for which water tagging results will be
                   calculated. 

  folderpath (Default = '')
   class: 'string', String of the path to the folder in which the output file will be placed.

  reg_name (Default = '')
   class: 'string', Name of region for which water tagging results are calculated, added to output file.

  filesuf (Default = '.xlsx')
   class: 'string', Suffix of output file.  

  '''
  #-----------------------------------------------------------------
  # Create multidimensional variable(s) [# of cases x lat x lon] 
  #-----------------------------------------------------------------

  # List of strings for month names
  months_str = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

  # Dataset for dimensions
  ds = xr.open_dataset(CASES[0])

  # Global variables [case,lat,lon]
  prect_global = xr.DataArray(None,dims=['case','months','lat','lon'],
                                   coords=dict(case=cases,months=months_str,
                                               lat=ds.lat,lon=ds.lon)).astype(float)
  d18Op_global = xr.DataArray(None,dims=['case','months','lat','lon'],
                                   coords=dict(case=cases,months=months_str,
                                               lat=ds.lat,lon=ds.lon)).astype(float)

  # Variables by tag [case,tag,lat,lon]
  Pi_by_tag = xr.DataArray(None,dims=['case','tag','months','lat','lon'],
                                coords=dict(case=cases,tag=tagcodes,months=months_str,
                                            lat=ds.lat,lon=ds.lon)).astype(float)
  d18Opsink_by_tag = xr.DataArray(None,dims=['case','tag','months','lat','lon'],
                                       coords=dict(case=cases,tag=tagcodes,months=months_str,
                                                   lat=ds.lat,lon=ds.lon)).astype(float)
  d18Opwt_by_tag = xr.DataArray(None,dims=['case','tag','months','lat','lon'],
                                       coords=dict(case=cases,tag=tagcodes,months=months_str,
                                                   lat=ds.lat,lon=ds.lon)).astype(float)
  # Average value for each tag [case,tag]
  prect = xr.DataArray(None,dims=['case','tag','months'],coords=dict(case=cases,tag=tagcodes,
                                                                     months=months_str))
  d18Op = xr.DataArray(None,dims=['case','tag','months'],coords=dict(case=cases,tag=tagcodes,
                                                                     months=months_str))

  # Global variable's average of the selected region
  prect_reg = xr.DataArray(None,dims=['case','months'],coords=dict(case=cases,months=months_str))
  d18Op_reg = xr.DataArray(None,dims=['case','months'],coords=dict(case=cases,months=months_str))

  # Sum of all tagged region values
  prect_sum = xr.DataArray(None,dims=['case','months'],coords=dict(case=cases,months=months_str))
  d18Op_sum = xr.DataArray(None,dims=['case','months'],coords=dict(case=cases,months=months_str))

  #-----------------------------------
  # Loop through each case and month
  #-----------------------------------
  
  for i in range(len(CASES)):
  
   print('Working on Excel for '+str(cases[i]))

   # Loop through each month

   for moni in range(12):

    mon = [moni]

    #-----------------------------------------------------------------
    # Global variables: prect and d18Op
    #-----------------------------------------------------------------

    # PRECT
    prect_global[i,moni,:,:] = seasavg.seasavg_prect_LL(path=CASES[i],begi=begi,endi=endi,months=mon,
                                                   wgt_mon=[1.])

    # d18Op
    d18Op_global[i,moni,:,:] = seasavg.seasavg_rainiso_LL(iso_type='d18O',path=CASES[i],begi=begi,endi=endi,
                                                         ptiny=1.E-18,months=mon,wgt_mon=[1.])

    #-----------------------------
    # Loop through each tag
    #-----------------------------

    for tag in range(len(tagnames)):

     Pi_by_tag[i,tag,moni,:,:], d18Opsink_by_tag[i,tag,moni,:,:] = seasavg_watertagging_vars(

          # Inputs for functions:
          tagcode=tagcodes[tag],months=mon,path=CASES[i],begi=begi,endi=endi,wgt_mon=[1.])

    #-----------------------------------
    # Perform calculation on variables
    #-----------------------------------

    # Initialize latitude weights
    lat = ds.lat
    lon = ds.lon
    lat_wgts = np.cos(np.deg2rad(lat))

    # Create list arrays of lat and lon values for regional averaging
    latarray, lonarray = lat_lon_index_array(lat=ds.lat,lon=ds.lon,slat=slat,nlat=nlat,
                                                                    wlon=wlon,elon=elon)

    # Loop through each tag to calculate prect and d18Op tag region values 
    for tag in range(len(tagnames)):

     # Use water tagging equation: d18Opwt@gc = d18Op_tag@gc * ( total16Op_tag@gc / prect_global@gc )
     d18Opwt_by_tag[i,tag,moni,:,:] = d18Opsink_by_tag[i,tag,moni,:,:] * ( 
                                        Pi_by_tag[i,tag,moni,:,:] / prect_global[i,moni,:,:] )

     # Area weight prect and d18Op for tagged region average
     prect_wgt         = Pi_by_tag[i,tag,moni,latarray,lonarray].weighted(lat_wgts[latarray])
     prect[i,tag,moni] = prect_wgt.mean(('lon','lat'))
     d18Op_wgt         = d18Opwt_by_tag[i,tag,moni,latarray,lonarray].weighted(lat_wgts[latarray])
     d18Op[i,tag,moni] = d18Op_wgt.mean(('lon','lat'))

    # Calculate tagged region average with global variables to check against tagged variables
    prect_global_wgt  = prect_global[i,moni,latarray,lonarray].weighted(lat_wgts[latarray])
    prect_reg[i,moni] = prect_global_wgt.mean(('lon','lat'))
    d18Op_global_wgt  = d18Op_global[i,moni,latarray,lonarray].weighted(lat_wgts[latarray])
    d18Op_reg[i,moni] = d18Op_global_wgt.mean(('lon','lat'))

    # Sum tag values together
    prect_sum[i,moni] = np.sum(prect[i,:,moni])
    d18Op_sum[i,moni] = np.sum(d18Op[i,:,moni])

  #--------------------------------------
  # Define region/point for output file 
  #--------------------------------------

  # Process coordinates
  latadj = abs(np.float64((lat[2]-lat[1])/2))
  lonadj = abs(np.float64((lon[2]-lon[1])/2))

  # Discover direction of latitude array
  if lat[0] > 0:
   lat_direction = 'N->S'
  elif lat[0] < 0:
   lat_direction = 'S->N'

  # Discover treatment of longitude array
  if any(i > 180. for i in lon):
   lon_treatment = 'degW = 180-360'
  elif any(i < 0. for i in lon):
   lon_treatment = 'degW = -180-0'

  # If longitude array treats degW as 180-360, modify lonw,lone
  if lon_treatment == 'degW = 180-360':
   if wlon < 0.:
    wlon = wlon + 360
   if elon < 0.:
    elon = elon + 360

  # Select coordinate indices
  lats = list(lat.values).index(lat.sel(lat=slat, method='nearest').lat)
  latn = list(lat.values).index(lat.sel(lat=nlat, method='nearest').lat)
  lonw = list(lon.values).index(lon.sel(lon=wlon, method='nearest').lon)
  lone = list(lon.values).index(lon.sel(lon=elon, method='nearest').lon)

  # Define perimeter
  rlats = np.float64(lat[lats]-latadj)
  rlatn = np.float64(lat[latn]+latadj)
  rlonw = np.float64(lon[lonw]-lonadj)
  rlone = np.float64(lon[lone]+lonadj)
  if rlonw > 180:
    rlonw = rlonw - 360
  if rlone > 180:
    rlone = rlone - 360

  #---------------------------------------
  # Construct Excel file from variables 
  #---------------------------------------

  # Ignore unnecessary warning about w.save() 
  simplefilter(action="ignore", category=FutureWarning)

  # Write Excel file name first
  w = pd.ExcelWriter('./'+str(folderpath)+'/Monthly_watertagging_values_for_'+str(reg_name)+   \
                                          '_lat'+str(round(rlats,2))+'to'+str(round(rlatn,2))+ \
                                          '_lon'+str(round(rlonw,2))+'to'+str(round(rlone,2))+ \
                                           str(filesuf))

  # Loop through each case and set values as a new sheet name
  for i in range(len(cases)):

   # Each tagged regions' values for precipitation and d18Op
   df1 = pd.DataFrame({'Tag':tagnames,'P_Jan':prect[i,:,0],'P_Feb':prect[i,:,1],'P_Mar':prect[i,:,2],
                                      'P_Apr':prect[i,:,3],'P_May':prect[i,:,4],'P_Jun':prect[i,:,5],
                                      'P_Jul':prect[i,:,6],'P_Aug':prect[i,:,7],'P_Sep':prect[i,:,8],
                                      'P_Oct':prect[i,:,9],'P_Nov':prect[i,:,10],'P_Dec':prect[i,:,11],
                                      'O_Jan':d18Op[i,:,0],'O_Feb':d18Op[i,:,1],'O_Mar':d18Op[i,:,2],
                                      'O_Apr':d18Op[i,:,3],'O_May':d18Op[i,:,4],'O_Jun':d18Op[i,:,5],
                                      'O_Jul':d18Op[i,:,6],'O_Aug':d18Op[i,:,7],'O_Sep':d18Op[i,:,8],
                                      'O_Oct':d18Op[i,:,9],'O_Nov':d18Op[i,:,10],'O_Dec':d18Op[i,:,11]})

   # Sums of all tagged regions together for each month
   df2 = pd.concat([df1,pd.DataFrame({'Tag':'Region summed',
                                      'P_Jan':np.float64([prect_sum[i,0]]), 'P_Feb':np.float64([prect_sum[i,1]]),
                                      'P_Mar':np.float64([prect_sum[i,2]]), 'P_Apr':np.float64([prect_sum[i,3]]), 
                                      'P_May':np.float64([prect_sum[i,4]]), 'P_Jun':np.float64([prect_sum[i,5]]),
                                      'P_Jul':np.float64([prect_sum[i,6]]), 'P_Aug':np.float64([prect_sum[i,7]]),
                                      'P_Sep':np.float64([prect_sum[i,8]]), 'P_Oct':np.float64([prect_sum[i,9]]),
                                      'P_Nov':np.float64([prect_sum[i,10]]),'P_Dec':np.float64([prect_sum[i,11]]),
                                      'O_Jan':np.float64([d18Op_sum[i,0]]), 'O_Feb':np.float64([d18Op_sum[i,1]]),
                                      'O_Mar':np.float64([d18Op_sum[i,2]]), 'O_Apr':np.float64([d18Op_sum[i,3]]),
                                      'O_May':np.float64([d18Op_sum[i,4]]), 'O_Jun':np.float64([d18Op_sum[i,5]]),
                                      'O_Jul':np.float64([d18Op_sum[i,6]]), 'O_Aug':np.float64([d18Op_sum[i,7]]),
                                      'O_Sep':np.float64([d18Op_sum[i,8]]), 'O_Oct':np.float64([d18Op_sum[i,9]]),
                                      'O_Nov':np.float64([d18Op_sum[i,10]]),'O_Dec':np.float64([d18Op_sum[i,11]])
                                      })]).reset_index(drop=True)

   # Average of specified region calculated from global variable (as opposed to individual tag variables)
   # This should be similar to the sumsm of all tagged regions together but will not match exactly
   df = pd.concat([df2,pd.DataFrame({'Tag':'Region avg with global variable',
                                     'P_Jan':[float(prect_reg[i,0])], 'P_Feb':[float(prect_reg[i,1])],
                                     'P_Mar':[float(prect_reg[i,2])], 'P_Apr':[float(prect_reg[i,3])],
                                     'P_May':[float(prect_reg[i,4])], 'P_Jun':[float(prect_reg[i,5])],
                                     'P_Jul':[float(prect_reg[i,6])], 'P_Aug':[float(prect_reg[i,7])],
                                     'P_Sep':[float(prect_reg[i,8])], 'P_Oct':[float(prect_reg[i,9])],
                                     'P_Nov':[float(prect_reg[i,10])],'P_Dec':[float(prect_reg[i,11])],
                                     'O_Jan':[float(d18Op_reg[i,0])], 'O_Feb':[float(d18Op_reg[i,1])],
                                     'O_Mar':[float(d18Op_reg[i,2])], 'O_Apr':[float(d18Op_reg[i,3])],
                                     'O_May':[float(d18Op_reg[i,4])], 'O_Jun':[float(d18Op_reg[i,5])],
                                     'O_Jul':[float(d18Op_reg[i,6])], 'O_Aug':[float(d18Op_reg[i,7])],
                                     'O_Sep':[float(d18Op_reg[i,8])], 'O_Oct':[float(d18Op_reg[i,9])],
                                     'O_Nov':[float(d18Op_reg[i,10])],'O_Dec':[float(d18Op_reg[i,11])]
                                     })]).reset_index(drop=True) 

   # Save each case as a new sheet name in the Excel file
   df.to_excel(w,sheet_name=cases[i],index=False)
   w.save()


