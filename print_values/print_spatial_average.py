#######################################################################################################
#
# These functions take an xr.DataArray global variable input and calculate a spatial average from
# that variable, printing the result (with metadata) to the command line.  
#
#######################################################################################################

import numpy as np
import xarray as xr
from climate_analysis.coordinate_functions import lat_lon_index_array

# Previously defined variables
time_dim_names = ['time','month','months','year','years']

####################
# Functions
####################

#---------------------
# Print global average
#---------------------

def print_global_average(case: str,
                         variable: xr.DataArray,
                         units: str):
    '''Prints spatially weighted global average of a global variable.

    Required parameters:
    ----------------------
    case: :class:'string'
          Case name or variable name of the global variable.

    variable: :class:'xarray.DataArray'
              Global variable from which to calculate global average. 
              NOTE: Must be of size [lat x lon] and dims=['lat','lon']

    units: :class:'string'
           String defining units of the variable.'''

    lat_wgts = np.cos(np.deg2rad(variable.lat))

    gwgt = variable.weighted(lat_wgts)
    print(str(case)+' global avg ('+str(units)+') = '+str(np.array(gwgt.mean(('lon','lat'))))) 

#---------------------
# Print region average
#---------------------

def print_region_average(case: str,
                         variable: xr.DataArray,
                         units: str,
                         slat: float,
                         nlat: float,
                         wlon: float,
                         elon: float):
                         
    '''Prints spatially weighted average over a specified region from the initial input of a 
    global variable. Bounds shown for region are grid cell edges. 

    Required parameters:
    ----------------------
    case: :class:'string'
          Case name or variable name of the global variable.

    variable: :class:'xarray.DataArray'
              Global variable from which to calculate global average. 
              NOTE: Must be of size [lat x lon] and dims=['lat','lon']

    units: :class:'string'
           String defining units of the variable.

    slat: :class:'float'
          Latitude value at the southern boundary of the desired region.
          Matches to nearest index in variable's coordinate array.
          NOTE: + value denote °N and - values denote °S

    nlat: :class:'float'
          Latitude value at the northern boundary of the desired region.
          Matches to nearest index in variable's coordinate array.
          NOTE: + value denote °N and - values denote °S

    wlon: :class:'float'
          Longitude value at the western boundary of the desired region.
          Matches to nearest index in variable's coordinate array.
          NOTE: + values denote °E and - values denote °W 

    elon: :class:'float'
          Longitude value at the eastern boundary of the desired region.
          Matches to nearest index in variable's coordinate array.
          NOTE: + values denote °E and - values denote °W               '''

    # Determine if variable has a time dimension, all times will be printed if True
    if any(x in time_dim_names for x in variable.dims):
     time_dim = True
    else:
     time_dim = False

    # Use function to calculate index arrays for lat and lon
    latarray, lonarray = lat_lon_index_array(lat=variable.lat,lon=variable.lon,slat=slat,nlat=nlat,wlon=wlon,elon=elon)

    # Spatially weight variable given region boundaries and print result (with metadata) to command line
    lat_wgts = np.cos(np.deg2rad(variable.lat))
    if time_dim == True:
     rwgt = variable[:,latarray,lonarray].weighted(lat_wgts[latarray])
    else:
     rwgt = variable[latarray,lonarray].weighted(lat_wgts[latarray])

    ### Print region boundaries (values = edge of grid cells)

    # Create spacing to show lat/lon values as edge of grid cells rather than midpoint
    latadj = abs(float((variable.lat[1] - variable.lat[2]) / 2))
    lonadj = abs(float((variable.lon[1] - variable.lon[2]) / 2))

    # Discover direction of latitude array
    if variable.lat[0] > 0:
     lat_direction = 'N->S'
    elif variable.lat[0] < 0:
     lat_direction = 'S->N'

    # Discover treatment of longitude array
    if any(i > 180. for i in variable.lon):
     lon_treatment = 'degW = 180-360' 
    elif any(i < 0. for i in variable.lon):
     lon_treatment = 'degW = -180-0' 

    # If longitude array treats degW as 180-360, modify lonw,lone
    if lon_treatment == 'degW = 180-360':
     if wlon < 0.:    
      wlon = wlon + 360
     if elon < 0.:
      elon = elon + 360

    # Select coordinate indices
    lats = list(variable.lat.values).index(variable.lat.sel(lat=slat, method='nearest').lat)
    latn = list(variable.lat.values).index(variable.lat.sel(lat=nlat, method='nearest').lat)
    lonw = list(variable.lon.values).index(variable.lon.sel(lon=wlon, method='nearest').lon)
    lone = list(variable.lon.values).index(variable.lon.sel(lon=elon, method='nearest').lon)

    # Define values for each coordinate to print
    lats_val = np.array(variable.lat[lats])
    latn_val = np.array(variable.lat[latn])
    lonw_val = np.array(variable.lon[lonw])
    lone_val = np.array(variable.lon[lone])

    # Switch lons back to - = °W before printing, if necessary
    if lon_treatment == 'degW = 180-360':
     if lonw_val > 180.:
      lonw_val = lonw_val - 360 
     if lone_val > 180.:
      lone_val = lone_val - 360 

    print(str(case)+' region avg '+str(lats_val-latadj)+'-'+str(latn_val+latadj)+'('+str(lats)+','+str(latn)+'),'+str(
              lonw_val-lonadj)+'-'+str(lone_val+lonadj)+'('+str(lonw)+','+str(lone)+') ('+str(units)+') = '+str(np.array(rwgt.mean(('lon','lat')))))
   
#---------------------
# Print point average
#---------------------

def print_point_average(case: str,
                        variable: xr.DataArray,
                        units: str,
                        latpt: float,
                        lonpt: float):

    '''Prints average at a specified point from the initial input of a 
    global variable. Bounds shown for point are grid cell edges. 

    Required parameters:
    ----------------------
    case: :class:'string'
          Case name or variable name of the global variable.

    variable: :class:'xarray.DataArray'
              Global variable from which to calculate global average. 
              NOTE: Must be of size [lat x lon] and dims=['lat','lon']

    units: :class:'string'
           String defining units of the variable.

    latpt: :class:'float'
          Latitude value for the desired point. Matches to nearest index of midpoint in variable's coordinate array.
          NOTE: + value denote °N and - values denote °S

    lonpt: :class:'float'
          Longitude value for the desired point. Matches to nearest index of midpoint in variable's coordinate array.
          NOTE: + values denote °E and - values denote °W                 ''' 

    # Determine if variable has a time dimension, all times will be printed if True
    if any(x in time_dim_names for x in variable.dims):
     time_dim = True
    else:
     time_dim = False

    # Create spacing to show lat/lon values as edge of grid cells rather than midpoint
    latadj = abs(float((variable.lat[1] - variable.lat[2]) / 2))
    lonadj = abs(float((variable.lon[1] - variable.lon[2]) / 2))

    # Discover treatment of longitude array
    if any(i > 180. for i in variable.lon):
     lon_treatment = 'degW = 180-360'
    elif any(i < 0. for i in variable.lon):
     lon_treatment = 'degW = -180-0'

    # If longitude array treats degW as 180-360, modify lonw,lone
    if lon_treatment == 'degW = 180-360':
     if lonpt < 0.:
      lonpt = lonpt + 360
   
    # Select coordinate indices
    lat_ind = list(variable.lat.values).index(variable.lat.sel(lat=latpt, method='nearest').lat)
    lon_ind = list(variable.lon.values).index(variable.lon.sel(lon=lonpt, method='nearest').lon)

    # Define values for each coordinate to print 
    lat_val = np.array(variable.lat[lat_ind])
    lon_val = np.array(variable.lon[lon_ind])

    # Switch lon back to - = °W before printing, if necessary
    if lon_treatment == 'degW = 180-360':
     if lon_val > 180.:
      lon_val = lon_val - 360

    if time_dim == True:
     print(str(case)+' point avg '+str(lat_val-latadj)+'-'+str(lat_val+latadj)+'('+str(lat_ind)+'),'+str(
                                      lon_val-lonadj)+'-'+str(lon_val+lonadj)+'('+str(lon_ind
                                  )+') ('+str(units)+') = '+str(np.array(variable[:,lat_ind,lon_ind])))
    else:
     print(str(case)+' point avg '+str(lat_val-latadj)+'-'+str(lat_val+latadj)+'('+str(lat_ind)+'),'+str(
                                      lon_val-lonadj)+'-'+str(lon_val+lonadj)+'('+str(lon_ind
                                  )+') ('+str(units)+') = '+str(np.array(variable[lat_ind,lon_ind])))






