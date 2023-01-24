#######################################################################################################
# 
# Functions for processing and indexing coordinates
#
#######################################################################################################

import numpy as np
import xarray as xr

#######################################################################################################
# Calculate index arrays for lat and lon based on given value ranges
#######################################################################################################

def lat_lon_index_array(lat = xr.DataArray,
                        lon = xr.DataArray,
                        slat = float,
                        nlat = float,
                        wlon = float,
                        elon = float):

   ''' Input lat/lon boundaries for a region (can be 1D along lat or lon) and return arrays of lat and
   lon indices corresponding to the coordinate range. This function can handle lat/lon coordinate 
   arrays of different sizes and sequential order. ''' 

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

   # Coordinate index arrays: latitude 
   if lat_direction == 'N->S':
    latarray = np.arange(latn,lats+1,1)
   elif lat_direction == 'S->N':
    latarray = np.arange(lats,latn+1,1)

   # Coordinate index arrays: longitude

   if lon_treatment == 'degW = 180-360':
    # if lonarray crosses prime meridian
    if lon[lonw] > lon[lone]:
     lonarray_west = np.arange(lonw,list(lon.values).index(lon.sel(lon=lon[-1]).lon)+1,1) # from index lonw to lon=0
     lonarray_east = np.arange(0,list(lon.values).index(lon.sel(lon=lon[lone]).lon)+1,1)  # from lon=0 to index lone
     lonarray = np.append(lonarray_west,lonarray_east)
    # if lonarray crosses international date line
    elif lon[lonw] < 180. and lon[lone] > 180.:
     lonarray_west = np.arange(lonw,list(lon.values).index(lon.sel(lon=180, method='nearest').lon),1) # from index lonw to lon=180
     lonarray_east = np.arange(list(lon.values).index(lon.sel(lon=180, method='nearest').lon),lone+1,1) # from lon=180 to index lone
     lonarray = np.append(lonarray_west,lonarray_east)
    elif lon[lonw] < lon[lone]:
     lonarray = np.arange(lonw,lone+1)
    elif lonw == lone:
     lonarray = lonw

   if lon_treatment == 'degW = -180-0':
    # if lonarray crosses prime meridian
    if lon[lonw] < 0. and lon[lone] > 0.:
     lonarray_west = np.arange(lonw,list(lon.values).index(lon.sel(lon=0, method='nearest').lon),1) # from index lonw to lon=0
     lonarray_east = np.arange(list(lon.values).index(lon.sel(lon=0, method='nearest').lon),lone+1,1) # from lon=0 to index lone      
     lonarray = np.append(lonarray_west,lonarray_east)
    # if lonarray crosses international date line
    elif lon[lonw] > 0. and lon[lone] < 0.:
     lonarray_west = np.arange(lonw,list(lon.values).index(lon.sel(lon=180, method='nearest').lon),1) # from index lonw to lon=180
     lonarray_east = np.arange(list(lon.values).index(lon.sel(lon=180, method='nearest').lon),lone+1,1) # from index lonw to lon=180
     lonarray = np.append(lonarray_west,lonarray_east)
    elif lonw < lone:
     lonarray = np.arange(lonw,lone+1)
    elif lonw == lone:
     lonarray = lonw

   return latarray, lonarray








