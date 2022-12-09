#########################################################################################################
#
# Plot seasonally averaged variable(s) from a number of climate model netCDF output files
#
# This script utilizes functions from map_avg_functions.py and map_plot_functions.py to create plots
#
#########################################################################################################

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cmaps
import geocat.viz.util as gv
import climate_analysis.seas_cycle_TLL as seascyc
import climate_analysis.seas_avg_LL as seasavg
from climate_analysis.mon_wgt_avg import mon_wgt_avg
import print_values.print_spatial_average as print_spatial_average
import plotting.plot_map_avg_functions as plotmap

#########################################################################################################
#
# Specifications for plot are made here 
#
#########################################################################################################

#------------------------------
# Specify data path variables
#------------------------------

# Component model to specify in case strings
model = 'cam' 

# File paths and names for each case
CASES = ['f.e12.F_1850_CAM5.wiso.f19.0ka.002.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21ka.modern.d18Osw.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21kaGHG.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21kaGlac.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21kaSL.002.watertags.'+model+'.h0.0006-0025.climo.nc']
cases = ['0ka',
         '21ka',
         '21kaGHG',
         '21kaGLAC',
         '21kaSL']

# Anything extra to add to output file name?
extra_name = 'last20yrs_wgtmon'

#--------------------------------
# Seasonal averaging variables
#--------------------------------

# Indices to define range of time dimension read in from files; if climatology file (Jan-Dec), specify begi=0 and endi=12
begi = 0
endi = 12

# Season to average over, indices corresponding to individual months, season string will automatically populate later 
MON = [0,1,2,3,4,5,6,7,8,9,10,11]
#MON = [5,6,7]

# Reference list for indices
# Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec
#   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11

#-----------------------------------
# Select variable(s) to plot
#-----------------------------------

# Description of categories 
# 'var': Regular variable (time x lat x lon), ex. TREFHT
# 'varlev': Variable with atm level dimension (time x lev x lat x lon), ex. T
# 'prect': Precipitation variable, automatically creates PRECC+PRECL as mm/day
# 'soilvar': Land variable that includes soil level dimension (time x levgrnd x lat x lon)
# 'rainiso': Precipitation isotope variable, automatically creates variable from specified parameters
# 'soiliso': Soil isotope variable that includes soil level dimension (time x levgrnd x lat x lon)
# 'vaporiso': Water vapor isotope variable, automatically creates variable at specified atmospheric level
# 'isoroot': Soil isotope variable weighted by rooting depth fraction

# Choose which variable function to specify
varfunc = 'prect'       # Options: 'var','varlev','prect','soilvar','rainiso','soiliso','vaporiso','isoroot'
#.......................................................................................................................
# If varfunc is 'var', 'varlev', or 'soilvar'; only 'UNITS' required for 'prect'
VAR1 = 'T'   
VAR2 = None   # options: None or '[name]'
MATH = None   # options: None or one of 'add'->VAR1+VAR2, 'sub'->VAR1-VAR2, 'mul'->VAR1*VAR2, 'div'->VAR1/VAR2, 'exp'->VAR1**VAR2
MULT = 1.0    # Scale variable(s) by scalar value
UNITS = 'mm/day'
LEVEL = 500            # Plot at a specific atmospheric pressure level? Specify value in hPa (ex. 850) 
SOIL_LEV = [0,1,2]     # If you need to plot a soil level, takes average of list, ex. [0] or [0,1,2] 
#......................................................................................................................
# If varfunc is an isotopic variable ('rainiso', 'soiliso', 'vaporiso', or 'isoroot')
ISO       = 'd18O'          # Options: 'd18O', 'dHDO', or 'dexcess'
IUNITS    = 'per mil'
ptiny     = 1.E-18
ILEVEL    = 700             # Plot at a specific atmospheric pressure level? Specify value in hPa (ex. 850) 
ISOIL_LEV = [0,1,2]         # Used in 'soiliso', takes average of list of soil level indices, ex. [0] or [0,1,2] 
ROOTWGT   = 'ROOTR_COLUMN'  # Used in 'isoroot', variable used for weighting by root depth fraction
#......................................................................................................................
#-------------------------------------------------------------------------------------
# Atmospheric pressure level variables, if necessary 
#-------------------------------------------------------------------------------------

# Define pressure levels with this array, ex. Pressure array goes from 0 hPa to 1000 hPa by 50 hPa 
plev = np.arange(0,1050,50)
lev_ind  = int(np.where(np.array(plev) == LEVEL)[0])   # lev_ind converts LEVEL into the corresponding index of plev
lev_name = str(LEVEL)+'hPa'

#---------------------------------------------------------------
# Specify vectors to overlay on plot, if necessary
#---------------------------------------------------------------

# Overlay a vector? If True, what type and level of the atmosphere?
overlay_vector = True  
overlay_type   = 'wind'    # 'wind','IVT'

# Variables for overlay_type == 'wind', uses 'plev' from above
WIND_LEVEL = 800   # Integer, in hPa
WIND_UNITS = 'm/s' # Text string

# Variables for overlay_type == 'IVT', uses 'plev' from above
top_lev = 50.   # in hPa
bot_lev = 1018. # in hPa

# Define file paths to read in for overlaying vector
if overlay_vector == True:
 if model != 'cam':
  CASES_vec = []
  for w in CASES:
   new_string = w.replace(model,"cam")
   CASES_vec.append(new_string)
 else:
  CASES_vec = CASES

#---------------------------------------------------------------------------------------
# Specify monthly weights, if necessary
#---------------------------------------------------------------------------------------

# Array with 12 months of weights 
wgt_by_mon = xr.DataArray(None,dims=['case','month'],coords=dict(case=cases,month=np.arange(1,13,1))).astype(float)

# Array with n (based on MON) months of weights, not specifying 'time' or 'month' coord so weighting with var_time will work
wgt_mon = xr.DataArray(np.zeros([len(cases),len(MON)]),dims=['case','time'],coords=dict(case=cases)).astype(float)

### Enter monthly weights for each case here
#...............................................................................................................................
# 0ka
wgt_by_mon[0,:] = np.array([ 0.08493151, 0.076712325, 0.08493151, 0.08219178, 0.08493151, 0.08219178, 0.08493151, 0.08493151, 0.08219178, 0.08493151, 0.08219178, 0.08493151 ])
# 21ka
wgt_by_mon[1:,:] = np.array([ 0.084931507, 0.076712329, 0.084931507, 0.082191781, 0.084931507, 0.082191781, 0.084931507, 0.084931507, 0.082191781, 0.084931507, 0.082191781, 0.084931507 ])

#...............................................................................................................................

# Modify depending if averaging over entire year or only a select few months
for i in range(len(cases)):
 if len(MON) == 12:
  wgt_mon[i,:] = wgt_by_mon[i,:]
 else:
  wgt_mon[i,:] = wgt_by_mon[i,MON]/np.sum(wgt_by_mon[i,MON])

#--------------------------------------------------------------
# Do you want to print out the average over a specified region?
#--------------------------------------------------------------

# Latitude:  +=north, -=south
# Longitude: 0-180=°E, -180-0=°W

# Do you want a global average?
Global = True 

# Do you want an average over a specific region?
Region = True 

southlat = 14.0
northlat = 18.0
westlon  = -92.5
eastlon  = -90.0

# Do you want an average at a specific grid point?
Point = True  

latpoint = 33.0
lonpoint = 255.0

#-----------------------------------------------------------
# Specify contour levels for contour and difference values
#-----------------------------------------------------------

# Manual Levels for contour plotting (True = set manual levels)
ManLevCntr = True  
c_hival = 10.0         # high value
c_loval = 0.0          # low value
c_spval = 0.5          # spacing
c_tkstd = 2.0          # tick/label stride
c_extnd = 'max'        # where to put triangles; 'both,'neither','max','min'

# Manual Levels for difference plotting (True = set manual levels)
ManLevDiff = True  
d_hival =  2.5         # high value
d_loval = -2.5         # low value
d_spval =  0.25        # spacing
d_tkstd =  1.0         # tick/label stride
d_extnd = 'both'       # where to put triangles; 'both,'neither','max','min'

#---------------------------------------------
# Zoom into any world region in particular?
#---------------------------------------------

# True=entire world, False=zoomed in to coordinate values in second block
World = False

if World == True:
 LatMin = -90
 LatMax = 90
 LonMin = -180
 LonMax = 180
else:
 LatMin =    0.0     # negative values = °S
 LatMax =   40.0     # positive values = °N
 LonMin = -120.0     # negative values = °W
 LonMax =  -60.0     # positive values = °E

#-------------------------------------------------
# Other plotting specifications are set here
#-------------------------------------------------

Ind_plots    = True                    # Include individual map plots in output file
proj         = ccrs.PlateCarree()      # Map projection
Contour_type = 'AreaFill'              # 'RasterFill' or 'AreaFill'
ColCntr      = cmaps.cmp_haxby_r       # Contour plot color table 
ColDiff      = cmaps.BlueYellowRed_r   # "_r" at end reverses color table 
folderpath   = 'pdfs'                  # folder to output file to
filesuf      = '.pdf'                  # type of output file

#########################################################################################################
#
# Read in and calculate seasonally averaged variables here 
#
#########################################################################################################

#-----------------------------------------------------------------
# Create multidimensional variable(s) [# of cases x lat x lon] 
#-----------------------------------------------------------------

datasize = xr.open_dataset(CASES[0])
var_avg_by_case = xr.DataArray(None,dims=['case','lat','lon'],coords=dict(case=cases,lat=datasize.lat,lon=datasize.lon)).astype(float)
if overlay_vector == True:
 vecsize = xr.open_dataset(CASES_vec[0])
 u_avg_by_case = xr.DataArray(None,dims=['case','lat','lon'],coords=dict(case=cases,lat=vecsize.lat,lon=vecsize.lon)).astype(float)
 v_avg_by_case = xr.DataArray(None,dims=['case','lat','lon'],coords=dict(case=cases,lat=vecsize.lat,lon=vecsize.lon)).astype(float)

#--------------------------------------------------------------------
# Process a few variables before starting to loop through each case
#--------------------------------------------------------------------

# Read in dimension variables
time = datasize.time
lat = datasize.lat
lon = datasize.lon

# Make season string automatically from season indices
month_by_letter = ['J','F','M','A','M','J','J','A','S','O','N','D']
if len(MON) == 12:   # annual average
 season = 'ANN'
else:
 season = ''.join([month_by_letter[i] for i in MON])

#------------------------------------------------------------
# Loop through each case and calculate seasonal average 
#------------------------------------------------------------

for i in range(len(cases)):

 # CONTOURS
 if varfunc not in ['rainiso','soiliso','vaporiso','isoroot']:

  if varfunc == 'var':
   var_avg_by_case[i,:,:] = seasavg.seasavg_var_LL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,math=MATH,mult=MULT,
                                                    months=MON,wgt_mon=wgt_mon[i,:])
  if varfunc == 'varlev':
   var_avg_by_case[i,:,:] = seasavg.seasavg_var_LL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,math=MATH,mult=MULT,
                                                    level=LEVEL,plev=plev,months=MON,wgt_mon=wgt_mon[i,:]) 
  if varfunc == 'prect':
   var_avg_by_case[i,:,:] = seasavg.seasavg_prect_LL(path=CASES[i],begi=begi,endi=endi,months=MON,wgt_mon=wgt_mon[i,:]) 

  if varfunc == 'soilvar':
   var_avg_by_case[i,:,:] = seasavg.seasavg_soilvar_LL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,math=MATH,
                                                        soillev=SOIL_LEV,mult=MULT,months=MON,wgt_mon=wgt_mon[i,:])

 else:

  # If isotopic variable, calculate seasonally averaged variable (due to weighting by precip, etc.)
  if varfunc == 'rainiso':
   var_avg_by_case[i,:,:] = seasavg.seasavg_rainiso_LL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,ptiny=ptiny,
                                                        months=MON,wgt_mon=wgt_mon[i,:])
  if varfunc == 'soiliso':
   var_avg_by_case[i,:,:] = seasavg.seasavg_soiliso_LL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,ptiny=ptiny,
                                                        soillev=SOIL_LEV,months=MON,wgt_mon=wgt_mon[i,:])
  if varfunc == 'vaporiso':
   var_avg_by_case[i,:,:] = seasavg.seasavg_vaporiso_LL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,ptiny=ptiny,
                                                         level=LEVEL,plev=plev,months=MON,wgt_mon=wgt_mon[i,:])
  if varfunc == 'isoroot':
   var_avg_by_case[i,:,:] = seasavg.seasavg_isoroot_LL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,ptiny=ptiny,
                                                         rootwgt=ROOTWGT,months=MON,wgt_mon=wgt_mon[i,:])

 print(var_avg_by_case[i,:,:])

 # VECTORS 
 if overlay_vector == True:
   
  if overlay_type == 'wind':
   u_avg_by_case[i,:,:], v_avg_by_case[i,:,:] = seasavg.seasavg_wind_vec_LL(path=CASES_vec[i],begi=begi,endi=endi,
                                                              level=WIND_LEVEL,plev=plev,months=MON,wgt_mon=wgt_mon[i,:])

#########################################################################################################
#
# Loop through cases and print spatially averaged values 
#
#########################################################################################################

#-----------------------------
# Determine units
#-----------------------------

if varfunc == 'var' or varfunc == 'varlev' or varfunc == 'prect' or varfunc == 'soilvar':
 UNITS = UNITS
else:
 UNITS = IUNITS

#-------------------
# Print values
#-------------------

if Global == True or Region == True or Point == True:
 print("**************************")
 print("SPATIALLY AVERAGED VALUES")
 print("**************************")

if Global == True:
 for i in range(len(cases)):
  print_spatial_average.print_global_average(case=cases[i],variable=var_avg_by_case[i,:,:],units=UNITS)

if Region == True:
 for i in range(len(cases)):
  print_spatial_average.print_region_average(case=cases[i],variable=var_avg_by_case[i,:,:],units=UNITS,
                                             slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon)

if Point == True:
 for i in range(len(cases)):
  print_spatial_average.print_point_average(case=cases[i],variable=var_avg_by_case[i,:,:],units=UNITS,
                                            latpt=latpoint,lonpt=lonpoint)

#########################################################################################################
#
# Plot the values 
#
#########################################################################################################

#--------------------------------------
# Determine var_name for output file
#--------------------------------------

if VAR2 == None:
 VAR2 = ''
 MATH = ''

if varfunc == 'varlev':
 var_name = str(VAR1)+str(MATH)+str(VAR2)+'$_{'+str(LEVEL)+'hPa}$'
elif varfunc == 'prect':  
 var_name = 'PRECT'
elif varfunc == 'rainiso':
 if ISO == 'd18O':
  var_name = '$\delta^{18}$O$_{P}$'
 elif ISO == 'dHDO':
  var_name = '$\delta$D$_{P}$'
elif varfunc == 'soiliso':
 if ISO == 'd18O':
  var_name = '$\delta^{18}$O$_{S}$'
 elif ISO == 'dHDO':
  var_name = '$\delta$D$_{S}$'
elif varfunc == 'vaporiso':
 if ISO == 'd18O':
  var_name = '$\delta^{18}$O$_{VAPOR_{'+str(ILEVEL)+'hPa}}$'
 elif ISO == 'dHDO':
  var_name = '$\delta$D$_{VAPOR_{'+str(ILEVEL)+'hPa}}$'
elif varfunc == 'isoroot':
 if ISO == 'd18O':
  var_name = '$\delta^{18}$O$_{S_{rootwgt}}$'
 elif ISO == 'dHDO':
  var_name = '$\delta$D$_{S_{rootwgt}}$'
else:
 var_name = str(VAR1)+str(MATH)+str(VAR2) 

#-------------------------------------
# Determine vec_name for output file 
#-------------------------------------

# Processing vec_name and other vector parameters
if overlay_vector == True:
 if overlay_type == 'wind':
  vec_name = str(WIND_LEVEL)+'hPa'+str(overlay_type) 
 elif overlay_type == 'IVT':
  vec_name = 'IVT'+str(ptop)+'-'+str(pbot)+'hPa'
else:
 u_avg_by_case = None
 v_avg_by_case = None
 vec_name = ''

#-----------------------------
# Determine contour levels 
#-----------------------------

if ManLevCntr == False:
   c_loval,c_hival,c_spval,c_tkstd = None,None,None,None
if ManLevDiff == False:
   d_loval,d_hival,d_spval,d_tkstd = None,None,None,None

#------------------
# Make the maps
#------------------

plotmap.plot_contour_map_avg(
                             # Required variables for plot
                              var=var_avg_by_case,cases=cases,seas=season,var_name=var_name,units=UNITS,
                             
                             # Vector variables for plot
                             overlay_vec=overlay_vector,u=u_avg_by_case,v=v_avg_by_case,vec_name=vec_name,vec_units=WIND_UNITS, 
                             
                             # Naming conventions for output file
                             folderpath=folderpath,filesuf=filesuf,extra_name=extra_name,
                             
                             # Mapping specifications                             
                             LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,
                             loval=c_loval,hival=c_hival,spval=c_spval,tkstd=c_tkstd,extnd=c_extnd,
                             proj=proj,cntr_type=Contour_type,colort=ColCntr,Ind_plots=Ind_plots,

                             # Plotting region and point on map
                             regbox=Region,regslat=southlat,regnlat=northlat,regwlon=westlon,regelon=eastlon,
                             point=Point,ptlat=latpoint,ptlon=lonpoint,

                             # Extras for making a specific plot that are not specified above (add as many or as few as needed)
                             figh=14.5,figw=9,lonlblsp=30.,xmin_mj=3,ymin_mj=3,lontksp=13,cbar_pad=0.06
                            
                             )

plotmap.plot_diff_contour_map_avg(
                                  # Required variables for plot
                                  var=var_avg_by_case,cases=cases,seas=season,var_name=var_name,units=UNITS,

                                  # Vector variables for plot
                                  overlay_vec=overlay_vector,u=u_avg_by_case,v=v_avg_by_case,vec_name=vec_name,vec_units=WIND_UNITS,

                                  # Naming conventions for output file
                                  folderpath=folderpath,filesuf=filesuf,extra_name=extra_name,

                                  # Mapping specifications                             
                                  LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,
                                  loval=d_loval,hival=d_hival,spval=d_spval,tkstd=d_tkstd,extnd=d_extnd,
                                  proj=proj,cntr_type=Contour_type,colort=ColDiff,Ind_plots=Ind_plots,

                                  # Plotting region and point on map
                                  regbox=Region,regslat=southlat,regnlat=northlat,regwlon=westlon,regelon=eastlon,
                                  point=Point,ptlat=latpoint,ptlon=lonpoint,

                                  # Extras for making a specific plot that are not specified above (add as many or as few as needed)
                                  figh=7.,figw=8.,ttlfts=10.,lonlblsp=30.,xmin_mj=3,ymin_mj=3,lontksp=13

                                  )
