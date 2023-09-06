#########################################################################################################
#
# Plot seasonally averaged variable(s) from a number of climate model netCDF output files
#
# This script utilizes functions from map_avg_functions.py and map_plot_functions.py to create plots
#
#########################################################################################################

import sys, os 
sys.path.append(os.path.dirname(os.getcwd())) # one dir back
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter
import cmaps
import geocat.viz.util as gv
from ESMplot.climate_analysis import seas_cycle_TLL as seascyc
from ESMplot.climate_analysis import seas_avg_LL as seasavg
from ESMplot.climate_analysis.mon_wgt_avg import mon_wgt_avg
from ESMplot.print_values import print_spatial_average 
from ESMplot.plotting import plot_map_avg_functions as plotmap

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
# 20yr water tagging experiments (cam only)
CASES = ['f.e12.F_1850_CAM5.wiso.f19.0ka.002.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21ka.modern.d18Osw.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21kaGHG.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21kaGlac.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
         'f.e12.F_1850_CAM5.wiso.f19.21kaSL.002.watertags.'+model+'.h0.0006-0025.climo.nc']
cases = ['0ka',
         '21ka',
         '21ka$_{GHG}$',
         '21ka$_{GLAC}$',
         '21ka$_{SL}$']

# Anything extra to add to output file name?
extra_name = 'last20yrs_wgtmon'

#--------------------------------
# Seasonal averaging variables
#--------------------------------

# Indices to define range of time dimension read in from files; 'beg' and 'end' specify entire file length 
begi = 'beg'  # 'beg' or index like 0
endi = 'end'  # 'end' or index like 12

# Season to average over, indices corresponding to individual months, season string will automatically populate later 
MON = [0,1,2,3,4,5,6,7,8,9,10,11]
#MON = [6,7,8,9]
#MON = [5,6,7]
#MON = [5]
#MON = [11,0,1]

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
# 'PminE': Precipitation-Evaporation, automatically reads in and adjusts P and E variables
# 'soilvar': Land variable that includes soil level dimension (time x levgrnd x lat x lon)
# 'rainiso': Precipitation isotope variable, automatically creates variable from specified parameters
# 'soiliso': Soil isotope variable that includes soil level dimension (time x levgrnd x lat x lon)
# 'vaporiso': Water vapor isotope variable, automatically creates variable at specified atmospheric level
# 'isoroot': Soil isotope variable weighted by rooting depth fraction

# Choose which variable function to specify
varfunc = 'prect'    # Options: 'var','varlev','prect','PminE','soilvar','rainiso','soiliso','vaporiso','isoroot'
#.......................................................................................................................
# If not an 'iso-' variable, specify units
UNITS = 'mm/day'  # String of units (ex. 'mm/day')
MULT  = 86400.0    # For 'var', 'varlev', and 'soilvar', scale variable(s) by scalar value
#.......................................................................................................................
# If varfunc is 'var', 'varlev', or 'soilvar'
VAR1  = 'QFLX'     
VAR2  = None   # options: None or '[name]'
MATH  = None   # options: None or one of 'add'->VAR1+VAR2, 'sub'->VAR1-VAR2, 'mul'->VAR1*VAR2, 'div'->VAR1/VAR2, 'exp'->VAR1**VAR2
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
lev_name = f'{LEVEL}hPa'

#---------------------------------------------------------------
# Specify vectors to overlay on plot, if necessary
#---------------------------------------------------------------

# Overlay a vector? If True, what type and level of the atmosphere?
overlay_vector = True  
overlay_type   = 'wind'    # 'wind','ivt'

# Variables for overlay_type == 'wind', uses 'plev' from above
WIND_LEVEL = 500   # Integer, in hPa
WIND_UNITS = 'm/s' # Text string
kwargs_WIND = dict(vec_name=f'{overlay_type}{WIND_LEVEL}hPa',
                   vec_units='m/s', vec_ref=20, vec_scale=100, vec_skip=6)

# Variables for overlay_type == 'ivt', uses 'plev' from above
ptop_lev = 50.   # in hPa
pbot_lev = 1018. # in hPa
kwargs_IVT = dict(vec_name=f'{overlay_type}{int(ptop_lev)}-{int(pbot_lev)}hPa',
                  vec_units='kg/(m*s)', vec_ref=250., vec_scale=2500.)

# Define file paths to read in for overlaying vector
CASES_vec = [w.replace(model, "cam") if overlay_vector == True and model != 'cam' else w for w in CASES]

# Define which kwargs to pass: wind or ivt
kwargs_vec = kwargs_WIND if overlay_type == 'wind' else kwargs_IVT if overlay_type == 'ivt' else None
if overlay_vector == False: u_avg_by_case, v_avg_by_case = None, None

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
wgt_by_mon[0,:] = np.array([0.08493151,0.076712325,0.08493151,0.08219178,0.08493151,0.08219178,
                            0.08493151,0.08493151, 0.08219178,0.08493151,0.08219178,0.08493151])
# 21ka
wgt_by_mon[1:,:] = np.array([0.084931507,0.076712329,0.084931507,0.082191781,0.084931507,0.082191781, 
                             0.084931507,0.084931507,0.082191781,0.084931507,0.082191781,0.084931507])

#...............................................................................................................................

# Modify depending if averaging over entire year or only a select few months
for i in range(len(cases)):
 wgt_mon[i,:] = wgt_by_mon[i,:] if len(MON) == 12 else wgt_by_mon[i,MON]/np.sum(wgt_by_mon[i,MON])

#--------------------------------------------------------------
# Do you want to print out the average over a specified region?
#--------------------------------------------------------------

# Latitude:  +=north, -=south
# Longitude: 0-180=°E, -180-0=°W

# Do you want a global average?
Global = True 

# Do you want an average over a specific region?
Region = False

southlat = 14.0
northlat = 18.0
westlon  = -92.5
eastlon  = -90.0

# Do you want an average at a specific grid point?
Point = False 

latpoint = 33.0
lonpoint = 255.0

#-----------------------------------------------------------
# Specify contour levels for contour and difference values
#-----------------------------------------------------------

# Manual Levels for contour plotting (True = set manual levels)
ManLevCntr = True  
c_hival = 8.0          # high value
c_loval = 0.0          # low value
c_spval = 0.5          # spacing
c_tkstd = 2.0          # tick/label stride
c_extnd = 'max'        # where to put triangles; 'both,'neither','max','min'

# Manual Levels for difference plotting (True = set manual levels)
ManLevDiff = True  
d_hival =  2.0         # high value
d_loval = -2.0         # low value
d_spval =  0.1         # spacing
d_tkstd =  0.5         # tick/label stride
d_extnd = 'both'       # where to put triangles; 'both,'neither','max','min'

#---------------------------------------------
# Zoom into any world region in particular?
#---------------------------------------------

# True=entire world, False=zoomed in to coordinate values in second block
World = True  

if World == True:
 LatMin = -90
 LatMax = 90
 LonMin = -180
 LonMax = 180
else:
 LatMin =  -60.0     # negative values = °S
 LatMax =   90.0     # positive values = °N
 LonMin = -180.0     # negative values = °W
 LonMax =   60.0     # positive values = °E

#-------------------------------------------------
# Plotting specifications are set here
#-------------------------------------------------

Ind_plots    = True                    # Include individual map plots in output file
proj         = ccrs.PlateCarree()      # Map projection
Contour_type = 'AreaFill'              # 'RasterFill' or 'AreaFill'
ColCntr      = cmaps.cmp_haxby_r       # Contour plot color table 
ColDiff      = cmaps.BlueYellowRed     # "_r" at end reverses color table 
folderpath   = 'pdfs'                  # folder to output file to
filesuf      = '.pdf'                  # type of output file

#-------------------------------------------------
# Other specifications are set here as kwargs
#-------------------------------------------------

contour_kwargs = dict(figh=14.5,figw=9.,lonlblsp=30.,xmin_mj=3,ymin_mj=3,lontksp=13,cbar_pad=0.06)

dffrnce_kwargs = dict(figh=7.,figw=8.,ttlfts=10.,lonlblsp=30.,xmin_mj=3,ymin_mj=3,lontksp=13)

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

 # CONTOURS (non-isotope variables)
 if varfunc == 'var':
  var_avg_by_case[i,:,:] = seasavg.seasavg_var_LL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,math=MATH,mult=MULT,
                                                   months=MON,wgt_mon=wgt_mon[i,:])
 if varfunc == 'varlev':
  var_avg_by_case[i,:,:] = seasavg.seasavg_var_LL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,math=MATH,mult=MULT,
                                                   level=LEVEL,plev=plev,months=MON,wgt_mon=wgt_mon[i,:]) 
 if varfunc == 'prect':
  var_avg_by_case[i,:,:] = seasavg.seasavg_prect_LL(path=CASES[i],begi=begi,endi=endi,months=MON,wgt_mon=wgt_mon[i,:]) 

 if varfunc == 'PminE':
  var_avg_by_case[i,:,:] = seasavg.seasavg_PminE_LL(path=CASES[i],begi=begi,endi=endi,months=MON,wgt_mon=wgt_mon[i,:])

 if varfunc == 'soilvar':
  var_avg_by_case[i,:,:] = seasavg.seasavg_soilvar_LL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,math=MATH,
                                                       soillev=SOIL_LEV,mult=MULT,months=MON,wgt_mon=wgt_mon[i,:])

 # CONTOURS (isotopic variables) 
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

  if overlay_type == 'ivt':
   u_avg_by_case[i,:,:], v_avg_by_case[i,:,:] = seasavg.seasavg_IVT_vec_LL(path=CASES_vec[i],begi=begi,endi=endi,ptop=ptop_lev,
                                                                       pbot=pbot_lev,plev=plev,months=MON,wgt_mon=wgt_mon[i,:])

#########################################################################################################
#
# Loop through cases and print spatially averaged values 
#
#########################################################################################################

#-------------------------------------------------------------
# Determine units based on isotopic or non-isotopic variable
#-------------------------------------------------------------

UNITS = UNITS if varfunc not in ['rainiso','soiliso','vaporiso','isoroot'] else IUNITS

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

#-------------------------------------------------------
# A few things to specify to make the code run smoothly
#-------------------------------------------------------

# Manage non-specified variables 
if VAR2 == None: VAR2, MATH = '', ''
if ManLevCntr == False: c_loval,c_hival,c_spval,c_tkstd = None,None,None,None
if ManLevDiff == False: d_loval,d_hival,d_spval,d_tkstd = None,None,None,None

# Specify name of variable
varnames = dict(var=f'{VAR1}{MATH}{VAR2}',varlev=f'{VAR1}{MATH}{VAR2}$_{{{LEVEL}hPa}}$',
                prect='PRECT',PminE='PminE',soilvar=f'{VAR1}{MATH}{VAR2}')
d18Onames = dict(rainiso='$\delta^{18}$O$_{P}$',soiliso='$\delta^{18}$O$_{S}$',
                 vaporiso=f'$\delta^{{18}}$O$_{{VAPOR_{{{ILEVEL}hPa}}}}$',isoroot=f'$\delta^{{18}}$O$_{{S_{{rootwgt}}}}$')
dHDOnames = dict(rainiso='$\delta$D$_{P}$',soiliso='$\delta$D$_{S}$',
                 vaporiso=f'$\delta$D$_{{VAPOR_{{{ILEVEL}hPa}}}}$',isoroot=f'$\delta$D$_{{S_{{rootwgt}}}}$')                
if varfunc in varnames.keys(): var_name = varnames[varfunc]
if ISO == 'd18O': 
 if varfunc in d18Onames.keys(): var_name = d18Onames[varfunc]
elif ISO == 'dHDO': 
 if varfunc in dHDOnames.keys(): var_name = dHDOnames[varfunc]

#------------------
# Make the maps
#------------------

# Contour map
plotmap.plot_contour_map_avg(
                             # Required variables for plot
                              var=var_avg_by_case,cases=cases,seas=season,var_name=var_name,units=UNITS,
                             
                             # Vector variables for plot
                             overlay_vec=overlay_vector,u=u_avg_by_case,v=v_avg_by_case,**kwargs_vec,
 
                             # Naming conventions for output file
                             folderpath=folderpath,filesuf=filesuf,extra_name=extra_name,
                             
                             # Mapping specifications                             
                             LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,   
                             loval=c_loval,hival=c_hival,spval=c_spval,tkstd=c_tkstd,extnd=c_extnd,
                             proj=proj,cntr_type=Contour_type,colort=ColCntr,Ind_plots=Ind_plots,

                             # Plotting region and point on map
                             regbox=Region,regslat=southlat,regnlat=northlat,regwlon=westlon,regelon=eastlon,
                             point=Point,ptlat=latpoint,ptlon=lonpoint,

                             # All other specifications must be contained within contour_kwargs
                             **contour_kwargs)

# Difference map
# Modify kwargs_vec to scale better with the contour map
if 'vec_ref' in kwargs_vec or 'vec_scale' in kwargs_vec:
 diff_kwargs_vec = {key: value/2 if key == 'vec_ref' else value/1.7 if key == 'vec_scale' else value \
                                                                for key, value in kwargs_WIND.items()}
plotmap.plot_diff_contour_map_avg(
                                  # Required variables for plot
                                  var=var_avg_by_case,cases=cases,seas=season,var_name=var_name,units=UNITS,

                                  # Vector variables for plot
                                  overlay_vec=overlay_vector,u=u_avg_by_case,v=v_avg_by_case,**diff_kwargs_vec,

                                  # Naming conventions for output file
                                  folderpath=folderpath,filesuf=filesuf,extra_name=extra_name,

                                  # Mapping specifications                             
                                  LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,
                                  loval=d_loval,hival=d_hival,spval=d_spval,tkstd=d_tkstd,extnd=d_extnd,
                                  proj=proj,cntr_type=Contour_type,colort=ColDiff,Ind_plots=Ind_plots,

                                  # Plotting region and point on map
                                  regbox=Region,regslat=southlat,regnlat=northlat,regwlon=westlon,regelon=eastlon,
                                  point=Point,ptlat=latpoint,ptlon=lonpoint,

                                  # All other specifications must be contained within dffrnce_kwargs
                                  **dffrnce_kwargs)

