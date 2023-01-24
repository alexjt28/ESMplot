#########################################################################################################
#
# Plot seasonal cycle of variable(s) from a number of climate model netCDF output files
#
# Final plots consist of a line plot and panel plots of each monthly average                        
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
import print_values.print_spatial_average as print_spatial_average
import plotting.plot_seascycle_functions as plotseas

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
extra_name = 'last20yrs'

#--------------------------------
# Range of dates 
#--------------------------------

# Indices to define range of time dimension read in from files; if climatology file (Jan-Dec), specify begi=0 and endi=12
begi = 0
endi = 12

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
varfunc = 'prect'      # Options: 'var','varlev','prect','PminE','soilvar','rainiso','soiliso','vaporiso','isoroot'
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
ptop_lev = 50.   # in hPa
pbot_lev = 1018. # in hPa

# Define file paths to read in for overlaying vector
if overlay_vector == True:
 if model != 'cam':
  CASES_vec = []
  for w in CASES:
   new_string = w.replace(model,"cam")
   CASES_vec.append(new_string)
 else:
  CASES_vec = CASES

#--------------------------------------------------------------
# Do you want to print out the average over a specified region?
#--------------------------------------------------------------

xybound = 'Region'   # Options: 'Global','Region','Point'

# Latitude:  +=north, -=south
# Longitude: 0-180=°E, -180-0=°W

# Do you want a global average?
Global = True

# Do you want an average over a specific region?
Region = True

southlat =  14.0
northlat =  18.0
westlon  = -92.5
eastlon  = -90.0

# Do you want an average at a specific grid point?
Point = True

latpoint = 33.0
lonpoint = 255.0

#-----------------------------------------------------------
# Specify line plot parameters 
#-----------------------------------------------------------

# Lists must be same size as CASES/cases
linestyle = ['-']*len(cases)
linewidth = ['3']*len(cases)
linecolor = ['k','tab:red','tab:orange','tab:blue','tab:cyan']
ManY      = False                 # Manual y-axis values
ybot      = 0.0                   # bottom of y-axis
ytop      = 8.0                   # top of y-axis
leg_title = 'Simulation'          # title of legend

#-----------------------------------------------------------
# Specify contour levels for contour and difference values
#-----------------------------------------------------------

# Manual Levels for contour plotting (True = set manual levels)
ManLevCntr = True  
c_hival = 5.0    # high value
c_loval = -5.0   # low value
c_spval = 0.5    # spacing
c_tkstd = 2.0    # tick/label stride
c_extnd = 'both' # where to put triangles; 'both,'neither','max','min'

# Manual Levels for difference plotting (True = set manual levels)
ManLevDiff = True  
d_hival =  2.0   # high value
d_loval = -2.0   # low value
d_spval =  0.2   # spacing
d_tkstd =  0.5   # tick/label stride
d_extnd = 'both' # where to put triangles; 'both,'neither','max','min'

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
 LatMin =    0.0     # negative values = °S
 LatMax =   40.0     # positive values = °N
 LonMin = -140.0     # negative values = °W
 LonMax =  -60.0     # positive values = °E

#-------------------------------------------------
# Other plotting specifications are set here
#-------------------------------------------------

lineplot     = True                    # Include line plot
mapplot      = True                    # Include map plot
mapdiff      = False                   # Are map plots difference plots? 
makegif      = True                    # GIF file of each case's seasonal cycle
proj         = ccrs.PlateCarree()      # Map projection
Contour_type = 'AreaFill'              # 'RasterFill' or 'AreaFill'
ColCntr      = cmaps.cmp_haxby_r       # Contour plot color table 
ColDiff      = cmaps.BlueYellowRed_r   # "_r" at end reverses color table 
folderpath   = 'pdfs'                  # folder to output file to
filesuf      = '.pdf'                  # type of output file

#########################################################################################################
#
# Read in and calculate seasonal cycle map variables here 
#
#########################################################################################################

#----------------------------------------------------------------------------
# Create multidimensional variable(s) [# of cases x 12 months x lat x lon] 
#----------------------------------------------------------------------------

datasize = xr.open_dataset(CASES[0])
var_avg_by_case = xr.DataArray(None,dims=['case','time','lat','lon'],
                               coords=dict(case=cases,time=np.arange(12),lat=datasize.lat,lon=datasize.lon)).astype(float)
if overlay_vector == True:
 vecsize = xr.open_dataset(CASES_vec[0])
 u_avg_by_case = xr.DataArray(None,dims=['case','time','lat','lon'],coords=dict(case=cases,time=np.arange(12),
                                                                                lat=vecsize.lat,lon=vecsize.lon)).astype(float)
 v_avg_by_case = xr.DataArray(None,dims=['case','time','lat','lon'],coords=dict(case=cases,time=np.arange(12),
                                                                                lat=vecsize.lat,lon=vecsize.lon)).astype(float)

#------------------------------------------------------------
# Loop through each case and calculate seasonal average 
#------------------------------------------------------------

for i in range(len(cases)):

 # First, calculate seasonal cycle
 if varfunc == 'var':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_var_TLL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,
                                                      var2=VAR2,math=MATH,mult=MULT)
 if varfunc == 'varlev':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_var_TLL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,
                                                    math=MATH,level=LEVEL,plev=plev,mult=MULT)
 if varfunc == 'prect':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_prect_TLL(path=CASES[i],begi=begi,endi=endi)

 if varfunc == 'PminE':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_PminE_TLL(path=CASES[i],begi=begi,endi=endi)

 if varfunc == 'soilvar':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_soilvar_TLL(var1=VAR1,path=CASES[i],begi=begi,endi=endi,var2=VAR2,
                                                             math=MATH,soillev=SOIL_LEV,mult=MULT)
 if varfunc == 'rainiso':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_rainiso_TLL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,ptiny=ptiny)

 if varfunc == 'soiliso':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_soiliso_TLL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,
                                                             ptiny=ptiny,soillev=ISOIL_LEV)
 if varfunc == 'vaporiso':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_vaporiso_TLL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,
                                              level=ILEVEL,plev=plev,ptiny=ptiny)
 if varfunc == 'isoroot':
  var_avg_by_case[i,:,:,:] = seascyc.seascyc_isoroot_TLL(iso_type=ISO,path=CASES[i],begi=begi,endi=endi,
                                                               ptiny=ptiny,rootwgt=ROOTWGT)

 if overlay_vector == True:

  if overlay_type == 'wind':
   u_avg_by_case[i,:,:,:], v_avg_by_case[i,:,:,:] = seascyc.seascyc_wind_vec_TLL(path=CASES_vec[i],begi=begi,endi=endi,
                                                                  level=WIND_LEVEL,plev=plev)
  if overlay_type == 'IVT':
   u_avg_by_case[i,:,:,:], v_avg_by_case[i,:,:,:] = seascyc.seascyc_IVT_vec_TLL(path=CASES_vec[i],begi=begi,endi=endi,
                                                                   ptop=ptop_lev,pbot=pbot_lev,plev=plev)

 print(var_avg_by_case[i,:,:,:])

#########################################################################################################
#
# Loop through cases and print spatially averaged values 
#
#########################################################################################################

#-----------------------------
# Determine units
#-----------------------------

if varfunc == 'var' or varfunc == 'varlev' or varfunc == 'prect' or varfunc == 'PminE' or varfunc == 'soilvar':
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
  print_spatial_average.print_global_average(case=cases[i],variable=var_avg_by_case[i,:,:,:],units=UNITS)

if Region == True:
 for i in range(len(cases)):
  print_spatial_average.print_region_average(case=cases[i],variable=var_avg_by_case[i,:,:,:],units=UNITS,
                                             slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon)

if Point == True:
 for i in range(len(cases)):
  print_spatial_average.print_point_average(case=cases[i],variable=var_avg_by_case[i,:,:,:],units=UNITS,
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
elif varfunc == 'PminE':
 var_name = 'PminE'
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

#------------------------------------------
# Determine vec parameters for output file
#------------------------------------------

# Processing vec_name and other vector parameters
if overlay_vector == True:
 if overlay_type == 'wind':
  vec_name  = str(overlay_type)+str(WIND_LEVEL)+'hPa'
  vec_units = 'm/s'
  vec_ref   = 10.
  vec_scale = 100.
 elif overlay_type == 'IVT':
  vec_name  = str(overlay_type)+str(int(ptop_lev))+'-'+str(int(pbot_lev))+'hPa'
  vec_units = 'kg/(m*s)'
  vec_ref   = 250.
  vec_scale = 2500.
else:
 u_avg_by_case = None
 v_avg_by_case = None
 vec_name = ''

#-----------------------------
# Determine line plot values 
#-----------------------------

if ManY == False:
   ybot,ytop = None,None

#-----------------------------
# Determine contour levels 
#-----------------------------

if ManLevCntr == False:
   c_loval,c_hival,c_spval,c_tkstd = None,None,None,None
if ManLevDiff == False:
   d_loval,d_hival,d_spval,d_tkstd = None,None,None,None

if mapdiff == True:
   colort = ColDiff
   loval,hival,spval,tkstd,extnd = d_loval,d_hival,d_spval,d_tkstd,d_extnd
elif mapdiff == False:
   colort = ColCntr
   loval,hival,spval,tkstd,extnd = c_loval,c_hival,c_spval,c_tkstd,c_extnd

#-------------------------------------
# Determine xy bounds for line plot       
#-------------------------------------

if xybound == 'Global':
 southlat = -90.
 northlat = 90.
 westlon  = -180.
 eastlon  = 180.
elif xybound == 'Region':
 southlat = southlat
 northlat = northlat
 westlon  = westlon
 eastlon  = eastlon
elif xybound == 'Point':
 southlat = latpoint 
 northlat = latpoint
 westlon  = lonpoint
 eastlon  = lonpoint

#------------------
# Make the plots
#------------------

plotseas.plot_seasonal_cycle(
                             # Required variables for plot
                             var=var_avg_by_case,cases=cases,var_name=var_name,units=UNITS,
                             slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon,

                             # Naming conventions for output file
                             folderpath=folderpath,filesuf=filesuf,extra_name=extra_name,

                             # Line plot specifications
                             lineplot=lineplot,linestyle=linestyle,linewidth=linewidth,linecolor=linecolor,
                             ybot=ybot,ytop=ytop,leg_title=leg_title,

                             # Extra line plot specs

                             # Vector variables for map plots
                             overlay_vec=overlay_vector,u=u_avg_by_case,v=v_avg_by_case,vec_name=vec_name,vec_units=vec_units,
                             vec_ref=vec_ref,vec_scale=vec_scale,
                           
                             # Mapping specifications                             
                             mapplot=mapplot,mapdiff=mapdiff,
                             LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,
                             loval=loval,hival=hival,spval=spval,tkstd=tkstd,extnd=extnd,
                             proj=proj,cntr_type=Contour_type,colort=colort, 

                             # Plotting region and point on map
                             regbox=Region,regslat=southlat,regnlat=northlat,regwlon=westlon,regelon=eastlon,
                             #point=Point,ptlat=latpoint,ptlon=lonpoint,

                             # GIF specifications
                             makegif=makegif,

                             # Extra map plot specs
                             lonlblsp=30.,xmin_mj=3,ymin_mj=3,lontksp=13 

                             )
