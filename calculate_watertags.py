####################################################################################################################
#
# Calculate and plot rainfall and d18Op results from water tagged iCESM simulation(s)
#
# *** What this script can do: ***
# 
# 1. Global maps with seasonally averaged values for each region (precipitation, precipitation percentage, d18Op)
# 2. Plot individual global maps for each tag region's seasonally averaged precipitation and d18Op
# 3. Output an excel file with each tag region's precipitation and d18Op values for each month.
#
# This script utilizes functions from ./watertagging/
#
# NOTE: If you are using this script to calculate your own set of water tagged simulations follow these 
#       instructions:
#
#       First, create your drawn land/ocean tags in ./watertagging/tagged_regions.py
#       Second, use this script and modify parameters below in the section called 
#               'Specifications for plots are made here'
#
###################################################################################################################

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cmaps
from matplotlib import colors
from watertagging.print_watertag_values import print_watertag_values,monthly_watertag_values_to_excel
from watertagging.watertag_plots import watertagging_values_on_map,plot_tagged_precip_and_d18Op
from watertagging.seas_avg_LL_watertags import seasavg_watertagging_vars
import climate_analysis.seas_avg_LL as seasavg
from climate_analysis.coordinate_functions import lat_lon_index_array

#########################################################
#
# Specifications for plots are made here
#
#########################################################

#------------------------------------------------
# Which plots to include in output...
#------------------------------------------------

# Global maps with values for precip, precip pct, and d18Op
TEXT_MAPS = True 

# Print values for each tag region to screen
PRINT_VAL = True  

# Individual global maps of...
IND_PRECIP = True 
IND_d18Op  = True  

# Excel sheet with monthly values for each tagged region by month
MAKE_EXCEL = True  

# For the above three outputs, should they be differences between cases?
DIFF = False  

#------------------------------------------------
# Specify data path variables 
#------------------------------------------------

# Component model to specify in case strings
model = 'cam'

# File paths and names for each case
# 20yr water tagging experiments (cam only)
CASES = ['/paleonas/ajthompson/postproc/f.e12.F_1850_CAM5.wiso.f19.0ka.002.watertags.2.'+model+'.h0.0006-0025.climo.nc',
        '/paleonas/ajthompson/postproc/f.e12.F_1850_CAM5.wiso.f19.21ka.modern.d18Osw.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
        '/paleonas/ajthompson/postproc/f.e12.F_1850_CAM5.wiso.f19.21kaGHG.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
        '/paleonas/ajthompson/postproc/f.e12.F_1850_CAM5.wiso.f19.21kaGlac.001.watertags.2.'+model+'.h0.0006-0025.climo.nc',
        '/paleonas/ajthompson/postproc/f.e12.F_1850_CAM5.wiso.f19.21kaSL.002.watertags.'+model+'.h0.0006-0025.climo.nc']
cases = ['0ka',
        '21ka',
        '21ka$_{GHG}$',
        '21ka$_{GLAC}$',
        '21ka$_{SL}$']

# Anything extra to add to output file name?
extra_name = 'last20yrs'

#--------------------------------
# Seasonal averaging variables
#--------------------------------

# Indices to define range of time dimension read in from files; if climatology file (Jan-Dec), specify begi=0 and endi=12
begi = 0
endi = 12

# Season to average over, indices corresponding to individual months, season string will automatically populate later 
MON = [0,1,2,3,4,5,6,7,8,9,10,11]
#MON = [6,7,8,9]

# Reference list for indices
# Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec
#   0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11

#---------------------------------------------------------------------------------------
# Specify monthly weights, if necessary
#---------------------------------------------------------------------------------------

# Array with 12 months of weights 
wgt_by_mon = xr.DataArray(None,dims=['case','month'],coords=dict(case=cases,month=np.arange(1,13,1))).astype(float)

# Array with n (based on MON) months of weights 
# NOTE: Here we are not specifying the 'time' or 'month' coord so weighting with var_time will work
wgt_mon = xr.DataArray(np.zeros([len(cases),len(MON)]),dims=['case','time'],coords=dict(case=cases)).astype(float)

### Enter monthly weights for each case here
#...............................................................................................................
# 0ka
wgt_by_mon[0,:] = np.array([0.08493151,0.076712325,0.08493151,0.08219178,0.08493151,0.08219178,
                            0.08493151,0.08493151 ,0.08219178,0.08493151,0.08219178,0.08493151])
# 21ka
wgt_by_mon[1:,:] = np.array([0.084931507,0.076712329,0.084931507,0.082191781,0.084931507,0.082191781,
                             0.084931507,0.084931507,0.082191781,0.084931507,0.082191781,0.084931507])
#...............................................................................................................

# Modify depending if averaging over entire year or only a select few months
for i in range(len(cases)):
 if len(MON) == 12:
  wgt_mon[i,:] = wgt_by_mon[i,:]
 else:
  wgt_mon[i,:] = wgt_by_mon[i,MON]/np.sum(wgt_by_mon[i,MON])

#----------------------------------------------------------------------------------
# Water tagging variables, plotting boundaries takes place in a separate function
#----------------------------------------------------------------------------------

# Long-form name of each tag (in order)
tagnames = ['Antarctica','North America/Greenland','South America (-Amazon)','Eurasia','Africa (-Congo)','Sundaland NW',
            'SundalandNE','Sundaland SW','Sundaland SE','Sahulland','Australia/Oceania','Amazon','Congo','North Pacific',
            'North Atlantic','North Barents/Arctic Sea','Tropical Pacific NE','Caribbean','Tropical Atlantic NW',
            'Tropical Atlantic NE','Mediterranean','Indian Ocean NW/Arabian Sea','Indian Ocean NE/Bay of Bengal',
            'Sundaland NW ocean','Sundaland NE ocean/South China Sea','Sundaland SW ocean','Sundaland SE ocean',
            'Tropical Pacific NW','Tropical Pacific North Central','Tropical Pacific SE','Tropical Atlantic SW',
            'Tropical Atlantic SE','Tropical Indian SW','Tropical Indian South Central','Tropical Indian SE',
            'Sahul region ocean','Tropical Pacific South Central','South Pacific','South Atlantic','South Indian']

# Code name of each tag (in order)
tagcodes = ['ANTA','NAMG','SAME','ERAS','AFRI','SLNW','SLNE','SLSW','SLSE','SAHL','AUST','AMAZ','CONG','NPAC','NATL',
            'ARCT','TPNE','CARB','TANW','TANE','MEDI','ARAB','BOFB','SONW','SONE','SOSW','SOSE','TPNW','TPNC','TPSE',
            'TASW','TASE','TISW','TISC','TISE','SAHO','TPSC','SPAC','SATL','SIND']

# Number of land vs. ocean tags
num_landtags  = 13
num_oceantags = 27

# Lat/lon values for plotting text in land and ocean tag region plots
landlat  = [-81,  38, -20, 58, 22, 11,  11, -2,  -2,  -7, -26,   0,  5]
landlon  = [  0,-100, -58, 55,  2, 99, 113, 99, 113, 135, 135, -61, 22]
oceanlat = [  36, 36,58,   8, 37, 18,18,44,10,25, 11, 11, -2, -2, 18,   8,  -8,-13,-13,-5,-10,-17,-20,  -8, -40,-40,-40]
oceanlon = [-150,-50,50,-120,-88,-45, 0,22,63,85,100,113,100,113,138,-165,-120,-25,  1,57, 83,107,134,-165,-135,-15, 90]

#------------------------------------------------------------------------------------------
# Define region for which water tagging results will be calculated  
#------------------------------------------------------------------------------------------

# Name the region
reg_name = 'Guatemala'

# Define bounds (for single grid cell, set lats as same value and lons as same value)
# negative values = °S and °W, positive values = °N and °E
southlat = 14.0
northlat = 18.0
westlon  = -92.5
eastlon  = -90.0

#-----------------------------------------------------------------
# Specify individual map plot contour levels for prect and d18Op 
#-----------------------------------------------------------------

# When DIFF == False, modify these values for contour levels
if DIFF == False: 

 # Precipitation contours are set manually to include tick for cutoff value
 p_hival   = 2.                          # high value
 p_loval   = 0.                          # low value
 p_spval   = 0.1                         # spacing
 p_mantick = [0.001,0.2,0.5,1.0,1.5,2.0] # manual colorbar ticks
 p_extnd   = 'max'                       # placement of triangles

 # d18Op contours are set based on the following parameters
 o_hival = 0.           # high value
 o_loval = -5.          # low value
 o_spval = 0.2          # spacing
 o_tkstd = 1.           # tick stride
 o_extnd = 'both'       # placement of triangles

# When DIFF == True, modify these values for contour levels
if DIFF == True:

 # Precipitation contours are set based on the following parameters
 p_hival   = 0.5                         # high value
 p_loval   = -0.5                        # low value
 p_spval   = 0.05                        # spacing
 p_mantick = [-0.5,-0.25,0.,0.25,0.5] # manual colorbar ticks
 p_extnd   = 'both'                      # placement of triangles

 # d18Op contours are set based on the following parameters
 o_hival = 5.           # high value
 o_loval = -5.          # low value
 o_spval = 0.2          # spacing
 o_tkstd = 1.           # tick stride
 o_extnd = 'both'       # placement of triangles

#------------------------------------------------------------
# For map plots, zoom into any world region in particular?
#------------------------------------------------------------

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
 LonMin = -120.0     # negative values = °W
 LonMax =  -60.0     # positive values = °E

#---------------------------------------------------------------
# Specify vectors to overlay on plot, if necessary
#---------------------------------------------------------------

# Overlay a vector? If True, what type and level of the atmosphere? 
overlay_vec  = True  
overlay_type = 'IVT'   # 'wind','IVT'

# Define pressure levels with this array, ex. Pressure array goes from 0 hPa to 1000 hPa by 50 hPa 
plev = np.arange(0,1050,50)

# Variables for overlay_type == 'wind', uses 'plev' from above
WIND_LEVEL = 700   # Integer, in hPa
WIND_UNITS = 'm/s' # Text string

# Variables for overlay_type == 'IVT', uses 'plev' from above
ptop_lev = 50.   # in hPa
pbot_lev = 1018. # in hPa

#-------------------------------------------------
# Other map plot specifications are set here
#-------------------------------------------------

# For map
proj         = ccrs.PlateCarree()      # Map projection
Contour_type = 'RasterFill'            # 'RasterFill' or 'AreaFill'

# When DIFF == False, modify these to set color tables
if DIFF == False:
 colorp       = cmaps.cmp_haxby_r       # Contour plot color table, "_r" at end reverses color table  
 coloro       = colors.LinearSegmentedColormap.from_list('name',cmaps.cmp_haxby(np.arange(50)))
                                       # Remove last color (white) from cmp_haxby

# When DIFF == True, modify these to set color tables
elif DIFF == True:
 colorp = cmaps.BlueYellowRed_r
 coloro = cmaps.BlueYellowRed

# For output file
folderpath   = 'pdfs'                  # folder to output file to
filesuf      = '.pdf'                  # type of output file

###############################################################################
#
# Loop through each case to create variables        
#
###############################################################################

#-------------------------------------------------------
# Make season string automatically from season indices
#-------------------------------------------------------

month_by_letter = ['J','F','M','A','M','J','J','A','S','O','N','D']
if len(MON) == 12:   # annual average
 season = 'ANN'
else:
 season = ''.join([month_by_letter[i] for i in MON])

#-----------------------------------------------------------------
# Create multidimensional variable(s) [# of cases x lat x lon] 
#-----------------------------------------------------------------

# Dataset for dimensions
ds = xr.open_dataset(CASES[0])

# Global variables [case,lat,lon]
prect_global = xr.DataArray(None,dims=['case','lat','lon'],
                                 coords=dict(case=cases,lat=ds.lat,lon=ds.lon)).astype(float)
d18Op_global = xr.DataArray(None,dims=['case','lat','lon'],
                                 coords=dict(case=cases,lat=ds.lat,lon=ds.lon)).astype(float)

# Variables by tag [case,tag,lat,lon]
Pi_by_tag = xr.DataArray(None,dims=['case','tag','lat','lon'],
                              coords=dict(case=cases,tag=tagcodes,lat=ds.lat,lon=ds.lon)).astype(float)
d18Opsink_by_tag = xr.DataArray(None,dims=['case','tag','lat','lon'],
                                     coords=dict(case=cases,tag=tagcodes,lat=ds.lat,lon=ds.lon)).astype(float)
d18Opwt_by_tag = xr.DataArray(None,dims=['case','tag','lat','lon'],
                                     coords=dict(case=cases,tag=tagcodes,lat=ds.lat,lon=ds.lon)).astype(float)
# Average value for each tag [case,tag]
prect = xr.DataArray(None,dims=['case','tag'],coords=dict(case=cases,tag=tagcodes))
d18Op = xr.DataArray(None,dims=['case','tag'],coords=dict(case=cases,tag=tagcodes))

# Global variable's average of the selected region
prect_reg = xr.DataArray(None,dims=['case'],coords=dict(case=cases)) 
d18Op_reg = xr.DataArray(None,dims=['case'],coords=dict(case=cases)) 

# Sum of all tagged region values
prect_sum = xr.DataArray(None,dims=['case'],coords=dict(case=cases))
d18Op_sum = xr.DataArray(None,dims=['case'],coords=dict(case=cases))

# Wind vectors if specified above
if overlay_vec == True:
 U = xr.DataArray(None,dims=['case','lat','lon'],coords=dict(case=cases,lat=ds.lat,lon=ds.lon)).astype(float)
 V = xr.DataArray(None,dims=['case','lat','lon'],coords=dict(case=cases,lat=ds.lat,lon=ds.lon)).astype(float)

#--------------------------
# Loop through each case
#--------------------------

for i in range(len(CASES)):

 print('Working on '+str(cases[i]))

 #-----------------------------------------------------------------
 # Global variables: prect and d18Op
 #-----------------------------------------------------------------

 # PRECT
 prect_global[i,:,:] = seasavg.seasavg_prect_LL(path=CASES[i],begi=begi,endi=endi,months=MON,
                                                                        wgt_mon=wgt_mon[i,:]) 
 # d18Op
 d18Op_global[i,:,:] = seasavg.seasavg_rainiso_LL(iso_type='d18O',path=CASES[i],begi=begi,endi=endi,
                                                        ptiny=1.E-18,months=MON,wgt_mon=wgt_mon[i,:]) 

 #-----------------------------
 # Loop through each tag
 #-----------------------------

 for tag in range(len(tagnames)):

  Pi_by_tag[i,tag,:,:], d18Opsink_by_tag[i,tag,:,:] = seasavg_watertagging_vars(

       # Inputs for functions:
       tagcode=tagcodes[tag],months=MON,path=CASES[i],begi=begi,endi=endi,wgt_mon=wgt_mon[i,:]) 

 #-----------------------------------
 # Perform calculation on variables
 #-----------------------------------

 # Initialize latitude weights
 lat = ds.lat
 lon = ds.lon
 lat_wgts = np.cos(np.deg2rad(lat))

 # Create list arrays of lat and lon values for regional averaging
 latarray, lonarray = lat_lon_index_array(lat=ds.lat,lon=ds.lon,slat=southlat,nlat=northlat,
                                                                 wlon=westlon,elon=eastlon)

 # Loop through each tag to calculate prect and d18Op tag region values 
 for tag in range(len(tagnames)):

  # Use water tagging equation: d18Opwt@gc = d18Op_tag@gc * ( total16Op_tag@gc / prect_global@gc )
  d18Opwt_by_tag[i,tag,:,:] = d18Opsink_by_tag[i,tag,:,:] * ( Pi_by_tag[i,tag,:,:] / prect_global[i,:,:] )

  # Area weight prect and d18Op for tagged region average
  prect_wgt    = Pi_by_tag[i,tag,latarray,lonarray].weighted(lat_wgts[latarray])
  prect[i,tag] = prect_wgt.mean(('lon','lat')) 
  d18Op_wgt    = d18Opwt_by_tag[i,tag,latarray,lonarray].weighted(lat_wgts[latarray])
  d18Op[i,tag] = d18Op_wgt.mean(('lon','lat'))

 # Calculate tagged region average with global variables to check against tagged variables
 prect_global_wgt = prect_global[i,latarray,lonarray].weighted(lat_wgts[latarray])
 prect_reg[i]      = prect_global_wgt.mean(('lon','lat'))
 d18Op_global_wgt = d18Op_global[i,latarray,lonarray].weighted(lat_wgts[latarray])
 d18Op_reg[i]      = d18Op_global_wgt.mean(('lon','lat'))

 # Sum tag values together
 prect_sum[i] = np.sum(prect[i,:])
 d18Op_sum[i] = np.sum(d18Op[i,:])

 # Define wind vector variables if necessary
 if overlay_vec == True:
  if overlay_type == 'wind':
   U[i,:,:], V[i,:,:] = seasavg.seasavg_wind_vec_LL(path=CASES[i],begi=begi,endi=endi,level=WIND_LEVEL,
                                                             plev=plev,months=MON,wgt_mon=wgt_mon[i,:])
  if overlay_type == 'IVT':
   U[i,:,:], V[i,:,:] = seasavg.seasavg_IVT_vec_LL(path=CASES[i],begi=begi,endi=endi,ptop=ptop_lev,
                                                    pbot=pbot_lev,plev=plev,months=MON,wgt_mon=wgt_mon[i,:])
 elif overlay_vec == False:
  U = None
  V = None

########################################################################################################
#
# Make map/text plots for each case (diffs if specified above) 
#
########################################################################################################

if TEXT_MAPS == True:
 print('Plotting text maps...')

 #------------------------------------------------------------------
 # Determine if difference plots or not then loop through each case
 #------------------------------------------------------------------
 
 if DIFF == False:
 
  for i in range(len(CASES)): # loop through all cases
 
   print(cases[i])
 
   # Print tag values for each region to frame
   if PRINT_VAL == True:
    if DIFF == False:
     print_watertag_values(precip=prect[i,:],d18Op=d18Op[i,:],precip_sum=prect_sum[i],d18Op_sum=d18Op_sum[i],
                           precip_reg_gbl=prect_reg[i],d18Op_reg_gbl=d18Op_reg[i],
                           lat=lat,lon=lon,case=cases[i],tagnames=tagnames,season=season,
                           slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon
                           )
 
 
   # Make plots of values for each tagged region
   if TEXT_MAPS == True:
    if DIFF == False:
     watertagging_values_on_map(precip=prect[i,:],d18Op=d18Op[i,:],case=cases[i],tagnames=tagnames,
                                num_landtags=num_landtags,num_oceantags=num_oceantags,path=CASES[i],
                                season=season,lat=lat,lon=lon,landlat=landlat,landlon=landlon,
                                oceanlat=oceanlat,oceanlon=oceanlon,slat=southlat,nlat=northlat,
                                wlon=westlon,elon=eastlon,folderpath=folderpath,filesuf=filesuf,
                                reg_name=reg_name,extra_name=extra_name,proj=proj,
                                LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax
                                )
 
 elif DIFF == True:
 
  for i in range(1,len(CASES)): # Loop through cases 1:end, take difference with first case
 
   print(cases[i]+'-'+cases[0])
 
   # Print tag values for each region to frame
   if PRINT_VAL == True:
    print_watertag_values(precip=prect[i,:]-prect[0,:],d18Op=d18Op[i,:]-d18Op[0,:],
                          precip_sum=prect_sum[i]-prect_sum[0],d18Op_sum=d18Op_sum[i]-d18Op_sum[0],
                          precip_reg_gbl=prect_reg[i]-prect_reg[0],d18Op_reg_gbl=d18Op_reg[i]-d18Op_reg[0],
                          lat=lat,lon=lon,case=str(cases[i]+'-'+cases[0]),tagnames=tagnames,season=season,
                          slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon
                          )
 
 
   # Make plots of values for each tagged region
   if TEXT_MAPS == True:
    watertagging_values_on_map(precip=prect[i,:],d18Op=d18Op[i,:],cntlp=prect[0,:],cntlo=d18Op[0,:],
                               diff=True,case=str(cases[i]+'-'+cases[0]),tagnames=tagnames,
                               num_landtags=num_landtags,num_oceantags=num_oceantags,path=CASES[i],
                               season=season,lat=lat,lon=lon,landlat=landlat,landlon=landlon,
                               oceanlat=oceanlat,oceanlon=oceanlon,slat=southlat,nlat=northlat,
                               wlon=westlon,elon=eastlon,folderpath=folderpath,filesuf=filesuf,
                               reg_name=reg_name,extra_name=extra_name,proj=proj,
                               LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax
                               )


#################################################################################################
#
# Make map plots for precipitation and d18Op for each tagged region (diffs if specified above) 
#
#################################################################################################

if IND_PRECIP == True or IND_d18Op == True:

 print('Plotting maps for each tagged region...')

 #------------------------------------------
 # Determine vec parameters for output file
 #------------------------------------------
 
 # Processing vec_name and other vector parameters
 if overlay_vec == True:
  if overlay_type == 'wind':
   vec_name  = str(overlay_type)+str(WIND_LEVEL)+'hPa'
   vec_units = 'm/s'
   vec_ref   = 10.
   vec_scale = 200.
  elif overlay_type == 'IVT':
   vec_name  = str(overlay_type)+str(int(ptop_lev))+'-'+str(int(pbot_lev))+'hPa'
   vec_units = 'kg/(m*s)'
   vec_ref   = 250.
   vec_scale = 3000.
 else:
  u_avg_by_case = None
  v_avg_by_case = None
  vec_name = ''

 if DIFF == False:

  for i in range(len(CASES)):

   print(cases[i])
 
   plot_tagged_precip_and_d18Op(P=IND_PRECIP,O=IND_d18Op,season=season,
                                prect=Pi_by_tag[i,:,:,:],d18Op=d18Opwt_by_tag[i,:,:,:],lat=lat,lon=lon,
                                num_landtags=num_landtags,num_oceantags=num_oceantags,
                                tagnames=tagnames,case=cases[i],colorp=colorp,coloro=coloro,proj=proj,
                                cntr_type=Contour_type,p_hival=p_hival,p_loval=p_loval,p_spval=p_spval,
                                p_mantick=p_mantick,p_extnd=p_extnd,o_hival=o_hival,o_loval=o_loval,
                                o_spval=o_spval,o_tkstd=o_tkstd,o_extnd=o_extnd,
                                slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon,
                                overlay_vec=overlay_vec,u=U[i,:,:],v=V[i,:,:],vec_units=vec_units,
                                vec_ref=vec_ref,vec_scale=vec_scale,vec_name=vec_name,
                                folderpath=folderpath,filesuf=filesuf,
                                reg_name=reg_name,extra_name=extra_name,
                                LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax
                                ) 

 elif DIFF == True:

  for i in range(1,len(CASES)): # Loop through cases 1:end, take difference with first case

   print(cases[i]+'-'+cases[0])

   plot_tagged_precip_and_d18Op(P=IND_PRECIP,O=IND_d18Op,season=season,cutoff=0.,
                                prect=Pi_by_tag[i,:,:,:]-Pi_by_tag[0,:,:,:],
                                d18Op=d18Opwt_by_tag[i,:,:,:]-d18Opwt_by_tag[0,:,:,:],
                                lat=lat,lon=lon,num_landtags=num_landtags,num_oceantags=num_oceantags,
                                tagnames=tagnames,case=str(cases[i]+'-'+cases[0]),
                                colorp=colorp,coloro=coloro,cntr_type=Contour_type,proj=proj,
                                p_hival=p_hival,p_loval=p_loval,p_spval=p_spval,
                                p_mantick=p_mantick,p_extnd=p_extnd,o_hival=o_hival,o_loval=o_loval,
                                o_spval=o_spval,o_tkstd=o_tkstd,o_extnd=o_extnd,
                                slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon,
                                overlay_vec=overlay_vec,u=U[i,:,:]-U[0,:,:],v=V[i,:,:]-V[0,:,:],
                                vec_units=vec_units,vec_ref=vec_ref,vec_scale=vec_scale,vec_name=vec_name,
                                folderpath=folderpath,filesuf=filesuf,
                                reg_name=reg_name,extra_name=extra_name,
                                LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax
                                )

#################################################################################################
#
# Make Excel file of monthly tagged values for each case 
#
#################################################################################################

if MAKE_EXCEL == True:

 monthly_watertag_values_to_excel(CASES=CASES,cases=cases,begi=begi,endi=endi,
                                  tagnames=tagnames,tagcodes=tagcodes,folderpath=folderpath,
                                  slat=southlat,nlat=northlat,wlon=westlon,elon=eastlon,
                                  reg_name=reg_name)






