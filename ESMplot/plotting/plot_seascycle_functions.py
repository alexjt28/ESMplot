#######################################################################################################
#
# These functions take an xr.DataArray input and generate seasonal cycle line and map plots for each
# case contained within the variable.
#
#######################################################################################################

import numpy as np
import xarray as xr
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.animation as animation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geocat.viz.util as gv
import cmaps
from ESMplot.plotting.plot_functions import save_multi_image,map_ticks_and_labels,draw_region_box
from ESMplot.climate_analysis.coordinate_functions import lat_lon_index_array

# Default variables
fig_let = ['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)','s)',
           't)','u)','v)','w)','x)','y)','z)']

#######################################################################################################
# Seasonal cycle plot, includes line and map plots 
#######################################################################################################

def plot_seasonal_cycle(
                        ### Required Parameters
                        var: xr.DataArray,
 
                        ### Optional Parameters
                        cases: list = [''], var_name: str = '', units: str = '', figw: float = 10., 
                        figh: float = 10., fdpi: float = 300., lineplot: bool = True, 
                        slat: float = None, nlat: float = None, wlon: float = None, elon: float = None, 
                        linestyle: list = ['-']*20, linewidth: list = ['3']*20, 
                        linecolor: list = ['k']*20, ybot: float = None, ytop: float = None,
                        xytitle: bool = True, xyttlloc: str = 'left', xyttlfts: float = 14.,
                        xyttlpad: float = 10., leg_frame: bool = True, leg_framea: float = 1.0,
                        leg_edgecol: str = 'k', leg_ncol: int = 1, leg_loc: str = 'upper left',
                        leg_fontsz: float = 12., leg_title: str = '', leg_ttlfs: float = 12.,
                        leg_bdrpad: float = 1., xy_lbl_sz: float = 15., xy_ypad: float = 10.,
                        xy_tk_sz: float = 14., gridw: float = 0.1, mapplot: bool = True, 
                        mapdiff: bool = False, proj: cartopy.crs.Projection = ccrs.PlateCarree(),
                        hspace: float = 0.14, wspace: float = 0.15, fig_let: list = fig_let,
                        ttlfts: float = 12., ttlhgt: float = 0.92, spttlloc: str = 'left',
                        spttlfts: float = 10., spttlpad: float = 10., LatMin: float = -90.,
                        LatMax: float = 90., LonMin: float = -180., LonMax: float = 180.,
                        coast: bool = True, coastlw: float = 0.75, lake: bool = False,
                        lakelw: float = 0.5, border: bool = True, borderlw: float = 0.75,
                        usstate: bool = False, usstatelw: float = 0.5, add_axes: bool = True,
                        latlblsp: float = 30., lonlblsp: float = 60., 
                        lattksp: int = 7, lontksp: int = 7, tklblsz: float = 9., 
                        tkmajlg: float = 4., tkminlg: float = 2., xmin_mj: int = 2, ymin_mj: int = 2,
                        xpads: float = 10., ypads: float = 10., 
                        glalpha: float = 0.0, glcolor: str = 'gray', cntr_type: str = 'AreaFill',
                        colort: colors.ListedColormap = cmaps.MPL_viridis,
                        loval: float = None, hival: float = None, spval: float = None,
                        tkstd: float = None, extnd: str = 'both', cbar_orient: str = 'horizontal',
                        cbar_pad: float = 0.06, cbar_shrk: float = 0.6, cbar_sp: str = 'proportional',
                        cbar_lblsz: float = 10., overlay_vec: bool = False,
                        u: xr.DataArray = None, v: xr.DataArray = None, vec_scale: float = 100.,
                        vec_wid: float = 0.002, vec_hal: float = 4., vec_hdl: float = 4.,
                        vec_hdw: float = 4., vec_skip: int = 3, vec_ref: float = 5.,
                        vec_units: str = '', vec_name: str = '', regbox: bool = False,
                        regcol: str = 'r', regline: str = '-', reglw: float = 1., 
                        regslat: float = None, regnlat: float = None,
                        regwlon: float = None, regelon: float = None,
                        point: bool = False, ptlat: float = None, ptlon: float = None,
                        shape_lats: list = [], shape_lons: list = [], shape_type: list = [],
                        shape_size: list = [], shape_col: list = [], 
                        extra_name: str = '', folderpath: str = '', filesuf: str = '.pdf',
                        makegif: bool = False, gif_itvl: int = 500
                        ):

  '''Reads in an xr.DataArray of 2+ cases of a seasonal cycle global variable and produces a line plot 
  and panel plots of each case in subsequent PDF pages.  

  Required Parameters                       
  --------------------------------------------------------------------------------------------------------
  var
   class 'xarray.DataArray', The variable containing the values to be plot as contours on a map. Must 
                             contain the following dimensions: [# of cases x lat x lon]. # of cases needs 
                             to be two (2) or greater for this function. 

  Optional Parameters                       
  --------------------------------------------------------------------------------------------------------
  cases: class 'list', List of strings containing the name of case(s) that will be plotted. These strings 
                       will be displayed in the title of each map plot. Default = ['']
  var_name: class 'string', Name of the variable that will be plotted. This string will appear in the 
                            title of each map plot and in the file output name. Default = ''
  units: class 'string', Specified units of variable. This string will be displayed underneath the 
                         colorbar. Default = ''
  figw, figh, fdpi: class 'float', Width (inches), height (inches), and DPI (dots per inch) of figure. 
                                   Default = 10., 10., 300.
  lineplot: class 'bool', If True, includes line plot of seasonal cycles in output file. Default = True
  slat, nlat, wlon, elon: class 'float', Southern/northern latitude and western/eastern longitude boundary
                                         values of the desired region. Matches to nearest index in 
                                         variable's coordinate array. Required if 'lineplot' = True.
                                         LAT: - = °S, + = °N, LON: - = °W, + = °E. 
                                         Default = None, None, None, None
  linestyle: class 'list', Line style for line plot lines as a list of strings. Default is solid lines.  
                           Default = list of str('-')
  linewidth: class 'list', Line width for line plot lines as a list of strings. Default is line width '3'
                           as ['3','3',...,'3'].
  linecolor: class 'list', Line color for line plot lines as a list of strings. Default is black lines as
                           ['k','k',...,'k'].       
  ybot, ytop: class 'float', Value limits at bottom and top of y-axis of line plot. If not specified, this 
                             gets set automatically based on variable. 
  xytitle: class 'bool', If True, sets title above line plot. Default = True
  xyttlloc: class 'str', Location of line plot title. Default = 'left' 
  xyttlfts: class 'float', Font size of line plot title. Default = 14.
  xyttlpad: class 'float', Padding of line plot title relative to top of line plot. Default = 10. 
  leg_frame: class 'bool', If True, sets frame around legend of line plot. Default = True
  leg_framea class 'float', Alpha (opacity) of frame around legend of line plot. Default = 1.
  leg_edgecol: class 'str', Edge color of frame around legend of line plot. Default = 'k'
  leg_ncol: class 'int', Number of text columns in legend of line plot. Default = 1            
  leg_loc: class 'str', Location of legend of line plot. Default = 'upper left' 
  leg_fontsz: class 'float', Font size for text of legend of line plot. Default = 12.
  leg_title: class 'str', Text to display for the title of the legend of line plot. Default = ''
  leg_ttlfs: class 'float', Font size for legend title of line plot. Default = 12.
  leg_bdrpad: class 'float', Fractional white space inside the legend border of line plot. Default = 1.
  xy_lbl_sz: class 'float', Size of line plot labels. Same value for x- and y-axis. Default = 15.
  xy_ypad: class 'float', Padding of y-axis labels of line plot. Default = 10.
  xy_tk_sz: class 'float', Size of line plot ticks. Same value for x- and y-axis. Default = 14.
  gridw: class 'float', Width of grid lines on line plot. gridw = 0.0 eliminates grid lines. Default = 0.1 
  mapplot: class 'bool', If True, includes map plots of each month in output file. Default = True
  mapdiff: class 'bool', If True, map plots are differences between two cases. Default = False
  proj: class 'cartopy.crs.Projection', Projection of map for each subplot. This parameter must be 
                                        cartopy.crs.PlateCarree() currently as cartopy's gridlines do not 
                                        support any other projection at this time. 
                                        Default = ccrs.PlateCarree()
  hspace, wspace: class 'float', Height and width of padding between subplots as a fraction of the average 
                                 Axes height/width, used in plt.subplots_adjust. If number of cases is an 
                                 odd number, this parameter is disabled and users should use figw and figh.
  fig_let: class 'list', List of letters to denote subplots. Spans a-z, but can be modified by specifying a 
                         new list. Default = ['a)','b)','c)', ... ,'z)']
  ttlfts: class 'float', Font size for image title above all monthly map plots. Default = 12.
  ttlhgt: class 'float', Vertical placement of image title. Default = 0.92 
  spttlloc: class 'str', Location of title for each monthly map plot. Default = 'left'
  spttlfts: class 'float', Font size of title for each monthly map plot. Default = 10.
  spttlpad: class 'float', Padding of title for each monthly map plot. Default = 10.
  LatMin, LatMax, LonMin, LonMax: class 'float', Southern/northern latitude and western/eastern longitude
                                                 that sets extent of map to be plot. Values: 
                                                 LAT: - = °S, + = °N, range = -90 to 90
                                                 LON: - = °W, + = °E, range = -180 to 180
                                                 Default = -90., 90., -180., 180.
  coast, coastlw: class 'bool','string', include coastlines on base map and if True, width of lines
                                         Default = True, 0.75
  lake, lakelw: class 'bool','string', include lakes on base map and if True, width of lines
                                       Default = False, 0.5  
  border, borderlw: class 'bool','string', include country borders on base map and if True, width of 
                                           lines. Default = True, 0.75
  usstate, usstatelw: class 'bool','float', include US state boundaries on base map and if True 
                                            width of lines. Default = False, 0.5
  add_axes: class 'bool', Add axes ticks and labels to each map. If True, the following tick and label
                          parameters are active. Default = True
  latlblsp, lonlblsp: class 'float', Spacing of latitude tick labels in degrees of latitude and 
                                     longitude. Default = 30., 60.
  lattksp, lontksp: class 'int', Tick spacing for major latitude and longitude ticks. Examples:
                                 lattksp = 7 -> -90 to 90 by 30 (7 total major ticks)
                                 lontksp = 7 -> -180 to 180 by 60 (7 total major ticks)
                                 Default = 7, 7
  tklblsz: class 'float', Font size of tick labels. Default = 9. 
  tkmajlg: class 'float', Size of major ticks around the outline of each map. Default = 4.
  tkminlg: class 'float', Size of minor ticks around the outline of each map. Default = 2.
  xmin_mj: class 'int', Number of minor ticks in between each major tick on the x-axis (longitude).
                        Default = 2
  ymin_mj: class 'int', Number of minor ticks in between each major tick on the y-axis (latitude). 
                        Default = 2
  xpads, ypads: class 'float', Padding of x-axis longitude and y-axis latitude tick labels from 
                               map outline. Default = 10., 10.
  glalpha: class 'float', Alpha (line opacity) for gridlines between lat/lon ticks. Default is no 
                          gridlines. Modify 1.0 >= alpha > 0.0 for visible gridlines.
                          Default = 0.0 (no gridlines)
  glcolor: class 'str', Color of gridlines if glalpha > 0.0. Default = 'gray'
  cntr_type: class 'str', Type of contour for plotting. Options: 'AreaFill' - interpolated, 
                          'RasterFill' - raster. Default = 'AreaFill'
  colort: class 'colors.ListedColormap', Color table used for plotting contours. Specify names from cmaps 
                                         (https://github.com/hhuangwx/cmaps) or Matplotlib 
                                         (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
                                         Default = cmaps.MPL_viridis
  loval: class 'float', Minimum value in contour levels. If not specified, will default to lowest value 
                        found in variable. 
  hival: class 'float', Maximum value in contour levels. If not specified, will default to highest value 
                        found in variable. 
  spval: class 'float', Stride (spacing value) in contour levels. If not specified, will default to 
                        (hival-loval)/10. 
  tkstd: class 'float', Stride of tick labels for color bar. If not specified, will default to spval. 
  extnd: :class:'str', Option to add triangles on either side of color bar that signify extended value 
                       range. Options: 'neither','max','min','both'. Default = 'both'
  cbar_orient: class 'str', Orientation of color bar. Options: 'horizontal','vertical'
                            Default = 'horizontal'
  cbar_pad: class 'float', Padding of color bar from bottom of map. Default = 0.06
  cbar_shrk: class 'float', Color bar shrink parameter. Modifies the sizing of the color bar. 
                            Default = 0.6
  cbar_sp: class 'str', Spacing of colors in color bar, If 'extnd' is not 'neither', this parameters
                        must be 'proportional'. Default = 'proportional'
  cbar_lblsz: class 'float', Label size for color bar tick labels. Default = 10.
  overlay_vec: class 'bool', Overlay vectors on map plots. If set to True, the following vector parameters 
                             will be active. See 
                             https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.quiver.html and 
                             https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiverkey.html 
                             for more information on the vector parameters listed below. Default = False
  u, v: class 'xr.DataArray', xarray.DataArray of the u- and v-component of the vectors to plot if 
                              overlay_vec = True. Like 'var', must have dimensions: [# of cases x lat x lon].
                              Default = None, None
  vec_scale: class 'float', Scales the length of the arrow inversely. Number of data units per arrow lenth 
                            unit. A smaller scale parameter makes the arrow longer. Scale_units are 
                            permanently set to 'height' for this function. Default = 100.
  vec_wid: class 'float', Shaft width of vector. Default = 0.002
  vec_hal: class 'float', Head axis length of vector. Default = 4.
  vec_hdl: class 'float', Head length of vector. Default = 4.
  vec_hdw: class 'float', Head width of vector. Default = 4.
  vec_skip: class 'int', Number of grid cells to skip when plotting vectors. Higher numbers correspond to 
                         less congestion between plotted vector arrows. Default = 4
  vec_ref: class 'float', Reference vector to be used as the length of the vector key. Default = 5.        
  vec_units: class 'str', Units of the vector that will be displayed in vector key following 'vec_ref'.
                          Default = ''
  vec_name: class 'str', Manual text to be included in output file name containing information about vector. 
                         For example, vec_name = '850hPawinds' for u and v wind plot at 850 hPa.
                         Default = ''
  regbox: class 'bool', Include box around specified region on each subplot. If True, the following region 
                        parameters will be active. Default = False
  regcol: class 'str', Color of region box. Default = 'r' (red) 
  regline: class 'str', Line style of line marking region boundary. Default = '-'
  reglw: class 'float', Width of line marking region boundary. Default = 1.
  regslat, regnlat, regwlon, regelon: class'float', Southern/northern/western/eastern border of region box. 
                                                    Function finds grid cell midpoint closest to specified 
                                                    value and plots box around edge of grid cell.
                                                    Latitude values: - = °S, + = °N.
                                                    Longitude values: - = °W, + = °E.
                                                    Default = None, None, None, None 
  point: class 'bool', Include point at specified latitude and longtiude value on each subplot. If true, the 
                       following point parameters will be active. Default = False
  ptlat, ptlon: class 'float', Latitude and longitude points of plotted point. Function finds grid cell 
                               midpoint closest to specified value and plots marker at the grid cell 
                               midpoint. Latitude values: - = °S, + = °N. Longitude values: - = °W, + = °E.
                               Default = None, None
  shape_lats, shape_lons: class 'list', List of latitude and longitude values for plotting shapes on each 
                                        subplot. Shapes will be plot at exact latitude and longitude
                                        specified. Latitude values: - = °S, + = °N. 
                                        Longitude values: - = °W, + = °E. List should consist of floats, 
                                        for example: [30., 35., 40.]. Default = [], []
  shape_type: class 'list', Marker type for each shape. Find examples at 
                            https://matplotlib.org/stable/api/markers_api.html. List should consist of 
                            strings, for example: ['o', 's', 'p']. Default = []
  shape_size: :class:'list', Marker size for each shape. List should consist of floats, for example: 
                             [10., 5., 10.]. Default = []
  shape_col: class 'list', Marker color for each shape. List should consist of strings, for example: 
                           ['k', 'k', 'k']. Default = []
  extra_name: class 'str', Optional string to be included at the end of the file output name. This is a 
                           place to include any additional information in the file output name. 
                           Default = ''
  folderpath: class 'str', String of the path to the folder with which to deposite the output file produced 
                           by this function. The file name is generated automatically based on input 
                           parameters. Do not include first and last forward slashes, for example: 'pdfs' 
                           or 'home/pdfs/untitledfolder'. Default = ''
  filesuf: class 'str', Type of output file specified as a string. Default = '.pdf'
  makegif: class 'bool', If True, includes gif of map plot of each month as a separate output file.
                         Default = False
  gif_itvl: class 'int', Delay between frames in GIF in milliseconds. Default = 500
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  Possible output: PDF file (line + map plots) and GIF file 
                   Output file(s) saved to a specified folder/file path 

  '''

  #*********************************
  # Set variables before plotting 
  #*********************************

  # Text for monthly plots
  month_txt = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

  # Make variables to plot cyclic
  varc = gv.xr_add_cyclic_longitudes(var, "lon")
  if overlay_vec == True:
   uc = gv.xr_add_cyclic_longitudes(u, "lon")
   vc = gv.xr_add_cyclic_longitudes(v, "lon")

  # Set names of difference cases
  if mapdiff == True:
   cases_diff = [None]*(len(cases)-1)
   for x in range(1,len(cases)):
    cases_diff[x-1] = str(cases[x])+'-'+str(cases[0])

  # Set variable for differences
  if mapdiff == True:
   diff = xr.DataArray(None,dims=['casediff','time','lat','lon'],
                       coords=dict(casediff=cases_diff,time=np.arange(12),lat=varc.lat,lon=varc.lon)).astype(float)
   for x in range(1,len(cases)):
    diff[x-1,:,:,:] = varc[x,:,:,:] - varc[0,:,:,:]
   if overlay_vec == True:
    udiff = xr.DataArray(None,dims=['casediff','time','lat','lon'],
                         coords=dict(casediff=cases_diff,time=np.arange(12),lat=varc.lat,lon=varc.lon)).astype(float)
    vdiff = xr.DataArray(None,dims=['casediff','time','lat','lon'],
                         coords=dict(casediff=cases_diff,time=np.arange(12),lat=varc.lat,lon=varc.lon)).astype(float)
    for x in range(1,len(cases)):
     udiff[x-1,:,:,:] = uc[x,:,:,:] - uc[0,:,:,:]
     vdiff[x-1,:,:,:] = vc[x,:,:,:] - vc[0,:,:,:]

  # Set contour level values if not specified by user
  if loval == None:
   loval = np.nanmin(var)
  if hival == None:
   hival = np.nanmax(var)
  if spval == None:
   spval = (np.nanmax(var)-np.nanmin(var))/10
  if tkstd == None:
   tkstd = (np.nanmax(var)-np.nanmin(var))/10

  # Setting file output name parameters
  if overlay_vec == False:
   vec_name = ''
  else:
   vec_name = ''.join(('_',vec_name))

  # Setting vector grid cell skipping parameters
  skip1D = (slice(None, None, vec_skip))
  skip2D = (slice(None, None, vec_skip), slice(None, None, vec_skip))

  #*********************************
  # 1. PLOT XY OF MONTHLY VALUES
  #*********************************
 
  if lineplot == True:

   print('Working on line plot...')

   #--------------------------------------------------
   # Make variable of line plot values from variable
   #--------------------------------------------------

   # Use function to calculate index arrays for lat and lon
   latarray, lonarray = lat_lon_index_array(lat=var.lat,lon=var.lon,slat=slat,nlat=nlat,wlon=wlon,elon=elon)

   # Spatially weight variable given region boundaries and print result (with metadata) to command line
   lat_wgts = np.cos(np.deg2rad(var.lat))

   # Initialize variable
   var_xy = xr.DataArray(None,dims=['case','months'],coords=dict(case=cases,months=np.arange(1,13,1)))

   # Load values for each case
   for i in range(len(cases)):
    rwgt = var[i,:,latarray,lonarray].weighted(lat_wgts[latarray])
    var_xy[i,:] = rwgt.mean(('lon','lat'))

   #-----------------------
   # Make values for title
   #-----------------------
 
   # Create spacing to show lat/lon values as edge of grid cells rather than midpoint
   latadj = abs(float((var.lat[1] - var.lat[2]) / 2))
   lonadj = abs(float((var.lon[1] - var.lon[2]) / 2))

   # Discover direction of latitude array
   if var.lat[0] > 0:
    lat_direction = 'N->S'
   elif var.lat[0] < 0:
    lat_direction = 'S->N'

   # Discover treatment of longitude array
   if any(i > 180. for i in var.lon):
    lon_treatment = 'degW = 180-360'
   elif any(i < 0. for i in var.lon):
    lon_treatment = 'degW = -180-0'

   # If longitude array treats degW as 180-360, modify lonw,lone
   if lon_treatment == 'degW = 180-360':
    if wlon < 0.:
     wlon = wlon + 360
    if elon < 0.:
     elon = elon + 360

   # Select coordinate indices
   lats = list(var.lat.values).index(var.lat.sel(lat=slat, method='nearest').lat)
   latn = list(var.lat.values).index(var.lat.sel(lat=nlat, method='nearest').lat)
   lonw = list(var.lon.values).index(var.lon.sel(lon=wlon, method='nearest').lon)
   lone = list(var.lon.values).index(var.lon.sel(lon=elon, method='nearest').lon)

   # Define values for each coordinate to print
   lats_val = np.array(var.lat[lats])
   latn_val = np.array(var.lat[latn])
   lonw_val = np.array(var.lon[lonw])
   lone_val = np.array(var.lon[lone])

   # Switch lons back to - = °W before printing, if necessary
   if lon_treatment == 'degW = 180-360':
    if lonw_val > 180.:
     lonw_val = lonw_val - 360
    if lone_val > 180.:
     lone_val = lone_val - 360

   #----------------
   # Define figure
   #----------------
   
   figxy = plt.figure(figsize=(figw,figh),dpi=fdpi)
   axy = figxy.add_subplot(111)
   
   #---------------------
   # Make xy plot
   #---------------------
 
   # Plot data
   for i in range(len(cases)):
    axy.plot(np.arange(1,13,1),var_xy[i,:],linestyle=linestyle[i],linewidth=linewidth[i],color=linecolor[i])

   # Set title
   #if xytitle == True: 
   # axy.set_title(str(var_name)+' ('+str((lats_val-latadj).round(2))+' to '+str((latn_val+latadj).round(2))+'°N, '+str(
   #                                      (lonw_val-lonadj).round(2))+' to '+str((lone_val+lonadj).round(2))+'°E)',
   #                                      loc=xyttlloc,fontsize=xyttlfts,pad=xyttlpad)  
 
   # Set legend
   axy.legend(cases,frameon=leg_frame,framealpha=leg_framea,edgecolor=leg_edgecol,ncol=leg_ncol,loc=leg_loc,
                    fontsize=leg_fontsz,title=leg_title,title_fontsize=leg_ttlfs,borderpad=leg_bdrpad)
   
   # Set X ticks and labels 
   axy.set_xlim([0.9,12.1])
   axy.set_xticks(np.arange(1,13,1))
   axy.set_xticklabels(month_txt,fontsize=xy_tk_sz)
   
   # Set Y ticks and labels
   axy.tick_params(axis='y',which='major',labelsize=xy_tk_sz)
   axy.set_ylabel(units,labelpad=xy_ypad,fontsize=xy_lbl_sz)
   if ybot != None and ytop != None:
    axy.set_ylim([ybot,ytop])
   
   # Grid  
   axy.grid(linewidth=gridw)

  #*********************************
  # 2. PLOT MAPS FOR EACH MONTH 
  #*********************************

  if mapplot == True:

   print('Working on map plots...')

   #------------------------------------
   # Determine how many plots to make
   #------------------------------------

   if mapdiff == False: 
    numpages = len(cases)
    varplot  = varc
    names    = cases 
    if overlay_vec == True:
     uplot = u
     vplot = v
   elif mapdiff == True:
    numpages = len(cases)-1
    varplot  = diff
    names    = cases_diff
    if overlay_vec == True:
     uplot = udiff
     vplot = vdiff

   #--------------------------------------------------------
   # Loop through each PDF page and make monthly map plots 
   #--------------------------------------------------------

   for i in range(numpages):  

    print(names[i]+'...')

    #------------------
    # Define figures
    #------------------

    fig, axes = plt.subplots(nrows=4,ncols=3,figsize=(figw,figh),dpi=fdpi,subplot_kw={'projection':proj})
    fig.subplots_adjust(hspace=hspace,wspace=wspace)

    #--------------
    # Set titles 
    #--------------

    # Main title
    fig.suptitle(str(var_name)+str(vec_name)+' '+names[i],fontsize=ttlfts,y=ttlhgt)

    # Subplot titles
    for r in range(4):
     for c in range(3):
      axes[r,c].set_title(month_txt[3*r+c],loc=spttlloc,fontsize=spttlfts,pad=spttlpad)

    #-------------------
    # Set map extents  
    #-------------------
   
    for r in range(4):
     for c in range(3):
      axes[r,c].set_extent([LonMin,LonMax,LatMin,LatMax],proj)

    #---------------------
    # Define figure map
    #---------------------

    for r in range(4):
     for c in range(3):
      if coast == True:
       axes[r,c].add_feature(cfeature.COASTLINE,linewidths=coastlw)
      if lake == True:
       axes[r,c].add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
      if border == True:
       axes[r,c].add_feature(cfeature.BORDERS,linewidths=borderlw)
      if usstate == True:
       axes[r,c].add_feature(cfeature.STATES,linewidths=usstatelw)

    #----------------------
    # Set axes for maps
    #----------------------

    if add_axes == True:

     # First three rows, first column (lat on left, no lon); [0,0] [1,0] [2,0]
     for r in range(3):
      map_ticks_and_labels(ax=axes[r,0],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=True,right=False)

     # First three rows, middle column (no lat or lon); [0,1] [1,1] [2,1]
     for r in range(3):
      map_ticks_and_labels(ax=axes[r,1],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=False,right=False)

     # First three rows, last column (lat on right, no lon); [0,2] [1,2] [2,2]
     for r in range(3):
      map_ticks_and_labels(ax=axes[r,2],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=False,right=True)

     # Last row, first column (lat on left, lon on bottom); [3,0] 
     map_ticks_and_labels(ax=axes[3,0],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                          xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                          glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                          top=False,bot=True,left=True,right=False)

     # Last row, middle column (no lat, lon on bottom); [3,1] 
     map_ticks_and_labels(ax=axes[3,1],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                          xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                          glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                          top=False,bot=True,left=False,right=False)
   
     # Last row, last column (lat on right, lon on bottom); [3,2] 
     map_ticks_and_labels(ax=axes[3,2],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                          xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                          glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                          top=False,bot=True,left=False,right=True)

    #----------------------
    # Plot variables 
    #----------------------

    for r in range(4):
     for c in range(3):
      if cntr_type == 'AreaFill':
       cntr = varplot[i,3*r+c,:,:].plot.contourf(ax=axes[r,c],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                                 add_colorbar=False,add_labels=False,extend=extnd)
      elif cntr_type == 'RasterFill':
       cntr = varplot[i,3*r+c,:,:].plot(ax=axes[r,c],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                        add_colorbar=False,add_labels=False,extend=extnd)
    
      if overlay_vec == True:
       vec = axes[r,c].quiver(uplot.lon[skip1D],uplot.lat[skip1D],uplot[i,3*r+c,:,:][skip2D],vplot[i,3*r+c,:,:][skip2D],scale_units='height',
                              pivot='mid',scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
       qkey = axes[r,c].quiverkey(vec,0.8,0.2,vec_ref,str(int(vec_ref))+str(vec_units),labelpos='N',coordinates='figure',color='k')

    #-------------------------------
    # Color bar for plots 
    #-------------------------------

    cbar = fig.colorbar(cntr,ax=axes,orientation=cbar_orient,ticks=MultipleLocator(tkstd),shrink=cbar_shrk,pad=cbar_pad,label=units,
                         spacing=cbar_sp)
    cbar.ax.tick_params(labelsize=cbar_lblsz)

    #------------------------------
    # Add region and point to map
    #------------------------------

    if regbox == True or point == True:

     # Adjusters to put box around entire grid cell
     latadj = abs(float((var.lat[1] - var.lat[2]) / 2))
     lonadj = abs(float((var.lon[1] - var.lon[2]) / 2))

    if regbox == True:

     # Handle values if °W = 180-360
     if regwlon < 0.:
      regwlon = regwlon + 360
     if regelon < 0.:
      regelon = regelon + 360

     # Pick best values for region
     slat = float(var.lat.sel(lat=regslat, method='nearest'))
     nlat = float(var.lat.sel(lat=regnlat, method='nearest'))
     wlon = float(var.lon.sel(lon=regwlon, method='nearest'))
     elon = float(var.lon.sel(lon=regelon, method='nearest'))

     # Handle values if °W = 180-360
     if wlon > 180.0:
      wlon = wlon - 360
     if elon > 180.0:
      elon = elon - 360

     for r in range(4):
      for c in range(3):
       draw_region_box(ax=axes[r,c],slat=slat-latadj,nlat=nlat+latadj,wlon=wlon-lonadj,elon=elon+lonadj,
                      linestyle=regline,facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)

    if point == True:

     # Pick best values for point
     latpt = float(var.lat.sel(lat=ptlat, method='nearest'))
     lonpt = float(var.lon.sel(lon=ptlon, method='nearest'))

     # Handle values if °W = 180-360
     if lonpt > 180.0:
      lonpt = lonpt - 360

     for r in range(4):
      for c in range(3):
       axes[r,c].plot(lonpt,latpt,'ro',markersize=5)

    #------------------------
    # Add shapes to map  
    #------------------------
 
    for r in range(4):
     for c in range(3):
      for s in range(len(shape_type)):
       axes[r,c].plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                       markeredgecolor=shape_col)
 
  #***********************
  # Save output to file
  #***********************

  if mapdiff == False:
   difftxt = ''
  elif mapdiff == True:
   difftxt = '_diff'
 
  save_multi_image('./'+str(folderpath)+'/seascycle'+difftxt+'_'+str(var_name)+str(vec_name)+'_'+str(extra_name)+str(filesuf))

  plt.close('all')

  ############################################################################################################
  # GIF animation
  ############################################################################################################

  if makegif == True:

   print('Working on animated gifs...')

   #------------------------------------
   # Determine how many plots to make
   #------------------------------------

   if mapdiff == False:
    numpages = len(cases)
    varplot  = varc
    names    = cases
    if overlay_vec == True:
     uplot = u
     vplot = v
   elif mapdiff == True:
    numpages = len(cases)-1
    varplot  = diff
    names    = cases_diff
    if overlay_vec == True:
     uplot = udiff
     vplot = vdiff

   # Make GIF for each case

   for i in range(numpages): 

    print(names[i]+'...')

    figg = plt.figure(figsize=(figw,figh),dpi=fdpi)
    axg  = figg.add_subplot(111, projection=proj)

    # Defining function for animation to loop through
    def animate_seas(x):
     axg.set_title(str(var_name)+str(vec_name)+' '+names[i]+' '+month_txt[x],loc=spttlloc,fontsize=spttlfts,pad=spttlpad)
     axg.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
     if coast == True:
      axg.add_feature(cfeature.COASTLINE,linewidths=coastlw)
     if lake == True:
      axg.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
     if border == True:
      axg.add_feature(cfeature.BORDERS,linewidths=borderlw)
     if usstate == True:
      axg.add_feature(cfeature.STATES,linewidths=usstatelw)
     if add_axes == True:
      map_ticks_and_labels(ax=axg,LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                           xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                           glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                           top=False,bot=True,left=True,right=True)

     if cntr_type == 'AreaFill':
       cntr = varplot[i,x,:,:].plot.contourf(ax=axg,transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                                 add_colorbar=False,add_labels=False,extend=extnd)
     elif cntr_type == 'RasterFill':
       cntr = varplot[i,x,:,:].plot(ax=axg,transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                        add_colorbar=False,add_labels=False,extend=extnd)
     if overlay_vec == True:
      vec = axg.quiver(uplot.lon[skip1D],uplot.lat[skip1D],uplot[i,x,:,:][skip2D],vplot[i,x,:,:][skip2D],scale_units='height',
                       pivot='mid',scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
      qkey = axg.quiverkey(vec,0.8,0.2,vec_ref,str(int(vec_ref))+str(vec_units),labelpos='N',coordinates='figure',color='k')

    # Create animation with self-defined function above
    anim = animation.FuncAnimation(figg,animate_seas,frames=12,interval=gif_itvl)
    # Color bar
    cbar = figg.colorbar(cntr,ax=axg,orientation=cbar_orient,ticks=MultipleLocator(tkstd),shrink=cbar_shrk,pad=cbar_pad,label=units,
                        spacing=cbar_sp)
    cbar.ax.tick_params(labelsize=cbar_lblsz)
    # Region/point drawn on map
    if regbox == True or point == True:
     # Adjusters to put box around entire grid cell
     latadj = abs(float((var.lat[1] - var.lat[2]) / 2))
     lonadj = abs(float((var.lon[1] - var.lon[2]) / 2))
    if regbox == True:
     # Handle values if °W = 180-360
     if regwlon < 0.:
      regwlon = regwlon + 360
     if regelon < 0.:
      regelon = regelon + 360
     # Pick best values for region
     slat = float(var.lat.sel(lat=regslat, method='nearest'))
     nlat = float(var.lat.sel(lat=regnlat, method='nearest'))
     wlon = float(var.lon.sel(lon=regwlon, method='nearest'))
     elon = float(var.lon.sel(lon=regelon, method='nearest'))
     # Handle values if °W = 180-360
     if wlon > 180.0:
      wlon = wlon - 360
     if elon > 180.0:
      elon = elon - 360
     draw_region_box(ax=axg,slat=slat-latadj,nlat=nlat+latadj,wlon=wlon-lonadj,elon=elon+lonadj,
                      linestyle=regline,facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)
    if point == True:
     # Pick best values for point
     latpt = float(var.lat.sel(lat=ptlat, method='nearest'))
     lonpt = float(var.lon.sel(lon=ptlon, method='nearest'))
     # Handle values if °W = 180-360
     if lonpt > 180.0:
      lonpt = lonpt - 360
     axg.plot(lonpt,latpt,'ro',markersize=5)

    # Save GIF as output file
    if mapdiff == False:
     difftxt = ''
    elif mapdiff == True:
     difftxt = '_diff'
    anim.save('./'+str(folderpath)+'/'+names[i]+'_seascycle'+difftxt+'_'+str(var_name)+str(vec_name)+'_'+str(extra_name)+'.gif',
              writer='pillow')

