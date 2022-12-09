#######################################################################################################
#
# These functions take an xr.DataArray input and generate seasonal average map plots for each case 
# contained within the variable. 
#
#######################################################################################################

import numpy as np
import xarray as xr
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.patches import Rectangle
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geocat.viz.util as gv
import cmaps
from plotting.plot_functions import save_multi_image,map_ticks_and_labels

# Default variables
fig_let = ['a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)','s)',
           't)','u)','v)','w)','x)','y)','z)']

#######################################################################################################
# Plot contours in panel and individual map plots 
#######################################################################################################

def plot_contour_map_avg(var: xr.DataArray,
                         cases: list = [''],
                         var_name: str = '',
                         seas: str = '',
                         units: str = '',
                         figw: float = 10.,
                         figh: float = 10.,
                         fdpi: float = 300.,
                         hspace: float = None,
                         wspace: float = None,
                         fig_let: list = fig_let,
                         ttlloc: str = 'left',
                         ttlfts: float = 12.,
                         ttlpad: float = 10.,
                         LatMin: float = -90.,
                         LatMax: float = 90.,
                         LonMin: float = -180.,
                         LonMax: float = 180.,
                         proj: ccrs.Projection = ccrs.PlateCarree(), 
                         coast: bool = True,
                         coastlw: float = 0.75,
                         border: bool = True,
                         borderlw: float = 0.75,
                         usstate: bool = False,
                         usstatelw: float = 0.5,
                         add_axes: bool = True,
                         latlblsp: float = 30.,
                         lonlblsp: float = 60.,
                         lattksp: int = 7,
                         lontksp: int = 7,
                         tklblsz: float = 10.,
                         tkmajlg: float = 4.,
                         tkminlg: float = 2.,
                         xmin_mj: int = 2,
                         ymin_mj: int = 2,
                         xpads: float = 15.,
                         ypads: float = 15.,
                         glalpha: float = 0.0,
                         glcolor: str = 'gray',
                         cntr_type: str = 'AreaFill',
                         colort: colors.ListedColormap = cmaps.MPL_viridis, 
                         loval: float = None,
                         hival: float = None,
                         spval: float = None,
                         tkstd: float = None,
                         extnd: str = 'both', 
                         cbar_orient: str = 'horizontal',
                         cbar_pad: float = 0.08,
                         cbar_shrk: float = 0.6,
                         cbar_sp: str = 'proportional', 
                         cbar_lblsz: float = 10.,
                         overlay_vec: bool = False, 
                         u: xr.DataArray = None,
                         v: xr.DataArray = None,
                         vec_scale: float = 100., 
                         vec_wid: float = 0.002,
			 vec_hal: float = 4.,
                         vec_hdl: float = 4.,
                         vec_hdw: float = 4.,
                         vec_skip: int = 2,
                         vec_ref: float = 5.,   
                         vec_units: str = '', 
                         vec_name: str = '',
                         regbox: bool = False,
                         regcol: str = 'r',
                         regline: str = '-',
                         regslat: float = None, 
                         regnlat: float = None,
                         regwlon: float = None,
                         regelon: float = None,
                         point: bool = False,
                         ptlat: float = None,
                         ptlon: float = None,
                         shape_lats: list = [], 
                         shape_lons: list = [], 
                         shape_type: list = [], 
                         shape_size: list = [],
                         shape_col: list = [], 
                         Ind_plots: bool = True,
                         extra_name: str = '',
                         folderpath: str = '',
                         filesuf: str = '.pdf'):

  '''Reads in an xr.DataArray of 2+ cases of a global variable and produces a panel plot and optional
  individual plots of each case in subsequent PDF pages.  

  Active issues with this script that need to be fixed:
  #1. Projection can't be changed from ccrs.PlateCarree() unless changing it in the function. This is
      an issue with geocat. 

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  var: :class:'xarray.DataArray'
        The variable containing the values to be plot as contours on a map. Must contain the following
        dimensions: [# of cases x lat x lon]. # of cases needs to be two (2) or greater for this function. 

  cases: :class:'string' - Optional, default = ''
         Name of case(s) that will be plotted. These strings will be displayed in the title of each map plot.

  var_name: :class:'string' - Optional, default = ''
            Name of the variable that will be plotted. This string will appear in the title of each map plot
            and in the file output name.

  seas: :class:'string' - Optional, default = ''
        Name of season for which variable is averaged over. This string will be displayed in the title of
        each map plot and in the file output name. If not explicitly called, this parameter will be ''.

  units: :class:'string' - Optional, default = ''
         Specified units of variable. This string will be displayed underneath the colorbar. If not 
         explicitly called, this parameter will be ''.

  figw: :class:'float' - Optional, default = 10.
        Width of figure in inches.

  figh: :class:'float' - Optional, default = 10.
        Height of figure in inches.
 
  fdpi: :class:'float' - Optional, default = 300.
        DPI (dots per inch) of figure. 

  hspace: :class:'float' - Optional, default = None
          Height of the padding between subplots as a fraction of the average Axes height, used in
          plt.subplots_adjust. If number of cases is an odd number, this parameter will not work.
          Use figw, figh instead. 
 
  wspace: :class:'float' - Optional, default = None
          Width of the padding between subplots as a fraction of the average Axes width, used in
          plt.subplots_adjust. If number of cases is an odd number, this parameter will not work.
          Use figw, figh instead. 

  fig_let: :class:'list' - Optional, default = ['a)','b)','c)', ... ,'z)']
           List of letters to denote subplots. Spans a-z, but can be modified by specifying a 
           new list.

  ttlloc: :class:'string' - Optional, default = 'left'
          Location of title above each subplot.
        
  ttlfts: :class:'float' - Optional, default = 12.
          Font size for title above each subplot.
 
  ttlpad: :class:'float' - Optional, default = 10.
          Padding for title above each subplot.

  LatMin: :class:'float' - Optional, default = -90.
          Southern latitude that sets extent of map to be plot. Values: - = °S, + = °N, range = -90 to 90

  LatMax: :class:'float' - Optional, default = 90.
          Northern latitude that sets extent of map to be plot. Values: - = °S, + = °N, range = -90 to 90

  LonMin: :class:'float' - Optional, default = -180.
          Western latitude that sets extent of map to be plot. Values: - = °W, + = °N, range = -180 to 180

  LonMax: :class:'float' - Optional, default = 180.
          Eastern latitude that sets extent of map to be plot. Values: - = °W, + = °N, range = -180 to 180
 
  proj: :class:'cartopy.crs.Projection' - Optional, default = ccrs.PlateCarree()
        Projection of map for each subplot. This parameter must be cartopy.crs.PlateCarree() currently as
        cartopy's gridlines do not support any other projection at this time. 

  coast: :class:'bool' - Optional, default = True
         Include coastlines on basemap for each subplot.

  coastlw: :class:'float' - Optional, default = 0.75
           Linewidth of coastline on basemap for each subplot. 

  border: :class:'bool' - Optional, default = True
          Include country borders on basemap for each subplot.

  borderlw: :class:'float' - Optional, default = 0.75
            Linewidth of country borders on basemap for each subplot. 
 
  usstate: :class:'bool' - Optional, default = False
           Include US state boundaries on basemap for each subplot.

  usstatelw: :class:'float' - Optional, default = 0.5
              Linewidth of US state boundaries on basemap for each subplot.

  add_axes: :class:'bool' - Optional, default = True
            Add axes ticks and labels to each subplot. If True, the following tick and label parameters
            are active.

  latlblsp: :class:'float' - Optional, default = 30.
            Spacing of latitude tick labels in degrees of latitude.

  lonlblsp: :class:'float' - Optional, default = 60.
            Spacing of longitude tick labels in degrees of longitude.

  lattksp: :class:'int' - Optional, default = 7
           Tick spacing for major latitude ticks. For example, 7 is -90 to 90 by 30 (7 total major ticks).

  lontksp: :class:'int' - Optional, default = 7
           Tick spacing for major longitude ticks. For example, 7 is -180 to 180 by 60 (7 total major ticks).

  tklblsz: :class:'float' - Optional, default = 10.
           Font size of tick labels.

  tkmajlg: :class:'float' - Optional, default = 4.
           Size of major ticks around the outline of each subplot.

  tkminlg: :class:'float' - Optional, default = 2.
           Size of minor ticks around the outline of each subplot.

  xmin_mj: :class:'int' - Optional, default = 2
           Number of minor ticks in between each major tick on the x-axis (longitude).

  ymin_mj: :class:'int' - Optional, default = 2
           Number of minor ticks in between each major tick on the y-axis (latitude). 

  xpads: :class:'float' - Optional, default = 15.
         Padding of x-axis longitude tick labels from subplot outline.

  ypads: :class:'float' - Optional, default = 15.
         Padding of y-axis latitude tick labels from subplot outline.                

  glalpha: :class:'float - Optional, default = 0.0 (no gridlines)
           Alpha (line opacity) for gridlines between lat/lon ticks. Default is no gridlines. Modify 
           1.0 >= alpha > 0.0 for visible gridlines.

  glcolor: :class:'str' - Optional, default = 'gray'
           Color of gridlines if glalpha > 0.0.

  cntr_type: :class:'str' - Optional, default = 'AreaFill'
             Type of contour for plotting. Options: 'AreaFill' - interpolated, 'RasterFill' - raster. 

  colort: :class:'colors.ListedColormap' - Optional, default = cmaps.MPL_viridis
          Color table used for plotting contours. Specify names from cmaps (https://github.com/hhuangwx/cmaps)
          or Matplotlib (https://matplotlib.org/stable/tutorials/colors/colormaps.html)

  loval: :class:'float' - Optional, default = None
         Minimum value in contour levels. If not specified, will default to lowest value found in variable.

  hival: :class:'float' - Optional, default = None
         Maximum value in contour levels. If not specified, will default to highest value found in variable.

  spval: :class:'float' - Optional, default = None
         Stride (spacing value) in contour levels. If not specified, will default to (hival-loval)/10. 

  tkstd: :class:'float' - Optional, default = None
         Stride of tick labels for color bar. If not specified, will default to spval. 

  extnd: :class:'str' - Optional, default = 'both'
         Option to add triangles on either side of color bar that signify extended value range. 
         Options: 'neither','max','min','both'

  cbar_orient: :class:'str' - Optional, default = 'horizontal'
               Orientation of color bar. Options: 'horizontal','vertical'

  cbar_pad: :class:'float' - Optional, default = 0.08
            Padding of color bar from bottom of subplots.

  cbar_shrk: :class:'float' - Optional, default = 0.6
             Color bar shrink parameter. Modifies the sizing of the color bar.

  cbar_sp: :class:'str' = 'proportional'
           Spacing of colors in color bar. If 'extend' is not 'neither', this parameter must be 'proportional'.

  cbar_lblsz: :class:'float' - Optional, default = 10.,
              Label size for color bar tick labels.

  overlay_vec: :class:'bool' - Optional, default = False
               Overlay vectors on map plots. If set to True, the following vector parameters will be active.
               See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.quiver.html and 
               https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiverkey.html for more 
               information on the vector parameters listed below.

  u: :class:'xr.DataArray' - Optional, default = None
     xarray.DataArray of the u-component of the vectors to plot if overlay_vec = True. Like 'var', must have
     dimensions: [# of cases x lat x lon].

  v: :class:'xr.DataArray' - Optional, default = None
     xarray.DataArray of the v-component of the vectors to plot if overlay_vec = True. Like 'var', must have
     dimensions: [# of cases x lat x lon].

  vec_scale: :class:'float' - Optional, default = 100.
             Scales the length of the arrow inversely. Number of data units per arrow lenth unit. A smaller
             scale parameter makes the arrow longer. Scale_units are permanently set to 'height' for this
             function. 

  vec_wid: :class:'float' - Optional, default = 0.002
           Shaft width of vector. 

  vec_hal: :class:'float' - Optional, default = 4.
           Head axis length of vector.

  vec_hdl: :class:'float' - Optional, default = 4.
           Head length of vector.

  vec_hdw: :class:'float' - Optional, default = 4.
           Head width of vector.

  vec_skip: :class:'int' - Optional, default = 2
            Number of grid cells to skip when plotting vectors. Higher numbers correspond to less 
            congestion between plotted vector arrows. 
 
  vec_ref: :class:'float' - Optional, default = 5.
           Reference vector to be used as the length of the vector key.           

  vec_units: :class:'str' - Optional, default = ''
             Units of the vector that will be displayed in vector key following 'vec_ref'.

  vec_name: :class:'str' - Optional, default = ''
            Manual text to be included in output file name containing information about vector. 
            For example, vec_name = '850hPawinds' for u and v wind plot at 850 hPa.

  regbox: :class:'bool' - Optional, default = False
          Include box around specified region on each subplot. If True, the following region parameters
          will be active.

  regcol: :class:'str' - Optional, default = 'r' (red)
          Color of region box. 

  regline: :class:'str' - Optional, default = '-'
          Line style of line marking region boundary.

  regslat: :class:'float' - Optional, default = None
           Southern border of region box. Latitude values: - = °S, + = °N. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  regnlat: :class:'float' - Optional, default = None
           Northern border of region box. Latitude values: - = °S, + = °N. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  regwlon: :class:'float' - Optional, default = None
           Western border of region box. Longitude values: - = °W, + = °E. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  regelon: :class:'float' - Optional, default = None
           Eastern border of region box. Longitude values: - = °W, + = °E. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  point: :class:'bool' - Optional, default = False
         Include point at specified latitude and longtiude value on each subplot. If true, the following
         point parameters will be active.

  ptlat: :class:'float' - Optional, default = None
         Latitude point of plotted point. Function finds grid cell midpoint closest to specified value and plots
         marker at the grid cell midpoint. Latitude values: - = °S, + = °N.

  ptlon: :class:'float' - Optional, default = None
         Longitude point of plotted point. Function finds grid cell midpoint closest to specified value and plots
         marker at the grid cell midpoint. Longitude values: - = °W, + = °E.

  shape_lats: :class:'list' - Optional, default = []
              List of latitude values for plotting shapes on each subplot. Shapes will be plot at exact latitude 
              specified. Latitude values: - = °S, + = °N. List should consist of floats, for example: 
              [30.0, 35.0, 40.0].

  shape_lons: :class:'list' - Optional, default = []
              List of longitude values for plotting shapes on each subplot. Shapes will be plot at exact longitude
              specified. Longitude values: - = °W, + = °E. List should consist of floats, for example:
              [0.0, 5.0, 10.0].

  shape_type: :class:'list' - Optional, default = []
              Marker type for each shape. Find examples at https://matplotlib.org/stable/api/markers_api.html.
              List should consist of strings, for example: ['o', 's', 'p'].

  shape_size: :class:'list' - Optional, default = []
              Marker size for each shape. List should consist of floats, for example: [10., 5., 10.].

  shape_col: :class:'list' - Optional, default = [] 
             Marker color for each shape. List should consist of strings, for example: ['k', 'k', 'k']. 

  Ind_plots: :class:'bool' - Optional, default = True
             Include individual map plots of each case as separate PDF pages appended to output file. 
             Individual plots will conform to all map parameters defined for panel plot.
 
  extra_name: :class:'str' - Optional, default = ''
              Optional string to be included at the end of the file output name. This is a place to include any
              additional information in the file output name.
 
  folderpath: :class:'str' - Optional, default = ''
              String of the path to the folder with which to deposite the output file produced by this function.
              The file name is generated automatically based on input parameters. Do not include first and last
              forward slashes, for example: 'pdfs' or 'home/pdfs/untitledfolder'.

  filesuf: :class:'str' - Optional, default = '.pdf'
           Type of output file specified as a string. 
 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: PDF file 
          Panel (+individuals) plot of the variable input(s) saved to a specified folder/file path  

  '''

  #-------------------------------------------------------------------------
  # Make variables to plot cyclic
  #-------------------------------------------------------------------------
  
  varc = gv.xr_add_cyclic_longitudes(var, "lon")
  if overlay_vec == True:
   uc = gv.xr_add_cyclic_longitudes(u, "lon")
   vc = gv.xr_add_cyclic_longitudes(v, "lon")

  #------------------------------
  # Set base figure parameters
  #------------------------------
  
  # Ceiling of half the number of cases, used for looping through axes
  rowceil = math.ceil(len(cases)/2)
  rowceil_min1 = math.ceil(len(cases)/2-1)

  # Determine if number of cases is even or odd, this determines plot dimensions
  if (len(cases) % 2) == 0: 
   numcases = 'even'
  elif (len(cases) % 2) == 1:
   numcases = 'odd'

  # Set contour level values if not specified by user
  if loval == None:
   loval = np.nanmin(var)
  if hival == None:
   hival = np.nanmax(var)
  if spval == None:
   spval = (np.nanmax(var)-np.nanmin(var))/10
  if tkstd == None:
   tkstd = (np.nanmax(var)-np.nanmin(var))/5

  # Setting file output name parameters
  if overlay_vec == False:
   vec_name = ''
  else: 
   vec_name = ''.join(('_',vec_name))

  # Setting vector grid cell skipping parameters
  skip1D = (slice(None, None, vec_skip)) 
  skip2D = (slice(None, None, vec_skip), slice(None, None, vec_skip))

  #----------------
  # Define figures
  #----------------

  print('Working on panel plot...')
 
  # Define fig and axes
  if numcases == 'even': # symmetrical plot
   fig, axes = plt.subplots(nrows=rowceil,ncols=2,figsize=(figw,figh),dpi=fdpi,subplot_kw={'projection':proj})
   fig.subplots_adjust(hspace=hspace,wspace=wspace)

  if numcases == 'odd': # use gridspec to place last axis in the middle
   fig = plt.figure(figsize=(figw,figh),dpi=fdpi)
   gs = fig.add_gridspec(nrows=100,ncols=100,hspace=hspace,wspace=wspace)  
   axes = []                   
   for r in range(rowceil_min1):
    for c in range(2):
     if c == 0:
      ax = fig.add_subplot(gs[r*int(100/rowceil):r*int(100/rowceil)+int(100/rowceil),0:48],projection=proj)
      axes.append(ax)
     if c == 1:
      ax = fig.add_subplot(gs[r*int(100/rowceil):r*int(100/rowceil)+int(100/rowceil),51:99],projection=proj)
      axes.append(ax)
   lastax = fig.add_subplot(gs[rowceil_min1*int(100/rowceil):-1,25:75],projection=proj)
   axes.append(lastax)
 
  # Reshape axes for plotting
  if numcases == 'even':  
   axes = axes.reshape(rowceil,2)

  #----------------------------------
  # Subplot map components
  #----------------------------------

  if numcases == 'even':
   for r in range(rowceil):
    for c in range(2):
     axes[r,c].set_title(fig_let[2*r+c]+' '+str(seas)+' '+str(var_name)+str(vec_name)+' '+cases[2*r+c],loc=ttlloc,
                                                                                        fontsize=ttlfts,pad=ttlpad) 
     axes[r,c].set_extent([LonMin,LonMax,LatMin,LatMax],proj)                                    
     if coast == True:
      axes[r,c].add_feature(cfeature.COASTLINE,linewidths=coastlw)
     if border == True:
      axes[r,c].add_feature(cfeature.BORDERS,linewidths=borderlw)
     if usstate == True:
      axes[r,c].add_feature(cfeature.STATES,linewidths=usstatelw)
 
  if numcases == 'odd':
   for i in range(len(cases)):
    axes[i].set_title(fig_let[i]+' '+str(seas)+' '+str(var_name)+str(vec_name)+' '+cases[i],loc=ttlloc,fontsize=ttlfts,
                                                                                                            pad=ttlpad) 
    axes[i].set_extent([LonMin,LonMax,LatMin,LatMax],proj)                            
    # Figure map
    if coast == True:
     axes[i].add_feature(cfeature.COASTLINE,linewidths=coastlw)
    if border == True:
     axes[i].add_feature(cfeature.BORDERS,linewidths=borderlw)
    if usstate == True:
     axes[i].add_feature(cfeature.STATES,linewidths=usstatelw) 

  #----------------------------------
  # Add map axes and ticks  
  #----------------------------------

  if add_axes == True:

   ### rows 1 to N-1
   if numcases == 'even':
    for r in range(rowceil_min1):
     for c in range(2):
      if c == 0:
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=True,right=False)
      if c == 1:
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=False,right=True)

   if numcases == 'odd':
    for i in range(len(cases)-1):
     if (i % 2) == 0: # all even indices are in left column
      map_ticks_and_labels(ax=axes[i],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                           xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                           glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                           top=False,bot=False,left=True,right=False)
     if (i % 2) == 1: # all odd indices are in left column
      map_ticks_and_labels(ax=axes[i],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                           xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                           glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                           top=False,bot=False,left=False,right=True)

   ### row N
   if numcases == 'even':
    for r in range(rowceil_min1,rowceil):
     for c in range(2):
      if c == 0:  # Left column
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=True,left=True,right=False)
      if c == 1:  # Right column
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=True,left=False,right=True)

   if numcases == 'odd':
    i = len(cases)-1
    map_ticks_and_labels(ax=axes[i],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                         xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                         glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                         top=False,bot=True,left=True,right=True)    

  #----------------------------------
  # Plot the variables 
  #----------------------------------

  if numcases == 'even':
   for r in range(rowceil):
    for c in range(2):
      if cntr_type == 'AreaFill':
       cntr = varc[2*r+c,:,:].plot.contourf(ax=axes[r,c],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                            add_colorbar=False,add_labels=False,extend=extnd)
      elif cntr_type == 'RasterFill':
       cntr = varc[2*r+c,:,:].plot(ax=axes[r,c],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                   add_colorbar=False,add_labels=False,extend=extnd)
      if overlay_vec == True:
       vec = axes[r,c].quiver(uc.lon[skip1D],uc.lat[skip1D],uc[2*r+c,:,:][skip2D],vc[2*r+c,:,:][skip2D],scale_units='height',pivot='mid',
                               scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
       qkey = axes[r,c].quiverkey(vec,0.8,0.2,vec_ref,str(vec_ref)+str(vec_units),labelpos='N',coordinates='figure',color='k')

  if numcases == 'odd':
   for i in range(len(cases)):
     if cntr_type == 'AreaFill':
      cntr = varc[i,:,:].plot.contourf(ax=axes[i],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                       add_colorbar=False,add_labels=False,extend=extnd)
     elif cntr_type == 'RasterFill':
      cntr = varc[i,:,:].plot(ax=axes[i],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                              add_colorbar=False,add_labels=False,extend=extnd)
     if overlay_vec == True:
      vec = axes[i].quiver(uc.lon[skip1D],uc.lat[skip1D],uc[i,:,:][skip2D],vc[i,:,:][skip2D],scale_units='height',pivot='mid',
                              scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
      qkey = axes[i].quiverkey(vec,0.8,0.2,vec_ref,str(vec_ref)+str(vec_units),labelpos='N',coordinates='figure',color='k')

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
   
   if numcases == 'even': 
    for r in range(rowceil):
     for c in range(2):
      axes[r,c].add_patch(Rectangle((wlon-lonadj,slat-latadj),abs(wlon-lonadj)-abs(elon+lonadj),(nlat+latadj)-(slat-latadj),
                                     linestyle=regline,facecolor='none',edgecolor=regcol,zorder=10))
   if numcases == 'odd':
    for i in range(len(cases)):
      axes[i].add_patch(Rectangle((wlon-lonadj,slat-latadj),abs(wlon-lonadj)-abs(elon+lonadj),(nlat+latadj)-(slat-latadj),
                                   linestyle=regline,facecolor='none',edgecolor=regcol,zorder=10))

  if point == True:
 
   # Pick best values for point
   latpt = float(var.lat.sel(lat=ptlat, method='nearest'))
   lonpt = float(var.lon.sel(lon=ptlon, method='nearest')) 

   # Handle values if °W = 180-360
   if lonpt > 180.0:
    lonpt = lonpt - 360

   if numcases == 'even':
    for r in range(rowceil):
     for c in range(2):
      axes[r,c].plot(lonpt,latpt,'ro',markersize=5)
   if numcases == 'odd':
    for i in range(len(cases)):
     axes[i].plot(lonpt,latpt,'ro',markersize=5)

  #------------------------
  # Add shapes to map  
  #------------------------
 
  if numcases == 'even':
   for r in range(rowceil):
    for c in range(2):
     for s in range(len(shape_type)):
      axes[r,c].plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                      markeredgecolor=shape_col)

  if numcases == 'odd':
   for i in range(len(cases)):
    for s in range(len(shape_type)):
     axes[i].plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                  markeredgecolor=shape_col[s])

  #-----------------------------------------------------------
  # Plot individual maps and append to pdf file if necessary
  #-----------------------------------------------------------

  if Ind_plots == True:

   print('Working on individual plots...')

   for i in range(len(cases)):

    figi = plt.figure(figsize=(figw,figh),dpi=fdpi)
    axi  = figi.add_subplot(111, projection=proj) 
    axi.set_title(fig_let[i]+' '+str(seas)+' '+str(var_name)+str(vec_name)+' '+cases[i],loc=ttlloc,fontsize=ttlfts,
                                                                                                   pad=ttlpad)  # title
    axi.set_extent([LonMin,LonMax,LatMin,LatMax],proj)                                                          # map extent
    if coast == True:
     axi.add_feature(cfeature.COASTLINE,linewidths=coastlw)
    if border == True:
     axi.add_feature(cfeature.BORDERS,linewidths=borderlw)
    if usstate == True:
     axi.add_feature(cfeature.STATES,linewidths=usstatelw)
    map_ticks_and_labels(ax=axi,LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                         xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                         glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                         top=False,bot=True,left=True,right=True)
    if cntr_type == 'AreaFill':
     cntr = varc[i,:,:].plot.contourf(ax=axi,transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                            add_colorbar=True,cbar_kwargs={'orientation':cbar_orient,'ticks':MultipleLocator(tkstd),
                                            'shrink':cbar_shrk,'pad':cbar_pad,'label':units,'spacing':cbar_sp},add_labels=False,extend=extnd)                       
    elif cntr_type == 'RasterFill':
     cntr = varc[i,:,:].plot(ax=axi,transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                    add_colorbar=True,cbar_kwargs={'orientation':cbar_orient,'ticks':MultipleLocator(tkstd),
                                    'shrink':cbar_shrk,'pad':cbar_pad,'label':units},add_labels=False,extend=extnd) 
    if overlay_vec == True:
     vec = axi.quiver(uc.lon[skip1D],uc.lat[skip1D],uc[i,:,:][skip2D],vc[i,:,:][skip2D],scale_units='height',pivot='mid',scale=vec_scale,
                                      width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
     qkey = axi.quiverkey(vec,0.8,0.2,vec_ref,str(vec_ref)+str(vec_units),labelpos='N',coordinates='figure',color='k')
    if regbox == True:
     axi.add_patch(Rectangle((wlon-lonadj,slat-latadj),abs(wlon-lonadj)-abs(elon+lonadj),(nlat+latadj)-(slat-latadj),
                                     linestyle=regline,facecolor='none',edgecolor=regcol,zorder=10))
    if point == True:
     axi.plot(lonpt,latpt,'ro',markersize=5)
    for s in range(len(shape_type)):
     axi.plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                      markeredgecolor=shape_col[s])

  #--------------------------
  # Save output PDF file
  #--------------------------

  save_multi_image('./'+str(folderpath)+'/'+str(seas)+'_'+str(var_name)+str(vec_name)+'_'+str(extra_name)+str(filesuf))
  plt.close('all')

#######################################################################################################
# Plot difference contours in panel and individual map plots 
#######################################################################################################

def plot_diff_contour_map_avg(var: xr.DataArray,
                              cases: list = [''],
                              var_name: str = '',
                              seas: str = '',
                              units: str = '',
                              figw: float = 10.,
                              figh: float = 10.,
                              fdpi: float = 300.,
                              hspace: float = None,
                              wspace: float = None,
                              fig_let: list = fig_let,
                              ttlloc: str = 'left',
                              ttlfts: float = 12.,
                              ttlpad: float = 10.,
                              LatMin: float = -90.,
                              LatMax: float = 90.,
                              LonMin: float = -180.,
                              LonMax: float = 180.,
                              proj: ccrs.Projection = ccrs.PlateCarree(),
                              coast: bool = True,
                              coastlw: float = 0.75,
                              border: bool = True,
                              borderlw: float = 0.75,
                              usstate: bool = False,
                              usstatelw: float = 0.5,
                              add_axes: bool = True,
                              latlblsp: float = 30.,
                              lonlblsp: float = 60.,
                              lattksp: int = 7,
                              lontksp: int = 7,
                              tklblsz: float = 10.,
                              tkmajlg: float = 4.,
                              tkminlg: float = 2.,
                              xmin_mj: int = 2,
                              ymin_mj: int = 2,
                              xpads: float = 15.,
                              ypads: float = 15.,
                              glalpha: float = 0.0,
                              glcolor: str = 'gray',
                              cntr_type: str = 'AreaFill',
                              colort: colors.ListedColormap = cmaps.ncl_default,
                              loval: float = None,
                              hival: float = None,
                              spval: float = None,
                              tkstd: float = None,
                              extnd: str = 'neither',
                              cbar_orient: str = 'horizontal',
                              cbar_pad: float = 0.08,
                              cbar_shrk: float = 0.6,
                              cbar_sp: str = 'proportional',
                              cbar_lblsz: float = 10.,
                              overlay_vec: bool = False,
                              u: xr.DataArray = None,
                              v: xr.DataArray = None,
                              vec_scale: float = 100.,
                              vec_wid: float = 0.002,
                              vec_hal: float = 4.,
                              vec_hdl: float = 4.,
                              vec_hdw: float = 4.,
                              vec_skip: int = 2,
                              vec_ref: float = 5.,
                              vec_units: str = '',
                              vec_name: str = '',
                              regbox: bool = False,
                              regcol: str = 'r',
                              regline: str = '-',
                              regslat: float = None,
                              regnlat: float = None,
                              regwlon: float = None,
                              regelon: float = None,
                              point: bool = False,
                              ptlat: float = None,
                              ptlon: float = None,
                              shape_lats: list = [],
                              shape_lons: list = [],
                              shape_type: list = [],
                              shape_size: list = [],
                              shape_col: list = [],
                              Ind_plots: bool = True,
                              extra_name: str = '',
                              folderpath: str = '',
                              filesuf: str = '.pdf'):

  '''Reads in an xr.DataArray of 3+ cases of a global variable and produces a panel difference plot and 
  optional individual difference plots of each case in subsequent PDF pages.  

  Active issues with this script that need to be fixed:
  #1. North/south limits of lat/lon aren't shown. Ex. Global plot excludes 90N and 90S.
  #2. Projection can't be changed from ccrs.PlateCarree() unless changing it in the function.
  #3. Default color table becomes re-registered if trying to specify the default table when calling
      the function. Ex. Cannot call cmaps.MPL_viridis when calling colort in function.  
  #4. hspace and wspace don't work when there are an odd number of plots in panel plot.

  Parameters                       
  --------------------------------------------------------------------------------------------------------
  var: :class:'xarray.DataArray'
        The variable containing the values to be plot as contours on a map. Must contain the following
        dimensions: [# of cases x lat x lon]. # of cases needs to be two (2) or greater for this function. 

  cases: :class:'string' - Optional, default = ''
         Name of case(s) that will be plotted. These strings will be displayed in the title of each map plot.

  var_name: :class:'string' - Optional, default = ''
            Name of the variable that will be plotted. This string will appear in the title of each map plot
            and in the file output name.

  seas: :class:'string' - Optional, default = ''
        Name of season for which variable is averaged over. This string will be displayed in the title of
        each map plot and in the file output name. If not explicitly called, this parameter will be ''.

  units: :class:'string' - Optional, default = ''
         Specified units of variable. This string will be displayed underneath the colorbar. If not 
         explicitly called, this parameter will be ''.

  figw: :class:'float' - Optional, default = 10.
        Width of figure in inches.

  figh: :class:'float' - Optional, default = 10.
        Height of figure in inches.
 
  fdpi: :class:'float' - Optional, default = 300.
        DPI (dots per inch) of figure. 

  hspace: :class:'float' - Optional, default = None
          Height of the padding between subplots as a fraction of the average Axes height, used in
          plt.subplots_adjust. If number of cases is an odd number, this parameter will not work.
          Use figw, figh instead. 
 
  wspace: :class:'float' - Optional, default = None
          Width of the padding between subplots as a fraction of the average Axes width, used in
          plt.subplots_adjust. If number of cases is an odd number, this parameter will not work.
          Use figw, figh instead. 

  fig_let: :class:'list' - Optional, default = ['a)','b)','c)', ... ,'z)']
           List of letters to denote subplots. Spans a-z, but can be modified by specifying a 
           new list.

  ttlloc: :class:'string' - Optional, default = 'left'
          Location of title above each subplot.
        
  ttlfts: :class:'float' - Optional, default = 12.
          Font size for title above each subplot.
 
  ttlpad: :class:'float' - Optional, default = 10.
          Padding for title above each subplot.

  LatMin: :class:'float' - Optional, default = -90.
          Southern latitude that sets extent of map to be plot. Values: - = °S, + = °N, range = -90 to 90

  LatMax: :class:'float' - Optional, default = 90.
          Northern latitude that sets extent of map to be plot. Values: - = °S, + = °N, range = -90 to 90

  LonMin: :class:'float' - Optional, default = -180.
          Western latitude that sets extent of map to be plot. Values: - = °W, + = °N, range = -180 to 180

  LonMax: :class:'float' - Optional, default = 180.
          Eastern latitude that sets extent of map to be plot. Values: - = °W, + = °N, range = -180 to 180
 
  proj: :class:'cartopy.crs.Projection' - Optional, default = ccrs.PlateCarree()
        Projection of map for each subplot. This parameter must be cartopy.crs.PlateCarree() currently as
        cartopy's gridlines do not support any other projection at this time. 

  coast: :class:'bool' - Optional, default = True
         Include coastlines on basemap for each subplot.

  coastlw: :class:'float' - Optional, default = 0.75
           Linewidth of coastline on basemap for each subplot. 

  border: :class:'bool' - Optional, default = True
          Include country borders on basemap for each subplot.

  borderlw: :class:'float' - Optional, default = 0.75
            Linewidth of country borders on basemap for each subplot. 
 
  usstate: :class:'bool' - Optional, default = False
           Include US state boundaries on basemap for each subplot.

  usstatelw: :class:'float' - Optional, default = 0.5
              Linewidth of US state boundaries on basemap for each subplot.

  add_axes: :class:'bool' - Optional, default = True
            Add axes ticks and labels to each subplot. If True, the following tick and label parameters
            are active.

  latlblsp: :class:'float' - Optional, default = 30.
            Spacing of latitude tick labels in degrees of latitude.

  lonlblsp: :class:'float' - Optional, default = 60.
            Spacing of longitude tick labels in degrees of longitude.

  lattksp: :class:'int' - Optional, default = 7
           Tick spacing for major latitude ticks. For example, 7 is -90 to 90 by 30 (7 total major ticks).

  lontksp: :class:'int' - Optional, default = 7
           Tick spacing for major longitude ticks. For example, 7 is -180 to 180 by 60 (7 total major ticks).

  tklblsz: :class:'float' - Optional, default = 10.
           Font size of tick labels.

  tkmajlg: :class:'float' - Optional, default = 4.
           Size of major ticks around the outline of each subplot.

  tkminlg: :class:'float' - Optional, default = 2.
           Size of minor ticks around the outline of each subplot.

  xmin_mj: :class:'int' - Optional, default = 2
           Number of minor ticks in between each major tick on the x-axis (longitude).

  ymin_mj: :class:'int' - Optional, default = 2
           Number of minor ticks in between each major tick on the y-axis (latitude). 

  xpads: :class:'float' - Optional, default = 15.
         Padding of x-axis longitude tick labels from subplot outline.

  ypads: :class:'float' - Optional, default = 15.
         Padding of y-axis latitude tick labels from subplot outline.                

  glalpha: :class:'float - Optional, default = 0.0 (no gridlines)
           Alpha (line opacity) for gridlines between lat/lon ticks. Default is no gridlines. Modify 
           1.0 >= alpha > 0.0 for visible gridlines.

  glcolor: :class:'str' - Optional, default = 'gray'
           Color of gridlines if glalpha > 0.0.

  cntr_type: :class:'str' - Optional, default = 'AreaFill'
             Type of contour for plotting. Options: 'AreaFill' - interpolated, 'RasterFill' - raster. 

  colort: :class:'colors.ListedColormap' - Optional, default = cmaps.MPL_viridis
          Color table used for plotting contours. Specify names from cmaps (https://github.com/hhuangwx/cmaps)
          or Matplotlib (https://matplotlib.org/stable/tutorials/colors/colormaps.html)

  loval: :class:'float' - Optional, default = None
         Minimum value in contour levels. If not specified, will default to lowest value found in variable.

  hival: :class:'float' - Optional, default = None
         Maximum value in contour levels. If not specified, will default to highest value found in variable.

  spval: :class:'float' - Optional, default = None
         Stride (spacing value) in contour levels. If not specified, will default to (hival-loval)/10. 

  tkstd: :class:'float' - Optional, default = None
         Stride of tick labels for color bar. If not specified, will default to spval. 

  extnd: :class:'str' - Optional, default = 'neither'
         Option to add triangles on either side of color bar that signify extended value range. 
         Options: 'neither','max','min','both'

  cbar_orient: :class:'str' - Optional, default = 'horizontal'
               Orientation of color bar. Options: 'horizontal','vertical'

  cbar_pad: :class:'float' - Optional, default = 0.08
            Padding of color bar from bottom of subplots.

  cbar_shrk: :class:'float' - Optional, default = 0.6
             Color bar shrink parameter. Modifies the sizing of the color bar.

  cbar_sp: :class:'str' = 'proportional'
           Spacing of colors in color bar. If 'extend' is not 'neither', this parameter must be 'proportional'.

  cbar_lblsz: :class:'float' - Optional, default = 10.,
              Label size for color bar tick labels.

  overlay_vec: :class:'bool' - Optional, default = False
               Overlay vectors on map plots. If set to True, the following vector parameters will be active.
               See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.quiver.html and 
               https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiverkey.html for more 
               information on the vector parameters listed below.

  u: :class:'xr.DataArray' - Optional, default = None
     xarray.DataArray of the u-component of the vectors to plot if overlay_vec = True. Like 'var', must have
     dimensions: [# of cases x lat x lon].

  v: :class:'xr.DataArray' - Optional, default = None
     xarray.DataArray of the v-component of the vectors to plot if overlay_vec = True. Like 'var', must have
     dimensions: [# of cases x lat x lon].

  vec_scale: :class:'float' - Optional, default = 100.
             Scales the length of the arrow inversely. Number of data units per arrow lenth unit. A smaller
             scale parameter makes the arrow longer. Scale_units are permanently set to 'height' for this
             function. 

  vec_wid: :class:'float' - Optional, default = 0.002
           Shaft width of vector. 

  vec_hal: :class:'float' - Optional, default = 4.
           Head axis length of vector.

  vec_hdl: :class:'float' - Optional, default = 4.
           Head length of vector.

  vec_hdw: :class:'float' - Optional, default = 4.
           Head width of vector.

  vec_skip: :class:'int' - Optional, default = 2
            Number of grid cells to skip when plotting vectors. Higher numbers correspond to less 
            congestion between plotted vector arrows. 
 
  vec_ref: :class:'float' - Optional, default = 5.
           Reference vector to be used as the length of the vector key.           

  vec_units: :class:'str' - Optional, default = ''
             Units of the vector that will be displayed in vector key following 'vec_ref'.

  vec_name: :class:'str' - Optional, default = ''
            Manual text to be included in output file name containing information about vector. 
            For example, vec_name = '850hPawinds' for u and v wind plot at 850 hPa.

  regbox: :class:'bool' - Optional, default = False
          Include box around specified region on each subplot. If True, the following region parameters
          will be active.

  regcol: :class:'str' - Optional, default = 'r' (red)
          Color of region box. 

  regline: :class:'str' - Optional, default = '-'
          Line style of line marking region boundary.

  regslat: :class:'float' - Optional, default = None
           Southern border of region box. Latitude values: - = °S, + = °N. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  regnlat: :class:'float' - Optional, default = None
           Northern border of region box. Latitude values: - = °S, + = °N. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  regwlon: :class:'float' - Optional, default = None
           Western border of region box. Longitude values: - = °W, + = °E. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  regelon: :class:'float' - Optional, default = None
           Eastern border of region box. Longitude values: - = °W, + = °E. Function finds grid cell midpoint
           closest to specified value and plots box around edge of grid cell.

  point: :class:'bool' - Optional, default = False
         Include point at specified latitude and longtiude value on each subplot. If true, the following
         point parameters will be active.

  ptlat: :class:'float' - Optional, default = None
         Latitude point of plotted point. Function finds grid cell midpoint closest to specified value and plots
         marker at the grid cell midpoint. Latitude values: - = °S, + = °N.

  ptlon: :class:'float' - Optional, default = None
         Longitude point of plotted point. Function finds grid cell midpoint closest to specified value and plots
         marker at the grid cell midpoint. Longitude values: - = °W, + = °E.

  shape_lats: :class:'list' - Optional, default = []
              List of latitude values for plotting shapes on each subplot. Shapes will be plot at exact latitude 
              specified. Latitude values: - = °S, + = °N. List should consist of floats, for example: 
              [30.0, 35.0, 40.0].

  shape_lons: :class:'list' - Optional, default = []
              List of longitude values for plotting shapes on each subplot. Shapes will be plot at exact longitude
              specified. Longitude values: - = °W, + = °E. List should consist of floats, for example:
              [0.0, 5.0, 10.0].

  shape_type: :class:'list' - Optional, default = []
              Marker type for each shape. Find examples at https://matplotlib.org/stable/api/markers_api.html.
              List should consist of strings, for example: ['o', 's', 'p'].

  shape_size: :class:'list' - Optional, default = []
              Marker size for each shape. List should consist of floats, for example: [10., 5., 10.].

  shape_col: :class:'list' - Optional, default = [] 
             Marker color for each shape. List should consist of strings, for example: ['k', 'k', 'k']. 

  Ind_plots: :class:'bool' - Optional, default = True
             Include individual map plots of each case as separate PDF pages appended to output file. 
             Individual plots will conform to all map parameters defined for panel plot.
 
  extra_name: :class:'str' - Optional, default = ''
              Optional string to be included at the end of the file output name. This is a place to include any
              additional information in the file output name.
 
  folderpath: :class:'str' - Optional, default = ''
              String of the path to the folder with which to deposite the output file produced by this function.
              The file name is generated automatically based on input parameters. Do not include first and last
              forward slashes, for example: 'pdfs' or 'home/pdfs/untitledfolder'.

  filesuf: :class:'str' - Optional, default = '.pdf'
           Type of output file specified as a string. 
 
  _______________________________________________________________________________________________________

  Returns
  ---------------------
  output: PDF file 
          Panel (+individuals) plot of the variable input(s) saved to a specified folder/file path  

  '''

  #-------------------------------------------------------------------------
  # Make variables to plot cyclic
  #-------------------------------------------------------------------------

  varc = gv.xr_add_cyclic_longitudes(var, "lon")
  if overlay_vec == True:
   uc = gv.xr_add_cyclic_longitudes(u, "lon")
   vc = gv.xr_add_cyclic_longitudes(v, "lon")

  #-------------------------------------------------------------------------
  # Make difference variables for plotting
  #-------------------------------------------------------------------------

  cases_diff = [None]*(len(cases)-1)
  for x in range(1,len(cases)):
   cases_diff[x-1] = str(cases[x])+'-'+str(cases[0])

  diff = xr.DataArray(None,dims=['casediff','lat','lon'],coords=dict(casediff=cases_diff,lat=varc.lat,lon=varc.lon)).astype(float)
  for x in range(1,len(cases)):
   diff[x-1,:,:] = varc[x,:,:] - varc[0,:,:] 

  if overlay_vec == True:
   udiff = xr.DataArray(None,dims=['casediff','lat','lon'],coords=dict(casediff=cases_diff,lat=uc.lat,lon=uc.lon)).astype(float)
   vdiff = xr.DataArray(None,dims=['casediff','lat','lon'],coords=dict(casediff=cases_diff,lat=vc.lat,lon=vc.lon)).astype(float)
   for x in range(1,len(cases)):
    udiff[x-1,:,:] = uc[x,:,:] - uc[0,:,:]
    vdiff[x-1,:,:] = vc[x,:,:] - vc[0,:,:]

  #------------------------------
  # Set base figure parameters
  #------------------------------

  # Ceiling of half the number of cases, used for looping through axes
  rowceil = math.ceil((len(cases)-1)/2)
  rowceil_min1 = math.ceil((len(cases)-1)/2-1)

  # Determine if number of cases is even or odd, this determines plot dimensions
  if ((len(cases)-1) % 2) == 0:
   numcases = 'even'
  elif ((len(cases)-1) % 2) == 1:
   numcases = 'odd'

  # Set contour level values if not specified by user
  if loval == None:
   loval = np.nanmin(diff)
  if hival == None:
   hival = -1.*loval 
  if spval == None:
   spval = (hival-loval)/10
  if tkstd == None:
   tkstd = (hival-loval)/5

  # Setting file output name parameters
  if overlay_vec == False:
   vec_name = ''
  else:
   vec_name = ''.join(('_',vec_name))

  # Setting vector grid cell skipping parameters
  skip1D = (slice(None, None, vec_skip))
  skip2D = (slice(None, None, vec_skip), slice(None, None, vec_skip))

  #----------------
  # Define figures
  #----------------

  print('Working on panel difference plot...')

  # Define fig and axes
  if numcases == 'even': # symmetrical plot
   fig, axes = plt.subplots(nrows=rowceil,ncols=2,figsize=(figw,figh),dpi=fdpi,subplot_kw={'projection':proj})
   fig.subplots_adjust(hspace=hspace,wspace=wspace)

  if numcases == 'odd': # use gridspec to place last axis in the middle
   fig = plt.figure(figsize=(figw,figh),dpi=fdpi)
   gs = fig.add_gridspec(nrows=100,ncols=100,hspace=hspace,wspace=wspace)
   axes = []
   for r in range(rowceil_min1):
    for c in range(2):
     if c == 0:
      ax = fig.add_subplot(gs[r*int(100/rowceil):r*int(100/rowceil)+int(100/rowceil),0:48],projection=proj)
      axes.append(ax)
     if c == 1:
      ax = fig.add_subplot(gs[r*int(100/rowceil):r*int(100/rowceil)+int(100/rowceil),51:99],projection=proj)
      axes.append(ax)
   lastax = fig.add_subplot(gs[rowceil_min1*int(100/rowceil):-1,25:75],projection=proj)
   axes.append(lastax)

  # Reshape axes for plotting
  if numcases == 'even':
   axes = axes.reshape(rowceil,2)

  #----------------------------------
  # Subplot map components
  #----------------------------------

  if numcases == 'even':
   for r in range(rowceil):
    for c in range(2):
     axes[r,c].set_title(fig_let[2*r+c]+' '+str(seas)+' '+str(var_name)+str(vec_name)+' '+cases_diff[2*r+c],loc=ttlloc,
                                                                                            fontsize=ttlfts,pad=ttlpad)
     axes[r,c].set_extent([LonMin,LonMax,LatMin,LatMax],proj)
     if coast == True:
      axes[r,c].add_feature(cfeature.COASTLINE,linewidths=coastlw)
     if border == True:
      axes[r,c].add_feature(cfeature.BORDERS,linewidths=borderlw)
     if usstate == True:
      axes[r,c].add_feature(cfeature.STATES,linewidths=usstatelw)

  if numcases == 'odd':
   for i in range(len(cases)-1):
    axes[i].set_title(fig_let[i]+' '+str(seas)+' '+str(var_name)+str(vec_name)+' '+cases_diff[i],loc=ttlloc,
                                                                                 fontsize=ttlfts,pad=ttlpad)
    axes[i].set_extent([LonMin,LonMax,LatMin,LatMax],proj)
    # Figure map
    if coast == True:
     axes[i].add_feature(cfeature.COASTLINE,linewidths=coastlw)
    if border == True:
     axes[i].add_feature(cfeature.BORDERS,linewidths=borderlw)
    if usstate == True:
     axes[i].add_feature(cfeature.STATES,linewidths=usstatelw)

  #----------------------------------
  # Add map axes and ticks  
  #----------------------------------

  if add_axes == True:

   ### rows 1 to N-1
   if numcases == 'even':
    for r in range(rowceil_min1):
     for c in range(2):
      if c == 0:  # Left column
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=True,right=False)
      if c == 1:  # Right column
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=False,right=True)

   if numcases == 'odd':
    for i in range(len(cases)-1):
      if (i % 2) == 0: # all even indices are in left column
       map_ticks_and_labels(ax=axes[i],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=True,right=False)
      if (i % 2) == 1: # all odd indices are in right column
       map_ticks_and_labels(ax=axes[i],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=False,left=False,right=True)

   ### row N
   if numcases == 'even':
    for r in range(rowceil_min1,rowceil):
     for c in range(2):
      if c == 0:  # Left column
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=True,left=True,right=False)
      if c == 1:  # Right column
       map_ticks_and_labels(ax=axes[r,c],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                            xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                            glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                            top=False,bot=True,left=False,right=True)

   if numcases == 'odd':
    i = (len(cases)-1)-1
    map_ticks_and_labels(ax=axes[i],LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                         xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                         glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                         top=False,bot=True,left=True,right=True)

  #----------------------------------
  # Plot the variables 
  #----------------------------------

  if numcases == 'even':
   for r in range(rowceil):
    for c in range(2):
      if cntr_type == 'AreaFill':
       cntr = diff[2*r+c,:,:].plot.contourf(ax=axes[r,c],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                            add_colorbar=False,add_labels=False,extend=extnd)
      elif cntr_type == 'RasterFill':
       cntr = diff[2*r+c,:,:].plot(ax=axes[r,c],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                   add_colorbar=False,add_labels=False,extend=extnd)
      if overlay_vec == True:
       vec = axes[r,c].quiver(uc.lon[skip1D],uc.lat[skip1D],udiff[2*r+c,:,:][skip2D],vdiff[2*r+c,:,:][skip2D],scale_units='height',pivot='mid',
                               scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
       qkey = axes[r,c].quiverkey(vec,0.8,0.2,vec_ref,str(vec_ref)+str(vec_units),labelpos='N',coordinates='figure',color='k')

  if numcases == 'odd':
   for i in range(len(cases)-1):
     if cntr_type == 'AreaFill':
      cntr = diff[i,:,:].plot.contourf(ax=axes[i],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                       add_colorbar=False,add_labels=False,extend=extnd)
     elif cntr_type == 'RasterFill':
      cntr = diff[i,:,:].plot(ax=axes[i],transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                              add_colorbar=False,add_labels=False,extend=extnd)
     if overlay_vec == True:
      vec = axes[i].quiver(uc.lon[skip1D],uc.lat[skip1D],udiff[i,:,:][skip2D],vdiff[i,:,:][skip2D],scale_units='height',pivot='mid',
                              scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
      qkey = axes[i].quiverkey(vec,0.8,0.2,vec_ref,str(vec_ref)+str(vec_units),labelpos='N',coordinates='figure',color='k')

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

   if numcases == 'even':
    for r in range(rowceil):
     for c in range(2):
      axes[r,c].add_patch(Rectangle((wlon-lonadj,slat-latadj),abs(wlon-lonadj)-abs(elon+lonadj),(nlat+latadj)-(slat-latadj),
                                     linestyle=regline,facecolor='none',edgecolor=regcol,zorder=10))
   if numcases == 'odd':
    for i in range(len(cases)-1):
      axes[i].add_patch(Rectangle((wlon-lonadj,slat-latadj),abs(wlon-lonadj)-abs(elon+lonadj),(nlat+latadj)-(slat-latadj),
                                   linestyle=regline,facecolor='none',edgecolor=regcol,zorder=10))

  if point == True:

   # Pick best values for point
   latpt = float(var.lat.sel(lat=ptlat, method='nearest'))
   lonpt = float(var.lon.sel(lon=ptlon, method='nearest'))

   # Handle values if °W = 180-360
   if lonpt > 180.0:
    lonpt = lonpt - 360

   if numcases == 'even':
    for r in range(rowceil):
     for c in range(2):
      axes[r,c].plot(lonpt,latpt,'ro',markersize=5)
   if numcases == 'odd':
    for i in range(len(cases)-1):
     axes[i].plot(lonpt,latpt,'ro',markersize=5)

  #------------------------
  # Add shapes to map  
  #------------------------

  if numcases == 'even':
   for r in range(rowceil):
    for c in range(2):
     for s in range(len(shape_type)):
      axes[r,c].plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                      markeredgecolor=shape_col)

  if numcases == 'odd':
   for i in range(len(cases)-1):
    for s in range(len(shape_type)):
     axes[i].plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                  markeredgecolor=shape_col[s])

  #-----------------------------------------------------------
  # Plot individual maps and append to pdf file if necessary
  #-----------------------------------------------------------

  if Ind_plots == True:

   print('Working on individual difference plots...')

   for i in range(len(cases)-1):

    figi = plt.figure(figsize=(figw,figh),dpi=fdpi)
    axi  = figi.add_subplot(111, projection=proj)
    axi.set_title(fig_let[i]+' '+str(seas)+' '+str(var_name)+str(vec_name)+' '+cases_diff[i],loc=ttlloc,fontsize=ttlfts,pad=ttlpad) 
    axi.set_extent([LonMin,LonMax,LatMin,LatMax],proj)                                   
    if coast == True:
     axi.add_feature(cfeature.COASTLINE,linewidths=coastlw)
    if border == True:
     axi.add_feature(cfeature.BORDERS,linewidths=borderlw)
    if usstate == True:
     axi.add_feature(cfeature.STATES,linewidths=usstatelw)
    map_ticks_and_labels(ax=axi,LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,lattksp=lattksp,lontksp=lontksp,
                         xmin_mj=xmin_mj,ymin_mj=ymin_mj,tkmajlg=tkmajlg,tkminlg=tkminlg,proj=proj,glalpha=glalpha,
                         glcolor=glcolor,lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,tklblsz=tklblsz,
                         top=False,bot=True,left=True,right=True)
    if cntr_type == 'AreaFill':
     cntr = diff[i,:,:].plot.contourf(ax=axi,transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                            add_colorbar=True,cbar_kwargs={'orientation':cbar_orient,'ticks':MultipleLocator(tkstd),
                                            'shrink':cbar_shrk,'pad':cbar_pad,'label':units,'spacing':cbar_sp},add_labels=False,extend=extnd)                     
    elif cntr_type == 'RasterFill':
     cntr = diff[i,:,:].plot(ax=axi,transform=proj,cmap=colort,levels=np.arange(loval,hival+spval,spval),
                                    add_colorbar=True,cbar_kwargs={'orientation':cbar_orient,'ticks':MultipleLocator(tkstd),
                                    'shrink':cbar_shrk,'pad':cbar_pad,'label':units},add_labels=False,extend=extnd)
    if overlay_vec == True:
     vec = axi.quiver(uc.lon[skip1D],uc.lat[skip1D],udiff[i,:,:][skip2D],vdiff[i,:,:][skip2D],scale_units='height',pivot='mid',scale=vec_scale,
                                      width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
     qkey = axi.quiverkey(vec,0.8,0.2,vec_ref,str(vec_ref)+str(vec_units),labelpos='N',coordinates='figure',color='k')
    if regbox == True:
     axi.add_patch(Rectangle((wlon-lonadj,slat-latadj),abs(wlon-lonadj)-abs(elon+lonadj),(nlat+latadj)-(slat-latadj),
                                     linestyle=regline,facecolor='none',edgecolor=regcol,zorder=10))
    if point == True:
     axi.plot(lonpt,latpt,'ro',markersize=5)
    for s in range(len(shape_type)):
     axi.plot(shape_lons[s],shape_lats[s],marker=shape_type[s],markersize=shape_size[s],markerfacecolor='none',
                      markeredgecolor=shape_col[s])

  #--------------------------
  # Save output PDF file
  #--------------------------

  save_multi_image('./'+str(folderpath)+'/'+str(seas)+'_diff_'+str(var_name)+str(vec_name)+'_'+str(extra_name)+str(filesuf))
  plt.close('all')

