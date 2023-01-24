#######################################################################################################
# 
# These functions create map plot visualizations of water tagged CESM results.
#
#######################################################################################################

import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geocat.viz.util as gv
import cmaps
from watertagging.tagged_regions import draw_land_tags, draw_ocean_tags
from plotting.plot_functions import save_multi_image,map_ticks_and_labels,draw_region_box
from warnings import simplefilter

# Default arrays
default_wgt_mon = xr.DataArray(np.ones(12),dims=['time']).astype(float)

#######################################################################################################
# Plot tag values for precipitation, precipitation percentage, and d18Op on a global map 
#######################################################################################################

def watertagging_values_on_map(precip: xr.DataArray, d18Op: xr.DataArray, case: str, tagnames: list, 
                               num_landtags: int, num_oceantags: int, path: str, season: str, 
                               lat: xr.DataArray, lon: xr.DataArray,
                               landlat: list, landlon: list, oceanlat: list, oceanlon: list, 
                               lo_cutoff: float = 0.5,
                               figw: float = 10.,figh: float = 3., fdpi: float = 300.,
                               wspace: float = -0.25, hspace: float = 0.05,
                               LatMin: float = -90., LatMax: float = 90.,
                               LonMin: float = -180., LonMax: float = 180.,
                               proj: ccrs.Projection = ccrs.PlateCarree(),
                               landcol: colors.ListedColormap = cmaps.MPL_gist_earth_r(np.arange(0,10,1)), 
                               oceancol: colors.ListedColormap = cmaps.WhiteBlue_r(np.arange(149,254,1)),
                               p_units: str = 'mm/day', o_units: str = 'per mil',
                               ttlfs: float = 8., ttlloc: str = 'left',
                               coast: bool = True, coastlw: float = 0.3, 
                               lake: bool = True, lakelw: float = 0.1,
                               border: bool = False, borderlw: float = 0.3,
                               usstate: bool = False, usstatelw: float = 0.3,
                               slat: float = None, nlat: float = None, 
                               wlon: float = None, elon: float = None,
                               regline: str = '-', regcol: str = 'r', reglw: float = 0.5, 
                               rnd_precip: int = 2, rnd_pctpr: int = 1, rnd_d18Op: int = 2, 
                               tag_fs: float = 3.,tag_lw: float = 0.5, 
                               tag_maj: str = '-', tag_min: str = '--', 
                               tag_col: str = 'k', tag_zorder: int = 5,
                               diff: bool = False, cntlp: xr.DataArray = None, cntlo: xr.DataArray = None,
                               folderpath: str = '.', reg_name: str = '', extra_name: str = '',
                               filesuf: str = '.pdf'):

  '''Creates three map plots showing values of 1) precipitation, 2) precipitation percentage, and
  3) d18Op for land and ocean water tag regions.

  Parameters
  --------------------------------------------------------------------------------------------------
  precip
   class: 'xarray.DataArray', variable containing averaged precipitation values for each tag region

  d18Op
   class: 'xarray.DataArray', variable containing averaged d18Op values for each tag region

  case
   class: 'string', name of the case being plotted

  tagnames
   class: 'list', list of strings indicating each tag region's name

  num_landtags, num_oceantags
   class: 'int', number of land and ocean tags

  path
   class: 'string', path to file that contains 'LANDFRAC' variable for coloring land/ocean

  season
   class: 'string', name of season for which the values are averaged over

  lat
   class: 'xarray.DataArray', latitude coordinate array

  lon
   class: 'xarray.DataArray', longitude coordinate array

  landlat, landlon
   class: 'list', list latitude and longitude values for plotting text on map for each land region

  oceanlat, oceanlon
   class: 'list', list latitude and longitude values for plotting text on map for each ocean region

  lo_cutoff (default = 0.5)
   class: 'float', cutoff proportion for plotting land vs. ocean, example: 0.5 = 50% land

  figw, figh, fdpi (default = 10., 3., 300.)
   class: 'float', figure width (in inches, height (in inches), and dpi 

  wspace, hspace (default = -0.25, 0.05)
   class: 'float', width and height of padding between subplots

  LatMin, LatMax, LonMin, LonMax (Default = -90., 90., -180., 180.)
   class: 'float', Coordinates defining visual extent of plot map. Values: - = °S or °W, + = °N or °E 

  proj (default = ccrs.PlateCarree())
   class: 'cartopy.crs.Projection', projection of map for each subplot

  landcol (default = cmaps.MPL_gist_earth_r(np.arange(0,10,1)))
   class: 'colors.ListedColormap', color table optimized to show land area

  oceancol (default = cmaps.WhiteBlue_r(np.arange(149,254,1)))
   class: 'colors.ListedColormap', color table optimized to show ocean area

  p_units (default = 'mm/day')
   class: 'string', units to be shown for precipitation 

  o_units (default = 'per mil')
   class: 'string', units to be shown for d18Op         

  ttlfs (default = 8.)
   class: 'float', font size for title above each subplot

  ttlloc (default = 'left')
   class: 'string', location of title above each subplot

  coast, coastlw (default = True, 0.3)
   class: 'bool','string', include coastlines on base map and if True, width of lines

  lake, lakelw (default = True, 0.1)
   class: 'bool','string', include lakes on base map and if True, width of lines

  border, borderlw (default = False, 0.3)
   class: 'bool','string', include country borders on base map and if True, width of lines

  usstate, usstatelw (default = False, 0.3)
   class: 'bool','string', include US state boundaries on base map and if True width of lines

  slat, nlat, wlon, elon (default = None, None, None, None)
   class: 'float', southern/northern latitude and western/eastern longitude defined region

  regline (default = '-')
   class: 'string', line style of line marking region boundary

  regcol (default = 'r')
   class: 'string', color of line marking region boundary

  reglw (default = 0.5)
   class: 'float', width of line marking region boundary

  rnd_precip, rnd_pctpr, rnd_d18Op (default = 2, 1, 2)
   class: 'int', number of significant digits for precip, precip%, and d18Op values shown
                 on global maps 

  tag_fs (default = 3.)
   class: 'float', font size of text indicating tag region values

  tag_lw (default = 0.5)
   class: 'float', width of lines delineating tag regions

  tag_maj (default = '-')
   class: 'string', line style of major lines delineating tag regions

  tag_min (default = '--')
   class: 'string', line style of minor lines delineating tag regions

  tag_col (default = 'k')
   class: 'string', color of lines delineating tag regions

  tag_zorder (default = 5)
   class: 'int', zorder of lines delineating tag regions

  diff (default = False)
   class: 'bool', toggle to True if plotting text values of a difference between two cases

  cntlp (default = None)
   class: 'xr.DataArray', if diff=True, this is the precip control variable, 'precip' variable
                          now becomes the experiment variable

  cntlo (default = None)
   class: 'xr.DataArray', if diff=True, this is the d18Op control variable, 'd18Op' variable
                          now becomes the experiment variable
                               
  folderpath (default = '.')
   class: 'string', path to the folder with which to deposit the output file
 
  reg_name (default = '')
   class: 'string', name of region for which water tagging results are calculated, added to 
                    output file

  extra_name (default = '')
   class: 'string', extra text to add at end of output file

  filesuf (default = '.pdf')
   class: 'string', type of output file 
  
  '''
  
  #--------------------------------------------
  # Handle variables if diff == True 
  #--------------------------------------------
  
  # If diff == True then we need an experiment case and a control case (cntlp, cntlo)   
  if diff == True:
   exptp = precip
   precip = precip - cntlp
   expto = d18Op
   d18Op = d18Op - cntlo

  #----------------------
  # Define plots
  #----------------------

  # Make 3 plots: precip, precip%, d18Op
  fig = plt.figure(figsize=(figw,figh),dpi=fdpi)
  ax1 = fig.add_subplot(231, projection=proj) # precip (land)
  ax2 = fig.add_subplot(234, projection=proj) # precip (ocean)
  ax3 = fig.add_subplot(232, projection=proj) # precip% (land)
  ax4 = fig.add_subplot(235, projection=proj) # precip% (ocean)
  ax5 = fig.add_subplot(233, projection=proj) # d18Op (land)
  ax6 = fig.add_subplot(236, projection=proj) # d18Op (ocean)
  plt.subplots_adjust(wspace=wspace,hspace=hspace)
  
  # Define land, set anything >[lo_cutoff] land as "land", else as "ocean"
  data_landfrac = xr.open_dataset(path)
  landfrac_all  = data_landfrac.LANDFRAC[0,:,:]
  landfrac_bin  = xr.where(landfrac_all > lo_cutoff, 1.0, 0.0)

  # Define land and ocean colors 
  LandColorTable  = colors.LinearSegmentedColormap.from_list('name',landcol) 
  OceanColorTable = colors.LinearSegmentedColormap.from_list('name',oceancol) 
 
  # Default is global plot, will change if LonMin, etc. are modified
  ax1.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
  ax2.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
  ax3.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
  ax4.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
  ax5.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
  ax6.set_extent([LonMin,LonMax,LatMin,LatMax],proj)
 
  #-------------------------------------
  # ax1 (land) and ax2 (ocean): precip 
  #-------------------------------------

  # Title for land+ocean plots
  ax1.set_title(str(case)+' precip ('+p_units+')',fontsize=ttlfs,loc=ttlloc)

  # Map features
  if coast == True:
   ax1.add_feature(cfeature.COASTLINE,linewidths=coastlw)
   ax2.add_feature(cfeature.COASTLINE,linewidths=coastlw)
  if lake == True:
   ax1.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
   ax2.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
  if border == True:
   ax1.add_feature(cfeature.BORDERS,linewidths=borderlw)
   ax2.add_feature(cfeature.BORDERS,linewidths=borderlw)
  if usstate == True:
   ax1.add_feature(cfeature.STATES,linewidths=usstatelw)
   ax2.add_feature(cfeature.STATES,linewidths=usstatelw)

  # Plot land and ocean colors
  landfrac_bin.plot(ax=ax1,transform=ccrs.PlateCarree(),cmap=LandColorTable,add_colorbar=False,add_labels=False)
  landfrac_bin.plot(ax=ax2,transform=ccrs.PlateCarree(),cmap=OceanColorTable,add_colorbar=False,add_labels=False)

  #-------------------------------------
  # ax3 (land) and ax4 (ocean): precip%
  #-------------------------------------

  # Title for land+ocean plots
  ax3.set_title(str(case)+' precip (%) ',fontsize=ttlfs,loc=ttlloc)

  # Map features
  if coast == True:
   ax3.add_feature(cfeature.COASTLINE,linewidths=coastlw)
   ax4.add_feature(cfeature.COASTLINE,linewidths=coastlw)
  if lake == True:
   ax3.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
   ax4.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
  if border == True:
   ax3.add_feature(cfeature.BORDERS,linewidths=borderlw)
   ax4.add_feature(cfeature.BORDERS,linewidths=borderlw)
  if usstate == True:
   ax3.add_feature(cfeature.STATES,linewidths=usstatelw)
   ax4.add_feature(cfeature.STATES,linewidths=usstatelw)

  # Plot land and ocean colors
  landfrac_bin.plot(ax=ax3,transform=ccrs.PlateCarree(),cmap=LandColorTable,add_colorbar=False,add_labels=False)
  landfrac_bin.plot(ax=ax4,transform=ccrs.PlateCarree(),cmap=OceanColorTable,add_colorbar=False,add_labels=False)

  #-------------------------------------
  # ax5 (land) and ax6 (ocean): d18Op  
  #-------------------------------------

  # Title for land+ocean plots
  ax5.set_title(str(case)+r' $\delta^{18}O_{P}$ ('+o_units+')',fontsize=ttlfs,loc=ttlloc)

  # Map features
  if coast == True:
   ax5.add_feature(cfeature.COASTLINE,linewidths=coastlw)
   ax6.add_feature(cfeature.COASTLINE,linewidths=coastlw)
  if lake == True:
   ax5.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
   ax6.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
  if border == True:
   ax5.add_feature(cfeature.BORDERS,linewidths=borderlw)
   ax6.add_feature(cfeature.BORDERS,linewidths=borderlw)
  if usstate == True:
   ax5.add_feature(cfeature.STATES,linewidths=usstatelw)
   ax6.add_feature(cfeature.STATES,linewidths=usstatelw)
  
  # Plot land and ocean colors
  landfrac_bin.plot(ax=ax5,transform=ccrs.PlateCarree(),cmap=LandColorTable,add_colorbar=False,add_labels=False)
  landfrac_bin.plot(ax=ax6,transform=ccrs.PlateCarree(),cmap=OceanColorTable,add_colorbar=False,add_labels=False)

  #----------------------------
  # Add region/point to plots 
  #----------------------------

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

  # Process coordinates
  latadj = abs(np.float64((lat[2]-lat[1])/2))
  lonadj = abs(np.float64((lon[2]-lon[1])/2))

  # Define perimeter
  rlats = np.float64(lat[lats]-latadj)
  rlatn = np.float64(lat[latn]+latadj)
  rlonw = np.float64(lon[lonw]-lonadj)
  rlone = np.float64(lon[lone]+lonadj)
  if rlonw > 180:
    rlonw = rlonw - 360
  if rlone > 180:
    rlone = rlone - 360

  # Plot rectangles
  draw_region_box(ax=ax1,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                  facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)
  draw_region_box(ax=ax2,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                  facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)
  draw_region_box(ax=ax3,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                  facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)
  draw_region_box(ax=ax4,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                  facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)
  draw_region_box(ax=ax5,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                  facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)
  draw_region_box(ax=ax6,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                  facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=10)

  #----------------------------
  # Add text for values
  #----------------------------

  # Precip text on map 
  text_precip = np.float64(precip)

  # Precip percentage text on map depends on diff
  if diff == False:
   text_pctpr = np.float64((precip/np.sum(precip))*100.)
  elif diff == True: # calculate each separately and subtract 
   text_pctpr = np.float64((exptp/np.sum(exptp))*100.) - np.float64((cntlp/np.sum(cntlp))*100.) 

  # d18Op text on map
  text_d18Op  = np.float64(d18Op)

  # Loop through land values
  for l in range(num_landtags):
   ax1.text(landlon[l],landlat[l],text_precip[l].round(rnd_precip),fontsize=tag_fs,ha='center',va='center')
   ax3.text(landlon[l],landlat[l],text_pctpr[l].round(rnd_pctpr),fontsize=tag_fs,ha='center',va='center')
   if text_d18Op[l] > 0:
    ax5.text(landlon[l],landlat[l],'+'+str(text_d18Op[l].round(rnd_d18Op)),color='r',fontsize=tag_fs,
                                                                                     ha='center',va='center')
   elif text_d18Op[l] < 0:
    ax5.text(landlon[l],landlat[l],text_d18Op[l].round(rnd_d18Op),color='b',fontsize=tag_fs,ha='center',va='center')
   else:
    ax5.text(landlon[l],landlat[l],text_d18Op[l].round(0),color='k',fontsize=tag_fs,ha='center',va='center')


  # Loop through ocean values
  for o in range(num_landtags,len(tagnames),1):
   ax2.text(oceanlon[o-num_landtags],oceanlat[o-num_landtags],text_precip[o].round(rnd_precip),fontsize=tag_fs,
                                                                                              ha='center',va='center')
   ax4.text(oceanlon[o-num_landtags],oceanlat[o-num_landtags],text_pctpr[o].round(rnd_pctpr),fontsize=tag_fs,
                                                                                             ha='center',va='center')
   if text_d18Op[o] > 0:
    ax6.text(oceanlon[o-num_landtags],oceanlat[o-num_landtags],'+'+str(text_d18Op[o].round(rnd_d18Op)),
                                                                   color='r',fontsize=tag_fs,ha='center',va='center')
   elif text_d18Op[o] < 0:
    ax6.text(oceanlon[o-num_landtags],oceanlat[o-num_landtags],text_d18Op[o].round(rnd_d18Op),color='b',fontsize=tag_fs,
                                                                                              ha='center',va='center')
   else:
    ax6.text(oceanlon[o-num_landtags],oceanlat[o-num_landtags],text_d18Op[o].round(0),color='k',fontsize=tag_fs,
                                                                                      ha='center',va='center')

  #----------------------------
  # Add tag regions to plots 
  #----------------------------

  for ltag in range(num_landtags):
   draw_land_tags(numtag=ltag,ax=ax1,lw=tag_lw,major=tag_maj,minor=tag_min,color=tag_col,zorder=tag_zorder)
   draw_land_tags(numtag=ltag,ax=ax3,lw=tag_lw,major=tag_maj,minor=tag_min,color=tag_col,zorder=tag_zorder)
   draw_land_tags(numtag=ltag,ax=ax5,lw=tag_lw,major=tag_maj,minor=tag_min,color=tag_col,zorder=tag_zorder)
  for otag in range(num_landtags,len(tagnames),1):
   draw_ocean_tags(numtag=otag,ax=ax2,lw=tag_lw,major=tag_maj,minor=tag_min,color=tag_col,zorder=tag_zorder)
   draw_ocean_tags(numtag=otag,ax=ax4,lw=tag_lw,major=tag_maj,minor=tag_min,color=tag_col,zorder=tag_zorder)
   draw_ocean_tags(numtag=otag,ax=ax6,lw=tag_lw,major=tag_maj,minor=tag_min,color=tag_col,zorder=tag_zorder)

  #---------------------
  # Save as pdf
  #---------------------

  # Setting file output name parameters
  if reg_name != '':
   reg_name = ''.join(('_',reg_name))
  if extra_name != '':
   extra_name = ''.join(('_',extra_name))

  # Name output file
  plt.savefig('./'+str(folderpath)+'/'+str(case)+'_'+str(season)+'_watertag_values'+reg_name+extra_name+str(filesuf),
              bbox_inches='tight')

  # Close figure to retain memory
  plt.close()

#######################################################################################################
# Plot precip and d18Op maps for each tagged region 
#######################################################################################################

def plot_tagged_precip_and_d18Op(P: bool, O: bool, prect: xr.DataArray, d18Op: xr.DataArray, 
                                 case: str, tagnames: list, num_landtags: int, num_oceantags: int,
                                 season: str, lat: xr.DataArray, lon: xr.DataArray,
                                 cutoff: float = 0.001, p_units: str = 'mm/day', o_units: str = 'per mil',
                                 figw: float = 10., figh: float = 10., fdpi: float = 300.,
                                 LatMin: float = -90., LatMax: float = 90.,
                                 LonMin: float = -180., LonMax: float = 180.,
                                 proj: ccrs.Projection = ccrs.PlateCarree(),
                                 coast: bool = True, coastlw: float = 0.5,
                                 lake: bool = True, lakelw: float = 0.3,
                                 border: bool = False, borderlw: float = 0.5,
                                 usstate: bool = False, usstatelw: float = 0.5,
                                 add_axes: bool = True, topax: bool = False, botax: bool = True,
                                 leftax: bool = True, rightax: bool = True,
                                 latlblsp: float = 30., lonlblsp: float = 60.,
                                 lattksp: int = 7, lontksp: int = 7, tklblsz: float = 10., 
                                 tkmajlg: float = 4., tkminlg: float = 2., xmin_mj: int = 2, 
                                 ymin_mj: int = 2, xpads: float = 15., ypads: float = 15., 
                                 glalpha: float = 0.0, glcolor: str = 'gray', 
                                 cntr_type: str = 'RasterFill', 
                                 colorp: colors.ListedColormap = cmaps.MPL_viridis, 
                                 coloro: colors.ListedColormap = cmaps.MPL_rainbow,
                                 p_loval: float = 0., p_hival: float = 2., p_spval: float = 0.1,
                                 p_mantick: list = [0.001,0.2,0.5,1.0,1.5,2.0], p_extnd: str = 'max',
                                 o_loval: float = -5., o_hival: float = 0., o_spval: float = 0.2,
                                 o_tkstd: float = 1., o_extnd: str = 'both', 
                                 cbar_orient: str = 'horizontal', cbar_pad: float = 0.08, 
                                 cbar_shrk: float = 0.6, cbar_lblsz: float = 10.,
                                 ttlfs: float = 12., ttlloc: str = 'left', ttlpad: str = 15.,
                                 regbox: bool = True, regcol: str = 'r', regline: str = '-', 
                                 reglw: float = 2., regzorder: float = 10.,
                                 slat: float = None, nlat: float = None, 
                                 wlon: float = None, elon: float = None,
                                 taglw: float = 1., tagmaj: str = '-', tagmin: str = '--',
                                 tagcol: str = 'k', tagzorder: int = 3, 
                                 overlay_vec: bool = False, u: xr.DataArray = None, v: xr.DataArray = None,
                                 vec_scale: float = 200., vec_wid: float = 0.001, vec_hal: float = 3.,
                                 vec_hdl: float = 3., vec_hdw: float = 4., vec_skip: int = 4,
                                 vec_ref: float = 10., vec_units: str = '', vec_name: str = '',
                                 diff: bool = False, folderpath: str = '.', reg_name: str = '', 
                                 extra_name: str = '',filesuf: str = '.pdf'
                                 ):

  '''Creates a multipage PDF file with a precipitation and d18Op map for each tagged region

  Parameters
  --------------------------------------------------------------------------------------------------
  P
   class: 'bool', If True, precipitation maps will be included in output file.

  O
   class: 'bool', If True, d18Op maps will be included in output file.

  prect 
   class: 'xarray.DataArray', variable containing averaged precipitation values for each tag region

  d18Op
   class: 'xarray.DataArray', variable containing averaged d18Op values for each tag region

  case
   class: 'string', name of the case being plotted

  tagnames
   class: 'list', list of strings indicating each tag region's name

  num_landtags, num_oceantags
   class: 'int', number of land and ocean tags

  season
   class: 'string', name of season for which the values are averaged over

  lat
   class: 'xarray.DataArray', latitude coordinate array

  lon
   class: 'xarray.DataArray', longitude coordinate array

  cutoff (Default = 0.001)
   class: 'float', On maps, locations where precipitation is below this value will not be shown.

  p_units (Default = 'mm/day')
   class: 'string', Units for precipitation.

  o_units (Default = 'per mil')
   class: 'string', Units for d18Op

  figw, figh, fdpi (default = 10., 10., 300.)
   class: 'float', figure width (in inches, height (in inches), and dpi

  LatMin, LatMax, LonMin, LonMax (Default = -90., 90., -180., 180.)
   class: 'float', Coordinates defining visual extent of plot map. Values: - = °S or °W, + = °N or °E 

  proj (Default = ccrs.PlateCarree())
   class:'cartopy.crs.Projection', Projection of map. This parameter must be cartopy.crs.PlateCarree() 
                                   currently as cartopy's gridlines do not support any other projection 
                                   at this time. 

  coast, coastlw (default = True, 0.5)
   class: 'bool','string', include coastlines on base map and if True, width of lines

  lake, lakelw (default = True, 0.3)
   class: 'bool','string', include lakes on base map and if True, width of lines

  border, borderlw (default = False, 0.5)
   class: 'bool','string', include country borders on base map and if True, width of lines

  usstate, usstatelw (default = False, 0.5)
   class: 'bool','string', include US state boundaries on base map and if True width of lines

  add_axes (Default = True)
   class: 'bool', Add axes ticks and labels to each map. 
                    
  topax, botax, leftax, rightax (Default = False, True, True, True)
   class: 'bool', If True, axes on that side of the map are active.  

  latlblsp, lonlblsp (Default = 30., 60.)
   class: 'float', Spacing of latitude tick labels in degrees of latitude and longitude.

  lattksp, lontksp (Default = 7, 7)
   class: 'int', Tick spacing for major latitude and longitude ticks. Examples:
                 lattksp = 7 -> -90 to 90 by 30 (7 total major ticks)
                 lontksp = 7 -> -180 to 180 by 60 (7 total major ticks)

  tklblsz (Default = 10.)
   class: 'float', Font size of tick labels.

  tkmajlg (Default = 4.)
   class: 'float', Size of major ticks around the outline of each map.       

  tkminlg (Default = 2.)
   class: 'float', Size of minor ticks around the outline of each map.    

  xmin_mj (Default = 2)
   class: 'int', Number of minor ticks in between each major tick on the x-axis (longitude).

  ymin_mj (Default = 2)
   class: 'int', Number of minor ticks in between each major tick on the y-axis (latitude). 

  xpads (Default = 15.)
   class: 'float', Padding of x-axis longitude tick labels from map outline.

  ypads (Default = 15.)
   class: 'float', Padding of y-axis latitude tick labels from map outline.                

  glalpha (Default = 0.0 (no gridlines))
   class: 'float', Alpha (line opacity) for gridlines between lat/lon ticks. Default is no gridlines. Modify 
                   1.0 >= alpha > 0.0 for visible gridlines.

  glcolor (Default = 'gray')
   class: 'str', Color of gridlines if glalpha > 0.0.

  cntr_type (Default = 'AreaFill')
   class: 'str', Type of contour for plotting. Options: 'AreaFill' - interpolated, 'RasterFill' - raster. 

  colorp (Default = cmaps.MPL_viridis)
   class: 'colors.ListedColormap', Color table for precipitation maps.

  coloro (Default = cmaps.MPL_rainbow) 
   class: 'colors.ListedColormap' Color table for d18Op maps. 

  p_loval, p_hival, p_spval (Default = 0., 2., 0.1)
   class: 'float', Low value, high value, and spacing value for precipitation contours.

  p_mantick (Default = [0.001,0.2,0.5,1.0,1.5,2.0])
   class: 'list', Manual color bar tick label values. These allow for manual setting of ticks in order to
                  include cutoff value as a color bar tick label. 

  p_extnd (Default = 'max')
   class: 'str', Location of triangles on precipitation color bar.

  o_loval, o_hival, o_spval (Default = -5., 0., 0.2)
   class: 'float', Low value, high value, and spacing value for d18Op contours.

  o_tkstd (Default = 1.)
   class: 'float', Tick stride for color bar for d18Op maps.

  o_extnd (Default = 'both')
   class: 'string', Location of triangles on d18Op color bar.

  cbar_orient (Default = 'horizontal')
   class: 'str', Orientation of color bar. Options: 'horizontal','vertical'

  cbar_pad (Default = 0.08)
   class: 'float', Padding of color bar from bottom of map.

  cbar_shrk (Default = 0.6)
   class: 'float', Color bar shrink parameter. Modifies the sizing of the color bar.

  cbar_lblsz (Default = 10.)
   class: 'float', Label size for color bar tick labels.

  ttlfs (Default = 12.)
   class: 'float', Font size for title above map. 

  ttlloc (Default = 'left')
   class: 'string', Location of title above map.

  ttlpad (Default = 15.)
   class: 'string', Padding of title above map.

  regbox (Default = True)
   class: 'bool', If True, box will be drawn over region for which water tagging results are 
                  calculated.

  regcol (Default = 'r')
   class: 'string', Color of line marking region boundary.

  regline (Default = '-')
   class: 'string', Line style of line marking region boundary.

  reglw (default = 2.)
   class: 'float', Width of line marking region boundary.

  regzorder (Default = 10.)
   class: 'float', Z-order of region boundary.

  slat, nlat, wlon, elon (Default = None, None, None, None)
   class: 'float', Southern/northern latitude and western/eastern longitude defined region.

  tag_lw (Default = 0.5)
   class: 'float', Width of lines delineating tag regions.

  tag_maj (Default = '-')
   class: 'string', Line style of major lines delineating tag regions.

  tag_min (Default = '--')
   class: 'string', Line style of minor lines delineating tag regions.

  tag_col (Default = 'k')
   class: 'string', Color of lines delineating tag regions.

  tag_zorder (Default = 5)
   class: 'int', Z-order of lines delineating tag regions.

  overlay_vec (Default = False)
   class: 'bool', If True, overlays vectors on map plots. 

  u, v (Default = None, None)
   class: 'xr.DataArray', If overlay_vec is True, these are xr.DataArray's of the u- and v-components.

  vec_scale (Default = 200.)
   class: 'float', Scales the length of the arrow inversely. Number of data units per arrow lenth unit. 
  
  vec_wid (Default = 0.001)
   class: 'float', Shaft width of vector. 

  vec_hal (Default = 3.)
   class: 'float', Head axis length of vector.

  vec_hdl (Default = 3.)
   class: 'float', Head length of vector.

  vec_hdw (Default = 4.)
   class: 'float', Head width of vector.

  vec_skip (Default = 4)
   class: 'int', Number of grid cells to skip when plotting vectors. Higher numbers correspond to less 
                 congestion between plotted vector arrows. 
 
  vec_ref (Default = 10.) 
   class: 'float', Reference vector to be used as the length of the vector key.           

  vec_units (Default = '')
   class: 'str', Units of the vector that will be displayed in vector key following 'vec_ref'.

  vec_name (Default = '')
   class: 'str',  Manual text to be included in output file name containing information about vector. 
                  For example, vec_name = '850hPawinds' for u and v wind plot at 850 hPa.

  diff (Default = False)
   class: 'bool', Toggle to True if plotting a difference between two cases. This option makes it easier
                  to switch between color tables for plotting differences.

  folderpath (Default = '.')
   class: 'string', Path to the folder with which to deposit the output file.

  reg_name (default = '')
   class: 'string', name of region for which water tagging results are calculated, added to 
                    output file

  extra_name (default = '')
   class: 'string', extra text to add at end of output file

  filesuf (Default = '.pdf')
   class: 'string', Type of output file.

  '''

  #---------------------------------------------------
  # Determine contour level spacing based on 'diff'
  #---------------------------------------------------

  # State of 'diff' determines how the color bar tick labels are specified

  if diff == False:
   # Add color bar tick label for cutoff value
   p_levels = np.insert(np.arange(p_loval,p_hival+p_spval,p_spval),1,cutoff)

  elif diff == True:
   # Normal color bar tick label spacing (low to high by spacing value)
   p_levels = np.arange(p_loval,p_hival+p_spval,p_spval)

  #------------------------------------------
  # Specify vector parameters if necessary
  #------------------------------------------

  # Setting vector grid cell skipping parameters
  skip1D = (slice(None, None, vec_skip))
  skip2D = (slice(None, None, vec_skip), slice(None, None, vec_skip))

  #----------------------------
  # Add region/point to plots 
  #----------------------------

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

  # Process coordinates
  latadj = abs(np.float64((lat[2]-lat[1])/2))
  lonadj = abs(np.float64((lon[2]-lon[1])/2))

  # Define perimeter
  rlats = np.float64(lat[lats]-latadj)
  rlatn = np.float64(lat[latn]+latadj)
  rlonw = np.float64(lon[lonw]-lonadj)
  rlone = np.float64(lon[lone]+lonadj)
  if rlonw > 180:
    rlonw = rlonw - 360
  if rlone > 180:
    rlone = rlone - 360
 
  #------------------------------------------------
  # Manage precipitation cutoff value if specified
  #------------------------------------------------

  if cutoff > 0.:
   d18Op = xr.where(prect < cutoff, float('NaN'), d18Op)

  #-----------------------------------------------------
  # Make variables to plot cyclic
  #-----------------------------------------------------
 
  # Adding cyclical longitude
  prect_c = gv.xr_add_cyclic_longitudes(prect,'lon')
  d18Op_c = gv.xr_add_cyclic_longitudes(d18Op,'lon')
  if overlay_vec == True:
   uc = gv.xr_add_cyclic_longitudes(u,'lon')
   vc = gv.xr_add_cyclic_longitudes(v,'lon')

  #--------------------------------------------------------
  # Make precipitation plots if specified
  #--------------------------------------------------------

  if P == True:

   for tag in range(len(tagnames)):

    # Ignore unnecessary warning from matplotlib about having multiple figures open
    simplefilter(action="ignore", category=RuntimeWarning)

    # Define figure and axis
    figp = plt.figure(figsize=(figw,figh),dpi=fdpi)
    axp = figp.add_subplot(111, projection=proj)

    # Default is global plot, will change if LonMin, etc. are modified
    axp.set_extent([LonMin,LonMax,LatMin,LatMax],proj)

    # Add features to map as indicated in parameters
    if coast == True:
     axp.add_feature(cfeature.COASTLINE,linewidths=coastlw) 
    if lake == True:
     axp.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k') 
    if border == True:
     axp.add_feature(cfeature.BORDERS,linewidths=borderlw) 
    if usstate == True:
     axp.add_feature(cfeature.STATES,linewidths=usstatelw) 

    # Add ticks and labels to map as indicated in parameters
    map_ticks_and_labels(ax=axp,LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,
                         lattksp=lattksp,lontksp=lontksp,xmin_mj=xmin_mj,ymin_mj=ymin_mj,
                         tkmajlg=tkmajlg,tkminlg=tkminlg,glalpha=glalpha,glcolor=glcolor,
                         lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,
                         tklblsz=tklblsz,proj=proj,top=topax,bot=botax,left=leftax,right=rightax)

    # Draw rectangle over region for which water tagging results will be calculated
    draw_region_box(ax=axp,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                    facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=regzorder)

    # Title above precipitation plots
    axp.set_title(str(case)+' Precip ('+p_units+') for '+str(tag+1)+'. '+tagnames[tag],
                  fontsize=ttlfs,loc=ttlloc,pad=ttlpad)

    # Make contour plots 
    if cntr_type == 'AreaFill': 
     prect_c[tag,:,:].plot.contourf(ax=axp,transform=ccrs.PlateCarree(),cmap=colorp,
                        levels=p_levels,
                        add_colorbar=True,add_labels=False,extend=p_extnd,
                        cbar_kwargs={'orientation':cbar_orient,'shrink':cbar_shrk,'pad':cbar_pad,
                        'label':p_units,'ticks':p_mantick,'format': '{}'.format})
    elif cntr_type == 'RasterFill':
     prect_c[tag,:,:].plot(ax=axp,transform=ccrs.PlateCarree(),cmap=colorp,
                        levels=p_levels,
                        add_colorbar=True,add_labels=False,extend=p_extnd,
                        cbar_kwargs={'orientation':cbar_orient,'shrink':cbar_shrk,'pad':cbar_pad,
                        'label':p_units,'ticks':p_mantick,'format': '{}'.format})

    # Overlay wind vectors if turned on
    if overlay_vec == True:
     vecp = axp.quiver(prect_c.lon[skip1D],prect_c.lat[skip1D],uc[skip2D],vc[skip2D],scale_units='height',pivot='mid',
               scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
     qkeyp = axp.quiverkey(vecp,0.8,0.2,vec_ref,str(int(vec_ref))+str(vec_units),labelpos='N',coordinates='figure',color='k')

    # Draw corresponding land or ocean tag region on plot
    if tag <= (num_landtags-1):
     draw_land_tags(numtag=tag,ax=axp,lw=taglw,major=tagmaj,minor=tagmin,color=tagcol,zorder=tagzorder)
    elif tag > (num_landtags-1):
     draw_ocean_tags(numtag=tag,ax=axp,lw=taglw,major=tagmaj,minor=tagmin,color=tagcol,zorder=tagzorder)

  #--------------------------------------------------------
  # Make d18Op plots if specified
  #--------------------------------------------------------

  if O == True:

   for tag in range(len(tagnames)):

    # Ignore unnecessary warning from matplotlib about having multiple figures open
    simplefilter(action="ignore", category=RuntimeWarning)

    # Define figure and axis
    figo = plt.figure(figsize=(figw,figh),dpi=fdpi)
    axo = figo.add_subplot(111, projection=proj)

    # Default is global plot, will change if LonMin, etc. are modified
    axo.set_extent([LonMin,LonMax,LatMin,LatMax],proj)

    # Add features to map as indicated in parameters
    if coast == True:
     axo.add_feature(cfeature.COASTLINE,linewidths=coastlw)
    if lake == True:
     axo.add_feature(cfeature.LAKES,linewidths=lakelw,facecolor='none',edgecolor='k')
    if border == True:
     axo.add_feature(cfeature.BORDERS,linewidths=borderlw)
    if usstate == True:
     axo.add_feature(cfeature.STATES,linewidths=usstatelw)

    # Add ticks and labels to map as indicated in parameters
    map_ticks_and_labels(ax=axo,LatMin=LatMin,LatMax=LatMax,LonMin=LonMin,LonMax=LonMax,
                         lattksp=lattksp,lontksp=lontksp,xmin_mj=xmin_mj,ymin_mj=ymin_mj,
                         tkmajlg=tkmajlg,tkminlg=tkminlg,glalpha=glalpha,glcolor=glcolor,
                         lonlblsp=lonlblsp,latlblsp=latlblsp,xpads=xpads,ypads=ypads,
                         tklblsz=tklblsz,proj=proj,top=topax,bot=botax,left=leftax,right=rightax)

    # Draw rectangle over region for which water tagging results will be calculated
    draw_region_box(ax=axo,slat=rlats,nlat=rlatn,wlon=rlonw,elon=rlone,linestyle=regline,
                    facecolor='none',edgecolor=regcol,linewidth=reglw,zorder=regzorder)

    # Title above d18Op plots
    axo.set_title(str(case)+r' $\delta^{18} O_{P}$ ('+o_units+') for '+str(tag+1)+'. '+tagnames[tag],
                  fontsize=ttlfs,loc=ttlloc,pad=ttlpad)

    # Make contour plots 
    if cntr_type == 'AreaFill':
     d18Op_c[tag,:,:].plot.contourf(ax=axo,transform=ccrs.PlateCarree(),cmap=coloro,
                        levels=np.arange(o_loval,o_hival+o_spval,o_spval),
                        add_colorbar=True,add_labels=False,extend=o_extnd,
                        cbar_kwargs={'orientation':cbar_orient,'shrink':cbar_shrk,'pad':cbar_pad,
                        'spacing':'proportional','ticks':MultipleLocator(o_tkstd),'label':o_units})
    elif cntr_type == 'RasterFill':
     d18Op_c[tag,:,:].plot(ax=axo,transform=ccrs.PlateCarree(),cmap=coloro,
                        levels=np.arange(o_loval,o_hival+o_spval,o_spval),
                        add_colorbar=True,add_labels=False,extend=o_extnd,
                        cbar_kwargs={'orientation':cbar_orient,'shrink':cbar_shrk,'pad':cbar_pad,
                        'spacing':'proportional','ticks':MultipleLocator(o_tkstd),'label':o_units})

    # Overlay wind vectors if turned on
    if overlay_vec == True:
     veco = axo.quiver(d18Op_c.lon[skip1D],d18Op_c.lat[skip1D],uc[skip2D],vc[skip2D],scale_units='height',pivot='mid',
                scale=vec_scale,width=vec_wid,headlength=vec_hdl,headwidth=vec_hdw,headaxislength=vec_hal)
     qkeyo = axo.quiverkey(veco,0.8,0.2,vec_ref,str(int(vec_ref))+str(vec_units),labelpos='N',coordinates='figure',color='k')

    # Draw corresponding land or ocean tag region on plot
    if tag <= (num_landtags-1):
     draw_land_tags(numtag=tag,ax=axo,lw=taglw,major=tagmaj,minor=tagmin,color=tagcol,zorder=tagzorder)
    elif tag > (num_landtags-1):
     draw_ocean_tags(numtag=tag,ax=axo,lw=taglw,major=tagmaj,minor=tagmin,color=tagcol,zorder=tagzorder)

   #--------------------------------
   # Make PDF output file
   #--------------------------------

   # Determine if prect/d18Op are included in output file name 
   Ptext = ''
   Otext = ''
   if P == True:
    Ptext = '_prect'
   if O == True: 
    Otext = '_d18Op'

   # Setting file output name parameters
   if overlay_vec == True:
    vec_name = ''.join(('_',vec_name))
   if reg_name != '':
    reg_name = ''.join(('_',reg_name))
   if extra_name != '':
    extra_name = ''.join(('_',extra_name))

   # Save multi image PDF file and close all figures when finished
   save_multi_image('./'+str(folderpath)+'/'+str(case)+'_'+str(season)+'_watertagged'+Ptext+Otext+      \
                    vec_name+reg_name+extra_name+str(filesuf))
   plt.close('all')
