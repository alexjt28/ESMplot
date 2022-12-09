#######################################################################################################
#
# These functions are generic functions to be used in other plotting functions 
#
#######################################################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
import geocat.viz.util as gv
import cartopy.crs as ccrs

#######################################################################################################
# Function to save multiple images within the same PDF file as separate pages 
#######################################################################################################

def save_multi_image(filename):
   pp = PdfPages(filename)
   fig_nums = plt.get_fignums()
   figs = [plt.figure(n) for n in fig_nums]
   for fig in figs:
       fig.savefig(pp, format='pdf',bbox_inches='tight')
   pp.close()

   ''' Function to save multiple images within the same PDF file as separate pages. Appends as many
   pages as are open in the file. This function was taken from the following source:
   https://www.tutorialspoint.com/saving-multiple-figures-to-one-pdf-file-in-matplotlib'''

#######################################################################################################
# Function to plot map ticks and coordinate labels 
#######################################################################################################

def map_ticks_and_labels(ax,LatMin,LatMax,LonMin,LonMax,lattksp,lontksp,xmin_mj,ymin_mj,tkmajlg,tkminlg,
                         glalpha,glcolor,lonlblsp,latlblsp,xpads,ypads,tklblsz,
                         proj=ccrs.PlateCarree(),top=False,bot=True,left=True,right=True):

   ''' Function to plot map ticks and labels for lat/lon coordinates. Uses geocat.viz and ax.gridlines
   functionality to allow for labels on any side of map plot. See plotting functions for more details
   regarding parameters. Parameters specific to this function are 'top', 'bot', 'left', and 'right',
   which are boolean determinators for whether labels will appear on that side of the map plot.''' 

   # Sets lat/lon limits for ticks, does not draw tick labels
   gv.set_axes_limits_and_ticks(ax,xlim=(LonMin,LonMax),ylim=(LatMin,LatMax),
                                   xticks=np.linspace(-180.,180.,lontksp),
                                   yticks=np.linspace(-90.,90.,lattksp),
                                   xticklabels=[],yticklabels=[])

   # Specifications for major and minor ticks
   gv.add_major_minor_ticks(ax, x_minor_per_major=xmin_mj, y_minor_per_major=ymin_mj)

   # Length of major and minor ticks
   ax.tick_params(which='major',length=tkmajlg)
   ax.tick_params(which='minor',length=tkminlg)

   # Draws lat/lon tick labels and specifies locations to draw, allows for gridlines to be drawn too
   gl = ax.gridlines(crs=proj,draw_labels=True,alpha=glalpha,color=glcolor,
                     xlim=(LonMin,LonMax),ylim=(LatMin,LatMax))
   gl.xlocator = MultipleLocator(lonlblsp)
   gl.ylocator = MultipleLocator(latlblsp)
   gl.top_labels = top
   gl.bottom_labels = bot  
   gl.left_labels = left
   gl.right_labels = right
   gl.xpadding = xpads
   gl.ypadding = ypads
   gl.xlabel_style = {'size': tklblsz}
   gl.ylabel_style = {'size': tklblsz}

