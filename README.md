# ESMplot

Welcome to ESMplot: the Earth System Model plotting package on Python! <br/>
Author: Alex Thompson (ajthompson@wustl.edu) <br/>
Latest Update: 2023-09-06 <br/>
<br/>
This is the beta version of ESMplot, a Python package designed for flexible visualization of Earth system model netCDF output. 

# Installation
For current prototype version, download entire directory 'ESMplot' (size: ~55MB) from GitHub and place in your own working directory.

# Conda environment
Use /conda_envs/environment_ESMplot_py310.yml to create a conda environment with python 3.10 to run this package.<br/>

With 'environment_ESMplot_py310.yml' in your working directory, type "conda env create -f environment_ESMplot_py310.yml"<br/>

Then type "conda activate ESMplot_py3.10" to activate the conda environment.<br/>

NOTE: ESMplot currently requires Python version 3.10 to run properly.<br/>

# Table of Contents

============================== <br/>
directory **conda_envs** <br/>
============================== <br/>

Contains 'environment_ESMplot_py310.yml' file for creating a conda environment from which to run ESMplot

============================== <br/>
directory **examples** <br/>
============================== <br/>

Example scripts that show how to make plots with ESMplot. Included are calculate_seasavg.py (calculate seasonal average), calculate_seascyc.py (calculate seasonal cycle), and calculate_watertags.py (calculate water tagging results).

============================== <br/>
directory **ESMplot** <br/>
============================== <br/>

Directory storing source code for ESMplot.

============================== <br/>
directory **pdfs** <br/>
============================== <br/>

Directory used for storing output files. 
