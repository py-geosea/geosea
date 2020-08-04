#
# Dependencies
#

sciPy: miscellaneous statistical functions

matplotlib: for plotting

pandas: data structures and data analysis tool

numPy: 1.7.1 or higher

obspy: for date and time types

multiprocessing



# General Folder Structure

RAW/    -> csv. files 

PYT/ -> python scripts

DAT/   -> output directory of converted station data and baselines

WOR/ -> working directory

PLT/  -> output directory of plots

# Import of GeoSEA Module

import geosea as geo

# Calculate all Possible Horizontal and Vertical Baselines in Network

geo.proc_bsl(*salinity*,*latitude*,*minmax*)

# Further Processing functions



