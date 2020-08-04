#!/usr/bin/env python

############################################################################
#
#  This file contains different functions used for the processing of Geodesy
#  data.
#
#  written by Florian Petersen and Katrin Hannemann
#
#  16 July 2020
#  
#  note that all *csv files have to be unique!
#
#  requires - pandas version 0.18.0 or larger!
#           - Python 3.X
#
#  Changelog:
#
#  21. Oct 2016
#  + use global variable GMT_DATEFORMAT
#  + change local variable dateformat in global variable IN_DATEFORMAT
#  + add documentation strings to all functions
#  + write function create_df() to create pandas.DataFrames and save data to
#       file
#  + include function create_df() in read()
#  24. Oct 2016
#  + include saving of baseline data to file and returning of baseline data
#       to read()
#  + change read() to use pandas.csv() to read from file
#  25. Oct 2016
#  + include function extract_df() which extracts new pandas.DataFrame from
#       existing pandas.DataFrame
#  + change read() to use extract_df() for extraction of baseline,
#       temperature, pressure, inclinometer, sound speed, battery and page
#       data
#  + include endtime in read() and make start and endtime optional arguments
#  26. Oct 2016
#  + make arguments index_list, column_list, label_list, dtype and date_pos
#       of function extract_df() optional to allow usage of function for
#       data slicing
#  + include function change2dateindex() which is used in extract_df()
#  + include function change_dtype() which is used in extract_df()
#  + include use of pandas.DataFrame in sort_bsl()
#  27. Oct 2016
#  + change outlier handling in sort_bsl() to work with pandas.DataFrame
#  + change function calc_baseline() to calc_baseline_amean() and include
#       new function calc_baseline_hmean() which uses the harmonic mean for
#       the calculation of the average velocity instead of the arithmetic
#       mean
#  + change function calc_baseline_theo() to calc_baseline_theo_amean() and
#       include new function calc_baseline_theo_hmean() which uses the
#       harmonic mean for the calculation of the average velocity instead of
#       the arithmetic mean
#  28. Oct 2016
#  + checked transformations of pressure used in the calculation of salinity
#	and sound speed and added references to function description
#  + for Leroy formula use pressure to depth conversion by Leroy and
#       Parthiot (1997)
#  3. Nov 2016
#  + change functions ata() and read_bsl() to use pandas.read_csv()
#  + change function range_SoundVelocity() to use pandas.DataFrame
#  4. Nov 2016
#  + change functions sal_wilson() and sal_medwin() to use pandas.DataFrame
#  8. Nov 2016
#  + include function search_df() to associate specific measurements to a
#       given measurement in time (e.g. prs to ssp)
#  + change function sort_bsl() to use function search_df()
#  + function to read in meta data like starttime etc.
#  21. Nov 2016
#  + include function lin_reg_0_intercept() to estimate slope of baseline data
#  24. Feb 2017
#  + include endtime in read_data() and start and endtime in read_bsl()
#  06. Mar 2017
#  + include function calc_bsl_baz() to calculate backazimuth of baselines
#      for either single ID pair or list of ID pairs and with option to return
#      either pandas.DataFrame() or simple list (single value in case of single
#      ID_pair)
#  20. Mar 2017
#  + include function calc_mov_average_std() which calculates moving average
#      and standard deviation for given time span and column list
#  + include function calc_fc_data() which calculates fractional change values
#      for given column list
#  + include function create_df4compare() which creates pandas.DataFrame() for
#      specific column of DataFrames in a list of DataFrames and performs a
#      linear interpolation in time to compare values of different DataFrames
#  + include function stationpair2nameext() which creates list with
#      extensions from list with ID pairs
#  + include function calc_diff_rel_col() which calculates relative changes
#      for columns to first moving average value and median of column
#  + include function stationlist2nameext() which creates list with extensions
#      from list with IDs
#
#  TO-DO LIST
#  + include function to use meta data?
#  + inculde automatice tide calculation with tpxo
#
############################################################################


################################Functions ######################################

import glob # Unix style pathname pattern expansion
import re # Regular expression operations
import numpy as np # fundamental package for scientific computing
from obspy.core import UTCDateTime # framework for processing seismological data
import datetime # Basic date and time types
import scipy.io # ecosystem of open-source software for mathematics, science, and engineering
import pandas as pd # open source, BSD-licensed library providing high-performance, easy-to-use data structures and data analysis tools
import matplotlib.pyplot as plt #
from obspy.geodetics.base import gps2dist_azimuth

import multiprocessing

global GMT_DATEFORMAT # Output date format
global IN_DATEFORMAT # Input date format
global PROJECTS # GeoSEA projects

GMT_DATEFORMAT = '%Y-%m-%dT%H:%M'
IN_DATEFORMAT = '%Y/%m/%d %H:%M:%S'

# MAR = MARSITE
# CHI = GeoSEA
# ETN = MARGOMET
PROJECTS = {'MAR' : '2014-11-16 00:00:00', 'CHI' : '2015-12-14 00:00:00', 'ETN' : '2016-04-15 00:00:00'}


#-------------------------------------------------------------------------------
#       ????
#-------------------------------------------------------------------------------

def corr_vert_bsl(ID, vertical_offset):

    H_new = []
    n = 0
    for id1 in ID:
        for id2 in ID:
            if id1 != id2:
                print(id1, '->', id2)
                prs1=(read_data(id1,'prs')-100)/10
                prs2=(read_data(id2,'prs')-100)/10
                G = prs2.mean()-prs1.mean()

                #print prs1.mean(), ' ', prs2.mean()

                print(G.values)

                H = read_bsl(id1,id2).bsl.rolling(freq='7d',window=1).median().dropna()

                first = H.first('1s').values
                last = H.last('1s').values

                off = (last - first)

				#print first

                alpha = np.arcsin(G.values/last)

				#print alpha

                A = np.cos(alpha)*last

                print(vertical_offset[n]/100)

                G_new = G.values - (vertical_offset[n]/100)

                print(G_new)

                H_new = np.sqrt(A**2 + G_new**2)

                if H_new == first:
                    H_corr = last - first
                else:
                    H_corr = H_new - first

                print(' ')
                print(off*100)
                print(H_corr*100)
                print(' ')
                print(H_corr*100-off*100)

                n = n+1

    return(H_corr)
# end def corr_vert_bsl(ID, vertical_offset):

def stationlist2nameext(IDs,ID_list):
    """Creates 1D list with extensions formed by index of IDs in ID_list.

    It needs:
    IDs ... list with single IDs which can include multiples
    ID_list ... list with single, unique IDs

    It returns:
    name_ext ... list with extensions '_<ID>'
    """

    name_ext = []
    for i in range(len(IDs)):
        st = str(ID_list.index(IDs[i])+1)
        name_ext.append('_'+st)
    # end for i in range(len(IDs))
    return(name_ext)
# end def stationlist2nameext(IDs,ID_list):

def stationpair2nameext(ID_pair,ID_list):
    """Creates 1D list with extensions for ID pairs formed by index of IDs in ID_list.

    It needs:
    ID_pair ... list with ID pairs (single entry: [ID1,ID2])
    ID_list ... list with single, unique IDs

    It returns:
    name_ext ... list with extensions '_<ID1>_<ID2>'
    """

    name_ext = []
    for i in range(len(ID_pair)):
        st1 = str(ID_list.index(ID_pair[i][0])+1)
        st2 = str(ID_list.index(ID_pair[i][1])+1)
        name_ext.append('_'+st1+'_'+st2)
    # end for i in range(len(ID_pair))
    return(name_ext)
# end def stationpair2nameext(ID_pair,ID_list):





############################### End of file ################################
