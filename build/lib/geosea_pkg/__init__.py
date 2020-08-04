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
import seawater as sw

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

def create_df(station,suffix,date,label_list,data_list):
    """Creates pandas.DataFrame from lists and writes data to file.

    It needs:
    station ... beacon ID
    suffix ... file suffix to specify data type
    date ... 1-dim list with times of measurement
    label_list ... 1-dim list with label(s) for data to be saved
    data_list ... 2-dim list with data to be saved (first dimension must
        correspond to number of labels, second dimension to length of date
        (if just one label, data_list is a 1-dim list)

    It returns:
    df ... pandas.DataFrame with stored data

    It additionally saves data to file.
    """

    # check if one or more data fields
    try:
        n = len(data_list[0])
    except TypeError:
	# just one label (m=1)
        n = len(data_list)
        m = 1
    else:
        # multiple labels
        m = len(data_list)
    # end try:

    if m == 1:
        # single label
        # dict for creating pandas data frame
        dict_data = {'date': date, label_list : data_list}
        # create list of labels for pandas.DataFrame
        column_list = ['date', label_list]
    else:
        # multiple label
        # dict for creating pandas data frame
        dict_data = {}
        dict_data['date'] = date
        # create list of labels for pandas.DataFrame
        column_list = ['date']
        # prepend 'date' to list
        for i in range(m):
            dict_data[label_list[i]] = data_list[i]
            column_list.append(label_list[i])
        # end for i in range(m):

    # create pandas.DataFrame
    df = pd.DataFrame(dict_data, columns = column_list)
    # transform date
    df['date'] = pd.to_datetime(df['date'])
    # set date as index of data frame
    df.index = df['date']
    # delete original date field
    del df['date']
    # drops entries which include a NaN, needs to be stored to df otherwise no
    # permanent change
    df = df.dropna()
    # writes data to file
    df.to_csv('../DATA/' + str(station) +'-'+ suffix+'.dat',sep='\t', header=False, date_format=GMT_DATEFORMAT)

    return(df)
# end def create_df(station,suffix,date,label_list,data_list)

def change2dateindex(old_df,date_pos):
    """Changes index of pandas.DataFrame to date.

    It needs:
    old_df ... pandas.DataFrame with at least one date column
    date_pos ... integer which gives position of date column in label_list

    It returns:
    df ... pandas.DataFrame with date as index
    """
    # copy old_df to df
    df = old_df.copy()

    # check if df.columns[date_pos] exists
    try:
        dummy = df.columns[date_pos]
    except IndexError:
        print("Data frame has no column {0:d}! No change of index!".format(date_pos))
    else:
        # transform date
        df[df.columns[date_pos]] = pd.to_datetime(df[df.columns[date_pos]])
        # set date as index of data frame
        df.index = df[df.columns[date_pos]]
        # delete original date field
        del df[df.columns[date_pos]]
    # end try

    return(df)
# end def change2dateindex(old_df,label_list,date_pos):

def change_dtype(old_df,dtype):
    """Changes datatype of columns in pandas.DataFrame.

    It needs:
    old_df ... pandas.DataFrame
    dtype ... 1-dim list of datatypes for columns in old_df (d - integer,
        f - float, s - string)

    It returns:
    df ... pandas.DataFrame with changed datatypes
    """

    # copy old_df to df
    df = old_df.copy()

    # check if dtype and df.columns have same length
    if len(dtype) == len(df.columns):
        # change datatype of columns
        for i in range(len(df.columns)):
            if dtype[i].upper() == 'D':
                # integer
                df[df.columns[i]] = df[df.columns[i]].astype(int)
            elif dtype[i].upper() == 'F':
                # float
                df[df.columns[i]] = df[df.columns[i]].astype(float)
            else:
                # if str, pass because type will remain object
                pass
            # end if dtype[i].upper() == 'D':
        # end for i in range(label_list):
    else:
        print("Length of dtype ({0:d}) and number of columns in data frame ({1:d}) do not match!\n Datatypes will not be changed!".format(len(dtype),len(df.columns)))
    # end if len(dtype) == len(df.columns):

    return(df)
# end def change_dtype(old_df,label_list,dtype)

def extract_df(old_df,index_list=None,column_list=None,label_list=None,dtype=None,date_pos=None):
    """Extracts pandas.DataFrame from old_df and sets datatype of columns.

    It needs:
    old_df ... pandas.DataFrame from which new pandas.DataFrame should be
        extracted
    index_list (optional) ... 1-dim list with index of old_df for selection
        of data, if not given all rows are chosen
    column_list (optional) ... 1-dim list with labels of columns to be
        extracted from old_df, if not given all columns are chosen
    label_list (optional) ... 1-dim list with label(s) for data to be saved
        (the number of entries and the order of column_list and label_list
        have to match)
    dtype (optional) ... 1-dim list of datatypes for chosen columns except
        date column (d - integer, f - float, s - string)
    date_pos (optional) ... integer which gives position of date column in
        column_list or label_list (or old_df if all columns are chosen) and
        is used to change the index to date

    It returns:
    df ... extracted pandas.DataFrame with stored data

    If a given index_list or column_list does not exist an error message
    will be printed and an empty pandas.DataFrame will be returned.
    """

    try:
        if index_list is not None and column_list is not None:
            # extract new pandas.DataFrame with index and column_list
            df = old_df.loc[index_list,column_list]
        elif index_list is None and column_list is not None:
            # extract new pandas.DataFrame with column_list (all rows)
            df = old_df.loc[:,column_list]
        elif index_list is not None and column_list is None:
            # extract new pandas.DataFrame with index_list (all columns)
            df = old_df.loc[index_list]
        elif index_list is None and column_list is None:
            # copy old_df to df, a copy is necessary otherwise old_df will be
            # changed every time df is changed
            df = old_df.copy()
        # end index_list is not None and column_list is not None:

    except KeyError:
        # print error message and create empty pandas.DataFrame
        print("\nGiven pandas.DataFrame does not contain index {0} and/ or column {1}!".format(index_list,column_list))
        print("Return empty pandas.DataFrame!\n")

        df = pd.DataFrame()
    else:
        if type(df) == pd.core.series.Series:
            if not (df.index == column_list).all:
                # if just one column has been extracted from old_df, df is a
                # series not a dataframe -> convert Series to DataFrame
                df = df.to_frame
            else:
                # if just one row has been extracted from old_df, df is a
                # series not a dataframe -> convert Series to DataFrame
                dic = {}
                for i in range(len(column_list)):
                   dic[column_list[i]] = df[column_list[i]]
                # end for i in range(len(column_list)):
                df = pd.DataFrame(dic,index=[0])
            # end if not (df.index == column_list).all:
        # end type(df) == pd.core.series.Series :
        if label_list is not None and not df.empty:
            if len(label_list) == len(df.columns):
                # change labels of columns in df
                df.columns = label_list
            else:
                print("Length of label_list ({0:d}) and number of columns in data frame ({1:d}) do not match!\n Labels will not be changed!".format(len(label_list),len(df.columns)))
            # end if len(label_list) == len(df.columns):
        # end if label_list is not None:

        if date_pos is not None and not df.empty:
            df = change2dateindex(df,date_pos)
        # end if date_pos is not None:

        if dtype is not None and not df.empty:
            df = change_dtype(df,dtype)
        # end if dtype is not None:

        # drops entries which include a NaN, needs to be stored to df
        # otherwise no permanent change
        df = df.dropna(how='all')
    # end try:

    return(df)
# end def extract_df(old_df,index,column_list,label_list,date_pos):

def read_ID(pathname=None):
    """ Reads the station Ids from *csv files.
        
        Note that the *csv files have to be unique for each station!
        
        pathname (optional) If the CSV file are not in the standard folder structre
        
        Returns:
        
        ID ... an 1-dim list with station IDs
        """
    
    if pathname is None:
        pathname = '../RAW/'
    # end if pathname is None:
    
    # returns list with file names which matches given wildcards (like ls)
        ifiles = glob.glob(pathname + 'Data_*_*_*.csv')
    
    ID = []
    for filename in ifiles:
        # extract beacon ID from file name
        filename_split = filename.split('_', 3)
        ID.append(filename_split[2])
    # end for filename in ifiles:
    
    # create unique ID list
    ID = list(set(ID))
    ID.sort()

    return(ID)

def read(starttime=None, endtime=None, pathname=None, writefile=True):
    """ Reads data from *csv files.

    Note that the *csv files have to be unique for each station!
    It needs:
    starttime (optional) ... no measurement before this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    endtime (optional) ... no measurement after this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    pathname (optional) ... location of input files (default ../RAW/)
    writefile (optional) ... if True files containing all read-in parameters
        will be created in current directory, if False data will just be
        returned (default True)

    It returns:
    ID ... an 1-dim list with station IDs
    st_series ... an 1-dim list with pandas.DataFrame with columns:
        temperature ('hrt'), pressure ('prs'), sound speed ('ssp'), temperature
        from pressure ('tpr'), inclinometer data ('pitch','roll'), battery ('bat','vlt')
        and pages ('pag') with corresponding times of measurement for each beacon
        (same order as items in ID)
    bsl_series ... an 1-dim list with pandas.DataFrame with baseline
        measurements: ID of other station ('range_ID'), traveltime ('range')
        and turn around time ('TAT') with corresponding times of measurement
        for each beacon (same order as items in ID)

    It further writes human readable files for pressure, inclinometer
    data, battery, and pages, respectively.
    """

    if pathname is None:
        pathname = '../RAW/'
    # end if pathname is None:

    # returns list with file names which matches given wildcards (like ls)
    ifiles = glob.glob(pathname + 'Data_*_*_*.csv')

    ID = []
    for filename in ifiles:
	# extract beacon ID from file name
        filename_split = filename.split('_', 3)
        ID.append(filename_split[2])
    # end for filename in ifiles:

    # create unique ID list
    ID = list(set(ID))
    ID.sort()
#-------------------------------------------------------------------------------
#       Open and merge all Raw Files
#-------------------------------------------------------------------------------
    st_series = []
    bsl_series = []

    # pre-define column names (needed because of different number of columns
    # per row in input files)
    my_cols = ["A","B","C","D","E","F","G","H","I","J"]
    
    print('-------------------------------------------------------------------------------\n')
    print('GeoSEA Python Module  v1.21   20 July 2020\n')
    print('GEOMAR Helmholtz Centre for Ocean Research Kiel')
    print('-------------------------------------------------------------------------------\n\n')
    
    for j,station in enumerate(ID):

        print('\nData Processing for Station: ' + station)
        print('-------------------------------------------------------------------------------')
        print('Open Files:')

        # create empty pandas.DataFrame for storing of data from file
        all_data = pd.DataFrame()
        for i,data in enumerate(ifiles):
            stationname = data.split('_', 3)
            if station in stationname:
                print(data)
                # reads data from csv-file using pandas
                # my_cols ... pre-defined column names
                # skiprow=13 ... skips first 13 rows
                # index_col=0 ... uses first column as index
                # final dataframe has columns 'B'-'J' (0-9) and index column
                # with ['PAG','BSL',...]
                curfile = pd.read_csv(data,names=my_cols,skiprows=13,index_col=0,low_memory=False)
                # append curfile to DataFrame, needs to be stored to all_data
                # otherwise no permanent change
                all_data = all_data.append(curfile)
        #print curfile
            # end if station in stationname:
        print('   ')
        # end for i,data in enumerate(ifiles):

        # remove duplicates, again has to be stored to all_data otherwise no
        # permanent change
        all_data = all_data.drop_duplicates()
        
        # transform date column 'B' to time format used by pandas
        all_data['B'] = pd.to_datetime(all_data['B'])
        
        if starttime is not None:
            # remove all entries before starttime
            all_data = all_data.loc[all_data['B'] >= starttime]
        # end if starttime is not None:

        if endtime is not None:
            # remove all entries after endtime
            all_data = all_data.loc[all_data['B'] <= endtime]
        # end if endtime is not None:
    
        ######## Sort Files to Sensor

        # position of date in following column_lists
        date_pos = 0

        sv_fr = 0
#-------------------------------------------------------------------------------
#       Travel Time measurement
#-------------------------------------------------------------------------------
        index = 'BSL'
        column_list = ['B','F','G','H']
        # columns contain date, ID of other station, traveltime measurement in
        # milliseconds and turn around time in milliseconds
        label_list = ['date','range_ID','range','TAT']
        # data types for columns except date
        dtype = ['d','f','f']
        df_bsl = extract_df(all_data,index,column_list,label_list,dtype,date_pos)
        
        if writefile:
            # writes data to file
            df_bsl.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            #df_bsl = read_data(str(station),'BSL')
            #df_bsl = df_bsl.reset_index().drop_duplicates(subset='date').set_index('date')
            #df_bsl.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
        # end if writefile:
        # df_bsl is not written to a file because first needs to be sorted
        
#-------------------------------------------------------------------------------
#       Sound speed and temperature for Fetch Stations
#-------------------------------------------------------------------------------
        if 'SVT' in all_data.index:
            index = 'SVT'
            print('SVT - Sound Speed and Temperature Sensor !')
            column_list = ['B','F']
            # columns contain date and sound speed measurement in metres per second
            label_list = ['date','hrt']
            # data types for columns except date
            dtype = ['f','f']
            df_svt = extract_df(all_data,index,column_list,label_list,dtype,date_pos)
        
            # removes sound speed measurements which are not in water
            #df_svt = df_svt.loc[df_svt['SSP']!=9996.]
            index = 'HRT'
            if writefile:
                # writes data to file
                df_svt.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
                df_svt = read_data(str(station),index)
                df_svt = df_svt.reset_index().drop_duplicates(subset='date').set_index('date')
                df_svt.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            # end if writefile:
        
#-------------------------------------------------------------------------------
#       Sound Speed
#-------------------------------------------------------------------------------
        if 'SSP' in all_data.index:
            index = 'SSP'
            column_list = ['B','E']
            # columns contain date and sound speed measurement in metres per second
            label_list = ['date','ssp']
            # data types for columns except date
            dtype = ['f']
            df_ssp = extract_df(all_data,index,column_list,label_list,dtype,date_pos)

            # removes sound speed measurements which are not in water
            df_ssp = df_ssp.loc[df_ssp['ssp']!=9996.]

            if writefile:
                # writes data to file
                df_ssp.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
                df_ssp = read_data(str(station),'SSP')
                df_ssp = df_ssp.reset_index().drop_duplicates(subset='date').set_index('date')
                df_ssp.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            # end if writefile:
            
#-------------------------------------------------------------------------------
#       Temperature
#-------------------------------------------------------------------------------
        if 'TMP' in all_data.index:
            index = 'TMP'
            # columns contain date and temperature in degree Celsius
            label_list = ['date','tmp']
            df_tp = extract_df(all_data,index,column_list,label_list,dtype,date_pos)
        
        if 'HRT' in all_data.index:
            index = 'HRT'
            column_list = ['B','E']
            # same label_list as 'TMP'
            label_list = ['date','hrt']
            df_hrt = extract_df(all_data,index,column_list,label_list,dtype,date_pos)

            #if index == 'TMP':# concatenat both temperature dataframes to one
            #df_hrt = pd.concat([df_tmp,df_hrt])

            if writefile:
            # writes data to file
                df_hrt.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
                df_hrt = read_data(str(station),'HRT')
                df_hrt = df_hrt.reset_index().drop_duplicates(subset='date').set_index('date')
                df_hrt.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            # end if writefile:

#-------------------------------------------------------------------------------
#       Pressure and Temperature data
#-------------------------------------------------------------------------------
        index = 'PRS'
        column_list = ['B','E']
        column_list2 = ['B','F']
        # columns contain date and pressure in kPa
        label_list = ['date','prs']
        label_list2 = ['date','tpr']
        df_tpr = extract_df(all_data,index,column_list2,label_list2,dtype,date_pos)
        df_prs = extract_df(all_data,index,column_list,label_list,dtype,date_pos)
        dtype = ['f']
        if writefile:
            # writes data to file
            df_prs.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            df_prs = read_data(str(station),'PRS')
            df_prs = df_prs.reset_index().drop_duplicates(subset='date').set_index('date')
            df_prs.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            # writes temperature from pressure sensor to file
            df_tpr.to_csv('../DATA/' + str(station) +'-'+ 'TPR'+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            df_tpr = read_data(str(station),'TPR')
            df_tpr = df_tpr.reset_index().drop_duplicates(subset='date').set_index('date')
            df_tpr.to_csv('../DATA/' + str(station) +'-'+ 'TPR'+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            # end if writefile:

#-------------------------------------------------------------------------------
#       Recorded pages in Bytes
#-------------------------------------------------------------------------------
        index = 'PAG'
        column_list = ['B','E']
        # columns contain date and page number
        label_list = ['date','pag']
        #  data types for columns except date
        dtype = ['d']
        df_pag = extract_df(all_data,index,column_list,label_list,dtype,date_pos)

        if writefile:
            # writes data to file
            df_pag.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            df_pag = read_data(str(station),'PAG')
            df_pag = df_pag.reset_index().drop_duplicates(subset='date').set_index('date')
            df_pag.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
        # end if writefile:

        # tranform page numbers to Bytes
        df_pag['size'] = df_pag['pag']*512/1000

        # total size of downloaded data in kB last entry in column 'pag'
        pag_size = df_pag['size'].iloc[-1]/1024

#-------------------------------------------------------------------------------
#       Battery Power
#-------------------------------------------------------------------------------
        index = 'BAT'
        column_list = ['B','E','F']
        # columns contain date, battery consumption in per cent and voltage in
        # volt
        label_list = ['date','bat','vlt']
        #  data types for columns except date
        dtype = ['d','f']
        df_bat = extract_df(all_data,index,column_list,label_list,dtype,date_pos)

        if writefile:
            # writes data to file
            df_bat.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            df_bat = read_data(str(station),'BAT')
            df_bat = df_bat.reset_index().drop_duplicates(subset='date').set_index('date')
            df_bat.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
        # end if writefile:

#-------------------------------------------------------------------------------
#       Inclinometer
#-------------------------------------------------------------------------------
        index = 'INC'
        # columns contain date, pitch and roll in radians
        label_list = ['date','pitch','roll']
        #  data types for columns except date
        dtype = ['f','f']
        df_inc = extract_df(all_data,index,column_list,label_list,dtype,date_pos)

        # transform radians to degrees
        df_inc['pitch'] = df_inc['pitch']*180/np.pi
        df_inc['roll'] = df_inc['roll']*180/np.pi

        if writefile:
            # writes data to file
            df_inc.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
            df_inc = read_data(str(station),'INC')
            df_inc = df_inc.reset_index().drop_duplicates(subset='date').set_index('date')
            df_inc.to_csv('../DATA/' + str(station) +'-'+ index+'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
        # end writefile:

        # Standard output
        print('Found: ' + str(len(df_bsl)) + '\t Baseline Records')
        print('Found: ' + str(len(df_prs)) + '\t Pressure Records')
        if sv_fr == 1:
            print('Found: ' + str(len(df_svt)) + '\t Sound Speed Records')
            print('Found: ' + str(len(df_svt)) + '\t Temperature Records')
        else:
            print('Found: ' + str(len(df_ssp)) + '\t Sound Speed Records')
            print('Found: ' + str(len(df_hrt)) + '\t HiRes Temperature Records')
        print('Found: ' + str(len(df_inc)) + '\t Inclination Records')
        print('Found: ' + str(len(df_bat)) + '\t Battery Records')
        print('Found: ' + str(pag_size) + '\t MB Data')


	# concatenate pandas data formats in one data format for temperature,
        # pressure, sound speed, inclinometer, battery, and pages
        if sv_fr == 1:
            df = pd.concat([df_svt, df_prs, df_inc, df_bat, df_pag], axis=1)
        else:
            df = pd.concat([df_ssp, df_prs, df_hrt, df_tpr, df_inc, df_bat, df_pag], axis=1)
	# append this to data formats of other stations
        st_series.append(df)

        # baseline data not included in pandas.concat as it holds multiple
        # entries per day for different baselines
        bsl_series.append(df_bsl)
    # end for j,station in enumerate(ID):

    if not writefile:
        print('\n')
        print('Data has not been stored in files!')
    # end if not writefile:


    return(ID,st_series,bsl_series)
# end def read(starttime, pathname=None):

def baseline_calc_amean ( SSP1, SSP2, BSL_range, TAT):
    """ Calculates measured baseline lengths (arithmetic mean sound speed).

    It needs:
    SSP1 ... sound speed at beacon 1 in metres per second
    SSP2 ... sound speed at beacon 2 in metres per second
    BSL_range ... measured traveltime in milliseconds
    TAT ... turn around time in milliseconds

    It returns:
    baseline length in metres
    """

    return(((SSP1+SSP2)/2)*(((BSL_range - TAT)/2)/1000))
# end def baseline_calc ( SSP1, SSP2, BSL_range, TAT):

def baseline_calc_theo_amean ( SSP1, SSP2, BSL_range):
    """ Calculates theoretical baseline length (arithmetic mean sound speed).

    It needs:
    SSP1 ... sound speed at beacon 1 in metres per second
    SSP2 ... sound speed at beacon 2 in metres per second
    BSL_range ... traveltime between beacons in milliseconds

    It returns:
    baseline length in metres
    """

    return(((SSP1+SSP2)/2)*(BSL_range))
# end def baseline_calc_theo ( SSP1, SSP2, BSL_range):

def baseline_calc_hmean ( SSP1, SSP2, BSL_range, TAT):
    """ Calculates measured baseline lengths (harmonic mean sound speed).

    It needs:
    SSP1 ... sound speed at beacon 1 in metres per second
    SSP2 ... sound speed at beacon 2 in metres per second
    BSL_range ... measured traveltime in milliseconds
    TAT ... turn around time in milliseconds

    It returns:
    baseline length in metres
    """

    return(((2*SSP1*SSP2)/(SSP2+SSP1))*(((BSL_range - TAT)/2)/1000))
# end def baseline_calc ( SSP1, SSP2, BSL_range, TAT):

def baseline_calc_theo_hmean ( SSP1, SSP2, BSL_range):
    """ Calculates theoretical baseline length (harmonic mean sound speed).

    It needs:
    SSP1 ... sound speed at beacon 1 in metres per second
    SSP2 ... sound speed at beacon 2 in metres per second
    BSL_range ... traveltime between beacons in milliseconds

    It returns:
    baseline length in metres
    """

    return(((2*SSP1*SSP2)/(SSP2+SSP1))*(BSL_range))
# end def baseline_calc_theo ( SSP1, SSP2, BSL_range):

def search_df(df,df2,minmax,def_val=None,new_column=None,column_num=None):
    """Searches closest (in time) entry of df2 in relation to df.

    It needs:
    df ... pandas.DataFrame with one or more columns and date as index
    df2 ... pandas.DataFrame with one column and date as index
    minmax ... maximum half time window length for searching for value in
        df2 relative to value in df
    def_val (optional) ... default value to be set to df2 if no match is
        found (default is NaN)
    new_column (optional) ... name of column for df_new (default is column
        name of df2)

    It returns:
    df_new ... pandas.DataFrame with two or more columns the first columns
        hold the values of df and the last holds the values extracted from
        df2
    err_count ... counter for number of failed trails for search in df2 for
        a matching value for df, if search failed value in df2 will be set
        to def_val
    """

    if column_num is None:
        column_num = 0
    
    if new_column is None:
        new_column = df2.columns[column_num]
   
    # end if new_column is None:

    if def_val is None:
        def_val = np.nan
    # end if def_val is None:

    # copy df
    df_new = df.copy()
    err_count = 0
    if not df.empty and not df2.empty:
        for i in range(len(df)):
            # set default value for value in df2
            val = def_val
            try:
                # search index of value in df2 for time of value in df with
                # minmax s tolerance
                # method = 'nearest' ensures that nearest measurement in time is
                # chosen
                ind = df2.index.get_loc(df.index[i],method='nearest',tolerance=pd.tseries.offsets.Second(minmax))
            except KeyError:
                # if date is the same in df or does not exist a
                # InvalidIndexError is raised
                # search index of value in df2 for time of value in df
                try:
                    val = df2.loc[df.index[i],df2.columns[column_num]]
                except KeyError:
                    # if date does not exist in df2 a KeyError is raised
                    err_count += 1
                # end try: (same date in df2?)
            else:
                # choose value based on index found
                val = df2.iloc[ind,column_num]
            # end try: (nearest date in df?)
            # store value in df_new
            df_new.at[df.index[i],new_column] = val

        # end for i in range(len(df)):
    # end if not df.empty and df2.empty:

    return(df_new,err_count)
# end def search_df(df,df2,minmax):

#-------------------------------------------------------------------------------
#       Replace of broken Sensor data
#-------------------------------------------------------------------------------
def replace(ID1,ID2, sensor, starttime=None, project=None):
    """ Replace incorrect Sensor data with a neighbour station
        
        It needs:
        ID1 ... ID of beacon with broken sensor
        ID2 ... ID of neighbour beacon
        
        sensor ... broken sensor
        starttime ... time before the sensor of ID1 is broken
        
        The replace data will be immediatly safed in ../DATA/
        
        The sensor data of ID2 will be interpolated with the index of ID1 by slinear method.
        The data needs to be shifted to avoid any leaps in the interpolated data.
        Check data afterwards!
        """
    
    if starttime is None:
        
        df_idata1 = read_data(ID1, str(sensor))
        df_idata2 = read_data(ID2, str(sensor))
        
        df_idata1_first = df_idata1 - df_idata1.first('1s').values.squeeze()
        df_idata2_first = df_idata2 - df_idata2.first('1s').values.squeeze()
        df_idata1_old_last = df_idata1 - df_idata1.first('1s').values.squeeze()
    
    elif project is None:
        
        df_idata1 = read_data(ID1, str(sensor), starttime=starttime)
        df_idata2 = read_data(ID2, str(sensor), starttime=starttime)
        df_idata1_old = read_data(ID1, str(sensor), endtime=starttime)
        
        df_idata1_first = df_idata1 - df_idata1.first('1s').values.squeeze()
        df_idata2_first = df_idata2 - df_idata2.first('1s').values.squeeze()
        
        df_idata1_old_last = df_idata1_old - df_idata1_old.last('1s').values.squeeze()

    else:
        
        df_idata1 = read_data(ID1, str(sensor), starttime=starttime)
        df_idata2 = read_data(ID2, str(sensor), starttime=starttime)
        
        # selection global PROJECTS starttime
        if project.upper() in PROJECTS:
            
            df_idata1_old = read_data(ID1, str(sensor),starttime=PROJECTS[project.upper()], endtime=starttime)

        # subtract the mean of both datasets
        df_idata1_first = df_idata1 - df_idata1.first('1s').values.squeeze()
        df_idata2_first = df_idata2 - df_idata2.first('1s').values.squeeze()

    df_idata1_old_last = df_idata1_old - df_idata1_old.last('1s').values.squeeze()

    # add Index of reference station to DataFrame and interpolate the sensor data in between
    df_data_newindex = df_idata2_first.reindex(df_idata1_first.index).append(df_idata2_first).sort_index()
    df_data_interpolate = df_data_newindex.interpolate(method ='slinear')
    df_odata1 = pd.DataFrame(df_data_interpolate, index = df_idata1.index)
    
    # Add back the mean of old sensor data
    df_odata1 = df_odata1 + df_idata1_old.last('1s').values.squeeze()
    
    if starttime is None:
        df_odata1.to_csv('../DATA/' + str(ID1) +'-'+ sensor.upper() +'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
    else:
        odata1_list = [df_idata1_old, df_odata1]
        df_odata1_all = pd.concat(odata1_list)
        
        # save the replaced data into the DATA directory
        df_odata1_all.to_csv('../DATA/' + str(ID1) +'-'+ sensor.upper() +'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)

    return(df_odata1_all)
    #end def replace (ID1,ID2, sensor, starttime=None, project=None):

    #-------------------------------------------------------------------------------
    #       Replace of HRT Data
    #-------------------------------------------------------------------------------

def replace_hrt(ID, starttime):
    """ Replace HRT Sensor data with the TMP data
    
        It needs:
        ID1 ... ID of beacon with broken sensor
    
        starttime
      
        The replace data will be immediatly safed in ../DATA/
    
        The sensor data of ID2 will be interpolated with the index of ID1 by slinear method.
        The data needs to be shifted to avoid any leaps in the interpolated data.
        Check data afterwards!
        """
    df_indata_hrt = read_data(ID, 'hrt', endtime=starttime)
    df_indata_tmp = read_data(ID, 'tmp', starttime=starttime)
    df_indata_tmp.columns = ['hrt']
    df_indata_hrt_last = df_indata_hrt - df_indata_hrt.last('1s').values.squeeze()
    df_indata_tmp_first = df_indata_tmp - df_indata_tmp.first('1s').values.squeeze()
    
    df_indata_hrt_concat = [df_indata_hrt_last, df_indata_tmp_first]
    df_outdata_hrt = pd.concat(df_indata_hrt_concat)
    
    df_outdata_hrt = df_outdata_hrt + df_indata_hrt.last('1s').values.squeeze()
        
    df_outdata_hrt.to_csv('../DATA/' + str(ID) +'-HRT.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
   
    return(df_outdata_hrt)

#-------------------------------------------------------------------------------
#       Complete Replace of Data
#-------------------------------------------------------------------------------

def replace_complete(ID1,ID2, sensor):
    """ Replace whole Sensor data with neighbour station
        
        It needs:
        ID1 ... ID of beacon with broken sensor
        ID2 ... ID of neighbour beacon
        
        sensor ... broken sensor
        
        The replace data will be immediatly safed in ../DATA/
        
        The sensor data of ID2 will be interpolated with the index of ID1 by slinear method.
        The data needs to be shifted to avoid any leaps in the interpolated data.
        Check data afterwards!
        """
  
    df_idata1 = read_data(ID1, str(sensor))
    df_idata2 = read_data(ID2, str(sensor))

    df_idata1_first = df_idata1
    df_idata2_first = df_idata2

    # add Index of reference station to DataFrame and interpolate the sensor data in between
    df_data_newindex = df_idata2_first.reindex(df_idata1_first.index).append(df_idata2_first).sort_index()
    df_data_interpolate = df_data_newindex.interpolate(method ='slinear')

    # Remove duplicate index in Dataframe
    df_data_interpolate = df_data_interpolate[~df_data_interpolate.index.duplicated(keep='first')]
    
    df_new = pd.DataFrame(df_data_interpolate.dropna(), index = df_idata1.index).dropna()
    
    # save the replaced data into the DATA directory
    df_new.to_csv('../DATA/' + str(ID1) +'-'+ sensor.upper() +'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
   
    return(df_new)
#end def replace_complete(ID1,ID2, sensor):

#-------------------------------------------------------------------------------
#       Calculate Baselines data
#-------------------------------------------------------------------------------

def sort_bsl(ID,bsl_all,st_series,minmax,outlier_flag=None,writefile=True):
    """Calculates baselines for all possible pairs.

    It needs:
    ID ... 1-dim array with all available beacon IDs
    bsl_all ... an 1-dim list with pandas.DataFrame with baseline
        measurements: ID of other station ('range_ID'), traveltime ('range')
        and turn around time ('TAT') with corresponding times of measurement
        for each beacon (same order as items in ID)
    ssp_all ... an 1-dim list with pandas.DataFrame with sound speed ('ssp')
        with corresponding times of measurement for each beacon (same order
        as items in ID)
    minmax ... half time window length for searching for sound speed at
        beacon 2
    outlier_flag (optional) ... if set to 1 all baselines with lengths
        +/-10m are removed
    suffix (optional) ... file suffix for saved baseline files (default - no
        suffix)
    writefile (optional) ... if True files containing all baseline parameters
        will be created in current directory, if False data will just be
        returned (default True)

    It returns:
    ID_pair ... a 2-dim list with IDs of beacon pairs
        ([[ID1,ID2],[ID1,ID3],...])
    final_bsls ... an 1-dim list with pandas.DataFrame with ID of beacon 1
        ('ID'), ID of beacon 2 ('range_ID'), calculated baseline lengths in
        metres ('bsl'), one way traveltime in seconds ('tt'), sound speed at
        beacon 1 ('ssp1') in metres per second, sound speed at beacon 2
        ('ssp2') in metres per second, measured traveltime in milliseconds
        ('range'), turn around time in milliseconds ('TAT') with
        corresponding times of measurement for each beacon pair (same order
        as list items in ID_pair)
    """

    cal_bsl_all = []
    pool = multiprocessing.Pool(processes=4)
    ID_pair = []
    st_svleroy = []
        
    for i, beacon_1 in enumerate(ID):
        for j, beacon_2 in enumerate(ID):
            if beacon_1 != beacon_2:

                print('Baseline Calculation for: ' + str(beacon_1) + ' <-> ' + str(beacon_2))
                print('-------------------------------------------------------------------------------')
                ID_pair.append([beacon_1,beacon_2])
                # create new pandas.DataFrame holding baseline measurements
                # between beacon_1 and beacon_2 which are not 0.0 milli seconds
                df_bsl = bsl_all[i].loc[(bsl_all[i]['range_ID']==int(beacon_2)) & (bsl_all[i]['range']!=0.0)]

                if not df_bsl.empty and not st_series[i].empty and not st_series[j].empty:
                    # set new column 'ID' of beacon 1
                    df_bsl['ID']=int(beacon_1)
                    # alternative version:
                    #df_bsl.loc[df_bsl.index,'ID'] = int(beacon_1)

                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'ssp1',0)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_hrt1',10)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_tpr1',11)
                    
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'prs1',1)
                    
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'hrt1',2)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'tpr1',3)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sal1',12)
                    
                    
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'ssp2',0)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_hrt2',10)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_tpr2',11)
                    
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'prs2',1)
                    
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'hrt2',2)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'tpr2',3)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sal2',12)
                    
                    # define criteria that ssp1 and ssp2 have to be not NaN
                    # for row selection allow also single sided baseline 
                    # calculation
    #---------------------------------------------------------------------------------------
    #            SSP
                    criteria_ssp = (df_bsl['ssp1']!=0.0) & (df_bsl['ssp2']!=0.0)
                    ssp2_check = (df_bsl['ssp2'] == 0.0) & (df_bsl['ssp1']!=0.0)
                    ssp1_check = (df_bsl['ssp1'] == 0.0) & (df_bsl['ssp2']!=0.0)
                    
                    df_bsl.loc[criteria_ssp,'tt'] = ((df_bsl[criteria_ssp]['range']-df_bsl[criteria_ssp]['TAT'])/2)/1000
                    df_bsl.loc[ssp2_check,'tt'] = ((df_bsl[ssp2_check]['range']-df_bsl[ssp2_check]['TAT'])/2)/1000
                    df_bsl.loc[ssp1_check,'tt'] = ((df_bsl[ssp1_check]['range']-df_bsl[ssp1_check]['TAT'])/2)/1000
                    
                    df_bsl.loc[criteria_ssp,'bsl'] = baseline_calc_hmean(df_bsl[criteria_ssp]['ssp1'],df_bsl[criteria_ssp]['ssp2'],df_bsl[criteria_ssp]['range'],df_bsl[criteria_ssp]['TAT'])
                    
                    df_bsl.loc[ssp1_check,'bsl'] = baseline_calc_hmean(df_bsl[ssp1_check]['ssp2'],df_bsl[ssp1_check]['ssp2'],df_bsl[ssp1_check]['range'],df_bsl[ssp1_check]['TAT'])

                    df_bsl.loc[ssp2_check,'bsl'] = baseline_calc_hmean(df_bsl[ssp2_check]['ssp1'],df_bsl[ssp2_check]['ssp1'],df_bsl[ssp2_check]['range'],df_bsl[ssp2_check]['TAT'])
                    
                    bsl_sucess = len(df_bsl.loc[pd.notnull(df_bsl['bsl'])])
    #---------------------------------------------------------------------------------------
    #            SV_HRT

                    criteria_sv_hrt = (df_bsl['sv_hrt1']!=0.0) & (df_bsl['sv_hrt2']!=0.0)
                    sv_hrt2_check = (df_bsl['sv_hrt2'] == 0.0) & (df_bsl['sv_hrt1']!=0.0)
                    sv_hrt1_check = (df_bsl['sv_hrt1'] == 0.0) & (df_bsl['sv_hrt2']!=0.0)
                    
                    df_bsl.loc[criteria_sv_hrt,'tt'] = ((df_bsl[criteria_sv_hrt]['range']-df_bsl[criteria_sv_hrt]['TAT'])/2)/1000
                    df_bsl.loc[sv_hrt2_check,'tt'] = ((df_bsl[sv_hrt2_check]['range']-df_bsl[sv_hrt2_check]['TAT'])/2)/1000
                    df_bsl.loc[sv_hrt1_check,'tt'] = ((df_bsl[sv_hrt1_check]['range']-df_bsl[sv_hrt1_check]['TAT'])/2)/1000
                    
                    df_bsl.loc[criteria_sv_hrt,'bsl_hrt'] = baseline_calc_hmean(df_bsl[criteria_sv_hrt]['sv_hrt1'],df_bsl[criteria_sv_hrt]['sv_hrt2'],df_bsl[criteria_sv_hrt]['range'],df_bsl[criteria_sv_hrt]['TAT'])
                    
                    df_bsl.loc[sv_hrt1_check,'bsl_hrt'] = baseline_calc_hmean(df_bsl[sv_hrt1_check]['sv_hrt2'],df_bsl[sv_hrt1_check]['sv_hrt2'],df_bsl[sv_hrt1_check]['range'],df_bsl[sv_hrt1_check]['TAT'])

                    df_bsl.loc[sv_hrt2_check,'bsl_hrt'] = baseline_calc_hmean(df_bsl[sv_hrt2_check]['sv_hrt1'],df_bsl[sv_hrt2_check]['sv_hrt1'],df_bsl[sv_hrt2_check]['range'],df_bsl[sv_hrt2_check]['TAT'])
                    
                    bsl_hrt_sucess = len(df_bsl.loc[pd.notnull(df_bsl['bsl_hrt'])])
    #---------------------------------------------------------------------------------------
    #            SV_TPR
    
                    criteria_sv_tpr = (df_bsl['sv_tpr1']!=0.0) & (df_bsl['sv_tpr2']!=0.0)
                    sv_tpr2_check = (df_bsl['sv_tpr2'] == 0.0) & (df_bsl['sv_tpr1']!=0.0)
                    sv_tpr1_check = (df_bsl['sv_tpr1'] == 0.0) & (df_bsl['sv_tpr2']!=0.0)
                    
                    df_bsl.loc[criteria_sv_tpr,'tt'] = ((df_bsl[criteria_sv_tpr]['range']-df_bsl[criteria_sv_tpr]['TAT'])/2)/1000
                    df_bsl.loc[sv_tpr2_check,'tt'] = ((df_bsl[sv_tpr2_check]['range']-df_bsl[sv_tpr2_check]['TAT'])/2)/1000
                    df_bsl.loc[sv_tpr1_check,'tt'] = ((df_bsl[sv_tpr1_check]['range']-df_bsl[sv_tpr1_check]['TAT'])/2)/1000
                    
                    df_bsl.loc[criteria_sv_hrt,'bsl_tpr'] = baseline_calc_hmean(df_bsl[criteria_sv_tpr]['sv_tpr1'],df_bsl[criteria_sv_tpr]['sv_tpr2'],df_bsl[criteria_sv_tpr]['range'],df_bsl[criteria_sv_tpr]['TAT'])
                    
                    df_bsl.loc[sv_tpr1_check,'bsl_tpr'] = baseline_calc_hmean(df_bsl[sv_tpr1_check]['sv_tpr2'],df_bsl[sv_tpr1_check]['sv_tpr2'],df_bsl[sv_tpr1_check]['range'],df_bsl[sv_tpr1_check]['TAT'])

                    df_bsl.loc[sv_tpr2_check,'bsl_tpr'] = baseline_calc_hmean(df_bsl[sv_tpr2_check]['sv_tpr1'],df_bsl[sv_tpr2_check]['sv_tpr1'],df_bsl[sv_tpr2_check]['range'],df_bsl[sv_tpr2_check]['TAT'])
                    
                    bsl_tpr_sucess = len(df_bsl.loc[pd.notnull(df_bsl['bsl_tpr'])])
                    
                    # calculate traveltime in seconds and store in new column of
                    # df_bsl

                    # calculate baseline length

                    
                    # count entries of df_bsl for which 'bsl' is not NaN

                else:
                    bsl_sucess = 0
                #end if not df_bsl.empty:

                df_bsl = extract_df(df_bsl,column_list=['ID','range_ID','range','TAT','tt','hrt1','hrt2','prs1','prs2','tpr1','tpr2','sal1','sal2','ssp1','ssp2','bsl','sv_hrt1','sv_hrt2','bsl_hrt','sv_tpr1','sv_tpr2','bsl_tpr'])
                if writefile:
                    df_bsl.to_csv('../DATA/' + beacon_1 +'-'+ beacon_2 +'-BSL.dat', header=True, date_format=GMT_DATEFORMAT)
                    # end if suffix == None:
                # end if writefile:

                print(str(len(df_bsl)) + '\t Ranges found')
                print(str(bsl_sucess) + '\t Successfull Calculated Baselines')
                if len(df_bsl) != 0:
                    print(str(SV_1_err_count) + '\t No SV Record in -> ' + str(beacon_1))
                    print(str(SV_2_err_count) + '\t No SV Record in -> ' + str(beacon_2))
                print(' \n')

                #  append to data formats of other stations
                cal_bsl_all.append(df_bsl)


            # end if beacon_1 != beacon_2:

        # end for j, beacon_2 in enumerate(ID):

    # end for i, beacon_1 in enumerate(ID):

    final_bsls = []
    for i in range(len(cal_bsl_all)):
        if outlier_flag == 1:
            ### cut off unrealistic Ranges and Baselines ###

            # This part removes all baselines with lengths +/-10m

            print('Cut Off unrealistic Ranges and Baselines')
            # calculate mean range exculding NaN value
            mean_bsl = cal_bsl_all[i].mean(axis=0,skipna=True)['bsl']
            # keep only those baselines within mean_bsl +/-10 m
            df_bsl = cal_bsl_all[i].loc[ (cal_bsl_all[i]['bsl']>mean_bsl-10) & (cal_bsl_all[i]['bsl']<mean_bsl+10)]
            print("Pair: {0:s} <-> {1:s}".format(ID_pair[i][0],ID_pair[i][1]))
            print("{0:d} baselines from {1:d} kept.".format(len(df_bsl),len(cal_bsl_all[i])))
        else:
            df_bsl = cal_bsl_all[i]
        # end if outlier_flag == 1:

        # re-arange order of columns
        df_bsl = extract_df(df_bsl,column_list=['ID','range_ID','bsl','tt','ssp1','ssp2','hrt1','hrt2','prs1','prs2','sal1','sal2','range','TAT'])
        final_bsls.append(df_bsl)

        if writefile:
            if suffix is None:
                final_bsls[i].to_csv('../DATA/' + ID_pair[i][0] +'-'+ ID_pair[i][1] +'.dat', header=True, date_format=GMT_DATEFORMAT)
            else:
                final_bsls[i].to_csv('../DATA/' + ID_pair[i][0] +'-'+ ID_pair[i][1] +'-'+ suffix + '.dat', header=True, date_format=GMT_DATEFORMAT)
            # end if suffix is None:
        # end if writefile:
    # end for i in range(len(cal_bsl_all):

    if not writefile:
        print('\n')
        print('Data has not been stored in files!')
    # end if not writefile:


    return(ID_pair,final_bsls)
# end def sort_bsl( ... ):

#-------------------------------------------------------------------------------
#       Calculate sound speed from Baseline
#-------------------------------------------------------------------------------

def range_SoundVelocity ( ID1, ID2, BSL ):
    """ Calculates sound speed from baseline measurements at one beacon.

    It needs:
    ID1 ... ID of beacon 1
    ID2 ... ID of beacon 2
    BSL ... pandas.DataFrame with baseline data at beacon 1 (ID of beacon 2
        ('range_ID'), measured traveltime in milliseconds ('range'), turn
        around time in milliseconds ('TAT'), sound speed at beacon 1
        ('ssp1') in metres per second, sound speed at beacon 2 ('ssp2') in
        metres per second one way traveltime in seconds ('tt'), calculated
        baseline lengths in metres ('bsl') with corresponding times of
        measurement)

    It returns:
    df_ssp_cal ... pandas.DataFrame with estimated sound speeds in metres
        per second with corresponding times

    If BSL contains baseline data from more than one ID_pair or if the IDs
    in BSL do not match the given IDs (ID1 and ID2), an error message will
    be printed and an empty pandas.DataFrame will be returned.
    """

    # create empty pandas.DataFrame for sound speed estimation
    df_ssp_cal = pd.DataFrame()

    dummy = BSL.loc[:,['ID','range_ID']].drop_duplicates()
    if (len(dummy) == 1) and (dummy.iloc[0][0] == int(ID1)) and  (dummy.iloc[0][1] == int(ID2)):
        # calculate mean range exculding NaN value
        mean_bsl = BSL.mean(axis=0,skipna=True)['bsl']
        # calculate sound speed from mean baseline lengths and travel time
        df_ssp_cal['ssp'] = mean_bsl/BSL.loc[:,'tt']
    else:
        print('Given baseline DateFrame contains either more than one ID pair orgiven ID pair does not match!')
        print('Empty DataFrame will be returned!')
    # end if (...):

    return(df_ssp_cal)
#end def range_SoundVelocity ( ID1, ID2, BSL ):

#-------------------------------------------------------------------------------
#       Complete processing
#-------------------------------------------------------------------------------

def read_proc_bsl (SAL,phi,minmax,outlier_flag=None,writefile=True):
    """ Complete Baseline processing of GeoSEA Raw data.

    It needs:
    SAL ... constant salinity value
    phi ... Latitude for Leroy formular
    minmax ... half time window length for searching for sound speed at
    beacon 2

    It returns:
    bsl ... list of pandas.DataFrame with calculated Baselines

    """
    ID,st_series,bsl_series = read()
    
    st_series_leroy = []
    
    for i, id in enumerate(ID):
        st_series_leroy.append(sv_leroy(st_series[i],SAL,phi))
        
    bsl = sort_bsl(ID,bsl_series,st_series_leroy,7600,outlier_flag,writefile)
    
    return(bsl)

#################### Seawater Module #########################

#-------------------------------------------------------------------------------
#       Salinity Wilson
#-------------------------------------------------------------------------------
def sal_wilson ( HRT, SSP, PRS ):
    """ Calculates salinity at one beacon using the Wilson formula.

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    SSP ... pandas.DataFrame with measured sound speed in metres per second
        at one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon

    It returns:
    df_sal ... pandas.DataFrame with calculated salinities in parts per
        thousand
    """

    # concatenate Data.Frames of temperature, sound speed and pressure
    df = pd.concat([HRT, SSP, PRS], axis=1)

    # transform pressure from kPa to kilogram-force per square centimetres
    # 1 kPa = 0.010197266 kg_f/cm^2
    df['prs'] = df['prs']*0.010197266

    df['S_t'] = -3.31986*df['hrt'] - 2.57236*10**(-3)*df['hrt']**2 + 2.32009*10**(-4)*df['hrt']**3 - 2.20892*10**(-6)*df['hrt']**4

    df['S_p'] = -1.14663*10**(-1)*df['prs'] - 1.00488*10**(-5)*df['prs']**2 - 2.14038 * 10**(-8)*df['prs']**3 + 2.23490 * 10**(-12) * df['prs']**4

    df['S_v'] = 7.21467 * 10**(-1) * (df['ssp'] - 1448.54) - 6.52575 * 10**(-4)* ( df['ssp'] - 1448.54)**2

    df['S_tpv'] = (df['ssp'] - 1448.54) * (1.17 * 10**(-2)*df['hrt'] + 1.04088 * 10**(-4) * df['prs'] + 1.29184 * 10**(-7) * df['prs']**2 + 4.43572 * 10**(-6) * df['prs'] * df['hrt'] - 9.63357 * 10**(-8) * df['prs'] * df['hrt']**2) + df['prs'] * (-1.45 * 10**(-3) * df['hrt'] - 2.5563 * 10**(-5) * df['hrt']**2 + 3.45997 * 10**(-7) * df['hrt']**3) + df['prs']**2 * (-7.29321 * 10**(-7) * df['hrt'] + 2.05515 * 10**(-8) * df['hrt']**2) + df['prs']**3 * (-1.19869 * 10**(-10) * df['hrt'])

    df['sal'] = 35 + df['S_t'] + df['S_p'] + df['S_v'] + df['S_tpv']

    # store salinity in data frame
    df_sal = extract_df(df,column_list=['sal'],dtype=['f'])

    return(df_sal)
# end def sal_wilson ( HRT, SSP, PRS ):

#-------------------------------------------------------------------------------
#       Salinity Medwin
#-------------------------------------------------------------------------------

def sal_medwin ( HRT, SSP, PRS ):
    """ Calculates salinity at one beacon using the Medwin formula.

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    SSP ... pandas.DataFrame with measured sound speed in metres per second
        at one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon

    It returns:
    df_sal ... pandas.DataFrame with calculated salinities in parts per
        thousand
    """

    # concatenate Data.Frames of temperature, sound speed and pressure
    df = pd.concat([HRT, SSP, PRS], axis=1)

    # pressure in depth (m)??
    df['prs'] = (df['prs']-100)/10

    df['sal'] = (((1449.2 + 4.6*df['hrt'] - 0.055*df['hrt']**2 + 0.00029*df['hrt']**3 + 0.016 * df['prs'] - df['ssp'])*(-1))/(1.34 - 0.010*df['hrt']))+35

    # store salinity in data frame
    df_sal = extract_df(df,column_list=['sal'],dtype=['f'])

    return(df_sal)

# end def sal_medwin ( ... ):

#-------------------------------------------------------------------------------
#       Sound Velocity Wilson
#-------------------------------------------------------------------------------

def sv_wilson ( HRT, PRS, SAL):
    """ Calculates sound speed using the Wilson formula.

    Wilson formula: Equation 1 in  Wayne D. Wilson, 1960. Speed of Sound in
        Sea Water as a Function of Temperature, Pressure, and Salinity, J.
        Acoust. Soc. Am. 32, 641, doi: 10.1121/1.1908167

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon
    SAL ... constant salinity in parts per thousand

    It returns:
    df_ssp ... pandas.DataFrame holding the estimated sound speeds
        ('ssp_w') in metres per second at one beacon
    """

    # concatenate Data.Frames of temperature and pressure
    if isinstance(SAL, pd.DataFrame):
        df = pd.concat([HRT, PRS, SAL], axis=1)
    else:
        # set constant salinity value
        df = pd.concat([HRT, PRS ], axis=1)
        df['sal']=SAL

    # transform pressure from kPa to kilogram-force per square centimetres
    # 1 kPa = 0.010197266 kg_f/cm^2
    df['prs'] = df['prs']*0.010197266

    df['V_t'] = 4.6233 * df['hrt'] - 5.4585 * 10**(-2) * df['hrt']**2 + 2.822 * 10**(-4) * df['hrt']**3 - 5.07 * 10**(-7) * df['hrt']**4

    df['V_p'] = 1.60518 * 10**(-1) * df['prs'] + 1.0279 * 10**(-5) * df['prs']**2 + 3.451 * 10**(-9) *df['prs']**3 - 3.503 * 10**(-12) * df['prs']**4

    df['V_s'] = 1.391 * (df['sal']-35) - 7.8 * 10**(-2) * (df['sal']-35)
    df['V_stp'] = (df['sal']-35) * (-1.197 * 10**(-2) * df['hrt'] + 2.61 * 10**(-4) * df['prs'] - 1.96* 10**(-7) * df['prs']**2 - 2.09 * 10**(-6) * df['prs'] * df['hrt']) + df['prs'] *(-2.796 * 10**(-4) * df['hrt'] + 1.3302 * 10**(-5) * df['hrt']**2 - 6.644*10**(-8) * df['hrt']**3) + df['prs']**2 * (-2.391 * 10**(-7) * df['hrt']+ 9.286 * 10**(-10) * df['hrt']**2) - 1.745 * 10**(-10) * df['prs']** 3 * df['hrt']

    df['ssp'] = 1449.22 + df['V_t'] + df['V_p'] + df['V_s'] + df['V_stp']

    # store sound speed in data frame
    df_ssp = extract_df(df,column_list=['ssp'],dtype=['f'])
    df_ssp = pd.concat([ df_ssp, HRT, PRS, df['sal']], axis=1)
    return(df_ssp)
# end def sv_wilson ( HRT, PRS, SAL):

#-------------------------------------------------------------------------------
#       Sound Velocity Leroy
#-------------------------------------------------------------------------------

#def sv_leroy (HRT, PRS, SAL, phi):
def sv_leroy (st_series, SAL, phi):
    """ Calculates sound speed using the Leroy formula.

    Leroy Formula: Equation 2 in Claude C. Leroy, Stephen P. Robinson and
        Mike J. Goldsmith, 2008. A new equation for the accurate calculation
        of sound speed in all oceans, J. Acoust. Soc. Am. 124, 2774,
        doi: 10.1121/1.2988296
    Conversion of pressure into depth: Equations 3 and 4 in C.C. Leroy and
        F. Parthiot, 1998. Depth-pressure relationships in the oceans and
        seas, J. Acoust. Soc. Am. 103, 1346, doi: 10.1121/1.421275 and C.C.
        Leroy, 2007. Erratum: "Depth-pressure relationships in the oceans
        and seas" [J. Acoust. Soc. Am.103(3), 1346-1352 (1998)], J. Acoust.
        Soc. Am. 121, 2447; doi: 10.1121/1.2534201 gives depth with an
        accuracy of +/- 1 metre

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon
    SAL ... constant salinity in parts per thousand
    phi ... latitude of working area in degrees

    It returns:
    df_ssp ... pandas.DataFrame holding the estimated sound speeds
        ('ssp_l') in metres per second
    """

    # concatenate Data.Frames of temperature and pressure

    if isinstance(SAL, pd.DataFrame):
        df = pd.concat([st_series, SAL], axis=1)
    else:
    # set constant salinity value
        df = pd.concat([st_series ], axis=1)
        df['sal']=SAL

    ##### transform pressure in kPa in depth in m #####

    # transform latitude in radians
    phi_rad = phi*np.pi/180
    # international formula for gravity (eq. 4 in Leroy and Parthiot, 1998)
    g = 9.780318*(1+5.2788*10**(-3)*(np.sin(phi_rad))**2-2.36*10**(-5)*(np.sin(phi_rad))**4)

    # transformation of pressure in kPa into MPa which is used in equation 3 in
    # Leroy and Parthiot (1998)
    df['prs'] = df['prs']/1000
    # transform pressure in MPa into depth in metres
    df['z'] = (9.72659*10**2*df['prs']-2.2512*10**(-1)*df['prs']**2+2.279*10**(-4)*df['prs']**3-1.82*10**(-7)*df['prs']**4)/(g+1.092*10**(-4)*df['prs'])

    # calculate sound speed in metres per second with eq.2 in Leroy et al.(2008)
    if 'hrt' in df:
        df['svl_hrt'] = 1402.5 + 5 * df['hrt'] - 5.44 * 10**(-2)*df['hrt']**2 + 2.1 * 10**(-4) * df['hrt']**3 + 1.33 * df['sal'] - 1.23 * 10**(-2) * df['sal'] * df['hrt'] + 8.7 * 10**(-5) * df['sal']* df['hrt']**2 + 1.56 * 10**(-2) * df['z'] + 2.55 * 10**(-7) * df['z']**2 - 7.3 * 10**(-12) * df['z']**3 + 1.2 * 10**(-6) * df['z'] * (phi - 45) - 9.5 * 10**(-13) * df['hrt'] * df['z']**3 + 3 * 10**(-7) * df['hrt']**2 * df['z'] + 1.43 * 10**(-5) * df['sal'] * df['z']
    if 'tpr' in df:
        df['svl_tpr'] = 1402.5 + 5 * df['tpr'] - 5.44 * 10**(-2)*df['tpr']**2 + 2.1 * 10**(-4) * df['tpr']**3 + 1.33 * df['sal'] - 1.23 * 10**(-2) * df['sal'] * df['tpr'] + 8.7 * 10**(-5) * df['sal']* df['tpr']**2 + 1.56 * 10**(-2) * df['z'] + 2.55 * 10**(-7) * df['z']**2 - 7.3 * 10**(-12) * df['z']**3 + 1.2 * 10**(-6) * df['z'] * (phi - 45) - 9.5 * 10**(-13) * df['tpr'] * df['z']**3 + 3 * 10**(-7) * df['tpr']**2 * df['z'] + 1.43 * 10**(-5) * df['sal'] * df['z']
      
    # store sound speed in data frame
    df_svl_hrt = extract_df(df,column_list=['svl_hrt'],dtype=['f'])
    df_svl_tpr = extract_df(df,column_list=['svl_tpr'],dtype=['f'])
    df_leroy = pd.concat([st_series, df_svl_hrt, df_svl_tpr, df['sal']], axis=1)

    return(df_leroy)
# end def sv_leroy (HRT, PRS, SAL, phi):

#-------------------------------------------------------------------------------
#       Sound Velocity Del Grosso
#-------------------------------------------------------------------------------

def sv_delgrosso ( HRT, PRS, SAL):
    """ Calculates sound speed using the Del Grosso formula.

    Del Grosso formula: V. A. Del Grosso, 1974. New equation for the speed
        of sound in natural waters (with comparisons to other equations), J.
        Acoust. Soc. Am. 56, 1084, doi: 10.1121/1.1903388

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degree Celsius at
        one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon
    SAL ... constant salinity in parts per thousand

    It returns:
    df_ssp ... pandas.DataFrame holding the estimated sound speeds in metres
        per second
    """

    # concatenate Data.Frames of temperature and pressure
    df = pd.concat([HRT, PRS], axis=1)
    # set constant salinity value
    df['sal']=SAL

    # transform pressure from kPa to kilogram-force per square centimetres
    # 1 kPa = 0.010197266 kg_f/cm^2
    df['prs'] = df['prs']*0.010197266

    df['Ct'] = 0.501209398873*10*df['hrt'] -0.550946843172*10**(-1)*df['hrt']**2 + 0.22153596924*10**(-3)*df['hrt']**3

    df['Cs'] = 0.132952290781*10*df['sal'] + 0.128955756844*10**(-3)*df['sal']**2

    df['Cp'] = 0.156059257041*df['prs'] + 0.244998688441*10**(-4)*df['prs']**2  -0.883392332513*10**(-8)*df['prs']**3

    df['Cstp'] = -0.127562783426*10**(-1)*df['hrt']*df['sal'] + 0.635191613389*10**(-2)*df['hrt']*df['prs'] +0.265484716608*10**(-7)*df['hrt']**2*df['prs']**2-0.159349479045*10**(-5)*df['hrt']*df['prs']**2 + 0.522116437235*10**(-9)*row['hrt']*df['prs']**3 -0.438031096213*10**(-6)*df['hrt']**3*df['prs'] - 0.161674495909*10**(-8)*df['sal']**2*df['prs']**2 + 0.96840315641*10**(-4)*df['hrt']**2*df['sal'] + 0.485639620015*10**(-5)*df['hrt']*df['sal']**2*df['prs'] -0.340597039004*10**(-3)*df['hrt']*df['sal']*df['prs']

    df['ssp_d'] = 1402.392 + df['Ct'] + df['Cs'] + df['Cp'] + df['Cstp']

    # store sound speed in data frame
    df_ssp = extract_df(df,column_list=['ssp_d'],dtype=['f'])

    return(df_ssp)
# end def sv_delgrosso ( HRT, PRS, SAL):

#-------------------------------------------------------------------------------
#       DATA Read of specific sensor
#-------------------------------------------------------------------------------

def read_data(ID, sensor, starttime=None, endtime=None, pathname=None,suffix=None):
    """ Reads sensor data from file created with read().

    It needs:
    ID ... beacon ID
    sensor ... sensor flag: pressure - 'prs', sound speed - 'ssp',
        temperature - 'hrt', pages - 'pag', inclinometer - 'inc', battery -
        'bat'
    starttime (optional) ... no measurement before this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    endtime (optional) ... no measurement after this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    pathname (optional) ... location of files created with read() (default
        is ../DATA/)

    It returns:
    df ... pandas.DataFrame with requested data

    If file given by <pathname><ID>-<sensor>.dat (<sensor> will be
    transformed to uppercase) does not exist, an error message will be
    printed and an empty pandas.DataFrame will be returned.
    """

    if pathname is None:
        pathname = '../DATA/'
    # end if pathname is None:

    # file suffix is in uppercase
    sensor = str(sensor)
    sensor = sensor.upper()

    # check if chosen sensor is valid
    if sensor in ['SSP','HRT','PRS','PAG','BAT','INC', 'SVT','BSL','SAL','TMP', 'TPR']:
        # create file name
        if suffix is None:
            filename = pathname + str(ID) + '-' + sensor + '.dat'
        else:
            filename = pathname + str(ID) + '-' + sensor + '-' + str(suffix) + '.dat'
            # end if suffix is None:

        # read data from file in pandas.DataFrame
        try:
            df_data = pd.read_csv(filename,sep='\t',index_col=0,parse_dates=True)

            if starttime is not None:
                # remove all entries before starttime
                df_data2 = df_data.loc[df_data.index >= starttime]
            else:
                df_data2 = df_data
            # if starttime is not None:
            if endtime is not None:
                # remove all entries after endtime
                df = df_data2.loc[df_data2.index <= endtime]
            else:
                df = df_data2
            # if endtime is not None:


        except IOError as err:
            print('{0}'.format(err))
            print('Return empty DataFrame.')
            df = pd.DataFrame()
        # end try (read file)
    else:
        print('No valid sensor: {0}!'.format(sensor))
        print('Return empty DataFrame.')
        df = pd.DataFrame()
    # end if sensor in ['ssp','hrt','prs','pag','bat','inc']:

    return (df)
# end def read_data(ID, sensor, pathname=None):

#-------------------------------------------------------------------------------
#       DATA Read of Baseline Data
#-------------------------------------------------------------------------------

def read_bsl(ID1, ID2, starttime=None, endtime=None, pathname=None,suffix=None):
    """ Read baseline data from file created with sort_bsl().

    It needs:
    ID1 ... ID of beacon 1
    ID2 ... ID of beacon 2
    starttime (optional) ... no measurement before this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    endtime (optional) ... no measurement after this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    pathname (optional) ... location of baseline files (default is ../DATA/)
    suffix (optional) ... file suffix for saved baseline files (default - no
        suffix)

    It returns:
    df ... pandas.DataFrame with requested baseline data

    If file given by <pathname><ID1>-<ID2>.dat (or
    <pathname><ID1>-<ID2>-<suffix>.dat) does not exist, an error message
    will be printed and an empty pandas.DataFrame will be returned.
    """

    if pathname is None:
        pathname = '../DATA/'
    # end if pathname is None:

    if suffix is None:
        filename = pathname + str(ID1) + '-' + str(ID2) + '.dat'
    else:
        filename = pathname + str(ID1) + '-' + str(ID2) + '-' + str(suffix) + '.dat'
    # end if suffix is None:

    try:
        df_tmp = pd.read_csv(filename,sep=',',index_col=0,header=0,parse_dates=True)
        if starttime is not None:
            # remove all entries before starttime
            df_tmp2 = df_tmp.loc[df_tmp.index >= starttime]
        else:
            df_tmp2 = df_tmp
        # if starttime is not None:
        if endtime is not None:
            # remove all entries after endtime
            df = df_tmp2.loc[df_tmp2.index <= endtime]
        else:
            df = df_tmp2
        # if endtime is not None:

    except IOError as err:
        print('{0}'.format(err))
        print('Return empty DataFrame.')
        df = pd.DataFrame()
    # end try (read file)

    return (df)
# def read_bsl(ID1, ID2, pathname=None,suffix=None):

#-------------------------------------------------------------------------------
#       Meta Data
#-------------------------------------------------------------------------------

def read_meta(filename,pathname=None,sep=',',**kwds):
    """Reads meta data from file.

    File contains meta data which is specified by header. Columns should be
    seperated by ',' (default other seperator needs to be specified). In 
    column one the station ID is expected this serves as index.

    It needs:
    filename ... name of input file with meta information
    pathname (optional)  ... location of file with meta information (default
        is '../INFO/')
    sep (optional) ... seperator in file (default is ',')
    kwds ... keywords, options to pass to pandas.read_csv()


    It returns:
    df ... pandas.DataFrame with station ID as index
    """

    if pathname is None:
        pathname = '../INFO/'
    # end if pathname is None:

    try:
        df = pd.read_csv(pathname+filename,sep=sep,index_col=0,header=0,**kwds)
    except IOError as err:
        print('{0}'.format(err))
        print('Return empty DataFrame.')
        df = pd.DataFrame()
    # end try (read file)

    return(df)
# end def read_meta(filename):

#-------------------------------------------------------------------------------
#       DATA Read of Tides
#-------------------------------------------------------------------------------

def read_tides(filename, pathname=None):
    """ Reads tide data from file created with tpxo.

    It needs:

    filename 'tides.txt'
    pathname (optional) ... location of files created with read() (default
        is ../DATA/)

    It returns:
    df ... pandas.DataFrame with requested data

    If file given by filename does not exist, an error message will be
    printed and an empty pandas.DataFrame will be returned.
    """

    if pathname is None:
        pathname = '../Info/'
    # end if pathname is None:

    # create file name
    filename = pathname + filename
    # read data from file in pandas.DataFrame
    try:
        df = pd.read_csv(filename,sep='\t',index_col=0,header=0,parse_dates=True)
        df_tide = pd.to_numeric(df.tide, errors = 'coerce')
    except IOError as err:
        print('{0}'.format(err))
        print('Return empty DataFrame.')
        df = pd.DataFrame()

    # end try (read file)
    d = {'tide': df_tide}
    df = pd.DataFrame(data=d, index=df.index)



    return (df)
# end def read_data(ID, sensor, pathname=None):

#-------------------------------------------------------------------------------
#       DATA Read of Air Pressure
#-------------------------------------------------------------------------------

def read_airpressure(filename, pathname=None):
    """ Reads Air pressure data from file to correct local tides.

    It needs:

    filename 'airpressure.txt'
    pathname (optional) ... location of files created with read() (default
        is ../DATA/)

    It returns:
    df ... pandas.DataFrame with requested data

    If file given by filename does not exist, an error message will be
    printed and an empty pandas.DataFrame will be returned.
    """

    if pathname is None:
        pathname = '../DATA/'
    # end if pathname is None:

    # create file name
    filename = pathname + filename
    # read data from file in pandas.DataFrame
    try:
        df = pd.read_csv(filename,sep='\t',index_col=0,parse_dates=True)
    except IOError as err:
        print('{0}'.format(err))
        print('Return empty DataFrame.')
        df = pd.DataFrame()
    # end try (read file)

    return (df)
# end def read_data(ID, sensor, pathname=None):

#-------------------------------------------------------------------------------
#       Calculate Pressure Differences of Network
#-------------------------------------------------------------------------------

def calc_vert_motion(ID, tidesfile=None, starttime=None, freq=None):
    """ Calculates vertical pressure differences by subtracting pressure each other.

    It needs:

    ID      list of station IDs
    tides  (optional) ... filename of global or regional Tidemodel. Format:
        YYYY-MM-DDTHH:MM   XXXXX.X)

    starttime (optional) ... no measurement before this time is used ( format
        'YYYY-MM-DD hh:mm:ss')

    freq (optional) ... frequancy of movingaverage in pressure data (Default = 7d)

    It returns:
    List of vertical motion differences in cm

    """
    # convert tide data from kPa to dbar and subtract the mean
    if tidesfile is not None:
        df_tides = read_tides(tidesfile)
        tides = df_tides
        tides_mean = tides - tides.mean()
	
    # read pressure data from file and convert Kpa to dbar
    prs_corr = []
    for id in ID:

        if starttime is None:
            prs_dummy = ((read_data(id, 'prs')-100)/10)

        else:
            prs_dummy = ((read_data(id, 'prs',starttime=starttime)-100)/10)

        # subtract mean and append to a list of pandas.DataFrame
        prs_dummy = prs_dummy - prs_dummy.mean()

        if tidesfile is None:
            prs_corr.append(prs_dummy)

        else:
	    # add Index of reference station to DataFrame and interpolate the pressure in between
            df_tide_newindex=tides_mean.reindex(prs_dummy.index).append(tides_mean).sort_index()
            print(df_tide_newindex)
            tide_interpolate = df_tide_newindex.interpolate(method ='slinear')
            prs_corr.append(pd.DataFrame({"prs": prs_dummy.prs.subtract(tide_interpolate.tide).dropna()}))
        offset = []

    if freq is None:
        rolling_freq = '7d'

    else:
        rolling_freq = str(freq)

    # Loop over all Statsions
    for i ,df_prs1 in enumerate(prs_corr):

        print(' ')
        print('Station: ', ID[i])
        print(' ')
        for j , df_prs2 in enumerate(prs_corr):

            if i != j:
                # add Index of reference station to DataFrame and interpolate the pressure in between
                df_prs1_newindex=df_prs1.reindex(df_prs2.index).append(df_prs1).sort_index()
                prs_interpolate = df_prs1_newindex.interpolate(method ='slinear')
                prs_diff = pd.DataFrame({"prs": df_prs2.prs.subtract(prs_interpolate.prs).dropna()})

                # calculation of rolling mean to smoth dataset ans remove high tide frequancies
                prs_diff_mean = prs_diff.rolling(freq=rolling_freq,window=1).median().dropna()
                prs_diff_mean.to_csv('../DATA/'+ str(ID[i]) + '-' + str(ID[j]) + '-PRS.dat', sep='\t', header=False, date_format='%Y-%m-%dT%H:%M')
                # calculation of pressure difference from last to first rolling median entry
                first = prs_diff_mean.first('1s').values
                last = prs_diff_mean.last('1s').values

                # append to List and convert to cm
                offset.append((last - first)*100)
                off = (last - first)*100

                print(ID[j], ' ', off)

    return(offset)
# end def calc_vert_motion(ID, tidesfile=None, starttime=None, freq=None):

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

#-------------------------------------------------------------------------------
#       Calculates harmonic mean of sound speed
#-------------------------------------------------------------------------------

def calc_hmssp_recp_v(bsl):
    """Calculates harmonic mean of sound speeds and its reciprocal.

    It needs:
    bsl ... pandas.Dataframe with ID of beacon 1 ('ID'), ID of beacon 2
        ('range_ID'), calculated baseline lengths in metres ('bsl'), one
        way traveltime in seconds ('tt'), sound speed at beacon 1 ('ssp1')
        in metres per second, sound speed at beacon 2 ('ssp2') in metres per
        second, measured traveltime in milliseconds ('range'), turn around
        time in milliseconds ('TAT') with corresponding times of measurement
        for beacon pair.

    It returns:
    bsl ... pandas.Dataframe with ID of beacon 1 ('ID'), ID of beacon 2
        ('range_ID'), calculated baseline lengths in metres ('bsl'), one
        way traveltime in seconds ('tt'), sound speed at beacon 1 ('ssp1')
        in metres per second, sound speed at beacon 2 ('ssp2') in metres per
        second, measured traveltime in milliseconds ('range'), turn around
        time in milliseconds ('TAT'), harmonic mean of 'ssp1' and 'ssp2'
        ('hmssp') and reciprocal of harmonic mean of 'ssp1' and 'ssp2'
        ('1/v') with corresponding times of measurement for beacon pair.
    """

    bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']!=0.),'hmssp'] = 2*bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']!=0.),'ssp1']*bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']!=0.),'ssp2']/(bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']!=0.),'ssp2']+bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']!=0.),'ssp1'])
    bsl.loc[(bsl['ssp1']==0.)&(bsl['ssp2']!=0.),'hmssp'] = bsl.loc[(bsl['ssp1']==0.)&(bsl['ssp2']!=0.),'ssp2']
    bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']==0.),'hmssp'] = bsl.loc[(bsl['ssp1']!=0.)&(bsl['ssp2']==0.),'ssp1']
    bsl['1/v'] = 1/bsl['hmssp']

    return(bsl)
# end def calc_hmssp_recp_v(bsl):

def est_const_bsl(bsl,starttime=None,endtime=None,intercept=False,val_tw=None):
    """Performs a linear regression (assuming the intercept at the origin).

    The corresponding formula is tt-S*1/v-c = 0 in which tt is the travel
    time of the acoustic signal in seconds and 1/v is the reciprocal of the
    harmonic mean of the sound speed. The slope S is equal to the constant
    baseline length and by default c is assumed to be 0, but can optionally
    also be determined (intercept=True).

    It needs:
    bsl ... pandas.Dataframe with ID of beacon 1 ('ID'), ID of beacon 2
        ('range_ID'), calculated baseline lengths in metres ('bsl'), one
        way traveltime in seconds ('tt'), sound speed at beacon 1 ('ssp1')
        in metres per second, sound speed at beacon 2 ('ssp2') in metres per
        second, measured traveltime in milliseconds ('range'), turn around
        time in milliseconds ('TAT')(eventually harmonic mean of 'ssp1' and
        'ssp2' ('hmssp') and reciprocal of harmonic mean of 'ssp1' and
        'ssp2' ('1/v'); if they do not exist, they will be calculated) with
        corresponding times of measurement for beacon pair.
    starttime (optional) ... string with starttime of time window for
        estimation of constant baseline length (format: 'YYYY-mm-dd
        HH:MM:SS', default: first entry in bsl)
    endtime (optional) ... string with endtime of time window for estimation
        of constant baseline length (format: 'YYYY-mm-dd HH:MM:SS', default:
        last entry in bsl)
    intercept (optional) ... specify whether intercept should be set to
        0 [False] or should be calculated [True] (default is False)
    val_tw (optional) ... specify time window for which estimated constant
        baseline length and standard deviation (as well as intercept) will be
        stored in returned pandas.Dataframe (format: ['YYYY-mm-dd HH:MM:SS',
        'YYYY-mm-dd HH:MM:SS'], default is starttime and endtime)

    It returns:
    bsl ... pandas.Dataframe with ID of beacon 1 ('ID'), ID of beacon 2
        ('range_ID'), calculated baseline lengths in metres ('bsl'), one
        way traveltime in seconds ('tt'), sound speed at beacon 1 ('ssp1')
        in metres per second, sound speed at beacon 2 ('ssp2') in metres per
        second, measured traveltime in milliseconds ('range'), turn around
        time in milliseconds ('TAT'), harmonic mean of 'ssp1' and 'ssp2'
        ('hmssp'), reciprocal of harmonic mean of 'ssp1' and 'ssp2' ('1/v'),
        constant baseline length ('bsl_const') in given time window and
        standard deviation of the measurements compared to the fitted line
        in seconds (sigma = sqrt(sum((tt-S*1/v)^2)/(len(1/v)-1)),
        'std_dev_tt') in given time window (and intercept ('intercept') )
        with corresponding times of measurement for beacon pair.
    """

    # check if columns 'hmssp' and '1/v' (harmonic mean of sound speeds and its
    # reciprocal already exist in bsl and if not then calculate them
    if not set(['hmssp','1/v']).issubset(bsl.columns):
        bsl = calc_hmssp_recp_v(bsl)
    # end if not set(['hmssp','1/v']).issubset(bsl.columns):

    # copy bsl to new pandas.Dataframe to cut it in time
    bsl_new = bsl.copy()
    # check if time window for estimation of constant baseline length is given
    if starttime is not None:
        bsl_new = bsl_new.loc[starttime:]
    else:
        # set startime to first index in bsl
        starttime = bsl_new.index[0]
    # end if starttime is not None:
    if endtime is not None:
        bsl_new = bsl_new.loc[:endtime]
    else:
        # set endtime to last index in bsl
        endtime = bsl_new.index[-1]
    # end if endtime is not None:

    # the numpy function numpy.linalg.lstsq() needs x as (M,N) matrix
    if not intercept:
        x = bsl_new['1/v'][:,np.newaxis]
    else:
        x = np.array(([[bsl_new['1/v'][j], 1] for j in range(len(bsl_new))]))
    # end if not intercept:
    S,residuals,_,_ = np.linalg.lstsq(x,bsl_new['tt'])
    sigma = np.sqrt(residuals/(len(x)-1))

    # set column 'bsl_const' for values between starttime and endtime to S and
    # column 'std_dev_tt' to estimated sigma in bsl
    if val_tw is not None:
        starttime = val_tw[0]
        endtime = val_tw[1]
    # end if val_tw is not None:
    if not intercept:
        bsl.loc[starttime:endtime,'bsl_const'] = S
    else:
        bsl.loc[starttime:endtime,'bsl_const'] = S[0]
        bsl.loc[starttime:endtime,'intercept'] = S[1]
    # end if not intercept:
    bsl.loc[starttime:endtime,'std_dev_tt'] = sigma

    return(bsl)
# end def est_const_bsl(bsl,starttime=None,endtime=None):

#-------------------------------------------------------------------------------
#       Calculation of Backazimuth
#-------------------------------------------------------------------------------

def calc_bsl_baz(ID_pair,metafile,out='list',lat_column='lat',lon_column='lon',pathname=None):
    """ Calculates backazimuth of given ID pair(s).

    It needs:
    ID_pair ... list with ID_pairs e.g. [[2201,2202],[2202,2201]] (integer)
    metafile ... filename of meta file with columns 'lat' and 'lon' holding
        the latitude and longitude information for each station will be
        read by read_meta() (see function read_meta() for format
        specification)
    out (optional) ... flag to give back either list ('list', returns single
        value for single ID pair) or pandas.DataFrame ('pandas') (default is
        'list')
    lat_column (optional) ... name of column holding latitude (default is
        'lat')
    lon_column (optional) ... name of column holding longitude (default is
        'lon')
    pathname (optional) ... location of file with meta information (default is
        '../INFO/')


    It returns either:
    baz ... backazimuth of ID pair or 1D-list with backazimuths of ID pairs
        (for out = 'list')
    or
    df ... pd.DataFrame with column 'baz' (backazimuth of ID pair) and
        index <ID>-<range_ID> (for out = 'pandas')
    """

    # read meta data
    meta = read_meta(metafile,pathname)
    if out == 'list':
        baz = []
    elif out == 'pandas':
        df = pd.DataFrame()
    else:
        print("No valid option for out chosen. Will return to default (simple list)!")
        baz = []
        out = 'list'
    # end if out == 'list':

    try:
        ID = ID_pair[0][0]
    except TypeError:
        flag = 'single'
    else:
        flag = 'multi'
    # end try:

    if flag == 'single':
        _,baz_s,_ = gps2dist_azimuth(meta.lat[ID_pair[0]],meta.lon[ID_pair[0]],meta.lat[ID_pair[1]],meta.lon[ID_pair[1]])
        if out == 'list':
            baz = baz_s
        else:
            df = df.append(pd.DataFrame(data={'baz':baz_s},index=[str(ID_pair[0])+'-'+str(ID_pair[1])]))
        # end if out == 'pandas':
    elif flag == 'multi':
        for i in range(len(ID_pair)):
           _,baz_s,_ = gps2dist_azimuth(meta.lat[ID_pair[i][0]],meta.lon[ID_pair[i][0]],meta.lat[ID_pair[i][1]],meta.lon[ID_pair[i][1]])
           if out == 'list':
               baz.append(baz_s)
           else:
               df = df.append(pd.DataFrame(data={'baz':baz_s},index=[str(ID_pair[i][0])+'-'+str(ID_pair[i][1])]))
           # end if out == 'pandas':
        # end for i in range(len(ID_pair)):
    else:
        print("This should never appear!")
    # end if flag is 'single':/ elif flag is multi:
    if out == 'list':
        return(baz)
    else:
        return(df)
# end def calc_bsl_baz():

#-------------------------------------------------------------------------------
#       Moving average and standard deviation
#-------------------------------------------------------------------------------

def calc_mov_average_std(df,timespan,log_period,column_list=None,method=None):
    """Moving average and standard deviation in given time span for all/
       given columns in df.

    It needs:
    df ... pandas.DataFrame with data for whichmoving average and standard
        deviation should be calculated.
    timespan ... Time period of each window. Each window will be a variable
        sized based on the observations included in the time-period.
    log_period ... logging period in minutes (i.e. time span between two
        samples in minutes)
    column_list (optional) ... List of columns for which the moving average
        and the standard deviation should be computed (default is None ->
        for all columns)
    method (optional) ... (from pandas documentation:) method to use for
        filling holes in reindexed DataFrame. Please note: this is only
        applicable to DataFrames/Series with a monotonically
        increasing/decreasing index.
        default/ None: do not fill gaps
        pad / ffill: propagate last valid observation forward to next valid
        backfill / bfill: use next valid observation to fill gap
        nearest: use nearest valid observations to fill gap

    It returns:
    df_new ... pandas.DataFrame with data passed by in df and additional
        columns for all/ given columns holding the mean and the standard
        deviation for all time windows (mean_<old column name> and
        std_<old column name>)
    """

    # copy data frame to new data frame
    df_new = df.copy()

    # if no column list is given, use all columns in df_new
    if column_list is None:
        column_list = df_new.columns
    # end if column_list is None:

    # calculate theoretical number of samples in one time window
    tw = pd.to_timedelta(pd.tseries.frequencies.to_offset(timespan))
    period = pd.to_timedelta(str(log_period)+'min')
    samples = int(tw/period)
    # half of samples to shift moving average and standard deviation to
    # the middle of the given time window
    sample_shift = int(samples/2)
    # reindex df_new to include samples every log_period
    date_index = pd.date_range(df_new.index[0],df_new.index[-1],freq=str(log_period)+'min')
    # if regular sampling the existing index will be
    df_new = df_new.reindex(date_index,method=method)

    # iterate over all column_list items
    ncols = []
    for item in column_list:
         if item in df_new.columns:
             ncol_mean = 'mean_' + item
             ncol_std = 'std_' + item
             df_new[ncol_mean] = df_new[item].rolling(window=timespan).mean()
             df_new[ncol_std] = df_new[item].rolling(window=timespan).std()
             ncols.append(ncol_mean)
             ncols.append(ncol_std)
         # end if item in df_new.columns:
    # for item in column_list:

    df_new.loc[:,ncols] = df_new.loc[:,ncols].shift(-sample_shift)
    df_new.set_value(df_new.index[:sample_shift],ncols,np.nan)

    return(df_new)
# end def calc_mov_average_std(df,timespan,log_period,column_list=None)

def calc_diff_rel_col(df_list,timespan,log_period,column_list,fac=None,method=None):
    """Calculates relative changes to first mean and median of column(s).

    It needs:
    df_list ... 1D list with pandas.DataFrames holding the baseline data
    timespan ... Time period of each window for calculation of moving average.
        Each window will be a variable sized based on the observations
        included in the time-period.
    log_period ... logging period in minutes (i.e. time span between two
        samples in minutes)
    column_list ... list of columns for which relative changes should be
        calculated
    fac ... list of factors the relative changes should be multiplied with
        should have the same length as column_list (default is None -> factor
        for each column will be set to 1, same will be done for different
        length of fac and column_list)
    method (optional) ... (from pandas documentation:) method to use for
        filling holes in reindexed DataFrame. Please note: this is only
        applicable to DataFrames/Series with a monotonically
        increasing/decreasing index.
        default/None: do not fill gaps
        pad / ffill: propagate last valid observation forward to next valid
        backfill / bfill: use next valid observation to fill gap
        nearest: use nearest valid observations to fill gap



    It returns:
    df_list ... 1D list with pandas.DataFrames holding the baseline data
         and extra columns for 'mean_<column>' (moving average of <column>),
         'std_<column>' (standard deviation of <column>), 'diff_<column>'
         (difference of moving average of <column> with mean of column * fac
         for <column>), and 'rel_<column>' (difference of moving average of
         <column> with first moving average of <column> * fac for <column>)
    """

    # copy df_list to keep old df_list untouched
    df_list_new = []
    for i in range(len(df_list)):
        df_list_new.append(df_list[i].copy())
    # end for i in range(len(df_list)):
    # First of all check whether fac is present and has same length as column
    # list
    if fac is None:
        fac = np.ones(len(column_list))
    elif len(fac) != len(column_list):
        print('Length of fac and colum_list do not match!')
        print('Will use factor 1 for all columns in column list!')
        fac = np.ones(len(column_list))
    # end if fac is None:

    for i in range(len(df_list_new)):
        df_list_new[i] = calc_mov_average_std(df_list_new[i],timespan,log_period,column_list,method)
        for j,col in enumerate(column_list):
            mean_col = 'mean_'+col
            diff_col = 'diff_'+col
            rel_col = 'rel_'+col
            df_list_new[i][diff_col] = (df_list_new[i][mean_col]-df_list_new[i][col].mean())*fac[j]
            df_list_new[i][rel_col] = (df_list_new[i][mean_col]-df_list_new[i][mean_col].bfill()[0])*fac[j]
        # end for j,col in enumerate(column_list):
    # end for i in range(len(bsl_list)):

    return(df_list_new)
# end def calc_diff_rel_col(df_list,timespan,log_period,column_list,fac=None):

def calc_fc_data(df,column_list=None,ref_val='first'):
    """Calculate fractional change values for all/ given columns.

    Fractional change values are calculated using a reference value (first
    moving average, median or mean):
    fc_<column> = (<column>-<reference value>)/<refernce value>*1e6

    If column_list contains '1/v' the corresponding values for 'hmssp' will
    be used for the reference value. Furthermore, it is also possible that
    the selected column is a moving average.

    It needs:
    df ... pandas.DataFrame for which fractional values should be calculated
    column_list (optional) ... List of columns for which the fractional value
        should be computed (default is None -> for all columns)
    ref_val (optional) ... reference value for calculation of fractional
        change values:
        'first' (default) --- first moving average
        'median' --- median of all baselines (if moving average does not
             exist this will be default)
        'mean' --- mean of all baselines
        '<column_name>' --- any column of df
        if none of these valid options are chosen df_new will be identical to
        df

    It returns:
    df_new ... pandas.DataFrame with data passed by in df and additional
        columns for all/ given columns holding the fractional change values
        (fc_<old column name>)
    """

    # copy data frame to new data frame
    df_new = df.copy()

    # if no column list is given, use all columns in df_new
    if column_list is None:
        column_list = df_new.columns
    # end if column_list is None:

    for item in column_list:
        if '1/v' in item:
            rcol = 'hmssp'
        else:
            rcol = item
        ncol = 'fc_' + item
        ref = np.nan
        if ref_val == 'first':
            if 'mean' in rcol:
                col = rcol
            else:
                col = 'mean_' + rcol
            if col in df_new.columns:
                ref = df_new[col].bfill()[0]
            else:
                # set to second default if moving average does not exist
                ref_val='median'
            # end if col in df_new.columns
        # end if ref_val == 'first':
        if ref_val == 'median':
            ref = df_new[rcol].median()
        elif ref_val == 'mean':
            ref = df_new[rcol].mean()
        elif ref_val in df_new.columns:
            ref = df_new[ref_val]
        elif np.isnan(ref):
            print('No valid option for ref_val chosen!')
            print('Return original data frame!')
            return(df_new)
        # end if ref_val == 'median':
        if 'mean_1/v' == item and 'mean_hmssp' in df_new.columns:
            df_new[ncol] = (1/df_new['mean_hmssp']-ref)/ref*1e6
        elif 'mean_1/v' != item:
            df_new[ncol] = (df_new[item]-ref)/ref*1e6
    #end for item in column_list:

    return(df_new)
# end def calc_fc_data(df,column_list=None,ref_val='first')

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


def create_df4compare(df_list,column,name_ext):
    """Creates DataFrame from same column of multiple DataFrames (df_list) and
    resample it linear in time.

    It needs:
    df_list ... 1D list with pandas.DataFrames which have the common column(s)
        <column(s)>
    column ... column or list of columns to be copied to new DataFrame
    name_ext ... (string) list with extensions for each entry in df_list for
         new column(s) name (<column(s)>_<name_ext>), must be of same size as
         df_list

    It returns:
    df_compare ... pandas.DataFrame which holds all extracted columns sampled
        at the same times and ready for comparison
    """

    # First of all compare length of df_list with length of name_ext
    if len(df_list) != len(name_ext):
        print('Length of list with pandas.DataFrames does not match length of list with extensions!')
        print('Empty DataFrame will be returned!')
        df_compare = pd.DataFrame()
    else:
        # create 1D list of wished column(s)
        compare = []
        for i in range(len(df_list)):
            if type(column) == str:
                if column in df_list[i].columns:
                    compare.append(df_list[i][column].to_frame())
                    compare[i].columns = [column+name_ext[i]]
                else:
                    # if column does not exist for df_list[i], append empty
                    # DataFrame
                    compare.append(pd.DataFrame())
                # end if column in df_list[i].columns:
            elif type(column) == list:
                exist_col = []
                name_col = []
                for col in column:
                    if col in df_list[i].columns:
                        exist_col.append(col)
                        name_col.append(col+name_ext[i])
                    # end if col in df_list[i].columns:
                # for col in column:
                if len(exist_col)>1:
                    compare.append(df_list[i][exist_col])
                else:
                    compare.append(df_list[i][exist_col].to_frame())
                # end if len(exist_col)>1:
                compare[i].columns = name_col
            else:
                compare.append(pd.DataFrame())
            # end if type(column) == str:
        # end for i in range(len(df_list)):

        df_compare = pd.concat(compare,axis=1)
        df_compare = df_compare.interpolate(method='time')
        # limit each column to time window in original DataFrame
        for i in range(len(df_list)):
            if type(column) == str:
                if column in df_list[i].columns:
                    df_compare.loc[df_list[i].apply(pd.Series.last_valid_index)[column]+pd.tseries.offsets.DateOffset(seconds=1):,column+name_ext[i]] = np.nan
                # end if column in df_list[i].columns:
            elif type(column) == list:
                for col in column:
                    if col in df_list[i].columns:
                        df_compare.loc[df_list[i].apply(pd.Series.last_valid_index)[col]+pd.tseries.offsets.DateOffset(seconds=1):,col+name_ext[i]] = np.nan
                    # end if col in df_list[i].columns:
                # for col in column:
            # if type(column) == str:
        # end for i in range(len(df_list)):
    # end

    return(df_compare)
# end def create_df4compare(df_list,column,name_ext):


############################### End of file ################################
