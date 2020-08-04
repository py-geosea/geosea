


import glob # Unix style pathname pattern expansion
import re # Regular expression operations
import numpy as np # fundamental package for scientific computing
from obspy.core import UTCDateTime # framework for processing seismological data
import datetime # Basic date and time types
import scipy.io # ecosystem of open-source software for mathematics, science, and engineering
import pandas as pd # open source, BSD-licensed library providing high-performance, easy-to-use data structures and data analysis tools
import matplotlib.pyplot as plt #
from obspy.geodetics.base import gps2dist_azimuth

### Import GeoSEA Modules ###
from .extract_df import *
from .change2dateindex import *
from .read_data import *
from .read_id import *

### Global Variables ###
GMT_DATEFORMAT = '%Y-%m-%dT%H:%M'

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
    
    ID = []

    if pathname is None:
        pathname = '../RAW/'
    
    ID = read_id(pathname)
    ifiles = glob.glob(pathname + 'Data_*_*_*.csv')
        
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
