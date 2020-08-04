#-------------------------------------------------------------------------------
#       Utils for Baseline calculation
#-------------------------------------------------------------------------------


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
