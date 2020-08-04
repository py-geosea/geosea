#-------------------------------------------------------------------------------
#       DATA Read of Baseline Data
#-------------------------------------------------------------------------------

import pandas as pd 

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

