
import pandas as pd
#from obspy.core import UTCDateTime
#import datetime # Basic date and time types
#import scipy.io # ecosystem of open-source software for mathematics, science, and engineering


### Import GeoSEA Modules ###
from .change2dateindex import *
from .change_dtype import *

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

