
import pandas as pd

def compare_df(df_list,column,name_ext):
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
