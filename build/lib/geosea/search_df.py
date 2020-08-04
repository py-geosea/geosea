#-------------------------------------------------------------------------------
#       Search for nearest entry of Sensor data
#-------------------------------------------------------------------------------

import pandas as pd

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

