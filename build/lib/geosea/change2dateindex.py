
import pandas as pd
import datetime # Basic date and time types

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
