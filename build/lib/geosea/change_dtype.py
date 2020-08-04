
import pandas as pd

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

