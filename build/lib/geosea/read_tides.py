#-------------------------------------------------------------------------------
#       DATA Read of Tides
#-------------------------------------------------------------------------------

import pandas as pd 

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
