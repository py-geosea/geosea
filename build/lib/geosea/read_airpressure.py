#-------------------------------------------------------------------------------
#       DATA Read of Air Pressure
#-------------------------------------------------------------------------------

import pandas as pd 

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
