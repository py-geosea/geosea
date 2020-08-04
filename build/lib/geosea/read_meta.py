#-------------------------------------------------------------------------------
#       Meta Data
#-------------------------------------------------------------------------------

import pandas as pd 

def read_meta(filename,pathname=None,sep=',',**kwds):
    """Reads meta data from file.

    File contains meta data which is specified by header. Columns should be
    seperated by ',' (default other seperator needs to be specified). In
    column one the station ID is expected this serves as index.

    It needs:
    filename ... name of input file with meta information
    pathname (optional)  ... location of file with meta information (default
        is '../INFO/')
    sep (optional) ... seperator in file (default is ',')
    kwds ... keywords, options to pass to pandas.read_csv()


    It returns:
    df ... pandas.DataFrame with station ID as index
    """

    if pathname is None:
        pathname = '../INFO/'
    # end if pathname is None:

    try:
        df = pd.read_csv(pathname+filename,sep=sep,index_col=0,header=0,**kwds)
    except IOError as err:
        print('{0}'.format(err))
        print('Return empty DataFrame.')
        df = pd.DataFrame()
    # end try (read file)

    return(df)
# end def read_meta(filename):
