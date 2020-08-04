#-------------------------------------------------------------------------------
#       Search for RAW Files to get Station IDs
#-------------------------------------------------------------------------------

import glob

def read_id(pathname=None):
    """ Reads the station Ids from *csv files.
    
    Note that the *csv files have to be unique for each station!
    
    pathname (optional) If the CSV file are not in the standard folder structre
    
    Returns:
    
    ID ... an 1-dim list with station IDs
    """

    if pathname is None:
        pathname = '../RAW/'
# end if pathname is None:

# returns list with file names which matches given wildcards (like ls)
    ifiles = glob.glob(pathname + 'Data_*_*_*.csv')

    ID = []
    for filename in ifiles:
    # extract beacon ID from file name
        filename_split = filename.split('_', 3)
        ID.append(filename_split[2])
    # end for filename in ifiles:

    # create unique ID list
    ID = list(set(ID))
    ID.sort()

    return(ID)
