#-------------------------------------------------------------------------------
#       Calculate Pressure Differences of Seafloor Geodetic Network
#-------------------------------------------------------------------------------

import pandas as pd

from .read_data import *
from .read_tides import *

def vert_bsl(ID, tidesfile=None, starttime=None, freq=None):
    """ Calculates vertical pressure differences by subtracting pressure from each other.

    It needs:

    ID      list of station IDs
    tides  (optional) ... filename of global or regional Tidemodel. Format:
        YYYY-MM-DDTHH:MM   XXXXX.X)

    starttime (optional) ... no measurement before this time is used ( format
        'YYYY-MM-DD hh:mm:ss')

    freq (optional) ... frequancy of movingaverage in pressure data (Default = 7d)

    It returns:
    List of vertical motion differences in cm

    """
    # convert tide data from kPa to dbar and subtract the mean
    if tidesfile is not None:
        df_tides = read_tides(tidesfile)
        tides = df_tides
        tides_mean = tides - tides.mean()
    
    # read pressure data from file and convert Kpa to dbar
    prs_corr = []
    for id in ID:

        if starttime is None:
            prs_dummy = ((read_data(id, 'prs')-100)/10)

        else:
            prs_dummy = ((read_data(id, 'prs',starttime=starttime)-100)/10)

        # subtract mean and append to a list of pandas.DataFrame
        prs_dummy = prs_dummy - prs_dummy.mean()

        if tidesfile is None:
            prs_corr.append(prs_dummy)

        else:
        # add Index of reference station to DataFrame and interpolate the pressure in between
            df_tide_newindex=tides_mean.reindex(prs_dummy.index).append(tides_mean).sort_index()
            print(df_tide_newindex)
            tide_interpolate = df_tide_newindex.interpolate(method ='slinear')
            prs_corr.append(pd.DataFrame({"prs": prs_dummy.prs.subtract(tide_interpolate.tide).dropna()}))
        offset = []

    if freq is None:
        rolling_freq = '7d'

    else:
        rolling_freq = str(freq)

    # Loop over all Statsions
    for i ,df_prs1 in enumerate(prs_corr):

        print(' ')
        print('Station: ', ID[i])
        print(' ')
        for j , df_prs2 in enumerate(prs_corr):

            if i != j:
                # add Index of reference station to DataFrame and interpolate the pressure in between
                df_prs1_newindex=df_prs1.reindex(df_prs2.index).append(df_prs1).sort_index()
                prs_interpolate = df_prs1_newindex.interpolate(method ='slinear')
                prs_diff = pd.DataFrame({"prs": df_prs2.prs.subtract(prs_interpolate.prs).dropna()})

                # calculation of rolling mean to smoth dataset ans remove high tide frequancies
                prs_diff_mean = prs_diff.rolling(freq=rolling_freq,window=1).median().dropna()
                prs_diff_mean.to_csv('../DATA/'+ str(ID[i]) + '-' + str(ID[j]) + '-PRS.dat', sep='\t', header=False, date_format='%Y-%m-%dT%H:%M')
                # calculation of pressure difference from last to first rolling median entry
                first = prs_diff_mean.first('1s').values
                last = prs_diff_mean.last('1s').values

                # append to List and convert to cm
                offset.append((last - first)*100)
                off = (last - first)*100

                print(ID[j], ' ', off)

    return(offset)
# end def calc_vert_motion(ID, tidesfile=None, starttime=None, freq=None):

