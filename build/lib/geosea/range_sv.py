#-------------------------------------------------------------------------------
#       Calculate sound speed from Baseline
#-------------------------------------------------------------------------------

import pandas as pd

def range_sv ( ID1, ID2, BSL ):
    """ Calculates sound speed from baseline measurements at one beacon.

    It needs:
    ID1 ... ID of beacon 1
    ID2 ... ID of beacon 2
    BSL ... pandas.DataFrame with baseline data at beacon 1 (ID of beacon 2
        ('range_ID'), measured traveltime in milliseconds ('range'), turn
        around time in milliseconds ('TAT'), sound speed at beacon 1
        ('ssp1') in metres per second, sound speed at beacon 2 ('ssp2') in
        metres per second one way traveltime in seconds ('tt'), calculated
        baseline lengths in metres ('bsl') with corresponding times of
        measurement)

    It returns:
    df_ssp_cal ... pandas.DataFrame with estimated sound speeds in metres
        per second with corresponding times

    If BSL contains baseline data from more than one ID_pair or if the IDs
    in BSL do not match the given IDs (ID1 and ID2), an error message will
    be printed and an empty pandas.DataFrame will be returned.
    """

    # create empty pandas.DataFrame for sound speed estimation
    df_ssp_cal = pd.DataFrame()

    dummy = BSL.loc[:,['ID','range_ID']].drop_duplicates()
    if (len(dummy) == 1) and (dummy.iloc[0][0] == int(ID1)) and  (dummy.iloc[0][1] == int(ID2)):
        # calculate mean range exculding NaN value
        mean_bsl = BSL.mean(axis=0,skipna=True)['bsl']
        # calculate sound speed from mean baseline lengths and travel time
        df_ssp_cal['ssp'] = mean_bsl/BSL.loc[:,'tt']
    else:
        print('Given baseline DateFrame contains either more than one ID pair orgiven ID pair does not match!')
        print('Empty DataFrame will be returned!')
    # end if (...):

    return(df_ssp_cal)
#end def range_SoundVelocity ( ID1, ID2, BSL ):
