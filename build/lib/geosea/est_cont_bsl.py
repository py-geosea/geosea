
import pandas as pd

from calc import *

def est_const_bsl(bsl,starttime=None,endtime=None,intercept=False,val_tw=None):
    """Performs a linear regression (assuming the intercept at the origin).

    The corresponding formula is tt-S*1/v-c = 0 in which tt is the travel
    time of the acoustic signal in seconds and 1/v is the reciprocal of the
    harmonic mean of the sound speed. The slope S is equal to the constant
    baseline length and by default c is assumed to be 0, but can optionally
    also be determined (intercept=True).

    It needs:
    bsl ... pandas.Dataframe with ID of beacon 1 ('ID'), ID of beacon 2
        ('range_ID'), calculated baseline lengths in metres ('bsl'), one
        way traveltime in seconds ('tt'), sound speed at beacon 1 ('ssp1')
        in metres per second, sound speed at beacon 2 ('ssp2') in metres per
        second, measured traveltime in milliseconds ('range'), turn around
        time in milliseconds ('TAT')(eventually harmonic mean of 'ssp1' and
        'ssp2' ('hmssp') and reciprocal of harmonic mean of 'ssp1' and
        'ssp2' ('1/v'); if they do not exist, they will be calculated) with
        corresponding times of measurement for beacon pair.
    starttime (optional) ... string with starttime of time window for
        estimation of constant baseline length (format: 'YYYY-mm-dd
        HH:MM:SS', default: first entry in bsl)
    endtime (optional) ... string with endtime of time window for estimation
        of constant baseline length (format: 'YYYY-mm-dd HH:MM:SS', default:
        last entry in bsl)
    intercept (optional) ... specify whether intercept should be set to
        0 [False] or should be calculated [True] (default is False)
    val_tw (optional) ... specify time window for which estimated constant
        baseline length and standard deviation (as well as intercept) will be
        stored in returned pandas.Dataframe (format: ['YYYY-mm-dd HH:MM:SS',
        'YYYY-mm-dd HH:MM:SS'], default is starttime and endtime)

    It returns:
    bsl ... pandas.Dataframe with ID of beacon 1 ('ID'), ID of beacon 2
        ('range_ID'), calculated baseline lengths in metres ('bsl'), one
        way traveltime in seconds ('tt'), sound speed at beacon 1 ('ssp1')
        in metres per second, sound speed at beacon 2 ('ssp2') in metres per
        second, measured traveltime in milliseconds ('range'), turn around
        time in milliseconds ('TAT'), harmonic mean of 'ssp1' and 'ssp2'
        ('hmssp'), reciprocal of harmonic mean of 'ssp1' and 'ssp2' ('1/v'),
        constant baseline length ('bsl_const') in given time window and
        standard deviation of the measurements compared to the fitted line
        in seconds (sigma = sqrt(sum((tt-S*1/v)^2)/(len(1/v)-1)),
        'std_dev_tt') in given time window (and intercept ('intercept') )
        with corresponding times of measurement for beacon pair.
    """

    # check if columns 'hmssp' and '1/v' (harmonic mean of sound speeds and its
    # reciprocal already exist in bsl and if not then calculate them
    if not set(['hmssp','1/v']).issubset(bsl.columns):
        bsl = calc_hmssp_recp_v(bsl)
    # end if not set(['hmssp','1/v']).issubset(bsl.columns):

    # copy bsl to new pandas.Dataframe to cut it in time
    bsl_new = bsl.copy()
    # check if time window for estimation of constant baseline length is given
    if starttime is not None:
        bsl_new = bsl_new.loc[starttime:]
    else:
        # set startime to first index in bsl
        starttime = bsl_new.index[0]
    # end if starttime is not None:
    if endtime is not None:
        bsl_new = bsl_new.loc[:endtime]
    else:
        # set endtime to last index in bsl
        endtime = bsl_new.index[-1]
    # end if endtime is not None:

    # the numpy function numpy.linalg.lstsq() needs x as (M,N) matrix
    if not intercept:
        x = bsl_new['1/v'][:,np.newaxis]
    else:
        x = np.array(([[bsl_new['1/v'][j], 1] for j in range(len(bsl_new))]))
    # end if not intercept:
    S,residuals,_,_ = np.linalg.lstsq(x,bsl_new['tt'])
    sigma = np.sqrt(residuals/(len(x)-1))

    # set column 'bsl_const' for values between starttime and endtime to S and
    # column 'std_dev_tt' to estimated sigma in bsl
    if val_tw is not None:
        starttime = val_tw[0]
        endtime = val_tw[1]
    # end if val_tw is not None:
    if not intercept:
        bsl.loc[starttime:endtime,'bsl_const'] = S
    else:
        bsl.loc[starttime:endtime,'bsl_const'] = S[0]
        bsl.loc[starttime:endtime,'intercept'] = S[1]
    # end if not intercept:
    bsl.loc[starttime:endtime,'std_dev_tt'] = sigma

    return(bsl)
# end def est_const_bsl(bsl,starttime=None,endtime=None):
