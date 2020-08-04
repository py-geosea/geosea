#-------------------------------------------------------------------------------
#       Complete processing
#-------------------------------------------------------------------------------

from .read import *
from .vert_bsl import *
from .hori_bsl import *

from .sw import *

def proc_bsl (SAL,phi,minmax,outlier_flag=None,writefile=True):
    """ Complete Baseline processing of GeoSEA Raw data.

    It needs:
    SAL ... constant salinity value
    phi ... Latitude for Leroy formular
    minmax ... half time window length for searching for sound speed at
    beacon 2

    It returns:
    bsl ... list of pandas.DataFrame with calculated Baselines

    """
    ID,st_series,bsl_series = read()
    
    st_series_leroy = []
    
    for i, id in enumerate(ID):
        st_series_leroy.append(sv_leroy(st_series[i],SAL,phi))
        
    bsl_horizonal = hori_bsl(ID,bsl_series,st_series_leroy,7600,outlier_flag,writefile)
    
    bsl_vertical = vert_bsl(ID)
    
    return(bsl)
