#-------------------------------------------------------------------------------
#       Calculate Baselines data
#-------------------------------------------------------------------------------

import pandas as pd
import multiprocessing


from .search_df import *
from .extract_df import *
from .calc import *

GMT_DATEFORMAT = '%Y-%m-%dT%H:%M'

def hori_bsl(ID,bsl_all,st_series,minmax,outlier_flag=None,writefile=True):
    """Calculates baselines for all possible pairs.

    It needs:
    ID ... 1-dim array with all available beacon IDs
    bsl_all ... an 1-dim list with pandas.DataFrame with baseline
        measurements: ID of other station ('range_ID'), traveltime ('range')
        and turn around time ('TAT') with corresponding times of measurement
        for each beacon (same order as items in ID)
    ssp_all ... an 1-dim list with pandas.DataFrame with sound speed ('ssp')
        with corresponding times of measurement for each beacon (same order
        as items in ID)
    minmax ... half time window length for searching for sound speed at
        beacon 2
    outlier_flag (optional) ... if set to 1 all baselines with lengths
        +/-10m are removed
    suffix (optional) ... file suffix for saved baseline files (default - no
        suffix)
    writefile (optional) ... if True files containing all baseline parameters
        will be created in current directory, if False data will just be
        returned (default True)

    It returns:
    ID_pair ... a 2-dim list with IDs of beacon pairs
        ([[ID1,ID2],[ID1,ID3],...])
    final_bsls ... an 1-dim list with pandas.DataFrame with ID of beacon 1
        ('ID'), ID of beacon 2 ('range_ID'), calculated baseline lengths in
        metres ('bsl'), one way traveltime in seconds ('tt'), sound speed at
        beacon 1 ('ssp1') in metres per second, sound speed at beacon 2
        ('ssp2') in metres per second, measured traveltime in milliseconds
        ('range'), turn around time in milliseconds ('TAT') with
        corresponding times of measurement for each beacon pair (same order
        as list items in ID_pair)
    """

    cal_bsl_all = []
    pool = multiprocessing.Pool(processes=4)
    ID_pair = []
    st_svleroy = []
        
    for i, beacon_1 in enumerate(ID):
        for j, beacon_2 in enumerate(ID):
            if beacon_1 != beacon_2:

                print('Baseline Calculation for: ' + str(beacon_1) + ' <-> ' + str(beacon_2))
                print('-------------------------------------------------------------------------------')
                ID_pair.append([beacon_1,beacon_2])
                # create new pandas.DataFrame holding baseline measurements
                # between beacon_1 and beacon_2 which are not 0.0 milli seconds
                df_bsl = bsl_all[i].loc[(bsl_all[i]['range_ID']==int(beacon_2)) & (bsl_all[i]['range']!=0.0)]

                if not df_bsl.empty and not st_series[i].empty and not st_series[j].empty:
                    # set new column 'ID' of beacon 1
                    df_bsl['ID']=int(beacon_1)
                    # alternative version:
                    #df_bsl.loc[df_bsl.index,'ID'] = int(beacon_1)

                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'ssp1',0)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_hrt1',10)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_tpr1',11)
                    
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'prs1',1)
                    
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'hrt1',2)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'tpr1',3)
                    df_bsl,SV_1_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sal1',12)
                    
                    
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'ssp2',0)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_hrt2',10)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sv_tpr2',11)
                    
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'prs2',1)
                    
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'hrt2',2)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'tpr2',3)
                    df_bsl,SV_2_err_count = search_df(df_bsl,st_series[i],minmax,0.0,'sal2',12)
                    
                    # define criteria that ssp1 and ssp2 have to be not NaN
                    # for row selection allow also single sided baseline
                    # calculation
    #---------------------------------------------------------------------------------------
    #            SSP
                    criteria_ssp = (df_bsl['ssp1']!=0.0) & (df_bsl['ssp2']!=0.0)
                    ssp2_check = (df_bsl['ssp2'] == 0.0) & (df_bsl['ssp1']!=0.0)
                    ssp1_check = (df_bsl['ssp1'] == 0.0) & (df_bsl['ssp2']!=0.0)
                    
                    df_bsl.loc[criteria_ssp,'tt'] = ((df_bsl[criteria_ssp]['range']-df_bsl[criteria_ssp]['TAT'])/2)/1000
                    df_bsl.loc[ssp2_check,'tt'] = ((df_bsl[ssp2_check]['range']-df_bsl[ssp2_check]['TAT'])/2)/1000
                    df_bsl.loc[ssp1_check,'tt'] = ((df_bsl[ssp1_check]['range']-df_bsl[ssp1_check]['TAT'])/2)/1000
                    
                    df_bsl.loc[criteria_ssp,'bsl'] = baseline_calc_hmean(df_bsl[criteria_ssp]['ssp1'],df_bsl[criteria_ssp]['ssp2'],df_bsl[criteria_ssp]['range'],df_bsl[criteria_ssp]['TAT'])
                    
                    df_bsl.loc[ssp1_check,'bsl'] = baseline_calc_hmean(df_bsl[ssp1_check]['ssp2'],df_bsl[ssp1_check]['ssp2'],df_bsl[ssp1_check]['range'],df_bsl[ssp1_check]['TAT'])

                    df_bsl.loc[ssp2_check,'bsl'] = baseline_calc_hmean(df_bsl[ssp2_check]['ssp1'],df_bsl[ssp2_check]['ssp1'],df_bsl[ssp2_check]['range'],df_bsl[ssp2_check]['TAT'])
                    
                    bsl_sucess = len(df_bsl.loc[pd.notnull(df_bsl['bsl'])])
    #---------------------------------------------------------------------------------------
    #            SV_HRT

                    criteria_sv_hrt = (df_bsl['sv_hrt1']!=0.0) & (df_bsl['sv_hrt2']!=0.0)
                    sv_hrt2_check = (df_bsl['sv_hrt2'] == 0.0) & (df_bsl['sv_hrt1']!=0.0)
                    sv_hrt1_check = (df_bsl['sv_hrt1'] == 0.0) & (df_bsl['sv_hrt2']!=0.0)
                    
                    df_bsl.loc[criteria_sv_hrt,'tt'] = ((df_bsl[criteria_sv_hrt]['range']-df_bsl[criteria_sv_hrt]['TAT'])/2)/1000
                    df_bsl.loc[sv_hrt2_check,'tt'] = ((df_bsl[sv_hrt2_check]['range']-df_bsl[sv_hrt2_check]['TAT'])/2)/1000
                    df_bsl.loc[sv_hrt1_check,'tt'] = ((df_bsl[sv_hrt1_check]['range']-df_bsl[sv_hrt1_check]['TAT'])/2)/1000
                    
                    df_bsl.loc[criteria_sv_hrt,'bsl_hrt'] = baseline_calc_hmean(df_bsl[criteria_sv_hrt]['sv_hrt1'],df_bsl[criteria_sv_hrt]['sv_hrt2'],df_bsl[criteria_sv_hrt]['range'],df_bsl[criteria_sv_hrt]['TAT'])
                    
                    df_bsl.loc[sv_hrt1_check,'bsl_hrt'] = baseline_calc_hmean(df_bsl[sv_hrt1_check]['sv_hrt2'],df_bsl[sv_hrt1_check]['sv_hrt2'],df_bsl[sv_hrt1_check]['range'],df_bsl[sv_hrt1_check]['TAT'])

                    df_bsl.loc[sv_hrt2_check,'bsl_hrt'] = baseline_calc_hmean(df_bsl[sv_hrt2_check]['sv_hrt1'],df_bsl[sv_hrt2_check]['sv_hrt1'],df_bsl[sv_hrt2_check]['range'],df_bsl[sv_hrt2_check]['TAT'])
                    
                    bsl_hrt_sucess = len(df_bsl.loc[pd.notnull(df_bsl['bsl_hrt'])])
    #---------------------------------------------------------------------------------------
    #            SV_TPR
    
                    criteria_sv_tpr = (df_bsl['sv_tpr1']!=0.0) & (df_bsl['sv_tpr2']!=0.0)
                    sv_tpr2_check = (df_bsl['sv_tpr2'] == 0.0) & (df_bsl['sv_tpr1']!=0.0)
                    sv_tpr1_check = (df_bsl['sv_tpr1'] == 0.0) & (df_bsl['sv_tpr2']!=0.0)
                    
                    df_bsl.loc[criteria_sv_tpr,'tt'] = ((df_bsl[criteria_sv_tpr]['range']-df_bsl[criteria_sv_tpr]['TAT'])/2)/1000
                    df_bsl.loc[sv_tpr2_check,'tt'] = ((df_bsl[sv_tpr2_check]['range']-df_bsl[sv_tpr2_check]['TAT'])/2)/1000
                    df_bsl.loc[sv_tpr1_check,'tt'] = ((df_bsl[sv_tpr1_check]['range']-df_bsl[sv_tpr1_check]['TAT'])/2)/1000
                    
                    df_bsl.loc[criteria_sv_hrt,'bsl_tpr'] = baseline_calc_hmean(df_bsl[criteria_sv_tpr]['sv_tpr1'],df_bsl[criteria_sv_tpr]['sv_tpr2'],df_bsl[criteria_sv_tpr]['range'],df_bsl[criteria_sv_tpr]['TAT'])
                    
                    df_bsl.loc[sv_tpr1_check,'bsl_tpr'] = baseline_calc_hmean(df_bsl[sv_tpr1_check]['sv_tpr2'],df_bsl[sv_tpr1_check]['sv_tpr2'],df_bsl[sv_tpr1_check]['range'],df_bsl[sv_tpr1_check]['TAT'])

                    df_bsl.loc[sv_tpr2_check,'bsl_tpr'] = baseline_calc_hmean(df_bsl[sv_tpr2_check]['sv_tpr1'],df_bsl[sv_tpr2_check]['sv_tpr1'],df_bsl[sv_tpr2_check]['range'],df_bsl[sv_tpr2_check]['TAT'])
                    
                    bsl_tpr_sucess = len(df_bsl.loc[pd.notnull(df_bsl['bsl_tpr'])])
                    
                    # calculate traveltime in seconds and store in new column of
                    # df_bsl

                    # calculate baseline length

                    
                    # count entries of df_bsl for which 'bsl' is not NaN

                else:
                    bsl_sucess = 0
                #end if not df_bsl.empty:

                df_bsl = extract_df(df_bsl,column_list=['ID','range_ID','range','TAT','tt','hrt1','hrt2','prs1','prs2','tpr1','tpr2','sal1','sal2','ssp1','ssp2','bsl','sv_hrt1','sv_hrt2','bsl_hrt','sv_tpr1','sv_tpr2','bsl_tpr'])
                if writefile:
                    df_bsl.to_csv('../DATA/' + beacon_1 +'-'+ beacon_2 +'-BSL.dat', header=True, date_format=GMT_DATEFORMAT)
                    # end if suffix == None:
                # end if writefile:

                print(str(len(df_bsl)) + '\t Ranges found')
                print(str(bsl_sucess) + '\t Successfull Calculated Baselines')
                if len(df_bsl) != 0:
                    print(str(SV_1_err_count) + '\t No SV Record in -> ' + str(beacon_1))
                    print(str(SV_2_err_count) + '\t No SV Record in -> ' + str(beacon_2))
                print(' \n')

                #  append to data formats of other stations
                cal_bsl_all.append(df_bsl)


            # end if beacon_1 != beacon_2:

        # end for j, beacon_2 in enumerate(ID):

    # end for i, beacon_1 in enumerate(ID):

    final_bsls = []
    for i in range(len(cal_bsl_all)):
        if outlier_flag == 1:
            ### cut off unrealistic Ranges and Baselines ###

            # This part removes all baselines with lengths +/-10m

            print('Cut Off unrealistic Ranges and Baselines')
            # calculate mean range exculding NaN value
            mean_bsl = cal_bsl_all[i].mean(axis=0,skipna=True)['bsl']
            # keep only those baselines within mean_bsl +/-10 m
            df_bsl = cal_bsl_all[i].loc[ (cal_bsl_all[i]['bsl']>mean_bsl-10) & (cal_bsl_all[i]['bsl']<mean_bsl+10)]
            print("Pair: {0:s} <-> {1:s}".format(ID_pair[i][0],ID_pair[i][1]))
            print("{0:d} baselines from {1:d} kept.".format(len(df_bsl),len(cal_bsl_all[i])))
        else:
            df_bsl = cal_bsl_all[i]
        # end if outlier_flag == 1:

        # re-arange order of columns
        df_bsl = extract_df(df_bsl,column_list=['ID','range_ID','bsl','tt','ssp1','ssp2','hrt1','hrt2','prs1','prs2','sal1','sal2','range','TAT'])
        final_bsls.append(df_bsl)

        if writefile:
            final_bsls[i].to_csv('../DATA/' + ID_pair[i][0] +'-'+ ID_pair[i][1] +'.dat', header=True, date_format=GMT_DATEFORMAT)
      
            # end if suffix is None:
        # end if writefile:
    # end for i in range(len(cal_bsl_all):

    if not writefile:
        print('\n')
        print('Data has not been stored in files!')
    # end if not writefile:


    return(ID_pair,final_bsls)
# end def sort_bsl( ... ):
