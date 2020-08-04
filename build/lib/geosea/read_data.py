#-------------------------------------------------------------------------------
#       DATA Read of specific sensor
#-------------------------------------------------------------------------------

import pandas as pd 


def read_data(ID, sensor, starttime=None, endtime=None, pathname=None,suffix=None):
    """ Reads sensor data from file created with read().

    It needs:
    ID ... beacon ID
    sensor ... sensor flag: pressure - 'prs', sound speed - 'ssp',
        temperature - 'hrt', pages - 'pag', inclinometer - 'inc', battery -
        'bat'
    starttime (optional) ... no measurement before this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    endtime (optional) ... no measurement after this time is used (format
        'YYYY-MM-DD hh:mm:ss')
    pathname (optional) ... location of files created with read() (default
        is ../DATA/)

    It returns:
    df ... pandas.DataFrame with requested data

    If file given by <pathname><ID>-<sensor>.dat (<sensor> will be
    transformed to uppercase) does not exist, an error message will be
    printed and an empty pandas.DataFrame will be returned.
    """

    if pathname is None:
        pathname = '../DATA/'
    # end if pathname is None:

    # file suffix is in uppercase
    sensor = str(sensor)
    sensor = sensor.upper()

    # check if chosen sensor is valid
    if sensor in ['SSP','HRT','PRS','PAG','BAT','INC', 'SVT','BSL','SAL','TMP', 'TPR']:
        # create file name
        if suffix is None:
            filename = pathname + str(ID) + '-' + sensor + '.dat'
        else:
            filename = pathname + str(ID) + '-' + sensor + '-' + str(suffix) + '.dat'
            # end if suffix is None:

        # read data from file in pandas.DataFrame
        try:
            df_data = pd.read_csv(filename,sep='\t',index_col=0,parse_dates=True)

            if starttime is not None:
                # remove all entries before starttime
                df_data2 = df_data.loc[df_data.index >= starttime]
            else:
                df_data2 = df_data
            # if starttime is not None:
            if endtime is not None:
                # remove all entries after endtime
                df = df_data2.loc[df_data2.index <= endtime]
            else:
                df = df_data2
            # if endtime is not None:


        except IOError as err:
            print('{0}'.format(err))
            print('Return empty DataFrame.')
            df = pd.DataFrame()
        # end try (read file)
    else:
        print('No valid sensor: {0}!'.format(sensor))
        print('Return empty DataFrame.')
        df = pd.DataFrame()
    # end if sensor in ['ssp','hrt','prs','pag','bat','inc']:

    return (df)
# end def read_data(ID, sensor, pathname=None):
