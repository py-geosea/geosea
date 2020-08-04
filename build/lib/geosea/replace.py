#-------------------------------------------------------------------------------
#       Replace of broken Sensor data
#-------------------------------------------------------------------------------

import pandas as pd
from .read_data import *

GMT_DATEFORMAT = '%Y-%m-%dT%H:%M'
IN_DATEFORMAT = '%Y/%m/%d %H:%M:%S'


def replace(ID1,ID2, sensor, starttime=None, project=None):
    """ Replace incorrect Sensor data with a neighbour station
        
        It needs:
        ID1 ... ID of beacon with broken sensor
        ID2 ... ID of neighbour beacon
        
        sensor ... broken sensor
        starttime ... time before the sensor of ID1 is broken
        
        The replace data will be immediatly safed in ../DATA/
        
        The sensor data of ID2 will be interpolated with the index of ID1 by slinear method.
        The data needs to be shifted to avoid any leaps in the interpolated data.
        Check data afterwards!
        """
    
    if starttime is None:
        
        df_idata1 = read_data(ID1, str(sensor))
        df_idata2 = read_data(ID2, str(sensor))
        
        df_idata1_first = df_idata1 - df_idata1.first('1s').values.squeeze()
        df_idata2_first = df_idata2 - df_idata2.first('1s').values.squeeze()
        df_idata1_old_last = df_idata1 - df_idata1.first('1s').values.squeeze()
    
    elif project is None:
        
        df_idata1 = read_data(ID1, str(sensor), starttime=starttime)
        df_idata2 = read_data(ID2, str(sensor), starttime=starttime)
        df_idata1_old = read_data(ID1, str(sensor), endtime=starttime)
        
        df_idata1_first = df_idata1 - df_idata1.first('1s').values.squeeze()
        df_idata2_first = df_idata2 - df_idata2.first('1s').values.squeeze()
        
        df_idata1_old_last = df_idata1_old - df_idata1_old.last('1s').values.squeeze()

    else:
        
        df_idata1 = read_data(ID1, str(sensor), starttime=starttime)
        df_idata2 = read_data(ID2, str(sensor), starttime=starttime)
        
        # selection global PROJECTS starttime
        if project.upper() in PROJECTS:
            
            df_idata1_old = read_data(ID1, str(sensor),starttime=PROJECTS[project.upper()], endtime=starttime)

        # subtract the mean of both datasets
        df_idata1_first = df_idata1 - df_idata1.first('1s').values.squeeze()
        df_idata2_first = df_idata2 - df_idata2.first('1s').values.squeeze()

    df_idata1_old_last = df_idata1_old - df_idata1_old.last('1s').values.squeeze()

    # add Index of reference station to DataFrame and interpolate the sensor data in between
    df_data_newindex = df_idata2_first.reindex(df_idata1_first.index).append(df_idata2_first).sort_index()
    df_data_interpolate = df_data_newindex.interpolate(method ='slinear')
    df_odata1 = pd.DataFrame(df_data_interpolate, index = df_idata1.index)
    
    # Add back the mean of old sensor data
    df_odata1 = df_odata1 + df_idata1_old.last('1s').values.squeeze()
    
    if starttime is None:
        df_odata1.to_csv('../DATA/' + str(ID1) +'-'+ sensor.upper() +'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
    else:
        odata1_list = [df_idata1_old, df_odata1]
        df_odata1_all = pd.concat(odata1_list)
        
        # save the replaced data into the DATA directory
        df_odata1_all.to_csv('../DATA/' + str(ID1) +'-'+ sensor.upper() +'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)

    return(df_odata1_all)
    #end def replace (ID1,ID2, sensor, starttime=None, project=None):

    #-------------------------------------------------------------------------------
    #       Replace of HRT Data
    #-------------------------------------------------------------------------------

def replace_hrt(ID, starttime):
    """ Replace HRT Sensor data with the TMP data
    
        It needs:
        ID1 ... ID of beacon with broken sensor
    
        starttime
      
        The replace data will be immediatly safed in ../DATA/
    
        The sensor data of ID2 will be interpolated with the index of ID1 by slinear method.
        The data needs to be shifted to avoid any leaps in the interpolated data.
        Check data afterwards!
        """
    df_indata_hrt = read_data(ID, 'hrt', endtime=starttime)
    df_indata_tmp = read_data(ID, 'tmp', starttime=starttime)
    df_indata_tmp.columns = ['hrt']
    df_indata_hrt_last = df_indata_hrt - df_indata_hrt.last('1s').values.squeeze()
    df_indata_tmp_first = df_indata_tmp - df_indata_tmp.first('1s').values.squeeze()
    
    df_indata_hrt_concat = [df_indata_hrt_last, df_indata_tmp_first]
    df_outdata_hrt = pd.concat(df_indata_hrt_concat)
    
    df_outdata_hrt = df_outdata_hrt + df_indata_hrt.last('1s').values.squeeze()
        
    df_outdata_hrt.to_csv('../DATA/' + str(ID) +'-HRT.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
   
    return(df_outdata_hrt)

#-------------------------------------------------------------------------------
#       Complete Replace of Data
#-------------------------------------------------------------------------------

def replace_complete(ID1,ID2, sensor):
    """ Replace whole Sensor data with neighbour station
        
        It needs:
        ID1 ... ID of beacon with broken sensor
        ID2 ... ID of neighbour beacon
        
        sensor ... broken sensor
        
        The replace data will be immediatly safed in ../DATA/
        
        The sensor data of ID2 will be interpolated with the index of ID1 by slinear method.
        The data needs to be shifted to avoid any leaps in the interpolated data.
        Check data afterwards!
        """
  
    df_idata1 = read_data(ID1, str(sensor))
    df_idata2 = read_data(ID2, str(sensor))

    df_idata1_first = df_idata1
    df_idata2_first = df_idata2

    # add Index of reference station to DataFrame and interpolate the sensor data in between
    df_data_newindex = df_idata2_first.reindex(df_idata1_first.index).append(df_idata2_first).sort_index()
    df_data_interpolate = df_data_newindex.interpolate(method ='slinear')

    # Remove duplicate index in Dataframe
    df_data_interpolate = df_data_interpolate[~df_data_interpolate.index.duplicated(keep='first')]
    
    df_new = pd.DataFrame(df_data_interpolate.dropna(), index = df_idata1.index).dropna()
    
    # save the replaced data into the DATA directory
    df_new.to_csv('../DATA/' + str(ID1) +'-'+ sensor.upper() +'.dat',sep='\t', header=True, date_format=GMT_DATEFORMAT)
   
    return(df_new)
#end def replace_complete(ID1,ID2, sensor):

