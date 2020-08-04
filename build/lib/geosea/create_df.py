#-------------------------------------------------------------------------------
#       Create Pandas DataFrame from List
#-------------------------------------------------------------------------------

import pandas as pd

GMT_DATEFORMAT = '%Y-%m-%dT%H:%M'

def create_df(station,suffix,date,label_list,data_list):
    """Creates pandas.DataFrame from lists and writes data to file.

    It needs:
    station ... beacon ID
    suffix ... file suffix to specify data type
    date ... 1-dim list with times of measurement
    label_list ... 1-dim list with label(s) for data to be saved
    data_list ... 2-dim list with data to be saved (first dimension must
        correspond to number of labels, second dimension to length of date
        (if just one label, data_list is a 1-dim list)

    It returns:
    df ... pandas.DataFrame with stored data

    It additionally saves data to file.
    """

    # check if one or more data fields
    try:
        n = len(data_list[0])
    except TypeError:
    # just one label (m=1)
        n = len(data_list)
        m = 1
    else:
        # multiple labels
        m = len(data_list)
    # end try:

    if m == 1:
        # single label
        # dict for creating pandas data frame
        dict_data = {'date': date, label_list : data_list}
        # create list of labels for pandas.DataFrame
        column_list = ['date', label_list]
    else:
        # multiple label
        # dict for creating pandas data frame
        dict_data = {}
        dict_data['date'] = date
        # create list of labels for pandas.DataFrame
        column_list = ['date']
        # prepend 'date' to list
        for i in range(m):
            dict_data[label_list[i]] = data_list[i]
            column_list.append(label_list[i])
        # end for i in range(m):

    # create pandas.DataFrame
    df = pd.DataFrame(dict_data, columns = column_list)
    # transform date
    df['date'] = pd.to_datetime(df['date'])
    # set date as index of data frame
    df.index = df['date']
    # delete original date field
    del df['date']
    # drops entries which include a NaN, needs to be stored to df otherwise no
    # permanent change
    df = df.dropna()
    # writes data to file
    df.to_csv('../DATA/' + str(station) +'-'+ suffix+'.dat',sep='\t', header=False, date_format=GMT_DATEFORMAT)

    return(df)
# end def create_df(station,suffix,date,label_list,data_list)

