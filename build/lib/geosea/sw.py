#################### Seawater Module #########################

#-------------------------------------------------------------------------------
#       Salinity Wilson
#-------------------------------------------------------------------------------

import pandas as pd
import numpy as np

from .extract_df import *

def sal_wilson ( HRT, SSP, PRS ):
    """ Calculates salinity at one beacon using the Wilson formula.

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    SSP ... pandas.DataFrame with measured sound speed in metres per second
        at one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon

    It returns:
    df_sal ... pandas.DataFrame with calculated salinities in parts per
        thousand
    """

    # concatenate Data.Frames of temperature, sound speed and pressure
    df = pd.concat([HRT, SSP, PRS], axis=1)

    # transform pressure from kPa to kilogram-force per square centimetres
    # 1 kPa = 0.010197266 kg_f/cm^2
    df['prs'] = df['prs']*0.010197266

    df['S_t'] = -3.31986*df['hrt'] - 2.57236*10**(-3)*df['hrt']**2 + 2.32009*10**(-4)*df['hrt']**3 - 2.20892*10**(-6)*df['hrt']**4

    df['S_p'] = -1.14663*10**(-1)*df['prs'] - 1.00488*10**(-5)*df['prs']**2 - 2.14038 * 10**(-8)*df['prs']**3 + 2.23490 * 10**(-12) * df['prs']**4

    df['S_v'] = 7.21467 * 10**(-1) * (df['ssp'] - 1448.54) - 6.52575 * 10**(-4)* ( df['ssp'] - 1448.54)**2

    df['S_tpv'] = (df['ssp'] - 1448.54) * (1.17 * 10**(-2)*df['hrt'] + 1.04088 * 10**(-4) * df['prs'] + 1.29184 * 10**(-7) * df['prs']**2 + 4.43572 * 10**(-6) * df['prs'] * df['hrt'] - 9.63357 * 10**(-8) * df['prs'] * df['hrt']**2) + df['prs'] * (-1.45 * 10**(-3) * df['hrt'] - 2.5563 * 10**(-5) * df['hrt']**2 + 3.45997 * 10**(-7) * df['hrt']**3) + df['prs']**2 * (-7.29321 * 10**(-7) * df['hrt'] + 2.05515 * 10**(-8) * df['hrt']**2) + df['prs']**3 * (-1.19869 * 10**(-10) * df['hrt'])

    df['sal'] = 35 + df['S_t'] + df['S_p'] + df['S_v'] + df['S_tpv']

    # store salinity in data frame
    df_sal = extract_df(df,column_list=['sal'],dtype=['f'])

    return(df_sal)
# end def sal_wilson ( HRT, SSP, PRS ):

#-------------------------------------------------------------------------------
#       Salinity Medwin
#-------------------------------------------------------------------------------

def sal_medwin ( HRT, SSP, PRS ):
    """ Calculates salinity at one beacon using the Medwin formula.

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    SSP ... pandas.DataFrame with measured sound speed in metres per second
        at one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon

    It returns:
    df_sal ... pandas.DataFrame with calculated salinities in parts per
        thousand
    """

    # concatenate Data.Frames of temperature, sound speed and pressure
    df = pd.concat([HRT, SSP, PRS], axis=1)

    # pressure in depth (m)??
    df['prs'] = (df['prs']-100)/10

    df['sal'] = (((1449.2 + 4.6*df['hrt'] - 0.055*df['hrt']**2 + 0.00029*df['hrt']**3 + 0.016 * df['prs'] - df['ssp'])*(-1))/(1.34 - 0.010*df['hrt']))+35

    # store salinity in data frame
    df_sal = extract_df(df,column_list=['sal'],dtype=['f'])

    return(df_sal)

# end def sal_medwin ( ... ):

#-------------------------------------------------------------------------------
#       Sound Velocity Wilson
#-------------------------------------------------------------------------------

def sv_wilson ( HRT, PRS, SAL):
    """ Calculates sound speed using the Wilson formula.

    Wilson formula: Equation 1 in  Wayne D. Wilson, 1960. Speed of Sound in
        Sea Water as a Function of Temperature, Pressure, and Salinity, J.
        Acoust. Soc. Am. 32, 641, doi: 10.1121/1.1908167

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon
    SAL ... constant salinity in parts per thousand

    It returns:
    df_ssp ... pandas.DataFrame holding the estimated sound speeds
        ('ssp_w') in metres per second at one beacon
    """

    # concatenate Data.Frames of temperature and pressure
    if isinstance(SAL, pd.DataFrame):
        df = pd.concat([HRT, PRS, SAL], axis=1)
    else:
        # set constant salinity value
        df = pd.concat([HRT, PRS ], axis=1)
        df['sal']=SAL

    # transform pressure from kPa to kilogram-force per square centimetres
    # 1 kPa = 0.010197266 kg_f/cm^2
    df['prs'] = df['prs']*0.010197266

    df['V_t'] = 4.6233 * df['hrt'] - 5.4585 * 10**(-2) * df['hrt']**2 + 2.822 * 10**(-4) * df['hrt']**3 - 5.07 * 10**(-7) * df['hrt']**4

    df['V_p'] = 1.60518 * 10**(-1) * df['prs'] + 1.0279 * 10**(-5) * df['prs']**2 + 3.451 * 10**(-9) *df['prs']**3 - 3.503 * 10**(-12) * df['prs']**4

    df['V_s'] = 1.391 * (df['sal']-35) - 7.8 * 10**(-2) * (df['sal']-35)
    df['V_stp'] = (df['sal']-35) * (-1.197 * 10**(-2) * df['hrt'] + 2.61 * 10**(-4) * df['prs'] - 1.96* 10**(-7) * df['prs']**2 - 2.09 * 10**(-6) * df['prs'] * df['hrt']) + df['prs'] *(-2.796 * 10**(-4) * df['hrt'] + 1.3302 * 10**(-5) * df['hrt']**2 - 6.644*10**(-8) * df['hrt']**3) + df['prs']**2 * (-2.391 * 10**(-7) * df['hrt']+ 9.286 * 10**(-10) * df['hrt']**2) - 1.745 * 10**(-10) * df['prs']** 3 * df['hrt']

    df['ssp'] = 1449.22 + df['V_t'] + df['V_p'] + df['V_s'] + df['V_stp']

    # store sound speed in data frame
    df_ssp = extract_df(df,column_list=['ssp'],dtype=['f'])
    df_ssp = pd.concat([ df_ssp, HRT, PRS, df['sal']], axis=1)
    return(df_ssp)
# end def sv_wilson ( HRT, PRS, SAL):

#-------------------------------------------------------------------------------
#       Sound Velocity Leroy
#-------------------------------------------------------------------------------

#def sv_leroy (HRT, PRS, SAL, phi):
def sv_leroy (st_series, SAL, phi):
    """ Calculates sound speed using the Leroy formula.

    Leroy Formula: Equation 2 in Claude C. Leroy, Stephen P. Robinson and
        Mike J. Goldsmith, 2008. A new equation for the accurate calculation
        of sound speed in all oceans, J. Acoust. Soc. Am. 124, 2774,
        doi: 10.1121/1.2988296
    Conversion of pressure into depth: Equations 3 and 4 in C.C. Leroy and
        F. Parthiot, 1998. Depth-pressure relationships in the oceans and
        seas, J. Acoust. Soc. Am. 103, 1346, doi: 10.1121/1.421275 and C.C.
        Leroy, 2007. Erratum: "Depth-pressure relationships in the oceans
        and seas" [J. Acoust. Soc. Am.103(3), 1346-1352 (1998)], J. Acoust.
        Soc. Am. 121, 2447; doi: 10.1121/1.2534201 gives depth with an
        accuracy of +/- 1 metre

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degrees Celsius at
        one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon
    SAL ... constant salinity in parts per thousand
    phi ... latitude of working area in degrees

    It returns:
    df_ssp ... pandas.DataFrame holding the estimated sound speeds
        ('ssp_l') in metres per second
    """

    # concatenate Data.Frames of temperature and pressure

    if isinstance(SAL, pd.DataFrame):
        df = pd.concat([st_series, SAL], axis=1)
    else:
    # set constant salinity value
        df = pd.concat([st_series ], axis=1)
        df['sal']=SAL

    ##### transform pressure in kPa in depth in m #####

    # transform latitude in radians
    phi_rad = phi*np.pi/180
    # international formula for gravity (eq. 4 in Leroy and Parthiot, 1998)
    g = 9.780318*(1+5.2788*10**(-3)*(np.sin(phi_rad))**2-2.36*10**(-5)*(np.sin(phi_rad))**4)

    # transformation of pressure in kPa into MPa which is used in equation 3 in
    # Leroy and Parthiot (1998)
    df['prs'] = df['prs']/1000
    # transform pressure in MPa into depth in metres
    df['z'] = (9.72659*10**2*df['prs']-2.2512*10**(-1)*df['prs']**2+2.279*10**(-4)*df['prs']**3-1.82*10**(-7)*df['prs']**4)/(g+1.092*10**(-4)*df['prs'])

    # calculate sound speed in metres per second with eq.2 in Leroy et al.(2008)
    if 'hrt' in df:
        df['svl_hrt'] = 1402.5 + 5 * df['hrt'] - 5.44 * 10**(-2)*df['hrt']**2 + 2.1 * 10**(-4) * df['hrt']**3 + 1.33 * df['sal'] - 1.23 * 10**(-2) * df['sal'] * df['hrt'] + 8.7 * 10**(-5) * df['sal']* df['hrt']**2 + 1.56 * 10**(-2) * df['z'] + 2.55 * 10**(-7) * df['z']**2 - 7.3 * 10**(-12) * df['z']**3 + 1.2 * 10**(-6) * df['z'] * (phi - 45) - 9.5 * 10**(-13) * df['hrt'] * df['z']**3 + 3 * 10**(-7) * df['hrt']**2 * df['z'] + 1.43 * 10**(-5) * df['sal'] * df['z']
    if 'tpr' in df:
        df['svl_tpr'] = 1402.5 + 5 * df['tpr'] - 5.44 * 10**(-2)*df['tpr']**2 + 2.1 * 10**(-4) * df['tpr']**3 + 1.33 * df['sal'] - 1.23 * 10**(-2) * df['sal'] * df['tpr'] + 8.7 * 10**(-5) * df['sal']* df['tpr']**2 + 1.56 * 10**(-2) * df['z'] + 2.55 * 10**(-7) * df['z']**2 - 7.3 * 10**(-12) * df['z']**3 + 1.2 * 10**(-6) * df['z'] * (phi - 45) - 9.5 * 10**(-13) * df['tpr'] * df['z']**3 + 3 * 10**(-7) * df['tpr']**2 * df['z'] + 1.43 * 10**(-5) * df['sal'] * df['z']
      
    # store sound speed in data frame
    df_svl_hrt = extract_df(df,column_list=['svl_hrt'],dtype=['f'])
    df_svl_tpr = extract_df(df,column_list=['svl_tpr'],dtype=['f'])
    df_leroy = pd.concat([st_series, df_svl_hrt, df_svl_tpr, df['sal']], axis=1)

    return(df_leroy)
# end def sv_leroy (HRT, PRS, SAL, phi):

#-------------------------------------------------------------------------------
#       Sound Velocity Del Grosso
#-------------------------------------------------------------------------------

def sv_delgrosso ( HRT, PRS, SAL):
    """ Calculates sound speed using the Del Grosso formula.

    Del Grosso formula: V. A. Del Grosso, 1974. New equation for the speed
        of sound in natural waters (with comparisons to other equations), J.
        Acoust. Soc. Am. 56, 1084, doi: 10.1121/1.1903388

    It needs:
    HRT ... pandas.DataFrame with measured temperature in degree Celsius at
        one beacon
    PRS ... pandas.DataFrame with measured pressure in kPa at one beacon
    SAL ... constant salinity in parts per thousand

    It returns:
    df_ssp ... pandas.DataFrame holding the estimated sound speeds in metres
        per second
    """

    # concatenate Data.Frames of temperature and pressure
    df = pd.concat([HRT, PRS], axis=1)
    # set constant salinity value
    df['sal']=SAL

    # transform pressure from kPa to kilogram-force per square centimetres
    # 1 kPa = 0.010197266 kg_f/cm^2
    df['prs'] = df['prs']*0.010197266

    df['Ct'] = 0.501209398873*10*df['hrt'] -0.550946843172*10**(-1)*df['hrt']**2 + 0.22153596924*10**(-3)*df['hrt']**3

    df['Cs'] = 0.132952290781*10*df['sal'] + 0.128955756844*10**(-3)*df['sal']**2

    df['Cp'] = 0.156059257041*df['prs'] + 0.244998688441*10**(-4)*df['prs']**2  -0.883392332513*10**(-8)*df['prs']**3

    df['Cstp'] = -0.127562783426*10**(-1)*df['hrt']*df['sal'] + 0.635191613389*10**(-2)*df['hrt']*df['prs'] +0.265484716608*10**(-7)*df['hrt']**2*df['prs']**2-0.159349479045*10**(-5)*df['hrt']*df['prs']**2 + 0.522116437235*10**(-9)*row['hrt']*df['prs']**3 -0.438031096213*10**(-6)*df['hrt']**3*df['prs'] - 0.161674495909*10**(-8)*df['sal']**2*df['prs']**2 + 0.96840315641*10**(-4)*df['hrt']**2*df['sal'] + 0.485639620015*10**(-5)*df['hrt']*df['sal']**2*df['prs'] -0.340597039004*10**(-3)*df['hrt']*df['sal']*df['prs']

    df['ssp_d'] = 1402.392 + df['Ct'] + df['Cs'] + df['Cp'] + df['Cstp']

    # store sound speed in data frame
    df_ssp = extract_df(df,column_list=['ssp_d'],dtype=['f'])

    return(df_ssp)
# end def sv_delgrosso ( HRT, PRS, SAL):

