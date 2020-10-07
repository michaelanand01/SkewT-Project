#!/usr/bin/env python3
C_to_K = 273.15
c_p_dry = 1005.7
c_v_dry = 719
eps = 0.622
k_dry = 0.2854

def sat_vapor_pressure(T):
    """Calcuates the saturation vapor pressure"""
    e_s(T) =  6.112*np.exp(17.67T/(T+243.5))
    return e_s(T)

def sat_vapor_temperature(e_s):
   """Calculates the saturation vapor temperature"""
    T = (243.5*np.log(e_s) - 440.8)/(19.48 - np.log(e_s)
    return T

def sat_mixing_ratio (p,T):
    """Calculates the saturation mixing ratio"""
    e_s = 6.112*np.exp(17.67T/(T+243.5))
    w_s = eps(e_s/(p-e_s))
    return w_s

def mixing_ratio_line(p, w_s):
    """ determines the position of the mixing ratio line"""
    return T

def RH(T, p, w):
    """Calculates the relative humidity of a parcel"""
    w_s = sat_mixing_ratio (p,T)
    RH = (w/w_s)*100
    return RH

def T_LCL(T,RH):
    """Finds the temperature of the LCL in Kelvin"""
    Denom_1 = 1/((T)-55)
    Denom_2 = (ln(RH/100))/2840
    T_LCL = 1/(Denom_1 - Denom_2)
    return T_LCL
