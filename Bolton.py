#!/usr/bin/env python3
C_to_K = 273.15
c_p_dry = 1005.7
c_v_dry = 719
eps = 0.622
k_dry = 0.2854
import numpy as np

def sat_vapor_pressure(T):
    """Calcuates the saturation vapor pressure"""
    e_s =  6.112*np.exp(17.67*T/(T+243.5))
    return e_s

def sat_vapor_temperature(e_s):
    """Calculates the saturation vapor temperature"""
    T = (243.5*np.log(e_s) - 440.8)/(19.48 - np.log(e_s))
    return T

def sat_mixing_ratio (p,T):
    """Calculates the saturation mixing ratio"""
    e_s = sat_vapor_pressure(T)
    w_s = eps*(e_s/(p-e_s))
    return w_s

def mixing_ratio_line(p, w_s):
    """ determines the position of the mixing ratio line"""
    e_s = (w_s*p)/(w_s+ eps)
    T = sat_vapor_temperature(e_s)
    return T

def RH(T, p, w):
    """Calculates the relative humidity of a parcel"""
    w_s = sat_mixing_ratio (p,T)
    RH = (w/w_s)*100
    return RH

def T_LCL(T,rh):
    """Finds the temperature of the LCL in Kelvin"""
    T_K = T + C_to_K
    Denom_1 = 1/((T_K)-55)
    Denom_2 = np.log(rh/100)/2840
    T_LCL =(1/(Denom_1 - Denom_2)) + 55
    return T_LCL

def theta_dry(theta, p, p_0=1000.0):
    """Calculates the dry adiabats"""
    theta_dry = theta*(p/p_0)**(k_dry)
    return theta_dry

def pseudoeq_potential_T(T, p, w, p_0=1000.0):
    """Calcuates the pseudoadiabatic equivalent potential temperature by taking account the
    moisture content in the air"""
    T_K = T + C_to_K
    rh = RH(T,p,w)
    Temp_LCL = T_LCL(T,rh)
    term1 = (T_K*(1000/p)**k_dry)
    term2 = (3.376/Temp_LCL)-0.00254
    term3 = w*(1+(0.81*10**-3)*w)
    term4 = np.exp(term2*term3)
    pseudo_adia_equiv_pot_T = term1 * term4
    return pseudo_adia_equiv_pot_T

def theta_ep_field (T, p, p_0=1000.0):
    """Determine and calculate the moist adiabat field"""
    w_s = sat_mixing_ratio (p,T)
    theta_ep = pseudoeq_potential_T(T, p, w_s, p_0=1000.0)
    return theta_ep

