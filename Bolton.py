#!/usr/bin/env python3

C_to_K = 273.15
c_p_dry = 1004
c_v_dry = 1875
eps = 0.622
k_dry = 0.2854

def sat_vapor_pressure(T):
    e_s(T) =  6.112*np.exp(17.67T/(T+243.5))
    return e_s(T)

def sat_vapor_temperature(e_s):
    T = (243.5*np.log(e_s) - 440.8)/(19.48 - np.log(e_s)
    return T
def sat_mixing_ratio (p,t):
    e_s = 6.112*np.exp(17.67T/(T+243.5))
    w = eps(e_s/(p-e_s))
    return w
