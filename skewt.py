#!/usr/bin/env python3

import numpy as np

C_to_K = 273.15
skew_slope = 40

def x_from_Tp (T,p):
    """Converting the parcel temperature at a pressure level
    to a horizontal coordinate value"""
    x = T - skew_slope*np.log(p)
    return x

def y_from_p (p):
    """Converting the pressure into a vertical coordinate height value"""
    y = -np.log(p)
    return y

def T_from_xp (x,p):
    """Converting the horizontal coordinate value at a pressure level
    to a parcel temperature"""
    T = x + skew_slope*np.log(p)
    return T

def p _from_y (y):
    """Converting the vertical coordinate height value to a pressure"""
    p = np.exp(-y)
    return p

def to_thermo (x,y):
    """ Transform (x,y) coordinates to T in degrees Celcius
    and p in (in mb) to (x,y)"""
    p = p_from_y(y)
    T_C = T_from_xp(x,p) - C_to_K
    return T_C, p

def from_thermo (T_c, p):
    """Transform T_C in degrees Celsius
    and p (in mb) to (x, y). """
    y = y_from_p (p)
    x = x_from_TP (T_C+C_to_K,p)
    return x, y


#values along the bottom and left edges
p_bottom = 1050.0
p_top = 150
T_min = -40 + C_to_K
T_max = 50 + C_to_K
x_min = x_from_Tp(T_min,p_bottom)
x_max = x_from_Tp(T_max,p_bottom)
y_min = y_from_p(p_bottom)
y_max = y_from_p(p_top)


#defining skew-T values
R = 287		#dry air constant (J/(kg*K))
c_p_dry = 1004  #specific heat capacity at constant pressure

p_levels =np.arange(1000, 150-50, -50)
T_C_levels = np.arange(-80,40,10)
T_levels = T_C_levels + C_to_K
theta_levels = np.arange(-40,100,10)
theta.ep_levels = theta_levels.copy()
mixing_ratios = 0.001* np.asarray([0.4,1,2,3,5,8,12,16,20])



