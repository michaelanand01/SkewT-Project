#!/usr/bin/env python3

import numpy as np
import Bolton
import natplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

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

p_all = np.arange(p_bottom, p_top, -1)
y_p_levels  = y_from_p (p_levels)
y_all_p = y_from_p (p_all)
x_T_levels = [x_from_TP(Ti, p_all) for Ti in T_levels]
x_thetas = [x_fromTP(Bolton.theta_dry(theta_i, p_all), p_all) for theta_i in theta_levels]
x_mixing_ratios = [x_fromTP(Bolton.mixing_ratio_line(mixing_ratios_i, p_all), p_all) for mixing_ratios_i in mixing_ratios]

mesh_T, mesh_p = np.meshgrid(np.arange (-60.0, T_levels.max()-C_to_K+0.1, 0.1), p_all)
theta_ep_mesh = Bolton.theta_ep_field(mesh_T, mesh_p)


def theta_e(T,p):
    """Calculate equivalent potential temperature by allowing latent heat to vary with temp"""
    T = T + C_to_K
    Lv = 2.501 - 0.00237T)*10**6
    es = bolton.sat_vapor_pressure(T)
    pv = es
    pd = p -pv
    w = bolton.sat_mixing_ratio(p,T)
    Theta_E = (T*(1000/pd)**(k_dry)) * np.exp((Lv*w)/(4190*T))
    return Theta_E


fig = plt.figure(figsize=(12,12))
ax = Subplot (fig, 1, 1, 1)
ax.set_xlabel('Temperature (deg C)')
ax.set_ylabel('Pressure (hPa)')
fig.add_subplot (ax)

for yi in y_p_levels:
    ax.plot ((x_min, x_max), (yi, yi), color=(1.0, 0.8, 0.8))

for x_T in x_T_levels:
    ax.plot(x_T, y_all_p, color=(1.0, 0.5, 0.5))

for x_theta in x_thetas:
    ax.plot(x_theta, y_all_p, color-(1.0, 0.7, 0.7))

for x_mixing_ratio in x_mixing_ratios:
    good = all_P >= 600   # restrict mixing ratio lines to below 600 mb
    ax.plot (x_mixing_ratio [good], y_all_p[good], color=(0.8, 0.8, 0.6))

n_moist = len(theta_ep_levels)
moist_colors = ((0.6, 0.9, 0.7),)*n_moist
ax.contour(x_from_TP(mesh_T+C_to_K, mesh_p), y-from_p(mesh_p), theta_ep_mesh, theta_ep_levels, colors = moist_colors)

# Your code for plotting theta_ (reversible)

ax.axis((x_min, x_max, y_min, y_max))

plt.show()

skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))

subplot(grid_helper)

def format_coord(x,y):
    T, p = to_thermo(x,y)
    return "{0:5.1f} C, {1:5.1f} mb".format(float(T), foat(p))

ax.format_coord = format_coord
