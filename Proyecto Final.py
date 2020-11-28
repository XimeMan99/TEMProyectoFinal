# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 14:03:15 2020

@author: sinap
"""

#H(x,z,t)=Hz+cos(377.14x)cos(40*109t-182.7z+)   k    [T]


import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d 
import numpy as np 
import matplotlib.pyplot as plt 
import math

#H = 100
E = 100
phi = 45
w = 40*np.pi*(10**9)
t = 0
xp = 2
yp = 2
zp = 2

#funcion que describe el plano, PLANO EN Z


def f(x,z,t):
    #fourier sierra
    #H= 2*(np.sin(377.25*x)-np.sin(377.25*x*2)/2+np.sin(377.25*x*3)/3-np.sin(377.25*x*4)/4+np.sin(377.25*x*5)/5-np.sin(377.25*x*6)/6+np.sin(377.25*x*7)/7-np.sin(377.25*x*8)/8+np.sin(377.25*x*9)/9-np.sin(377.25*x*10)/10) 
    #fourier cuadrada
    #H=(4/np.pi)*(np.sin(377.25*x)+np.sin(377.25*x*3)/3+np.sin(377.25*x*5)/5+np.sin(377.25*x*7)/7+np.sin(377.25*x*9)/9+np.sin(377.25*x*11)/11+np.sin(377.25*x*13)/13+np.sin(377.25*x*15)/15+np.sin(377.25*x*17)/17+np.sin(377.25*x*19)/19)
    #fourier triangular
    H=(8/np.pi**2)*(np.sin(377.25*x)-np.sin(377.25*x*3)/9+np.sin(377.25*x*5)/25-np.sin(377.25*x*7)/49+np.sin(377.25*x*9)/81-np.sin(377.25*x*11)/121+np.sin(377.25*x*13)/169-np.sin(377.25*x*15)/225+np.sin(377.25*x*17)/289-np.sin(377.25*x*19)/361)
    """
    #H_z
    H1 = H*(np.cos(377.25*x)*np.cos(7.4*np.pi*(10**9)*t)*np.exp(-369.2*z))
    H2 = H*(np.cos(377.25*x)*np.cos(14.8*np.pi*(10**9)*t)*np.exp(-343.95*z))
    H3 = H*(np.cos(377.25*x)*np.cos(22.2*np.pi*(10**9)*t)*np.exp(-297.13*z))
    H4 = H*(np.cos(377.25*x)*np.cos(29.6*np.pi*(10**9)*t)*np.exp(-215.07*z))
    H5 = H*(np.cos(377.25*x)*np.cos(37*np.pi*t*(10**9) - 88.19*z))
    H6 = H*(np.cos(377.25*x)*np.cos(44.4*np.pi*t*(10**9) - 271.69*z))
    H7 = H*(np.cos(377.25*x)*np.cos(51.8*np.pi*t*(10**9) - 389.7*z))
    H8 = H*(np.cos(377.25*x)*np.cos(59.2*np.pi*t*(10**9) - 491.86*z))
    H9 = H*(np.cos(377.25*x)*np.cos(66.6*np.pi*t*(10**9) - 486.51*z))
    H10 = H*(np.cos(377.25*x)*np.cos(74*np.pi*t*(10**9) - 676.8*z))
    """
    """
    #estos son los E_y pero me da hueva hacer otra funcion para estos
    H1 = 77.46 *H*(np.sin(377.25*x)*np.sin(7.4*np.pi*(10**9)*t)*np.exp(-369.2*z))
    H2 = 154.92*H*(np.sin(377.25*x)*np.sin(14.8*np.pi*(10**9)*t)*np.exp(-343.95*z))
    H3 = 232.39*H*(np.sin(377.25*x)*np.sin(22.2*np.pi*(10**9)*t)*np.exp(-297.13*z))
    H4 = 309.85*H*(np.sin(377.25*x)*np.sin(29.6*np.pi*(10**9)*t)*np.exp(-215.07*z))
    H5 = 387.31*H*(np.sin(377.25*x)*np.sin(37*np.pi*t*(10**9) - 88.19*z))
    H6 = 464.77*H*(np.sin(377.25*x)*np.sin(44.4*np.pi*t*(10**9) - 271.69*z))
    H7 = 542.23*H*(np.sin(377.25*x)*np.sin(51.8*np.pi*t*(10**9) - 389.7*z))
    H8 = 619.19*H*(np.sin(377.25*x)*np.sin(59.2*np.pi*t*(10**9) - 491.86*z))
    H9 = 697.16*H*(np.sin(377.25*x)*np.sin(66.6*np.pi*t*(10**9) - 486.51*z))
    H10 = 774.62*H*(np.sin(377.25*x)*np.sin(74*np.pi*t*(10**9) - 676.8*z))
    """
    
    #H_x
    H1 = 0.98 *H*(np.sin(377.25*x)*np.cos(7.4*np.pi*(10**9)*t)*np.exp(-369.2*z))
    H2 = 0.912*H*(np.sin(377.25*x)*np.cos(14.8*np.pi*(10**9)*t)*np.exp(-343.95*z))
    H3 = 0.79*H*(np.sin(377.25*x)*np.cos(22.2*np.pi*(10**9)*t)*np.exp(-297.13*z))
    H4 = 0.57*H*(np.sin(377.25*x)*np.cos(29.6*np.pi*(10**9)*t)*np.exp(-215.07*z))
    H5 = 0.234*H*(np.sin(377.25*x)*np.sin(37*np.pi*t*(10**9) - 88.19*z))
    H6 = -0.72*H*(np.sin(377.25*x)*np.sin(44.4*np.pi*t*(10**9) - 271.69*z))
    H7 = -1.033*H*(np.sin(377.25*x)*np.sin(51.8*np.pi*t*(10**9) - 389.7*z))
    H8 = -1.304*H*(np.sin(377.25*x)*np.sin(59.2*np.pi*t*(10**9) - 491.86*z))
    H9 = -1.56*H*(np.sin(377.25*x)*np.sin(66.6*np.pi*t*(10**9) - 486.51*z))
    H10 = -1.79*H*(np.sin(377.25*x)*np.sin(74*np.pi*t*(10**9) - 676.8*z))
    
    
    return H1+H2+H3+H4+H5+H6+H7+H8+H9+H10

#definimos los otros planos de los cuales se utilizan las variables
x = np.linspace(0, 30, 30)
z = np.linspace(0, 30, 30)
t = np.linspace(0, 30, 30)

#hacemos meshgrid de x y y
X, Z = np.meshgrid(x, z)
#hacemos meshgrid para la funcion de superficie Z
H = f(X,Z,t) #=Z

#ploteamos y embellecemos grafica
ax = plt.axes(projection='3d')
ax.plot_surface(X, H, Z, rstride=1, cstride=1,cmap='viridis', edgecolor='none')
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
ax.set_title('$H_z(x,y,z,t)$')
ax.view_init(50, 70)




"""
MODO TE21 
-para onda f>fc
#H_z = H*np.cos(337.44*x)*np.cos(168.7*y)*np.cos(w*t-182.67*z+phi)
#H_x = -0.433*H*np.sin(337.44*x)*np.cos(168.7*y)*np.sin(w*t-182.67*z+phi)
#H_y = -0.217*H*np.cos(337.44*x)*np.sin(168.7*y)*np.sin(w*t-182.67*z+phi)
#E_x = -187.2*H*np.cos(337.44*x)*np.sin(168.7*y)*np.sin(w*t-182.67*z+phi)
#E_y = -374.4*H*np.sin(337.44*x)*np.cos(168.7*y)*np.sin(w*t-182.67*z+phi)

MODO TM11
-para onda f>fc
#E_z = E*np.sin(266.7*x)*np.sin(266.7*y)*np.cos(w*t-182.9*z+phi)
#E_x = 0.34*E*np.cos(266.7*x)*np.sin(266.7*y)*np.sin(w*t-182.9*z+phi)
#E_y = 0.34*E*np.sin(266.7*x)*np.cos(266.7*y)*np.sin(w*t-182.9*z+phi)
#H_x = -(2*10**-3)*E*np.sin(266.7*x)*np.cos(266.7*y)*np.sin(w*t-182.9*z+phi)
#H_y = (2*10**-3)*E*np.cos(266.7*x)*np.sin(266.7*y)*np.sin(w*t-182.9*z+phi)

MODO TE10 

-para onda f>fc
#H_z = H*np.cos(377.14*x) * np.cos(w*t - 182.7*zp + phi)
#E_y = 418.7*H*np.sin(377.14*x) * np.sin(w*t - 182.7*z + phi)
#H_x = -0.485*H*np.sin(377.14*x)*np.sin(w*t-182.7*z+phi)

-para onda f<fc
#H_z = H*np.cos(377.14*x) * np.cos(w*t) * math.exp(-313.68*zp)
"""




