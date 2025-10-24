"""
Author : Julien Dandoy
Date : 18/10/2025
Description : Standard atmosphere model
"""

import numpy as np

def stdatm(z):
    """
    Standard atmosphere model from sea-leavel up to 20 km

    z : Altitude [m]
    """

    # Constants
    g  = 9.81     # [m/s^2]
    R  = 287.0529 # [J/kg/K] 

    # Model parameter
    lbda = -6.5e-3 # [K/m]

    # Model outputs
    if z >= 0 and z < 11000:

        # Reference conditions
        p0   = 101325 # [Pa]
        T0   = 288.15 # [K]
        rho0 = 1.225  # [kg/m^3]

        T   = (1 + lbda*z/T0) * T0
        p   = pow(T/T0, -g/(R*lbda)) * p0
        rho = pow(T/T0, -g/(R*lbda) -1) * rho0
    
    elif z >= 11000 and z <= 20000:
        
        # Reference conditions at 11 km
        p11   = 22632.1 # [Pa]
        T11   = 216.65  # [K]
        rho11 = 0.36391 # [kg/m^3]

        T   = T11
        p   = np.exp(-g*(z - 11000)/(R*T11)) * p11
        rho = np.exp(-g*(z - 11000)/(R*T11)) * rho11

    else:
        raise ValueError("Altitude z is out of range [0, 20000] m")

    return p, T, rho
