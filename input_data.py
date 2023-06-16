"""
%--------------------------------------------------------------------------%
%           Aeroelastic Tailoring of Helicopter Main Rotor Blade           %
%                                                                          %
%           @author: Francisco Campos Moreira Rodrigues                    %
%           @date: October, 2022                                           %
%                                                                          %
%                       SCRIPT: input_data.py                              % 
%--------------------------------------------------------------------------%
"""

import numpy as np
from math import sin, cos, atan
# from gekko import GEKKO

def data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air):
    
    y = (R_min+L*i) #radial position [mm]
    T = W #[N]
    V_inf = miu_air*V_tip #free-stream velocity [mm/s]
    V_C = 0 # Climb velocity [mm/s]
    
    U_T = rotor_speed*y + miu_air*V_tip*sin(azimut) #in-plane velocity [mm/s] [mm/s]
    if miu_air == 0:
        v_i = np.sqrt(T/(2*rho*A)) #induced velocity (HF) [mm/s]
    else:
        v_i = T/(2*rho*A*V_inf) #induced velocity (FF) [mm/s]
    U_P = V_C+v_i #[mm/s] #out-of-plane velocity [mm/s]
    
    #  Total velocity
    U = np.sqrt(U_P**2+U_T**2) #[mm/s]
    
    # Induced angle
    phi = np.arctan(U_P/U_T)*180/np.pi #[degrees]
    
    # Reynolds & Mach numbers
    Reynolds = (U*c)/kin_viscosity
    # Ma = U/sound
    Ma = 0.1
    
    return(Reynolds, Ma, U, phi)
    
    

