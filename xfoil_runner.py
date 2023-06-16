"""
%--------------------------------------------------------------------------%
%           Aeroelastic Tailoring of Helicopter Main Rotor Blade           %
%                                                                          %
%           @author: Francisco Campos Moreira Rodrigues                    %
%           @date: October, 2022                                           %
%                                                                          %
%                       SCRIPT: xfoil_runner.py                            % 
%--------------------------------------------------------------------------%
"""

import os
import subprocess
import numpy as np
import math
import file_writer as fwrite
# from matplotlib import pyplot as plt

#%% XFOIL runner

def Xfoil (polar_data, nodal_data, nodes_U, nodes_L):
    subprocess.call("xfoil.exe < input_file.txt", shell=True)

    #%% DATA
    
    # Store Pressure coefficient
    nodal = np.loadtxt("{0}".format(nodal_data), skiprows=2)
    
    # # Store origin (0.0, 0.0) index
    # for i in range(len(nodal)):
    #     if nodal[i,0]==0.0:
    #         zero=i
    #         break

    # # Upper Surface
    # Cp_U = nodal[0:nodes_U-1,2]
    # Cp_U = np.append(Cp_U, nodal[zero,2])
            
    # # Lower Surface
    # Cp_L = nodal[-nodes_L:-1,2]
    # Cp_L = np.insert(Cp_L, 0, nodal[zero,2])

    # Upper Surface
    Cp_U = nodal[0:nodes_U,2]

    # Lower Surface
    Cp_L = nodal[-(nodes_L+1):-1,2]
           
    # Store Lift, Drag, Dp and Momentum coefficients
    coeff = np.loadtxt("{0}".format(polar_data), skiprows=12)
    if len(coeff) == 0:
        Cl = 'None'
        Cd = 'None'
    else:          
        Cl = coeff[1]
        Cd = coeff[2]
    return (Cp_U, Cp_L, Cl, Cd)