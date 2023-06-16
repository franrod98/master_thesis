"""
%--------------------------------------------------------------------------%
%           Aeroelastic Tailoring of Helicopter Main Rotor Blade           %
%                                                                          %
%           @author: Francisco Campos Moreira Rodrigues                    %
%           @date: October, 2022                                           %
%                                                                          %
%                            MAIN SCRIPT                                   % 
%--------------------------------------------------------------------------%
"""

import os
import sys
import time
import tkinter as tk
from tkinter import filedialog
import numpy as np
from math import atan
import input_data
import organize_txt_TESTE
import file_writer as fwrite
import xfoil_runner as xfrun
import transformation

#%% CLASSES

# Material
class material:
    def __init__(self, name, rho, E, poisson):
        self.name = name
        self.rho = rho
        self.E = E
        self.poisson = poisson

# Angles of interest
class angles:
    def __init__(self, AoA0, AoA1, initial):
        self.zero = AoA0
        self.one = AoA1
        self.initial = initial

# Convergence
class convergence:
    def __init__(self, residual, max_res, iteration):
        self.res = residual
        self.max_res = max_res
        self.iter = iteration
        
#%% INPUT DATA

# Blade data
m = 8322.4 #[kg]
W = m*9.81 #[N]
rotor_speed = 27 #[rad/s]
Nb = 4
c = 527.3 #[mm]
Elem = 50
R_max = 8177.8 #[mm]
R_min = 0.31*R_max #[mm]
L = (R_max-R_min)/Elem #[mm]
A = np.pi*(R_max**2) #[mm**2]
pitch0 = 10.4360 #forward flight [degrees]
azimut = 0*np.pi/180 #forward flight [rad]

# Flow conditions
kin_viscosity = 15.16 #[(mm**2)/s]
miu_air = 0.17 # advance ratio
sound = 343E3 #[mm/s]
rho = 1.204E-12 #[tonne/(mm**3)]
V_tip = rotor_speed*R_max
magnitude = -0.5*rho
P_atm = -0.101325 #[MPa]

# Material properties
ny = material ('Nylon', 1.18E-9, 3E3, 0.41) #Nylon
ti = material('Titanium', 4.43E-9, 114E3, 0.34) #Ti-6A1-4V

# Material lists

mat_structures = [ny, ti, ti, ti, ti, ti, ti, ti, ti, ti, ti, ti, ti, ti, ti]
# mat_structures = [ny]
mat_airfoil = [ti]

# Xfoil
iter_xf = 500

# FEA simulation
step = 'quad_sparX'
job0 = 'CONVG'
my_dir = 'C:/Users/mario/Documents/TÉCNICO/5º ANO- 2º Semestre/TESE/CÓDIGO/ANACONDA/abaqus_testing'
OPEN = 1
structures = len(mat_structures)
seeds = 15
interactions = True

root = tk.Tk()
root.withdraw()

if interactions == True:
    manual=input('Manual definition of contact pairs? (Y/N) ')
    if manual.islower() == True:
        manual = manual.swapcase()
    if manual!='Y' and manual!='N':
        while manual !='Y' and manual !='N':
            manual=input('Error: Choose "Y" for manual selection or "N" for automatic selection. ')
            if manual.islower() == True:
                manual = manual.swapcase()
    if manual=='N':
        input('\nChoose file for surface selection  ')
        surf_path = filedialog.askopenfilename()
        surf_path = surf_path.replace(my_dir+'/','')
        
start = time.time()

#%% Delete all existing files

for fname in os.listdir(my_dir):
    if fname.startswith('COORDS') or fname.startswith('EDGE') or fname.startswith('MAP'):
        os.remove(os.path.join(my_dir, fname))
        print('\nFile {0} deleted.'.format(fname))

#%% Initialize project variables

convg = convergence (1, 1, 1)
pitch = angles ([], [], pitch0) #degrees (null alfa_tip)
alpha = angles ([], [], [])
b = angles ([], [], [])
ind_angle = angles ([], [], [])
pitch_diff = []
theta_T = [] # Torsion deflection
theta_T_norm = []

variable_U = []
variable_L = []
U = []

CL0 = []
CD0 = []
Lift0 = 0
Drag0 = 0
CL1 = []
CD1 = []
Lift1 = 0
Drag1 = 0

#================================================ CREATE MODEL ========================================================

#%% Initialize .cae file

os.system('cmd /k "abaqus cae noGui=NACA0012_start.py -- {0} {1} {2} {3} {4}"\n'.format(step, structures, R_min, L, Elem))
    
#%% Define materials

# Parts excluding blade skin
for i in range(0, structures):
    os.system('cmd /k "abaqus cae noGui=NACA0012_structures_mat.py -- {0} {1} {2} {3} {4}"\n'.format(i, mat_structures[i].name, mat_structures[i].rho, mat_structures[i].E, mat_structures[i].poisson))

# Blade skin
os.system('cmd /k "abaqus cae noGui=NACA0012_airfoil_mat.py -- {0} {1} {2} {3}"\n'.format(mat_airfoil[0].name, mat_airfoil[0].rho, mat_airfoil[0].E, mat_airfoil[0].poisson))

#%% Select surface pairs

if interactions is True:
    if manual == 'N':
        os.system('cmd /k "abaqus cae noGui={0}"'.format(surf_path))
    else:
        print('\nSelecting contact pairs...')
        input('<Press enter to continue> ')

#%% Create blade elements

# partitions.blade_elements (R_min, Elem, L, c)
os.system('cmd /k "abaqus cae noGui=NACA0012_elements.py -- {0} {1} {2} {3}"\n'.format(R_min, L, Elem, c))

#%% Define interactions

os.system('cmd /k "abaqus cae noGui=NACA0012_interactions.py')

#%% Mesh part (to extract airfoil profile)

mesh0 = time.time()
os.system('cmd /k "abaqus cae noGui=NACA0012_mesh.py -- {0} {1} {2} {3}"\n'.format(structures, Elem, seeds, c))
mesh1 = time.time()
ty_res = time.gmtime(mesh1-mesh0)
res = time.strftime("%H:%M:%S",ty_res)
print('\nMeshing time: {0}\n'.format(res))

#==================================================== ANALYSIS ============================================================

#%% Initialize files 

mapped_U=open('MAPU.txt', 'w')
mapped_L=open('MAPL.txt', 'w')

index_U = open('INDEX_U.txt', 'w')
text_U = 'NODE INDEXES (UPPER SURFACE)'
index_L = open('INDEX_L.txt', 'w')
text_L = 'NODE INDEXES (LOWER SURFACE)'

for i in range(0, Elem+1):

    #%% Organize data in accordance with XFoil
    
    coord_U = 'EDGE{0}_U'.format(i)
    coord_L = 'EDGE{0}_L'.format(i)
    
    # Sort coordinate vectores and store original indexes
    sort_INDEX_U, sort_INDEX_L, x_U, y_U, z_U, x_L, y_L, z_L = organize_txt_TESTE.organize_data (coord_U, coord_L)
    x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted = organize_txt_TESTE.sort_vectors (sort_INDEX_U, sort_INDEX_L, x_U, y_U, z_U, x_L, y_L, z_L)
    
    # Define b0
    if i == 0:
        b0 = atan(z_U_sorted[0]/x_U_sorted[0])
        b.initial = b0

    # b0 = 0
    # b.initial = b0    

    coords_U_sorted = np.vstack((x_U_sorted, y_U_sorted, z_U_sorted)).T
    coords_L_sorted = np.vstack((x_L_sorted, y_L_sorted, z_L_sorted)).T
    
    # Upper Surface
    if any(n<0 for n in coords_U_sorted[:,0]):
        count = 0
        for j in coords_U_sorted:
            if round(j[0],6) == 0:
                extract = coords_U_sorted[np.r_[count:len(coords_U_sorted)],:]
                extract_index = sort_INDEX_U[np.r_[count:len(sort_INDEX_U)]]
                order = np.argsort(extract[:,2])[::-1]
                for k in order:
                    coords_U_sorted[count,:] = extract[k,:]
                    sort_INDEX_U[count] = extract_index[k]
                    count = count + 1
                break
            count = count + 1
            
    # Lower Surface
    if any(n<0 for n in coords_L_sorted[:,0]):
        count = 0
        pos = 0
        for j in coords_L_sorted:
            if round(j[0],6) == 0:
                extract = coords_L_sorted[np.r_[0:count+1],:]
                extract_index = sort_INDEX_L[np.r_[0:count+1]]
                order = np.argsort(extract[:,2])[::-1]
                for k in order:
                    coords_L_sorted[pos,:] = extract[k,:]
                    sort_INDEX_L[pos] = extract_index[k]
                    pos = pos + 1
                break
            count = count + 1   
            
    coords_U_local = np.zeros((len(coords_U_sorted),3))
    coords_L_local = np.zeros((len(coords_L_sorted),3))
    for j in range(0,len(coords_U_sorted)):
        coords_U_local[j,:] = np.add(coords_U_sorted[j,:], -coords_U_sorted[-1,:])       
    for j in range(0,len(coords_L_sorted)):
        coords_L_local[j,:] = np.add(coords_L_sorted[j,:], -coords_U_sorted[-1,:])
    _, bi, _, _ = transformation.find_angles (coords_U_local)
    if bi<b.initial:
        Ry = transformation.T_matrix(-(b.initial-bi))
        theta = pitch.initial+(b.initial-bi)*180/np.pi          
    else:
        Ry = transformation.T_matrix(bi-b.initial)
        theta = pitch.initial-(bi-b.initial)*180/np.pi
    b.zero = np.append(b.zero, bi)
    pitch.zero = np.append(pitch.zero, theta)
    coords_U_local = transformation.local_coords (Ry, coords_U_local)
    coords_L_local = transformation.local_coords (Ry, coords_L_local)
    
    X_U_optim = coords_U_local[:,0]
    Y_U_optim = coords_U_local[:,2]
    X_L_optim = coords_L_local[:,0]
    Y_L_optim = coords_L_local[:,2]
    
    # Interpolate airfoil
    if OPEN == 1:
        X_L_final, Y_L_final, _, _, _  = organize_txt_TESTE.close_airfoil (X_U_optim, Y_U_optim, X_L_optim, Y_L_optim)
    else:
        X_L_final = X_L_optim
        Y_L_final = Y_L_optim
        
    nodes_U = len(X_U_optim)
    nodes_L = len(X_L_optim)
     
    # Write file containing node coordinates
    organize_txt_TESTE.XFOIL_file (X_U_optim, Y_U_optim, X_L_final, Y_L_final, 'COORDS{0}'.format(i)) 

    #%% Evaluate pressure (Cp) and lift (Cl) coefficients from Xfoil
    
    # Input data for XFoil
    Re, Ma, Ui, phi = input_data.data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air)
    U = np.append(U, Ui)
    AoA = theta-phi
    alpha.zero = np.append(alpha.zero, AoA)
    ind_angle.zero = np.append(ind_angle.zero, phi)
    
    # Write input file for XFoil
    airfoil = 'COORDS{0}'.format(i)
    polar_data, nodal_data = fwrite.write (round(AoA,5), Re, Ma, airfoil, iter_xf)
    
    # Run XFoil  
    Cp_U, Cp_L, Cl, Cd = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
    
    # Correct non-converging profiles
    if Cl == 'None':
        Cl0 = Cl
        inc = 0
        AoA0 = AoA
        while Cl0 == 'None':
            # Lower limit
            inc = inc-0.01
            AoA = AoA0+inc
            polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
            Cp_U_lower, Cp_L_lower, Cl0, Cd0 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)            
        print('\nCOORDS{0}: Alfa = {1} --> Cl = {2}'.format(i, AoA, Cl0))
        # Upper limit
        Cl1 = 'None'
        inc = 0
        # AoA = AoA+inc*-1
        AoA = AoA0
        while Cl1 == 'None':
            inc = inc+0.01
            AoA = AoA0+inc
            polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
            Cp_U_upper, Cp_L_upper, Cl1, Cd1 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
        print('             Alfa = {0} --> Cl = {1}'.format(AoA, Cl1))
        
        # Obtain average value for Cp and Cl
        for j in range(0,len(Cp_U)):
            Cp_U[j] = (Cp_U_upper[j]+Cp_U_lower[j])/2
        for j in range(0,len(Cp_L)):
            Cp_L[j] = (Cp_L_upper[j]+Cp_L_lower[j])/2
        Cl = (Cl0+Cl1)/2
        Cd = (Cd0+Cd1)/2
        print('\nCl average = {0}'.format(Cl))
    
    # Store pressure coefficiente and fluid velocity
    variable_U[:] = Cp_U*(Ui**2)*magnitude+P_atm
    variable_L[:] = Cp_L*(Ui**2)*magnitude+P_atm
    
    # Store Lift and Drag coefficients
    CL0 = np.append(CL0, Cl)
    CD0 = np.append(CD0, Cd)
    
    #%% Create mapped field for distributed pressure
    
    L0 = L*i # position of blade section plane
    
    # Upper Surface
    for j in range(0, nodes_U):       
        Cp_U_0 = variable_U[j]     
        X_U_0 = coords_U_sorted[j,0]*c
        Y_U_0 = coords_U_sorted[j,2]*c
        mapped_U.write('{0}, {1}, {2}, {3}\n'.format(X_U_0, L0, Y_U_0, Cp_U_0))
        
    # Lower Surface
    for j in range(0, nodes_L):       
        Cp_L_0 = variable_L[j]      
        X_L_0 = coords_L_sorted[j,0]*c
        Y_L_0 = coords_L_sorted[j,2]*c    
        mapped_L.write('{0}, {1}, {2}, {3}\n'.format(X_L_0, L0, Y_L_0, Cp_L_0))
        
    #%% Store original indexes
    
    # Upper Surface
    output = '{0}'.format(i)+' '
    for j in sort_INDEX_U:
        output = output+str(j)+' '
    text_U ='\n'.join((text_U, output))
          
    # Lower Surface
    output = '{0}'.format(i)+' '
    for j in sort_INDEX_L:
        output = output+str(j)+' '
    text_L ='\n'.join((text_L, output))

#%% Close files

# Index files
index_U.write(text_U)
index_U.close()
index_L.write(text_L)
index_L.close()

# Mapped fields
mapped_U.close()
mapped_L.close()
 
#%% Apply distributed pressure

os.system('cmd /k "abaqus cae noGui=NACA0012_Pressure_test.py -- {0} {1} {2} {3} {4} {5} {6}"\n'.format(rotor_speed, P_atm, R_min, convg.iter, 'NACA0012_analysis', c, magnitude))

#%% Execute job

job = job0 + '-{0}'.format(convg.iter)
for fname in os.listdir(my_dir):
    if fname.startswith(job) and fname.endswith('.step'):
        continue
    elif fname.startswith(job):
        os.remove(os.path.join(my_dir, fname))
        print('\nFile {0} deleted.'.format(fname))

jobTime0 = time.time()
os.system('cmd /k "abaqus cae noGui=NACA0012_analysis_TEST.py -- {0} {1}"\n'.format('NACA0012_analysis',job))
if os.path.exists("{0}.odb".format(job)):
    print('\n\nABAQUS/CAE odb file created for iteration {0}.\n'.format(convg.iter))
else:
    print('\n\nERROR: ABAQUS/CAE odb file not created for iteration {0}.\n'.format(convg.iter))
    sys.exit()
jobTime1 = time.time()
ty_res = time.gmtime(jobTime1-jobTime0)
res = time.strftime("%H:%M:%S",ty_res)
print('\nRunning job time: {0}\n'.format(res))

#%% Update files containing edge nodes

os.system('cmd /k "abaqus cae noGui=NACA0012_post_processing.py -- {0} {1} {2} {3}"\n'.format(job, convg.iter, Elem, c))

#%% Organize data

for i in range(0,Elem+1):
    coord_U = 'EDGE{0}_U'.format(i)
    coord_L = 'EDGE{0}_L'.format(i)
    index_U, index_L = organize_txt_TESTE.get_index (i)    
    _, _, x_U, y_U, z_U, x_L, y_L, z_L = organize_txt_TESTE.organize_data (coord_U, coord_L)
    x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted = organize_txt_TESTE.sort_vectors (index_U, index_L, x_U, y_U, z_U, x_L, y_L, z_L)
    coords_U_sorted = np.vstack((x_U_sorted, y_U_sorted, z_U_sorted)).T
    coords_L_sorted = np.vstack((x_L_sorted, y_L_sorted, z_L_sorted)).T
    coords_U_local = np.zeros((len(coords_U_sorted),3))
    coords_L_local = np.zeros((len(coords_L_sorted),3))
    for j in range(0,len(coords_U_sorted)):
        coords_U_local[j,:] = np.add(coords_U_sorted[j,:], -coords_U_sorted[-1,:])       
    for j in range(0,len(coords_L_sorted)):
        coords_L_local[j,:] = np.add(coords_L_sorted[j,:], -coords_U_sorted[-1,:])
    _, bi, _, _ = transformation.find_angles (coords_U_local)
    if bi<b.initial:
        Ry = transformation.T_matrix(-(b.initial-bi))
        theta = pitch.initial+(b.initial-bi)*180/np.pi          
    else:
        Ry = transformation.T_matrix(bi-b.initial)
        theta = pitch.initial-(bi-b.initial)*180/np.pi
    b.one = np.append(b.one, bi)    
    pitch.one = np.append(pitch.one, theta)

    # Input data for XFoil
    Re, Ma, Ui, phi = input_data.data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air)
    U = np.append(U, Ui)
    AoA = theta-phi
    alpha.one = np.append(alpha.one, AoA)
    ind_angle.one = np.append(ind_angle.one, phi)

# Evaluate residual
for i in range(0, Elem+1):
    pitch_diff = np.append(pitch_diff, pitch.one[i]-pitch.zero[i])
convg.res = max(pitch_diff, key=abs)
theta_T = np.append(theta_T, convg.res)

while convg.res > 0.0035:
        
    alpha.zero = []
    alpha.one= []
    pitch.zero = []
    pitch.one = []
    ind_angle.zero = []
    ind_angle.one = []
    
    pitch_diff = []
    convg.res = []
    
    variable_U = []
    variable_L = []
    U = []
    
    CL0 = []
    CD0 = []
    Lift0 = 0
    Drag0= 0
    CL1 = []
    CD1 = []
    Lift1 = 0
    Drag1 = 0
    
           
    #%% Initialize files for mapped analytical fields
    
    mapped_U=open('MAPU.txt', 'w')
    mapped_L=open('MAPL.txt', 'w')
    
    #%% Organize data
    
    for i in range(0,Elem+1):
        
        for fname in os.listdir(my_dir):
            if fname.startswith('COORDS'):
                os.remove(os.path.join(my_dir, fname))
                print('\nFile {0} deleted.'.format(fname))
        
        coord_U = 'EDGE{0}_U'.format(i)
        coord_L = 'EDGE{0}_L'.format(i)
        index_U, index_L = organize_txt_TESTE.get_index (i)    
        _, _, x_U, y_U, z_U, x_L, y_L, z_L = organize_txt_TESTE.organize_data (coord_U, coord_L)
        x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted = organize_txt_TESTE.sort_vectors (index_U, index_L, x_U, y_U, z_U, x_L, y_L, z_L)
        
            # Define b0
        if i == 0:
            b0 = atan(z_U_sorted[0]/x_U_sorted[0])
            b.initial = b0
        
        coords_U_sorted = np.vstack((x_U_sorted, y_U_sorted, z_U_sorted)).T
        coords_L_sorted = np.vstack((x_L_sorted, y_L_sorted, z_L_sorted)).T
        coords_U_local = np.zeros((len(coords_U_sorted),3))
        coords_L_local = np.zeros((len(coords_L_sorted),3))
        for j in range(0,len(coords_U_sorted)):
            coords_U_local[j,:] = np.add(coords_U_sorted[j,:], -coords_U_sorted[-1,:])       
        for j in range(0,len(coords_L_sorted)):
            coords_L_local[j,:] = np.add(coords_L_sorted[j,:], -coords_U_sorted[-1,:])
        _, bi, _, _ = transformation.find_angles (coords_U_local)
        if bi<b.initial:
            Ry = transformation.T_matrix(-(b.initial-bi))
            theta = pitch.initial+(b.initial-bi)*180/np.pi           
        else:
            Ry = transformation.T_matrix(bi-b.initial)
            theta = pitch.initial-(bi-b.initial)*180/np.pi
            
        b.zero = np.append(b.zero, bi)
        pitch.zero = np.append(pitch.zero, theta)
        coords_U_local = transformation.local_coords (Ry, coords_U_local)
        coords_L_local = transformation.local_coords (Ry, coords_L_local)
        X_U_optim = coords_U_local[:,0]
        Y_U_optim = coords_U_local[:,2]
        X_L_optim = coords_L_local[:,0]
        Y_L_optim = coords_L_local[:,2]   
        
        # Interpolate airfoil
        if OPEN == 1:
            X_L_final, Y_L_final, _, _, _  = organize_txt_TESTE.close_airfoil (X_U_optim, Y_U_optim, X_L_optim, Y_L_optim)
        else:
            X_L_final = X_L_optim
            Y_L_final = Y_L_optim
            
        nodes_U = len(X_U_optim)
        nodes_L = len(X_L_optim)
        
        organize_txt_TESTE.XFOIL_file (X_U_optim, Y_U_optim, X_L_final, Y_L_final, 'COORDS{0}'.format(i))
        
        # Input data for XFoil
        Re, Ma, Ui, phi = input_data.data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air)
        U = np.append(U, Ui)
        AoA = theta-phi
        alpha.zero = np.append(alpha.zero, AoA)
        ind_angle.zero = np.append(ind_angle.zero, phi)
        
        # Write input file for XFoil
        airfoil = 'COORDS{0}'.format(i)
        polar_data, nodal_data = fwrite.write (round(AoA,5), Re, Ma, airfoil, iter_xf)
        
        # Run XFoil  
        Cp_U, Cp_L, Cl, Cd = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
        
        # Correct non-converging profiles
        if Cl == 'None':
            Cl0 = Cl
            inc = 0
            AoA0 = AoA
            while Cl0 == 'None':
                # Lower limit
                inc = inc-0.01
                AoA = AoA0+inc
                polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                Cp_U_lower, Cp_L_lower, Cl0, Cd0 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)            
            print('\nCOORDS{0}: Alfa = {1} --> Cl = {2}'.format(i, AoA, Cl0))
            # Upper limit
            Cl1 = 'None'
            inc = 0
            # AoA = AoA+inc*-1
            AoA = AoA0
            while Cl1 == 'None':
                inc = inc+0.01
                AoA = AoA0+inc
                polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                Cp_U_upper, Cp_L_upper, Cl1, Cd1 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
            print('             Alfa = {0} --> Cl = {1}'.format(AoA, Cl1))
            
            # Obtain average value for Cp and Cl
            for j in range(0,len(Cp_U)):
                Cp_U[j] = (Cp_U_upper[j]+Cp_U_lower[j])/2
            for j in range(0,len(Cp_L)):
                Cp_L[j] = (Cp_L_upper[j]+Cp_L_lower[j])/2
            Cl = (Cl0+Cl1)/2
            Cd = (Cd0+Cd1)/2
            print('\nCl average = {0}'.format(Cl))
        
        # Store pressure coefficiente and fluid velocity
        variable_U[:] = Cp_U*(Ui**2)*magnitude+P_atm
        variable_L[:] = Cp_L*(Ui**2)*magnitude+P_atm
        
        # Store Lift coefficients
        CL1 = np.append(CL1, Cl)
        CD1 = np.append(CD1, Cd)
        
        #%% Create mapped field for distributed pressure
        
        L0 = L*i
        
        # Upper Surface
        for j in range(0, nodes_U):       
            Cp_U_0 = variable_U[j]     
            X_U_0 = coords_U_sorted[j,0]*c
            Y_U_0 = coords_U_sorted[j,2]*c
            mapped_U.write('{0}, {1}, {2}, {3}\n'.format(X_U_0, L0, Y_U_0, Cp_U_0))
            
        # Lower Surface
        for j in range(0, nodes_L):       
            Cp_L_0 = variable_L[j]      
            X_L_0 = coords_L_sorted[j,0]*c
            Y_L_0 = coords_L_sorted[j,2]*c    
            mapped_L.write('{0}, {1}, {2}, {3}\n'.format(X_L_0, L0, Y_L_0, Cp_L_0))
    
    #%% Close files
    
    mapped_U.close()
    mapped_L.close()
        
    #%% Restart analysis
    
    convg.iter = convg.iter+1
    job = job0 + '-{0}'.format(convg.iter)
    for fname in os.listdir(my_dir):
        if fname.startswith(job):
            os.remove(os.path.join(my_dir, fname))
            print('\nFile {0} deleted.'.format(fname))
    os.system('cmd /k "abaqus cae noGui=NACA0012_restart_new.py -- {0} {1} {2} {3} {4} {5} {6}"\n'.format(convg.iter, job0, job, P_atm, rotor_speed, R_min, c))
    
    if os.path.exists("{0}.odb".format(job)):
        print('\n\nABAQUS/CAE odb file created for iteration {0}.\n'.format(convg.iter))
    else:
        print('\n\nERROR: ABAQUS/CAE odb file not created for iteration {0}.\n'.format(convg.iter))
        sys.exit()
    
    #%% Update files containing edge nodes
    
    os.system('cmd /k "abaqus cae noGui=NACA0012_post_processing.py -- {0} {1} {2} {3}"\n'.format(job, convg.iter, Elem, c))
    
    #%% Organize data
        
    for i in range(0,Elem+1):
        
        coord_U = 'EDGE{0}_U'.format(i)
        coord_L = 'EDGE{0}_L'.format(i)
        index_U, index_L = organize_txt_TESTE.get_index (i)    
        _, _, x_U, y_U, z_U, x_L, y_L, z_L = organize_txt_TESTE.organize_data (coord_U, coord_L)
        x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted = organize_txt_TESTE.sort_vectors (index_U, index_L, x_U, y_U, z_U, x_L, y_L, z_L)
        coords_U_sorted = np.vstack((x_U_sorted, y_U_sorted, z_U_sorted)).T
        coords_L_sorted = np.vstack((x_L_sorted, y_L_sorted, z_L_sorted)).T
        coords_U_local = np.zeros((len(coords_U_sorted),3))
        coords_L_local = np.zeros((len(coords_L_sorted),3))
        for j in range(0,len(coords_U_sorted)):
            coords_U_local[j,:] = np.add(coords_U_sorted[j,:], -coords_U_sorted[-1,:])       
        for j in range(0,len(coords_L_sorted)):
            coords_L_local[j,:] = np.add(coords_L_sorted[j,:], -coords_U_sorted[-1,:])
        _, bi, _, _ = transformation.find_angles (coords_U_local)
        if bi<b.initial:
            Ry = transformation.T_matrix(-(b.initial-bi))
            theta = pitch.initial+(b.initial-bi)*180/np.pi           
        else:
            Ry = transformation.T_matrix(bi-b.initial)
            theta = pitch.initial-(bi-b.initial)*180/np.pi

        b.one = np.append(b.one, bi)
        pitch.one = np.append(pitch.one, theta)
        coords_U_local = transformation.local_coords (Ry, coords_U_local)
        coords_L_local = transformation.local_coords (Ry, coords_L_local)
        X_U_optim = coords_U_local[:,0]
        Y_U_optim = coords_U_local[:,2]
        X_L_optim = coords_L_local[:,0]
        Y_L_optim = coords_L_local[:,2]   
        
        # Interpolate airfoil
        if OPEN == 1:
            X_L_final, Y_L_final, _, _, _  = organize_txt_TESTE.close_airfoil (X_U_optim, Y_U_optim, X_L_optim, Y_L_optim)
        else:
            X_L_final = X_L_optim
            Y_L_final = Y_L_optim
            
        nodes_U = len(X_U_optim)
        nodes_L = len(X_L_optim)
        
        organize_txt_TESTE.XFOIL_file (X_U_optim, Y_U_optim, X_L_final, Y_L_final, 'COORDS{0}'.format(i))
        
        # Input data for XFoil
        Re, Ma, Ui, phi = input_data.data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air)
        U = np.append(U, Ui)
        AoA = theta-phi
        alpha.one = np.append(alpha.one, AoA)
        ind_angle.one = np.append(ind_angle.one, phi)
        
        # Write input file for XFoil
        airfoil = 'COORDS{0}'.format(i)
        polar_data, nodal_data = fwrite.write (round(AoA,5), Re, Ma, airfoil, iter_xf)
        
        # Run XFoil  
        Cp_U, Cp_L, Cl, Cd = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
        
        # Correct non-converging profiles
        if Cl == 'None':
            Cl0 = Cl
            inc = 0
            AoA0 = AoA
            while Cl0 == 'None':
                # Lower limit
                inc = inc-0.01
                AoA = AoA0+inc
                polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                Cp_U_lower, Cp_L_lower, Cl0, Cd0 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)            
            print('\nCOORDS{0}: Alfa = {1} --> Cl = {2}'.format(i, AoA, Cl0))
            # Upper limit
            Cl1 = 'None'
            inc = 0
            # AoA = AoA+inc*-1
            AoA = AoA0
            while Cl1 == 'None':
                inc = inc+0.01
                AoA = AoA0+inc
                polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                Cp_U_upper, Cp_L_upper, Cl1, Cd1 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
            print('             Alfa = {0} --> Cl = {1}'.format(AoA, Cl1))
            
            # Obtain average value for Cp and Cl
            for j in range(0,len(Cp_U)):
                Cp_U[j] = (Cp_U_upper[j]+Cp_U_lower[j])/2
            for j in range(0,len(Cp_L)):
                Cp_L[j] = (Cp_L_upper[j]+Cp_L_lower[j])/2
            Cl = (Cl0+Cl1)/2
            Cd = (Cd0+Cd1)/2
            print('\nCl average = {0}'.format(Cl))
        
    # Evaluate residual
    for i in range(0, Elem+1):
        pitch_diff = np.append(pitch_diff, pitch.one[i]-pitch.zero[i])
    convg.res = max(pitch_diff, key=abs)
    theta_T = np.append(theta_T, convg.res)
    
#=============================================== HOVERING FLIGHT ==============================================================

#%% Initialize variables of interest

miu_air = 0 # hovering flight
magnitude = 0.5*rho
P_atm = 0.101325 #[MPa]
init_guess = 9.9951
pitch.initial = init_guess
Lift1 = 10E6

while abs(1-Lift1/(W/Nb)) > 0.05:
    
    theta_T = []
    theta_T_norm = []
    
    while convg.res > 0.035:
        
        alpha.zero = []
        alpha.one= []
        pitch.zero = []
        pitch.one = []
        ind_angle.zero = []
        ind_angle.one = []
        
        pitch_diff = []
        convg.res = []
        
        variable_U = []
        variable_L = []
        U = []
        
        CL0 = []
        CD0 = []
        Lift0 = 0
        Drag0= 0
        CL1 = []
        CD1 = []
        Lift1 = 0
        Drag1 = 0
        
        #%% Initialize files for mapped analytical fields
        
        mapped_U=open('MAPU.txt', 'w')
        mapped_L=open('MAPL.txt', 'w')
        
        #%% Organize data
        
        for i in range(0,Elem+1):
            
            coord_U = 'EDGE{0}_U'.format(i)
            coord_L = 'EDGE{0}_L'.format(i)
            index_U, index_L = organize_txt_TESTE.get_index (i)    
            _, _, x_U, y_U, z_U, x_L, y_L, z_L = organize_txt_TESTE.organize_data (coord_U, coord_L)
            x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted = organize_txt_TESTE.sort_vectors (index_U, index_L, x_U, y_U, z_U, x_L, y_L, z_L)
                
                # Define b0
            if i == 0:
                b0 = atan(z_U_sorted[0]/x_U_sorted[0])
                b.initial = b0            
            
            coords_U_sorted = np.vstack((x_U_sorted, y_U_sorted, z_U_sorted)).T
            coords_L_sorted = np.vstack((x_L_sorted, y_L_sorted, z_L_sorted)).T
            coords_U_local = np.zeros((len(coords_U_sorted),3))
            coords_L_local = np.zeros((len(coords_L_sorted),3))
            for j in range(0,len(coords_U_sorted)):
                coords_U_local[j,:] = np.add(coords_U_sorted[j,:], -coords_U_sorted[-1,:])       
            for j in range(0,len(coords_L_sorted)):
                coords_L_local[j,:] = np.add(coords_L_sorted[j,:], -coords_U_sorted[-1,:])
            _, bi, _, _ = transformation.find_angles (coords_U_local)
            if bi<b.initial:
                Ry = transformation.T_matrix(-(b.initial-bi))
                theta = pitch.initial+(b.initial-bi)*180/np.pi           
            else:
                Ry = transformation.T_matrix(bi-b.initial)
                theta = pitch.initial-(bi-b.initial)*180/np.pi

            b.zero = np.append(b.zero, bi)
            pitch.zero = np.append(pitch.zero, theta)
            coords_U_local = transformation.local_coords (Ry, coords_U_local)
            coords_L_local = transformation.local_coords (Ry, coords_L_local)
            X_U_optim = coords_U_local[:,0]
            Y_U_optim = coords_U_local[:,2]
            X_L_optim = coords_L_local[:,0]
            Y_L_optim = coords_L_local[:,2]   
            
            # Interpolate airfoil
            if OPEN == 1:
                X_L_final, Y_L_final, _, _, _  = organize_txt_TESTE.close_airfoil (X_U_optim, Y_U_optim, X_L_optim, Y_L_optim)
            else:
                X_L_final = X_L_optim
                Y_L_final = Y_L_optim
                
            nodes_U = len(X_U_optim)
            nodes_L = len(X_L_optim)
            
            organize_txt_TESTE.XFOIL_file (X_U_optim, Y_U_optim, X_L_final, Y_L_final, 'COORDS{0}'.format(i))
            
            # Input data for XFoil
            Re, Ma, Ui, phi = input_data.data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air)
            U = np.append(U, Ui)
            AoA = theta-phi
            alpha.zero = np.append(alpha.zero, AoA)
            ind_angle.zero = np.append(ind_angle.zero, phi)
            
            # Write input file for XFoil
            airfoil = 'COORDS{0}'.format(i)
            polar_data, nodal_data = fwrite.write (round(AoA,5), Re, Ma, airfoil, iter_xf)
            
            # Run XFoil  
            Cp_U, Cp_L, Cl, Cd = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
            
            # Correct non-converging profiles
            if Cl == 'None':
                Cl0 = Cl
                inc = 0
                AoA0 = AoA
                while Cl0 == 'None':
                    # Lower limit
                    inc = inc-0.01
                    AoA = AoA0+inc
                    polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                    Cp_U_lower, Cp_L_lower, Cl0, Cd0 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)            
                print('\nCOORDS{0}: Alfa = {1} --> Cl = {2}'.format(i, AoA, Cl0))
                # Upper limit
                Cl1 = 'None'
                inc = 0
                # AoA = AoA+inc*-1
                AoA = AoA0
                while Cl1 == 'None':
                    inc = inc+0.01
                    AoA = AoA0+inc
                    polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                    Cp_U_upper, Cp_L_upper, Cl1, Cd1 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
                print('             Alfa = {0} --> Cl = {1}'.format(AoA, Cl1))
                
                # Obtain average value for Cp and Cl
                for j in range(0,len(Cp_U)):
                    Cp_U[j] = (Cp_U_upper[j]+Cp_U_lower[j])/2
                for j in range(0,len(Cp_L)):
                    Cp_L[j] = (Cp_L_upper[j]+Cp_L_lower[j])/2
                Cl = (Cl0+Cl1)/2
                Cd = (Cd0+Cd1)/2
                print('\nCl average = {0}'.format(Cl))
            
            # Store pressure coefficiente and fluid velocity
            variable_U[:] = Cp_U*(Ui**2)*magnitude+P_atm
            variable_L[:] = Cp_L*(Ui**2)*magnitude+P_atm
            
            # Store Lift coefficients
            CL0 = np.append(CL0, Cl)
            CD0 = np.append(CD0, Cd)
            
            #%% Create mapped field for distributed pressure
            
            L0 = L*i
            
            # Upper Surface
            for j in range(0, nodes_U):       
                Cp_U_0 = variable_U[j]     
                X_U_0 = coords_U_sorted[j,0]*c
                Y_U_0 = coords_U_sorted[j,2]*c
                mapped_U.write('{0}, {1}, {2}, {3}\n'.format(X_U_0, L0, Y_U_0, Cp_U_0))
                
            # Lower Surface
            for j in range(0, nodes_L):       
                Cp_L_0 = variable_L[j]      
                X_L_0 = coords_L_sorted[j,0]*c
                Y_L_0 = coords_L_sorted[j,2]*c    
                mapped_L.write('{0}, {1}, {2}, {3}\n'.format(X_L_0, L0, Y_L_0, Cp_L_0))
        
        #%% Close files
        
        mapped_U.close()
        mapped_L.close()
            
        #%% Restart analysis
        
        convg.iter = convg.iter+1
        job = job0 + '-{0}'.format(convg.iter)
        for fname in os.listdir(my_dir):
            if fname.startswith(job):
                os.remove(os.path.join(my_dir, fname))
                print('\nFile {0} deleted.'.format(fname))
        os.system('cmd /k "abaqus cae noGui=NACA0012_restart_new.py -- {0} {1} {2} {3} {4} {5} {6}"\n'.format(convg.iter, job0, job, P_atm, rotor_speed, R_min, c))
        
        if os.path.exists("{0}.odb".format(job)):
            print('\n\nABAQUS/CAE odb file created for iteration {0}.\n'.format(convg.iter))
        else:
            print('\n\nERROR: ABAQUS/CAE odb file not created for iteration {0}.\n'.format(convg.iter))
            sys.exit()
        
        #%% Update files containing edge nodes
        
        os.system('cmd /k "abaqus cae noGui=NACA0012_post_processing.py -- {0} {1} {2} {3}"\n'.format(job, convg.iter, Elem, c))
        
        #%% Organize data
        
        for i in range(0,Elem+1):
            
            coord_U = 'EDGE{0}_U'.format(i)
            coord_L = 'EDGE{0}_L'.format(i)
            index_U, index_L = organize_txt_TESTE.get_index (i)    
            _, _, x_U, y_U, z_U, x_L, y_L, z_L = organize_txt_TESTE.organize_data (coord_U, coord_L)
            x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted = organize_txt_TESTE.sort_vectors (index_U, index_L, x_U, y_U, z_U, x_L, y_L, z_L)
            coords_U_sorted = np.vstack((x_U_sorted, y_U_sorted, z_U_sorted)).T
            coords_L_sorted = np.vstack((x_L_sorted, y_L_sorted, z_L_sorted)).T
            coords_U_local = np.zeros((len(coords_U_sorted),3))
            coords_L_local = np.zeros((len(coords_L_sorted),3))
            for j in range(0,len(coords_U_sorted)):
                coords_U_local[j,:] = np.add(coords_U_sorted[j,:], -coords_U_sorted[-1,:])       
            for j in range(0,len(coords_L_sorted)):
                coords_L_local[j,:] = np.add(coords_L_sorted[j,:], -coords_U_sorted[-1,:])
            _, bi, _, _ = transformation.find_angles (coords_U_local)
            if bi<b.initial:
                Ry = transformation.T_matrix(-(b.initial-bi))
                theta = pitch.initial+(b.initial-bi)*180/np.pi           
            else:
                Ry = transformation.T_matrix(bi-b.initial)
                theta = pitch.initial-(bi-b.initial)*180/np.pi
            
            b.one = np.append(b.one, bi)
            pitch.one = np.append(pitch.one, theta)
            coords_U_local = transformation.local_coords (Ry, coords_U_local)
            coords_L_local = transformation.local_coords (Ry, coords_L_local)
            X_U_optim = coords_U_local[:,0]
            Y_U_optim = coords_U_local[:,2]
            X_L_optim = coords_L_local[:,0]
            Y_L_optim = coords_L_local[:,2]   
            
            # Interpolate airfoil
            if OPEN == 1:
                X_L_final, Y_L_final, _, _, _  = organize_txt_TESTE.close_airfoil (X_U_optim, Y_U_optim, X_L_optim, Y_L_optim)
            else:
                X_L_final = X_L_optim
                Y_L_final = Y_L_optim
                
            nodes_U = len(X_U_optim)
            nodes_L = len(X_L_optim)
            
            organize_txt_TESTE.XFOIL_file (X_U_optim, Y_U_optim, X_L_final, Y_L_final, 'COORDS{0}'.format(i))
            
            # Input data for XFoil
            Re, Ma, Ui, phi = input_data.data (R_min, L, i, W, rotor_speed, rho, A, c, kin_viscosity, sound, V_tip, azimut, miu_air)
            U = np.append(U, Ui)
            AoA = theta-phi
            alpha.one = np.append(alpha.one, AoA)
            ind_angle.one = np.append(ind_angle.one, phi)
            
            # Write input file for XFoil
            airfoil = 'COORDS{0}'.format(i)
            polar_data, nodal_data = fwrite.write (round(AoA,5), Re, Ma, airfoil, iter_xf)
            
            # Run XFoil  
            Cp_U, Cp_L, Cl, Cd = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
            
            # Correct non-converging profiles
            if Cl == 'None':
                Cl0 = Cl
                inc = 0
                AoA0 = AoA
                while Cl0 == 'None':
                    # Lower limit
                    inc = inc-0.01
                    AoA = AoA0+inc
                    polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                    Cp_U_lower, Cp_L_lower, Cl0, Cd0 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)            
                print('\nCOORDS{0}: Alfa = {1} --> Cl = {2}'.format(i, AoA, Cl0))
                # Upper limit
                Cl1 = 'None'
                inc = 0
                # AoA = AoA+inc*-1
                AoA = AoA0
                while Cl1 == 'None':
                    inc = inc+0.01
                    AoA = AoA0+inc
                    polar_data, nodal_data = fwrite.write (AoA, Re, Ma, airfoil, iter_xf)
                    Cp_U_upper, Cp_L_upper, Cl1, Cd1 = xfrun.Xfoil (polar_data, nodal_data, nodes_U, nodes_L)
                print('             Alfa = {0} --> Cl = {1}'.format(AoA, Cl1))
                
                # Obtain average value for Cp and Cl
                for j in range(0,len(Cp_U)):
                    Cp_U[j] = (Cp_U_upper[j]+Cp_U_lower[j])/2
                for j in range(0,len(Cp_L)):
                    Cp_L[j] = (Cp_L_upper[j]+Cp_L_lower[j])/2
                Cl = (Cl0+Cl1)/2
                Cd = (Cd0+Cd1)/2
                print('\nCl average = {0}'.format(Cl))
            
            # Store Lift coefficients
            CL1 = np.append(CL1, Cl)
            CD1 = np.append(CD1, Cd)
        
        for i in range(1, Elem+1):
            dL = 0.5*rho*(U[i]**2)*c*CL1[i]*L
            Lift1 = Lift1+dL
            dD = 0.5*rho*(U[i]**2)*c*CD1[i]*L
            Drag1 = Drag1+dD
        
        # Evaluate residual
        for i in range(0, Elem+1):
            pitch_diff = np.append(pitch_diff, pitch.one[i]-pitch.zero[i])
        convg.res = max(pitch_diff, key=abs)
        theta_T = np.append(theta_T, convg.res)
    
    if (1-Lift1/(W/Nb)) > 0:
        pitch.initial = pitch.initial + 0.1
    else:
        pitch.initial = pitch.initial - 0.1
        

#%% Finish Program

end = time.time()
ty_res = time.gmtime(end-start)
res = time.strftime("%H:%M:%S",ty_res)
print('\nProgram execution time: {0}\n'.format(res))    