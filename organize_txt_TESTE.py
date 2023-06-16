"""
%--------------------------------------------------------------------------%
%           Aeroelastic Tailoring of Helicopter Main Rotor Blade           %
%                                                                          %
%           @author: Francisco Campos Moreira Rodrigues                    %
%           @date: October, 2022                                           %
%                                                                          %
%                       SCRIPT: organize_txt.py                            % 
%--------------------------------------------------------------------------%
"""

import numpy as np

def organize_data (coord_U, coord_L):
    
    #%% Read .txt files containing coordinates
    
    # (U)pper surface
    x_U = []
    y_U = []
    z_U = []
    COORD_U = np.loadtxt("{0}.txt".format(coord_U), skiprows=1)
    x_U = COORD_U[:,0]
    y_U = COORD_U[:,1]
    z_U = COORD_U[:,2]
    
    # (L)ower surface
    x_L = []
    y_L = []
    z_L = []
    COORD_L = np.loadtxt("{0}.txt".format(coord_L), skiprows=1)
    x_L = COORD_L[:,0]
    y_L = COORD_L[:,1]
    z_L = COORD_L[:,2]
    
    #%% Store indexes of the original vectors
    sort_INDEX_U = np.argsort(x_U)[::-1]    
    sort_INDEX_L = np.argsort(x_L)

    return (sort_INDEX_U, sort_INDEX_L, x_U, y_U, z_U, x_L, y_L, z_L)

def sort_vectors (sort_INDEX_U, sort_INDEX_L, x_U, y_U, z_U, x_L, y_L, z_L):
    
    # Sort vectors with the stored indexes of X_U
    numNodes_U = len(x_U)
    x_U_sorted = np.zeros(numNodes_U)
    y_U_sorted = np.zeros(numNodes_U)
    z_U_sorted = np.zeros(numNodes_U)
    
    # Sort vectors with the stored indexes of X_L
    numNodes_L = len(x_L)
    x_L_sorted = np.zeros(numNodes_L)
    y_L_sorted = np.zeros(numNodes_L)    
    z_L_sorted = np.zeros(numNodes_L)
    
    # Adjust vectors
    for i in range(0, numNodes_U):
        x_U_sorted[i] = x_U[int(sort_INDEX_U[i])]
        z_U_sorted[i] = z_U[int(sort_INDEX_U[i])]
        if str(y_U) != 'none':
            y_U_sorted[i] = y_U[int(sort_INDEX_U[i])]
    # Adjust vectors
    for i in range(0, numNodes_L):
        x_L_sorted[i] = x_L[int(sort_INDEX_L[i])]
        z_L_sorted[i] = z_L[int(sort_INDEX_L[i])]
        if str(y_L) != 'none':
            y_L_sorted[i] = y_L[int(sort_INDEX_L[i])]
    
    return (x_U_sorted, y_U_sorted, z_U_sorted, x_L_sorted, y_L_sorted, z_L_sorted)

def get_index (i):
    
    file = open('INDEX_U.txt', 'r')
    lines=file.readlines()[1:]
    for line in lines:
        data = line.split(' ')
        values = [int(float(x)) for x in data[0:-1]]
        if values[0] == i:
            index_U = values[1:]
    file = open('INDEX_L.txt', 'r')
    lines=file.readlines()[1:]
    for line in lines:
        data = line.split(' ')
        values = [int(float(x)) for x in data[0:-1]]
        if values[0] == i:
            index_L = values[1:]     
    
    return(index_U, index_L)

def optimize_computation (X_U_final, Y_U_final, X_L_final, Y_L_final):
    
    #%% Optimize computational time
    
    # Upper Surface 
    extend = np.loadtxt("NACA0012_U.txt", usecols=(0,1))
    for i in range(0, len(extend)):
        if (X_U_final[-1]<extend[i,0]) and (X_U_final[-2]>extend[i,0]):
            if (Y_U_final[-2]>extend[i,1]): # ADDED!!!
                X_U_final = np.insert(X_U_final, -1, extend[i,0])
                Y_U_final = np.insert(Y_U_final, -1, extend[i,1])
    
    # Lower Surface
    plus = 0
    position = 1
    extend = np.loadtxt("NACA0012_L.txt", usecols=(0,1))
    for i in range(0, len(extend)):
        if (X_L_final[0+plus]<extend[i,0]) and (X_L_final[1+plus]>extend[i,0]):
            if (Y_L_final[1+plus]<extend[i,1]): # ADDED!!!
                plus = plus+1
                X_L_final = np.insert(X_L_final, position, extend[i,0])
                Y_L_final = np.insert(Y_L_final, position, extend[i,1])
                position = position+1
    
    return (X_U_final, Y_U_final, X_L_final, Y_L_final)

def close_airfoil (X_U_sorted, Y_U_sorted, X_L_sorted, Y_L_sorted):
        
    #%% Interpolation
    
    import interpolation as inter
    
    nodes_U = len(X_U_sorted)
    nodes_L = len(X_L_sorted)
    nodes_inter = int(nodes_L/2)
    X_inter = np.zeros(nodes_inter)
    Y_inter = np.zeros(nodes_inter)
    
    for i in range(1, nodes_inter+1):
        X_inter[-i] = X_L_sorted[-i]
        Y_inter[-i] = Y_L_sorted[-i]
        
    X_trail, Y_trail = inter.polynomial (X_inter, Y_inter, X_U_sorted)
    X_L_final = np.append(X_L_sorted, X_trail)
    Y_L_final = np.append(Y_L_sorted, Y_trail)
    
    trail = abs(len(X_L_sorted)-len(X_L_final))
    
    return (X_L_final, Y_L_final, nodes_U, nodes_L, trail)

def XFOIL_file (X_U_final, Y_U_final, X_L_final, Y_L_final, file_name):

    #%% Write final COORD.txt
    
    COORD = open('{0}.txt'.format(file_name), 'w')
    for i in range (len(X_U_final)):
        if (Y_U_final[i] < 0):
            COORD.write('     {:.6f}  {:.6f}\n'.format(X_U_final[i],Y_U_final[i]))
        else:
            COORD.write('     {:.6f}   {:.6f}\n'.format(X_U_final[i],Y_U_final[i]))    
    for i in range (len(X_L_final)):
        if (Y_L_final[i] < 0):
            COORD.write('     {:.6f}  {:.6f}\n'.format(X_L_final[i],Y_L_final[i]))
        else:
            COORD.write('     {:.6f}   {:.6f}\n'.format(X_L_final[i],Y_L_final[i]))
    COORD.close()

    return ()