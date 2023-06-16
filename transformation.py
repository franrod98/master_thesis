import numpy as np
# from math import sin, cos
# from sympy import sin, cos
# from gekko import GEKKO
# import sympy as sym
from scipy.optimize import fsolve

def find_angles (coords_U_sorted):
    
    # Points of interest
    
    origin = coords_U_sorted[-1,:]
    x_axis = np.add(coords_U_sorted[0,:],-origin)
    xy_plane = np.add(coords_U_sorted[2,:],-origin)
    
    p1 = np.array([0, 0, 0])
    p2 = np.array(x_axis)
    p3 = np.array(xy_plane)
    
    # Define in-plane vectors
    v1 = p3 - p1
    v2 = p2 - p1
    
    # Cross product to define normal vector
    c_prod = np.cross(v1, v2)
    y_axis = c_prod

    #%% Find rotations along all axis
    
    # X axis
    vector1 = (1,0)
    vector2 = (y_axis[1], y_axis[2])
    
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    
    a0 = np.arccos(dot_product)
    
    # Y axis
    vector1 = (1,0)
    vector2 = (x_axis[0], x_axis[2])
    
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    
    b0 = np.arccos(dot_product)
    
    if x_axis[2]<0:
        b0 = -b0
    
    # Z axis  
    vector1 = (1,0)
    vector2 = (x_axis[0], x_axis[1])
    
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    
    g0 = np.arccos(dot_product)
    
    return (a0, b0, g0, origin)

# def find_angles (coords_U_sorted, X_AXIS):
       
#     origin = coords_U_sorted[-1,:]
#     x_axis = coords_U_sorted[0,:]

#     # Translation
#     x0 = -origin[0]
#     y0 = -origin[1]
#     z0 = -origin[2]

# #%% GEEKO
    
#     m = GEKKO()
#     cos = m.cos
#     sin = m.sin
#     a, b, g = [m.Var() for i in range(3)]
    
    # # Matrix (2)
    # m.Equations([cos(b)*cos(g)*x_axis[0]+(sin(a)*sin(g)+cos(a)*sin(b)*cos(g))*x_axis[1]+(sin(a)*sin(g)+cos(a)*sin(b)*cos(g))*x_axis[2]+x0==X_AXIS, \
    #               cos(b)*sin(g)*x_axis[0]+(cos(a)*cos(g)+sin(a)*sin(b)*sin(g))*x_axis[1]+(-sin(a)*cos(g)+cos(a)*sin(b)*sin(g))*x_axis[2]+y0==0, \
    #               -sin(b)*x_axis[0]+sin(a)*cos(b)*x_axis[1]+cos(a)*cos(b)*x_axis[2]+z0==0])
    
    # # Matrix (1)
    # m.Equations([cos(b)*cos(g)*x_axis[0]+(sin(a)*sin(b)*cos(g)+cos(a)*sin(g))*x_axis[1]+(sin(a)*sin(g)-cos(a)*sin(b)*cos(g))*x_axis[2]+x0==1, \
    #               -cos(b)*sin(g)*x_axis[0]+(cos(a)*cos(g)-sin(a)*sin(b)*sin(g))*x_axis[1]+(sin(a)*cos(g)+cos(a)*sin(b)*sin(g))*x_axis[2]+y0==0, \
    #               sin(b)*x_axis[0]-sin(a)*cos(b)*x_axis[1]+cos(a)*cos(b)*x_axis[2]+z0==0])
    
    # m.solve(disp=False)
    
    # return(a.value[0], b.value[0], g.value[0], x0, y0, z0)

#%% SYMPY
    
    # a,b,g = sym.symbols('a,b,g') 
    # T1 = sym.Eq(cos(b)*cos(g)*x_axis[0]+(sin(a)*sin(b)*cos(g)+cos(a)*sin(g))*x_axis[1]+(sin(a)*sin(g)-cos(a)*sin(b)*cos(g))*x_axis[2]+x0, 1) 
    # T2 = sym.Eq(-cos(b)*sin(g)*x_axis[0]+(cos(a)*cos(g)-sin(a)*sin(b)*sin(g))*x_axis[1]+(sin(a)*cos(g)+cos(a)*sin(b)*sin(g))*x_axis[2]+y0, 0)
    # T3 = sym.Eq(sin(b)*x_axis[0]-sin(a)*cos(b)*x_axis[1]+cos(a)*cos(b)*x_axis[2]+z0, 0)
    # sym.solve([T1, T2, T3], (a, b, g))
    
    # return(a, b, g, x0, y0, z0)
    
#%% SCIPY

    # guess = np.array([a0, b0, g0])
    # T = fsolve(non_lin, guess, args=(x_axis, x0, y0, z0, X_AXIS))
    # a = T[0]
    # b = T[1]
    # g = T[2]
    
    # return(a, b, g, x0, y0, z0)
    
def T_matrix (b):

    from math import sin, cos
    # # Transformation matrix (1)
    # T = np.array([[ cos(b)*cos(g), (sin(a)*sin(b)*cos(g) + cos(a)*sin(g)), (sin(a)*sin(g) - cos(a)*sin(b)*cos(g)), x0],
    #               [-cos(b)*sin(g), (cos(a)*cos(g)-sin(a)*sin(b)*sin(g)), (sin(a)*cos(g) + cos(a)*sin(b)*sin(g)), y0],
    #               [        sin(b), -sin(a)*cos(b), cos(a)*cos(b), z0],
    #               [ 0, 0, 0, 1]])
    
    # # Transformation matrix (2)
    # T = np.array([[ cos(b)*cos(g), (sin(a)*sin(b)*cos(g) - cos(a)*sin(g)), (sin(a)*sin(g) + cos(a)*sin(b)*cos(g)), x0],
    #               [cos(b)*sin(g), (cos(a)*cos(g)+sin(a)*sin(b)*sin(g)), (-sin(a)*cos(g) + cos(a)*sin(b)*sin(g)), y0],
    #               [        -sin(b), sin(a)*cos(b), cos(a)*cos(b), z0],
    #               [ 0, 0, 0, 1]])
    Ry = np.array([[ cos(b), 0, sin(b)],
                   [ 0      , 1,      0],
                   [-sin(b), 0, cos(b)]])
 
    return (Ry)

def local_coords (T, coordinates):
    # Transform coordinates
    new_coords = np.empty((0,3))
    for j in coordinates:
        coord = np.array([j[0], j[1], j[2]])
        transform = np.transpose(T.dot(coord))
        transformed = np.array([[transform[0], transform[1], transform[2]]])
        new_coords = np.append(new_coords, transformed, axis=0)
        
    return (new_coords)