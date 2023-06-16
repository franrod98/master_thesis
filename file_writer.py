""" Write main Xfoil input file """

import os

def write(AoA, Re, Ma, airfoil, iter_xf):
    # File specifications
    xfoil_file = 'input_file.txt'
    # airfoil = 'COORD.txt'
    polar_data = 'polar_data.txt'
    nodal_data = 'nodal_data.txt'
    PSAV = 'saved_airfoil.txt'
    
    # Delete all previous files
    if os.path.exists("{0}".format(polar_data)):
        os.remove("{0}".format(polar_data))
    if os.path.exists("{0}".format(PSAV)):
        os.remove("{0}".format(PSAV))
    if os.path.exists("{0}".format(nodal_data)):
        os.remove("{0}".format(nodal_data))
    if os.path.exists("{0}".format(xfoil_file)):
        os.remove("{0}".format(xfoil_file))
  
    # Write input file
    input_file = open("{0}".format(xfoil_file), 'w')
    input_file.write("LOAD {0}.txt\n".format(airfoil))
    input_file.write("\n")
    input_file.write("PSAV {0}\n".format(PSAV))
    # input_file.write("PPAR\nN {0}\n\n\n".format(ppar))
    input_file.write("OPER\n")
    input_file.write("iter\n{0}\n".format(iter_xf))
    input_file.write("visc\n")
    input_file.write("{0}\n".format(Re))
    input_file.write("Mach {0}\n".format(Ma))
    input_file.write("PACC\n")
    input_file.write("{0}\n\n".format(polar_data))
    input_file.write("Alfa {0}\n".format(AoA))
    input_file.write("CPWR {0}\n".format(nodal_data))
    input_file.write("\n")
    input_file.write("quit\n")
    input_file.close()
    
    return(polar_data, nodal_data)