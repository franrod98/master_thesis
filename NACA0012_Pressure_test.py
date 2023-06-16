# -*- coding: utf-8 -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
#======================================================== INPUT DATA =========================================================================
rotor_speed=float(sys.argv[-7])
P_atm=float(sys.argv[-6])
R_min=float(sys.argv[-5])
iteration=int(float(sys.argv[-4]))
cae_file=sys.argv[-3]
c=float(sys.argv[-2])
#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/{0}.cae'.format(cae_file))
#===================================================== ANALYTICAL FIELD ======================================================================
mdb.models['Model-1'].rootAssembly.regenerate()
# Upper Surface
path="C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/MAPU.txt"
datalist=[]
with open(path,"rb") as fp:
    for row in fp.readlines():
        tmp=row.split(",")
        try:
            datalist.append((float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])))
        except:pass
mdb.models['Model-1'].MappedField(name='UPPER SURFACES', description='', 
        regionType=POINT, partLevelData=False, localCsys=None, 
        pointDataFormat=XYZ, fieldDataType=SCALAR, xyzPointData=datalist)
mdb.models['Model-1'].rootAssembly.regenerate()
# Lower Surface
path="C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/MAPL.txt"
datalist=[]
with open(path,"rb") as fp:
    for row in fp.readlines():
        tmp=row.split(",")
        try:
            datalist.append((float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])))
        except:pass
mdb.models['Model-1'].MappedField(name='LOWER SURFACES', description='', 
        regionType=POINT, partLevelData=False, localCsys=None, 
        pointDataFormat=XYZ, fieldDataType=SCALAR, xyzPointData=datalist)
#===================================================== DISTRIBUTED PRESSURE ==================================================================
mdb.models['Model-1'].rootAssembly.regenerate()
# Upper Surface
assembly='NACA0012_U-1'
mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
    distributionType=FIELD, field='UPPER SURFACES', magnitude=1, name='UPPER SURFACES', 
    region=mdb.models['Model-1'].rootAssembly.instances[assembly].surfaces['UPPER SURFACES'])
mdb.models['Model-1'].rootAssembly.regenerate()
# Lower Surface
assembly='NACA0012_L-1'
mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
    distributionType=FIELD, field='LOWER SURFACES', magnitude=1, name='LOWER SURFACES', 
    region=mdb.models['Model-1'].rootAssembly.instances[assembly].surfaces['LOWER SURFACES'])
mdb.models['Model-1'].rootAssembly.regenerate()
#===================================================== ATMOSPHERIC PRESSURE ==================================================================
# # Upper Surface
# assembly='NACA0012_U-1'
# mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
#     distributionType=UNIFORM, field='', magnitude=P_atm, name='P_atm -- U', 
#     region=mdb.models['Model-1'].rootAssembly.instances[assembly].surfaces['UPPER SURFACES'])
# mdb.models['Model-1'].rootAssembly.regenerate()
# # Lower Surface
# assembly='NACA0012_L-1'
# mdb.models['Model-1'].Pressure(amplitude=UNSET, createStepName='Step-1', 
#     distributionType=UNIFORM, field='', magnitude=P_atm, name='P_atm -- L', 
#     region=mdb.models['Model-1'].rootAssembly.instances[assembly].surfaces['LOWER SURFACES'])
# mdb.models['Model-1'].rootAssembly.regenerate()
#==================================================== CENTRIFUGAL FORCE =======================================================================
# Upper Surface
assembly='NACA0012_U-1'
mdb.models['Model-1'].RotationalBodyForce(centrifugal=ON, createStepName=
    'Step-1', magnitude=rotor_speed, name='Centrifugal_U', point1=(c/2, -R_min, 0.0), 
    point2=(c/2, -R_min, 10.0), region=Region(
    cells=mdb.models['Model-1'].rootAssembly.instances[assembly].cells), rotaryAcceleration=OFF) ### CHANGE!!!
# Lower Surface
assembly='NACA0012_L-1'
mdb.models['Model-1'].RotationalBodyForce(centrifugal=ON, createStepName=
    'Step-1', magnitude=rotor_speed, name='Centrifugal_L', point1=(c/2, -R_min, 0.0), 
    point2=(c/2, -R_min, 10.0), region=Region(
    cells=mdb.models['Model-1'].rootAssembly.instances[assembly].cells), rotaryAcceleration=OFF)
#====================================================== SAVE .cae FILE =======================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/{0}.cae'.format(cae_file))