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
iteration=int(float(sys.argv[-7]))
job0=sys.argv[-6]
job=sys.argv[-5]
P_atm=float(sys.argv[-4])
rotor_speed=float(sys.argv[-3])
R_min=float(sys.argv[-2])
c=float(sys.argv[-1])
#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')
#======================================================== COPY MODEL =========================================================================
job0=job0+'-{0}'.format(iteration-1)
model0='Model-{0}'.format(iteration-1)
model='Model-{0}'.format(iteration)
step0='Step-{0}'.format(iteration-1)
new_step='Step-{0}'.format(iteration)
mdb.Model(name=model, objectToCopy=mdb.models[model0])
mdb.models[model].setValues(restartJob=job0, restartStep=step0)
#========================================================= NEW STEP ==========================================================================
mdb.models[model].StaticStep(initialInc=0.01, minInc=1e-06, name=new_step, nlgeom=ON, previous=step0)
mdb.models[model].fieldOutputRequests['F-Output-1'].setValues(frequency=2, 
    variables=('S', 'E', 'EE', 'LE', 'U', 'UT', 'UR', 'RF', 'CF', 'P', 'COORD'), position=NODES)
#====================================================== RESTART REQUEST ======================================================================
mdb.models[model].steps[new_step].Restart(frequency=2, numberIntervals=0, 
    overlay=ON, timeMarks=OFF)
#===================================================== ANALYTICAL FIELD ======================================================================
del mdb.models[model].loads['UPPER SURFACES']
del mdb.models[model].loads['LOWER SURFACES']
del mdb.models[model].loads['Centrifugal_U']
del mdb.models[model].loads['Centrifugal_L']
# assembly='NACA 0012-1'
mdb.models[model].rootAssembly.regenerate()
# Upper Surface
path="C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/MAPU.txt"
datalist=[]
with open(path,"rb") as fp:
    for row in fp.readlines():
        tmp=row.split(",")
        try:
            datalist.append((float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])))
        except:pass
mdb.models[model].MappedField(name='UPPER SURFACES', description='', 
        regionType=POINT, partLevelData=False, localCsys=None, 
        pointDataFormat=XYZ, fieldDataType=SCALAR, xyzPointData=datalist)
mdb.models[model].rootAssembly.regenerate()
# Lower Surface
path="C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/MAPL.txt"
datalist=[]
with open(path,"rb") as fp:
    for row in fp.readlines():
        tmp=row.split(",")
        try:
            datalist.append((float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])))
        except:pass
mdb.models[model].MappedField(name='LOWER SURFACES', description='', 
        regionType=POINT, partLevelData=False, localCsys=None, 
        pointDataFormat=XYZ, fieldDataType=SCALAR, xyzPointData=datalist)
#===================================================== ATMOSPHERIC PRESSURE ==================================================================
# # Upper Surface
# assembly='NACA0012_U-1'
# mdb.models[model].Pressure(amplitude=UNSET, createStepName=new_step, 
#     distributionType=UNIFORM, field='', magnitude=P_atm, name='P_atm -- U', 
#     region=mdb.models[model].rootAssembly.instances[assembly].surfaces['UPPER SURFACES'])
# mdb.models[model].rootAssembly.regenerate()
# # Lower Surface
# assembly='NACA0012_L-1'
# mdb.models[model].Pressure(amplitude=UNSET, createStepName=new_step, 
#     distributionType=UNIFORM, field='', magnitude=P_atm, name='P_atm -- L', 
#     region=mdb.models[model].rootAssembly.instances[assembly].surfaces['LOWER SURFACES'])
# mdb.models[model].rootAssembly.regenerate()
#===================================================== DISTRIBUTED PRESSURE ==================================================================
mdb.models[model].rootAssembly.regenerate()
# Upper Surface
assembly='NACA0012_U-1'
mdb.models[model].Pressure(amplitude=UNSET, createStepName=new_step, 
    distributionType=FIELD, field='UPPER SURFACES', magnitude=1, name='UPPER SURFACES', 
    region=mdb.models[model].rootAssembly.instances[assembly].surfaces['UPPER SURFACES'])
mdb.models[model].rootAssembly.regenerate()
# Lower Surface
assembly='NACA0012_L-1'
mdb.models[model].Pressure(amplitude=UNSET, createStepName=new_step, 
    distributionType=FIELD, field='LOWER SURFACES', magnitude=1, name='LOWER SURFACES', 
    region=mdb.models[model].rootAssembly.instances[assembly].surfaces['LOWER SURFACES'])
mdb.models[model].rootAssembly.regenerate()
#==================================================== CENTRIFUGAL FORCE =======================================================================
# Upper Surface
assembly='NACA0012_U-1'
mdb.models[model].RotationalBodyForce(centrifugal=ON, createStepName=
    new_step, magnitude=rotor_speed, name='Centrifugal_U', point1=(c/2, -R_min, 0.0), 
    point2=(c/2, -R_min, 10.0), region=Region(
    cells=mdb.models[model].rootAssembly.instances[assembly].cells), rotaryAcceleration=OFF) ### CHANGE!!!
# Lower Surface
assembly='NACA0012_L-1'
mdb.models[model].RotationalBodyForce(centrifugal=ON, createStepName=
    new_step, magnitude=rotor_speed, name='Centrifugal_L', point1=(c/2, -R_min, 0.0), 
    point2=(c/2, -R_min, 10.0), region=Region(
    cells=mdb.models[model].rootAssembly.instances[assembly].cells), rotaryAcceleration=OFF)
#======================================================== SUBMIT JOB ===========================================================================
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=100, memoryUnits=PERCENTAGE, model=model, modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=job, nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    RESTART, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs[job].submit(consistencyChecking=OFF)
#====================================================== SAVE .cae FILE =======================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')