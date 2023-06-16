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
import numpy as np
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')
# #======================================================== INPUT DATA =========================================================================
# interactions=int(float(sys.argv[-1]))
#======================================================= INTERACTIONS ========================================================================
all_instances=mdb.models['Model-1'].rootAssembly.instances.keys()
all_master=[]
# Count number of master surfaces
for instance in all_instances:
	all_surfaces=mdb.models['Model-1'].rootAssembly.instances[instance].surfaces.keys()
	for surface in all_surfaces:
		if surface.startswith('Master'):
			emp_str=''
			for i in surface:
				if i.isdigit():
					emp_str = emp_str+i
			all_master.append(emp_str)
all_master=np.array(all_master)
all_master = all_master.astype(np.int)
####
# mdb.models['Model-1'].ContactProperty('IntProp-1')
# mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
#     formulation=FRICTIONLESS)
####
for i in all_master:
	for instance in all_instances:
		all_surfaces=mdb.models['Model-1'].rootAssembly.instances[instance].surfaces.keys()
		if "Master{0}".format(i) in all_surfaces:
			master=mdb.models['Model-1'].rootAssembly.instances[instance].surfaces['Master{0}'.format(i)]
		elif "Slave{0}".format(i) in all_surfaces:
			slave=mdb.models['Model-1'].rootAssembly.instances[instance].surfaces["Slave{0}".format(i)]
	####
	# mdb.models['Model-1'].SurfaceToSurfaceContactStd(adjustMethod=NONE, 
 #    clearanceRegion=None, createStepName='Step-1', datumAxis=None, 
 #    initialClearance=0.0, interactionProperty='IntProp-1', master=master,
 #    name='Int-{0}'.format(i), slave=slave, sliding=SMALL, thickness=ON)
    #### 
	mdb.models['Model-1'].Tie(adjust=ON, master=master,
		name='Constraint-{0}'.format(i), positionToleranceMethod=COMPUTED, 
		slave=slave, thickness=ON, tieRotations=ON)
#====================================================== SAVE .cae FILE ==========================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')