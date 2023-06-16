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
j=int(float(sys.argv[-5]))
material=sys.argv[-4]
material_rho=float(sys.argv[-3])
E=float(sys.argv[-2])
poisson=float(sys.argv[-1])
#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')
#====================================================== CREATE MATERIAL ========================================================================
create=True
sections=mdb.models['Model-1'].sections.keys()
for i in sections:
	if material==i:
		create=False
if create==True:
	mdb.models['Model-1'].Material(name=material)
	mdb.models['Model-1'].materials[material].Density(table=((material_rho, ), ))
	mdb.models['Model-1'].materials[material].Elastic(table=((E, poisson), ))
	mdb.models['Model-1'].HomogeneousSolidSection(material=material, name=
	    material, thickness=None)
#====================================================== ASSIGN MATERIAL ========================================================================
mdb.models['Model-1'].parts['Structure{0}'.format(j)].SectionAssignment(offset=0.0, 
    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
    cells=mdb.models['Model-1'].parts['Structure{0}'.format(j)].cells), sectionName=material, 
    thicknessAssignment=FROM_SECTION)
#====================================================== SAVE .cae FILE =======================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')

