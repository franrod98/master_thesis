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
step=sys.argv[-5]
structures=int(float(sys.argv[-4]))
R_min=float(sys.argv[-3])
L=float(sys.argv[-2])
Elem=int(float(sys.argv[-1]))
#======================================================== IMPORT PARTS ========================================================================
mdb.openStep('C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/{0}.step'.format(step), scaleFromFile=OFF)
# Structures
for i in range(0,structures):
    mdb.models['Model-1'].PartFromGeometryFile(bodyNum=i+1, combine=False, 
        dimensionality=THREE_D, geometryFile=mdb.acis, name=
        'Structure{0}'.format(i), type=DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Structure{0}'.format(i)].Set(cells=mdb.models['Model-1'].parts['Structure{0}'.format(i)].cells,
        name='STR{0}_CELLS'.format(i))
    mdb.models['Model-1'].parts['Structure{0}'.format(i)].Set(edges=mdb.models['Model-1'].parts['Structure{0}'.format(i)].edges,
        name='STR{0}_EDGES'.format(i))
# Upper Part
upperPart=mdb.models['Model-1'].PartFromGeometryFile(bodyNum=structures+1, combine=False, 
    dimensionality=THREE_D, geometryFile=mdb.acis, name='NACA0012_U', type=DEFORMABLE_BODY)
# Lower Part
lowerPart=mdb.models['Model-1'].PartFromGeometryFile(bodyNum=structures+2, combine=False, 
    dimensionality=THREE_D, geometryFile=mdb.acis, name='NACA0012_L', type=DEFORMABLE_BODY)
part_U=upperPart.name
part_L=lowerPart.name
#===================================================== RECOMPUTE GEOMETRY =====================================================================
all_parts=mdb.models['Model-1'].parts.keys()
for i in all_parts:
    mdb.models['Model-1'].parts[i].ConvertToPrecise(method=RECOMPUTE_GEOMETRY)
for i in all_parts:
    mdb.models['Model-1'].parts[i].checkGeometry()
mdb.models['Model-1'].rootAssembly.regenerate()
#===================================================== SEPARATE FROM HUB =====================================================================
# for i in all_parts:
#     datum_plane=mdb.models['Model-1'].parts[i].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
#     plane_id=datum_plane.id
#     mdb.models['Model-1'].parts[i].PartitionCellByDatumPlane(cells=mdb.models['Model-1'].parts[i].cells,
#         datumPlane=mdb.models['Model-1'].parts[i].datums[plane_id])
# Upper Part
datum_plane=mdb.models['Model-1'].parts[part_U].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
plane_id=datum_plane.id
mdb.models['Model-1'].parts[part_U].PartitionCellByDatumPlane(cells=mdb.models['Model-1'].parts[part_U].cells,
    datumPlane=mdb.models['Model-1'].parts[part_U].datums[plane_id])
# Lower Part
datum_plane=mdb.models['Model-1'].parts[part_L].DatumPlaneByPrincipalPlane(offset=0.0, principalPlane=XZPLANE)
plane_id=datum_plane.id
mdb.models['Model-1'].parts[part_L].PartitionCellByDatumPlane(cells=mdb.models['Model-1'].parts[part_L].cells,
    datumPlane=mdb.models['Model-1'].parts[part_L].datums[plane_id])
#======================================== UPPER AND LOWER FACE (NACA0012 + HUB) SETS =========================================================
# Upper Surface
part=part_U
name='UPPER SURFACES'
HUB='HUB_U'
edge_num=len(mdb.models['Model-1'].parts[part].edges)
for i in range(0,edge_num):
    check=mdb.models['Model-1'].parts[part].edges[i].pointOn
    if round(check[0][0],3)==0 and round(check[0][2],4)==0:
        edge=mdb.models['Model-1'].parts[part].edges[i].getSize()
        # NACA 0012
        if round(edge,0)==round(L*Elem,0):
            faces_id=mdb.models['Model-1'].parts[part].edges[i].getFaces()
            for j in faces_id:
                normal=mdb.models['Model-1'].parts[part].faces[j].getNormal()
                if normal[2]>0:
                    mdb.models['Model-1'].parts[part].Set(faces=mdb.models['Model-1'].parts[part].faces[j:j+1], name=name)
        # HUB
        elif round(edge,0)==round(R_min,0):
            faces_id=mdb.models['Model-1'].parts[part].edges[i].getFaces()
            for j in faces_id:
                normal=mdb.models['Model-1'].parts[part].faces[j].getNormal()
                if normal[2]>0:
                    mdb.models['Model-1'].parts[part].Set(faces=mdb.models['Model-1'].parts[part].faces[j:j+1], name=HUB)
# Lower Surface
part=part_L
name='LOWER SURFACES'
HUB='HUB_L'
edge_num=len(mdb.models['Model-1'].parts[part].edges)
for i in range(0,edge_num):
    check=mdb.models['Model-1'].parts[part].edges[i].pointOn
    if round(check[0][0],3)==0 and round(check[0][2],4)==0:
        edge=mdb.models['Model-1'].parts[part].edges[i].getSize()
        # NACA 0012
        if round(edge,0)==round(L*Elem,0):
            faces_id=mdb.models['Model-1'].parts[part].edges[i].getFaces()
            for j in faces_id:
                normal=mdb.models['Model-1'].parts[part].faces[j].getNormal()
                if normal[2]<0:
                    mdb.models['Model-1'].parts[part].Set(faces=mdb.models['Model-1'].parts[part].faces[j:j+1], name=name)
        # HUB
        elif round(edge,0)==round(R_min,0):
            faces_id=mdb.models['Model-1'].parts[part].edges[i].getFaces()
            for j in faces_id:
                normal=mdb.models['Model-1'].parts[part].faces[j].getNormal()
                if normal[2]<0:
                    mdb.models['Model-1'].parts[part].Set(faces=mdb.models['Model-1'].parts[part].faces[j:j+1], name=HUB)
#=============================================== UPPER AND LOWER SURFACE SETS ================================================================
# Upper Surface
mdb.models['Model-1'].parts[part_U].Surface(side1Faces=
    mdb.models['Model-1'].parts[part_U].sets['UPPER SURFACES'].faces, name='UPPER SURFACES')
# Lower Surface
mdb.models['Model-1'].parts[part_L].Surface(side1Faces=
    mdb.models['Model-1'].parts[part_L].sets['LOWER SURFACES'].faces, name='LOWER SURFACES')
#================================================ UPPER AND LOWER CELL SETS ==================================================================
# Upper Surface
face=mdb.models['Model-1'].parts[part_U].sets['UPPER SURFACES'].faces[0]
find_cell=face.getCells()
cell=find_cell[0]
mdb.models['Model-1'].parts[part_U].Set(cells=
    mdb.models['Model-1'].parts[part_U].cells[cell:cell+1], name='UPPER CELLS')
# Lower
face=mdb.models['Model-1'].parts[part_L].sets['LOWER SURFACES'].faces[0]
find_cell=face.getCells()
cell=find_cell[0]
mdb.models['Model-1'].parts[part_L].Set(cells=
    mdb.models['Model-1'].parts[part_L].cells[cell:cell+1], name='LOWER CELLS')
#===================================================== PART INSTANCES ========================================================================
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
all_parts=mdb.models['Model-1'].parts.keys()
for i in all_parts:
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name=i+'-1', part=mdb.models['Model-1'].parts[i])
#============================================================ STEP ===========================================================================

# mdb.models['Model-1'].FrequencyStep(acousticCoupling=AC_OFF, eigensolver=
#     SUBSPACE, maxIterations=250, name='Step-1', normalization=MASS, numEigen=10, 
#     previous='Initial', simLinearDynamics=OFF, vectors=18)

mdb.models['Model-1'].StaticStep(initialInc=0.01, minInc=1e-06, name='Step-1', nlgeom=ON, previous='Initial')
# mdb.models['Model-1'].steps['Step-1'].setValues(adaptiveDampingRatio=0.05, 
#     continueDampingFactors=False, stabilizationMagnitude=0.0002, 
#     stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(frequency=2, 
    variables=('S', 'E', 'EE', 'LE', 'U', 'UT', 'UR', 'RF', 'CF', 'P', 'COORD'), position=NODES)
#====================================================== RESTART REQUEST ======================================================================
mdb.models['Model-1'].steps['Step-1'].Restart(frequency=2, numberIntervals=0, 
    overlay=ON, timeMarks=OFF)
#====================================================== HUB CONSTRAINTS ======================================================================
constraint='Fix HUB'
mdb.models['Model-1'].EncastreBC(createStepName='Step-1', localCsys=None, name=constraint+'_U',
    region=mdb.models['Model-1'].rootAssembly.instances['NACA0012_U-1'].sets['HUB_U'])
mdb.models['Model-1'].EncastreBC(createStepName='Step-1', localCsys=None, name=constraint+'_L',
    region=mdb.models['Model-1'].rootAssembly.instances['NACA0012_L-1'].sets['HUB_L'])
#====================================================== SAVE .cae FILE =======================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')