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
#======================================================== INPUT DATA =========================================================================
structures=int(float(sys.argv[-4]))
Elem=int(float(sys.argv[-3]))
seeds=int(float(sys.argv[-2]))
c=int(float(sys.argv[-1]))
#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')
#======================================================= GENARATE MESH =======================================================================
mdb.models['Model-1'].rootAssembly.regenerate()
###### Mesh refinement
# # Upper surface
# part='NACA0012_U'
# datum_plane=mdb.models['Model-1'].parts[part].DatumPlaneByPrincipalPlane(
#     offset=0.05*c, principalPlane=YZPLANE) #20/500=0.04
# plane_id=datum_plane.id
# mdb.models['Model-1'].parts[part].PartitionFaceByDatumPlane(
#     faces=mdb.models['Model-1'].parts[part].sets['UPPER SURFACES'].faces, datumPlane=mdb.models['Model-1'].parts[part].datums[plane_id])
# mdb.models['Model-1'].parts[part].PartitionFaceByDatumPlane(
#     faces=mdb.models['Model-1'].parts[part].sets['HUB_U'].faces, datumPlane=mdb.models['Model-1'].parts[part].datums[plane_id])
# # Lower surface
# part='NACA0012_L'
# datum_plane=mdb.models['Model-1'].parts[part].DatumPlaneByPrincipalPlane(
#     offset=0.05*c, principalPlane=YZPLANE) #20/500=0.04
# plane_id=datum_plane.id
# mdb.models['Model-1'].parts[part].PartitionFaceByDatumPlane(
#     faces=mdb.models['Model-1'].parts[part].sets['LOWER SURFACES'].faces, datumPlane=mdb.models['Model-1'].parts[part].datums[plane_id])
# mdb.models['Model-1'].parts[part].PartitionFaceByDatumPlane(
#     faces=mdb.models['Model-1'].parts[part].sets['HUB_L'].faces, datumPlane=mdb.models['Model-1'].parts[part].datums[plane_id])
# for i in range(0, Elem+1):
#     # Upper part
#     part='NACA0012_U'
#     edge_num=len(mdb.models['Model-1'].parts[part].sets['Edge{0}_U'.format(i)].edges)
#     for j in range(0, edge_num):
#         check=mdb.models['Model-1'].parts[part].sets['Edge{0}_U'.format(i)].edges[j].pointOn
#         if check[0][0]<0.04*c:
#             # mdb.models['Model-1'].parts[part].seedEdgeByNumber(edges=
#             #     (mdb.models['Model-1'].parts[part].sets['Edge{0}_U'.format(i)].edges[j],), number=5)
#             mdb.models['Model-1'].parts[part].seedEdgeBySize(edges=
#                 (mdb.models['Model-1'].parts[part].sets['Edge{0}_U'.format(i)].edges[j],), size=0.2*seeds)   
#     # Lower part
#     part='NACA0012_L'
#     edge_num=len(mdb.models['Model-1'].parts[part].sets['Edge{0}_L'.format(i)].edges)
#     for j in range(0, edge_num):
#         check=mdb.models['Model-1'].parts[part].sets['Edge{0}_L'.format(i)].edges[j].pointOn
#         if check[0][0]<0.04*c:    
#             # mdb.models['Model-1'].parts[part].seedEdgeByNumber(edges=
#             #     (mdb.models['Model-1'].parts[part].sets['Edge{0}_L'.format(i)].edges[j],), number=5)
#             mdb.models['Model-1'].parts[part].seedEdgeBySize(edges=
#                 (mdb.models['Model-1'].parts[part].sets['Edge{0}_L'.format(i)].edges[j],), size=0.2*seeds)
# #######
# # Mesh parts
# all_parts=mdb.models['Model-1'].parts.keys()
# for part in all_parts:
#     unmeshable=False
#     cell_num=len(mdb.models['Model-1'].parts[part].cells)
#     for i in range(0, cell_num):
#         control=mdb.models['Model-1'].parts[part].getMeshControl(region=
#             mdb.models['Model-1'].parts[part].cells[i], attribute=TECHNIQUE)  
#         if control==UNMESHABLE:
#             unmeshable=True
#     if unmeshable==True:
#         mdb.models['Model-1'].parts[part].setMeshControls(elemShape=
#             TET, regions=mdb.models['Model-1'].parts[part].cells, technique=FREE)
#         mdb.models['Model-1'].parts[part].setElementType(elemTypes=(ElemType(
# 		    elemCode=C3D10, elemLibrary=STANDARD), ), 
# 		    regions=Region(mdb.models['Model-1'].parts[part].cells))
#     else:
#         mdb.models['Model-1'].parts[part].setMeshControls(elemShape=
#             HEX, regions=mdb.models['Model-1'].parts[part].cells, technique=SWEEP)
#         mdb.models['Model-1'].parts[part].setElementType(elemTypes=(ElemType(
# 		    elemCode=C3D20, elemLibrary=STANDARD), ), 
# 		    regions=Region(mdb.models['Model-1'].parts[part].cells))
#     if part.startswith('Structure'):
#     ### DELETE ####
#       #   mdb.models['Model-1'].parts[part].setMeshControls(elemShape=
#       #       TET, regions=mdb.models['Model-1'].parts[part].cells, technique=FREE)
#       #   mdb.models['Model-1'].parts[part].setElementType(elemTypes=(ElemType(
# 		    # elemCode=C3D10, elemLibrary=STANDARD), ), 
# 		    # regions=Region(mdb.models['Model-1'].parts[part].cells))
#     ###############
#         mdb.models['Model-1'].parts[part].seedPart(deviationFactor=0.1, 
#             minSizeFactor=0.1, size=seeds) # CHANGE!!!
#         mdb.models['Model-1'].parts[part].generateMesh(meshTechniqueOverride=ON)
#     else:
#     ### DELETE ####
#         mdb.models['Model-1'].parts[part].setMeshControls(elemShape=
#             TET, regions=mdb.models['Model-1'].parts[part].cells, technique=FREE)
#         mdb.models['Model-1'].parts[part].setElementType(elemTypes=(ElemType(
# 		    elemCode=C3D10, elemLibrary=STANDARD), ), 
# 		    regions=Region(mdb.models['Model-1'].parts[part].cells))
#     ###############
#         mdb.models['Model-1'].parts[part].seedPart(deviationFactor=0.1, 
#             minSizeFactor=0.1, size=seeds)
#     	mdb.models['Model-1'].parts[part].generateMesh()
#### CONVERGENCE (DELETE) ##########################
all_parts=mdb.models['Model-1'].parts.keys()
for part in all_parts:
    mdb.models['Model-1'].parts[part].seedPart(deviationFactor=0.1, 
    	minSizeFactor=0.1, size=seeds)
    mdb.models['Model-1'].parts[part].setMeshControls(elemShape=
        TET, regions=mdb.models['Model-1'].parts[part].cells, technique=FREE)
    mdb.models['Model-1'].parts[part].setElementType(elemTypes=(ElemType(
	    elemCode=C3D10, elemLibrary=STANDARD), ), 
	    regions=Region(mdb.models['Model-1'].parts[part].cells))
    mdb.models['Model-1'].parts[part].generateMesh()
############
#====================================================== STORE EDGE DATA ======================================================================
for i in range(0,Elem+1):
    # Upper Surface
    part='NACA0012_U'
    edge_U='Edge{0}_U'.format(i)
    edge_file='EDGE{0}_U.txt'.format(i)
    nr_nodes=len(mdb.models['Model-1'].parts[part].sets[edge_U].nodes)
    text='This is a headliner'
    for j in range(0, nr_nodes):
        coords=mdb.models['Model-1'].parts[part].sets[edge_U].nodes[j].coordinates
        label=mdb.models['Model-1'].parts[part].sets[edge_U].nodes[j].label
        output=str(coords[0]/c)+'\t'+str(coords[1]/c)+'\t'+str(coords[2]/c)
        text='\n'.join((text,output))
    store=file(edge_file, 'w')
    store.write(text)
    store.close()
    # Lower Surface
    part='NACA0012_L'
    edge_L='Edge{0}_L'.format(i)
    edge_file='EDGE{0}_L.txt'.format(i)
    nr_nodes=len(mdb.models['Model-1'].parts[part].sets[edge_L].nodes)
    text='This is a headliner'
    for j in range(0, nr_nodes):
        coords=mdb.models['Model-1'].parts[part].sets[edge_L].nodes[j].coordinates
        label=mdb.models['Model-1'].parts[part].sets[edge_L].nodes[j].label
        output=str(coords[0]/c)+'\t'+str(coords[1]/c)+'\t'+str(coords[2]/c)
        text='\n'.join((text,output))
    store=file(edge_file, 'w')
    store.write(text)
    store.close()
#====================================================== SAVE .cae FILE ==========================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')    