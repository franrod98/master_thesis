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
R_min=float(sys.argv[-4])
L=float(sys.argv[-3])
Elem=int(float(sys.argv[-2]))
c=float(sys.argv[-1])
#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')
#====================================================== CROSS-SECTIONS =======================================================================
mdb.models['Model-1'].rootAssembly.regenerate()
cellSet_U='UPPER CELLS'
cellSet_L='LOWER CELLS'
faceSet_U='UPPER SURFACES'
faceSet_L='LOWER SURFACES'
for i in range(1,Elem+1):
    if (i!=Elem):
        # Upper Part
        part='NACA0012_U'
        datum_plane=mdb.models['Model-1'].parts[part].DatumPlaneByPrincipalPlane(
            offset=L*i, principalPlane=XZPLANE)
        plane_id=datum_plane.id
        mdb.models['Model-1'].parts[part].PartitionCellByDatumPlane(
            cells=mdb.models['Model-1'].parts[part].sets[cellSet_U].cells, datumPlane=mdb.models['Model-1'].parts[part].datums[plane_id])
        # Lower Part
        part='NACA0012_L'
        datum_plane=mdb.models['Model-1'].parts[part].DatumPlaneByPrincipalPlane(
            offset=L*i, principalPlane=XZPLANE)
        plane_id=datum_plane.id
        mdb.models['Model-1'].parts[part].PartitionCellByDatumPlane(
            cells=mdb.models['Model-1'].parts[part].sets[cellSet_L].cells, datumPlane=mdb.models['Model-1'].parts[part].datums[plane_id])
#======================================================== EDGE SETS ==========================================================================
    # Upper Part
    part='NACA0012_U'
    nr_faces=len(mdb.models['Model-1'].parts[part].sets[faceSet_U].faces)
    for j in range(0,nr_faces):
        face0=mdb.models['Model-1'].parts[part].sets[faceSet_U].faces[j]
        check=face0.pointOn
        upper_element='E{0}_U'.format(i)
        if (check[0][1]<L*i) and (check[0][1]>L*(i-1)):
            mdb.models['Model-1'].parts[part].Set(faces=
                mdb.models['Model-1'].parts[part].sets[faceSet_U].faces[j:j+1], name=upper_element)
            edges=mdb.models['Model-1'].parts[part].sets[upper_element].faces[0].getEdges()
            for k in edges:
                find_edge=edges=mdb.models['Model-1'].parts[part].edges[k].pointOn
                if round(find_edge[0][1],5)==round(L*(i-1),5):
                    mdb.models['Model-1'].parts[part].Set(edges=mdb.models['Model-1'].parts[part].edges[k:k+1],
                    name='Edge{0}_U'.format(i-1))
                elif round(find_edge[0][1],5)==round(L*i,5):
                    mdb.models['Model-1'].parts[part].Set(edges=mdb.models['Model-1'].parts[part].edges[k:k+1],
                    name='Edge{0}_U'.format(i))      
    # Lower Part
    part='NACA0012_L'
    nr_faces=len(mdb.models['Model-1'].parts[part].sets[faceSet_L].faces)
    for j in range(0,nr_faces):
        face0=mdb.models['Model-1'].parts[part].sets[faceSet_L].faces[j]
        check=face0.pointOn
        upper_element='E{0}_L'.format(i)
        if (check[0][1]<L*i) and (check[0][1]>L*(i-1)):
            mdb.models['Model-1'].parts[part].Set(faces=
                mdb.models['Model-1'].parts[part].sets[faceSet_L].faces[j:j+1], name=upper_element)
            edges=mdb.models['Model-1'].parts[part].sets[upper_element].faces[0].getEdges()
            for k in edges:
                find_edge=edges=mdb.models['Model-1'].parts[part].edges[k].pointOn
                if round(find_edge[0][1],5)==round(L*(i-1),5):
                    mdb.models['Model-1'].parts[part].Set(edges=mdb.models['Model-1'].parts[part].edges[k:k+1],
                    name='Edge{0}_L'.format(i-1))
                elif round(find_edge[0][1],5)==round(L*i,5):
                    mdb.models['Model-1'].parts[part].Set(edges=mdb.models['Model-1'].parts[part].edges[k:k+1],
                    name='Edge{0}_L'.format(i))  
mdb.models['Model-1'].rootAssembly.regenerate()
#====================================================== SAVE .cae FILE =======================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')