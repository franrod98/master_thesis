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
from odbAccess import *
import numpy as np
#======================================================== INPUT DATA =========================================================================
job=sys.argv[-4]
iteration=int(float(sys.argv[-3]))
Elem=int(float(sys.argv[-2]))
c=float(sys.argv[-1])
#====================================================== OPEN .odb FILE =======================================================================
odb = openOdb(path='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/{0}.odb'.format(job), readOnly=OFF)
# #====================================================== OPEN .cae FILE =======================================================================
# openMdb(
#     pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/NACA0012_analysis.cae')
#======================================================= POST-PROCESS ========================================================================
step='Step-{0}'.format(iteration)
step2read=odb.steps[step]
lastFrame=step2read.frames[-1]
displacement=lastFrame.fieldOutputs['U']
coord=lastFrame.fieldOutputs['COORD']
stress=lastFrame.fieldOutputs['S']
# Upper Surface
assembly='NACA0012_U-1'
for i in range(0,Elem+1):
	text='This is a headliner'
	name='EDGE{0}_U'.format(i)
	# edge=displacement.getSubset(region=odb.rootAssembly.instances[assembly].nodeSets[name])
	edge=coord.getSubset(region=odb.rootAssembly.instances[assembly].nodeSets[name])
	# node=0
	for e in edge.values:
		# coord_x=odb.rootAssembly.instances[assembly].nodeSets[name].nodes[node].coordinates[0]
		# coord_y=odb.rootAssembly.instances[assembly].nodeSets[name].nodes[node].coordinates[1]
		# coord_z=odb.rootAssembly.instances[assembly].nodeSets[name].nodes[node].coordinates[2]
		# output =str((coord_x+e.data[0])/c)+'\t'+str((coord_y+e.data[1])/c)+'\t'+str((coord_z+e.data[2])/c)
		output =str(e.data[0]/c)+'\t'+str(e.data[1]/c)+'\t'+str(e.data[2]/c)
		text='\n'.join((text,output))
		# node=node+1
	data=file('{0}.txt'.format(name), 'w')
	data.write(text)
	data.close()
# Lower Surface
assembly='NACA0012_L-1'
for i in range(0,Elem+1):
	text='This is a headliner'
	name='EDGE{0}_L'.format(i)
	# edge=displacement.getSubset(region=odb.rootAssembly.instances[assembly].nodeSets[name])
	edge=coord.getSubset(region=odb.rootAssembly.instances[assembly].nodeSets[name])
	# node=0
	for e in edge.values:
		# coord_x=odb.rootAssembly.instances[assembly].nodeSets[name].nodes[node].coordinates[0]
		# coord_y=odb.rootAssembly.instances[assembly].nodeSets[name].nodes[node].coordinates[1]
		# coord_z=odb.rootAssembly.instances[assembly].nodeSets[name].nodes[node].coordinates[2]
		# output =str((coord_x+e.data[0])/c)+'\t'+str((coord_y+e.data[1])/c)+'\t'+str((coord_z+e.data[2])/c)
		output =str(e.data[0]/c)+'\t'+str(e.data[1]/c)+'\t'+str(e.data[2]/c)
		text='\n'.join((text,output))
		# node=node+1
	data=file('{0}.txt'.format(name), 'w')
	data.write(text)
	data.close()