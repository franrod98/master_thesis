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
cae_file=sys.argv[-2]
job=sys.argv[-1]

#====================================================== OPEN .cae FILE =======================================================================
openMdb(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/{0}.cae'.format(cae_file))
#======================================================== SUBMIT JOB ===========================================================================
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=100, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=job, nodalOutputPrecision=SINGLE, 
    numCpus=2, numDomains=2, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs[job].submit(consistencyChecking=OFF)
#====================================================== SAVE .cae FILE ==========================================================================
mdb.saveAs(
    pathName='C:/Users/mario/Documents/T\xc9CNICO/5\xba ANO- 2\xba Semestre/TESE/C\xd3DIGO/ANACONDA/abaqus_testing/{0}.cae'.format(cae_file))