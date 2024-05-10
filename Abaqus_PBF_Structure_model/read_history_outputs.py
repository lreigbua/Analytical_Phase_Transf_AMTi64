# Read history data | Authors: Luis Reig
# Abaqus Python 2.7

# This script was developed to read history output variables from a given element and save them to a .csv file. 
# The job name and the desired elements must be specified. 

from asyncore import write
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
import random
from array import *
from odbAccess import openOdb
from odbAccess import *
from abaqus import * 
from abaqusConstants import *
import odbAccess
import math as m
import numpy as np    
import os        # Operating system
import shutil    # copying or moving files
from caeModules import *
from driverUtils import executeOnCaeStartup
import shutil

# SPECIFY JOB NAME! without .odb
Job_name='structure_model'

odbname=Job_name
path='./'                    # set odb path here (if in working dir no need to change!)
myodbpath=path+odbname+'.odb'    
odb=openOdb(path=myodbpath)

#Specify elements to study
#such as 'Element "instance name"."element number" Int Point "int point number"'
#Instance name must be in caps!
element='Element PART-1-1.1 Int Point 1'

step_name = 'Step-1'

step=odb.steps[step_name]
region=step.historyRegions[element]

#extract data
TEMP=region.historyOutputs['TEMP'].data
alpha=region.historyOutputs['SDV_fPhase_Alpha'].data
beta=region.historyOutputs['SDV_fPhase_Beta'].data
mart=region.historyOutputs['SDV_fPhase_Alphaprime'].data

##print to csv

#temp
np.savetxt('TEMP.csv', np.array(TEMP) ,delimiter=',', fmt = '%f')

#fractions
np.savetxt('alpha.csv', np.array(alpha),delimiter=',', fmt = '%f')
np.savetxt('Beta.csv', np.array(beta),delimiter=',', fmt = '%f')
np.savetxt('mart.csv', np.array(mart),delimiter=',', fmt = '%f')

odb.close()

