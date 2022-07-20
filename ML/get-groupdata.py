"""
Created on Wed Jul  6 22:48:58 2022
This file was actually used to create the group catalogue file of all halos.

@author: clara
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from pathlib import Path
import sys
sys.path.insert(1, '/home/clilje/summer_project')
from illustris_python import groupcat,lhalotree,snapshot,sublink,util

size = 50
res = 1
dark = False
#z = 0
snapnum = 99

if dark == True:
    basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
else:
   basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen','SubhaloMassType']


subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift
#print(a)
#print(z)

#print(groupcat.loadHeader(basePath, snapnum))
print(subhalos['SubhaloMassType'])

with open('50-1-subhalo-info.csv', 'w', encoding='UTF8', newline='') as subfile:
    header = ['SubhaloIndex','SubhaloPosX','SubhaloPosY','SubhaloPosZ','SubhaloHalfmassRad',
              'SubhaloMass','SubhaloLen', 'SubhaloGasMass','SubhaloDMMass','SubhaloStarMass',
              'SubhaloBHMass','SubhaloSpin','SubhaloVelDisp','SubhaloVmax']
    fwriter = csv.writer(subfile, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],
                      subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], 
                      subhalos['SubhaloMass'], subhalos['SubhaloLen'], 
                      subhalos['SubhaloMassType'][:, 0], subhalos['SubhaloMassType'][:, 1], 
                      subhalos['SubhaloMassType'][:, 4], subhalos['SubhaloMassType'][:, 5],
                      subhalos['SubhaloSpin'],subhalos['SubhaloVelDisp'],subhalos['SubhaloVmax']]).transpose()
    fwriter.writerows(data)
