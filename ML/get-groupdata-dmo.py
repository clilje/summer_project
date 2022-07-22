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

def distancefromcentre(cx, cy, cz, x, y, z, ):
    """
    

    Parameters
    ----------
    cx : Halo centre x-coord float
    cy : Halo centre y-coord float
    cz : Halo centre z-coord float
    x : Particle Position x-coord
    y : Particle Position y-coord
    z : Particle Position z-coord

    Returns
    -------
    Distance between particle and halo centre

    """
    return (np.sqrt((np.power(np.subtract(x,cx), 2)+ np.power(np.subtract(y,cy), 2) + np.power(np.subtract(z,cz), 2)))) # distance between the centre and given point



size = 50
res = 1
dark = True
#z = 0
snapnum = 99

if dark == True:
    basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
else:
   basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen','SubhaloMassType','SubhaloSpin',
          'SubhaloVelDisp','SubhaloVmax', 'SubhaloSFR','SubhaloGrNr']
foffields = ['GroupPos','GroupMass']

subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
fof_halos = groupcat.loadHalos(basePath,99,fields=foffields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift
#print(a)
#print(z)

radial_distance = distancefromcentre(fof_halos['GroupPos'][subhalos['SubhaloGrNr']][0], 
                                     fof_halos['GroupPos'][subhalos['SubhaloGrNr']][1], 
                                     fof_halos['GroupPos'][subhalos['SubhaloGrNr']][2], 
                                     subhalos['SubhaloPos'][:, 0], subhalos['SubhaloPos'][:, 1], 
                                     subhalos['SubhaloPos'][:, 2])


#print(groupcat.loadHeader(basePath, snapnum))
print(subhalos['SubhaloMassType'])

with open('50-1-subhalo-info-dark.csv', 'w', encoding='UTF8', newline='') as subfile:
    header = ['SubhaloIndex','SubhaloPosX','SubhaloPosY','SubhaloPosZ','SubhaloHalfmassRad','SubhaloMass',
              'SubhaloLen', 'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp','SubhaloVmax',
              'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter']
    fwriter = csv.writer(subfile, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],
                      subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], 
                      subhalos['SubhaloMass'], subhalos['SubhaloLen'], subhalos['SubhaloMassType'][:, 1], 
                      subhalos['SubhaloSpin'][:, 0],subhalos['SubhaloSpin'][:, 1],subhalos['SubhaloSpin'][:, 2],
                      subhalos['SubhaloVelDisp'],subhalos['SubhaloVmax'],subhalos['SubhaloSFR'], 
                      fof_halos['GroupMass'][subhalos['SubhaloGrNr']],radial_distance]).transpose()
    fwriter.writerows(data)
