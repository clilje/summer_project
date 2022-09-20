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
dark = False
#z = 0
snapnum = np.arange(9,100,10)
for i in snapnum:
    if dark == True:
        basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
    else:
       basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 
    fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen','SubhaloMassType',
              'SubhaloSpin','SubhaloVelDisp','SubhaloVmax','SubhaloBHMdot', 'SubhaloSFR','SubhaloGrNr']
    foffields = ['GroupPos','GroupMass']
    
    subhalos = groupcat.loadSubhalos(basePath,i,fields=fields)
    fof_halos = groupcat.loadHalos(basePath,i,fields=foffields)
    num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))
    
    a = (groupcat.loadHeader(basePath, i)['Time']).astype('int')  #scalefactor
    z = (groupcat.loadHeader(basePath, i)['Redshift']).astype('float')  #redshift
    print(i)
    print(z)