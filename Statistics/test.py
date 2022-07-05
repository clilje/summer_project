# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 10:56:06 2022

@author: clara
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from pathlib import Path
import pandas as pd
import os
import sys
sys.path.insert(1, '/home/clilje/summer_project')
from illustris_python import groupcat,lhalotree,snapshot,sublink,util
import time

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




x = 49
pdheader = ['ID','Type','x','y','z','mass','vx','vy','vz']
filename = 'HaloParticles50-1-pd/snap_99_halo_'+str(x)
snapnum = 99
basePath = '/disk01/rmcg/downloaded/tng/tng50-1'
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen']


subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))
#data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], subhalos['SubhaloMass']]).transpose()
halo_50 = [subhalos['SubhaloPos'][x, 0],subhalos['SubhaloPos'][x, 1],subhalos['SubhaloPos'][x, 2], subhalos['SubhaloHalfmassRad'][x], subhalos['SubhaloMass'][x], subhalos['SubhaloLen'][x]]
print(halo_50)

gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
starparts = snapshot.loadHalo(basePath, snapnum, x, 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
bhparts = snapshot.loadHalo(basePath, snapnum, x, 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
dmparts = snapshot.loadHalo(basePath, snapnum, x, 'dm', fields=['Coordinates','ParticleIDs','Velocities'])

partx = np.append(gasparts['Coordinates'][:,0],[starparts['Coordinates'][:,0],bhparts['Coordinates'][:,0],dmparts['Coordinates'][:,0]])
party = np.append(gasparts['Coordinates'][:,1],[starparts['Coordinates'][:,1],bhparts['Coordinates'][:,1],dmparts['Coordinates'][:,1]])
partz = np.append(gasparts['Coordinates'][:,2],[starparts['Coordinates'][:,2],bhparts['Coordinates'][:,2],dmparts['Coordinates'][:,2]])

print(len(partx))
print(type(partx))
print(partx.shape())
print(type(halo_50[0]))

x = 50
pdheader = ['ID','Type','x','y','z','mass','vx','vy','vz']
filename = 'HaloParticles50-1-pd/snap_99_halo_'+str(x)
snapnum = 99
basePath = '/disk01/rmcg/downloaded/tng/tng50-1'
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen']


subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))
#data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], subhalos['SubhaloMass']]).transpose()
halo_50 = [subhalos['SubhaloPos'][x, 0],subhalos['SubhaloPos'][x, 1],subhalos['SubhaloPos'][x, 2], subhalos['SubhaloHalfmassRad'][x], subhalos['SubhaloMass'][x], subhalos['SubhaloLen'][x]]
print(halo_50)

gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
starparts = snapshot.loadHalo(basePath, snapnum, x, 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
bhparts = snapshot.loadHalo(basePath, snapnum, x, 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
dmparts = snapshot.loadHalo(basePath, snapnum, x, 'dm', fields=['Coordinates','ParticleIDs','Velocities'])

partx = np.append(gasparts['Coordinates'][:,0],[starparts['Coordinates'][:,0],bhparts['Coordinates'][:,0],dmparts['Coordinates'][:,0]])
party = np.append(gasparts['Coordinates'][:,1],[starparts['Coordinates'][:,1],bhparts['Coordinates'][:,1],dmparts['Coordinates'][:,1]])
partz = np.append(gasparts['Coordinates'][:,2],[starparts['Coordinates'][:,2],bhparts['Coordinates'][:,2],dmparts['Coordinates'][:,2]])

print(len(partx))
#dis = distancefromcentre(halo_50[0], halo_50[1], halo_50[2], partx, party, partz)
#print(dis)
#print(np.min(dis))



x = 51
pdheader = ['ID','Type','x','y','z','mass','vx','vy','vz']
filename = 'HaloParticles50-1-pd/snap_99_halo_'+str(x)
snapnum = 99
basePath = '/disk01/rmcg/downloaded/tng/tng50-1'
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen']


subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))
#data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], subhalos['SubhaloMass']]).transpose()
halo_50 = [subhalos['SubhaloPos'][x, 0],subhalos['SubhaloPos'][x, 1],subhalos['SubhaloPos'][x, 2], subhalos['SubhaloHalfmassRad'][x], subhalos['SubhaloMass'][x], subhalos['SubhaloLen'][x]]
print(halo_50)

gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
starparts = snapshot.loadHalo(basePath, snapnum, x, 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
bhparts = snapshot.loadHalo(basePath, snapnum, x, 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
dmparts = snapshot.loadHalo(basePath, snapnum, x, 'dm', fields=['Coordinates','ParticleIDs','Velocities'])

partx = np.append(gasparts['Coordinates'][:,0],[starparts['Coordinates'][:,0],bhparts['Coordinates'][:,0],dmparts['Coordinates'][:,0]])
party = np.append(gasparts['Coordinates'][:,1],[starparts['Coordinates'][:,1],bhparts['Coordinates'][:,1],dmparts['Coordinates'][:,1]])
partz = np.append(gasparts['Coordinates'][:,2],[starparts['Coordinates'][:,2],bhparts['Coordinates'][:,2],dmparts['Coordinates'][:,2]])

print(len(partx))
#dis = distancefromcentre(halo_50[0], halo_50[1], halo_50[2], partx, party, partz)
#print(dis)
#print(np.min(dis))