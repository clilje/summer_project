# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 09:40:18 2022

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


#Read in basic Subhalo Group Info
subhalo_info = pd.read_csv('../50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy()



xx = [2,53,201,1000]
for g in xx:
    
    #read in particle data
    chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['x'],dtype={'x':object})
    chunk = chunk[chunk['x'] != 'x']
    partx = chunk['x'].to_numpy().astype(float)
    #print(partx)
    chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['y'],dtype={'y':object})
    chunk = chunk[chunk['y'] != 'y']
    party = chunk['y'].to_numpy().astype(float)
    chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['z'],dtype={'z':object})
    chunk = chunk[chunk['z'] != 'z']
    partz = chunk['z'].to_numpy().astype(float)
    chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['mass'],dtype={'mass':object})
    chunk = chunk[chunk['mass'] != 'mass']
    mass = chunk['mass'].to_numpy().astype(float)
    
    
    
    fig = plt.figure()
    ax = plt.axes(projection ='3d')
    xyz = np.arange(len(partx))
    index = np.random.choice(xyz,200)
    ax.scatter(partx[index], party[index], partz[index], marker='+',color='blue',alpha=0.1)
    #ax.scatter(halo_50[0], halo_50[1], halo_50[2], marker='+',color='red')
    ax.scatter(positionsX[g],positionsY[g], positionsZ[g],marker='x', color='black')
    #ax.scatter(CM[0], CM[1], CM[2], marker='+',color='pink')
    
    ax.set_xlabel('x [ckpc/h]')
    
    ax.set_ylabel('y [ckpc/h]')
    ax.set_zlabel('z [ckpc/h]')
    fig.savefig('halocomp-npdarkcentre-'+str(g))