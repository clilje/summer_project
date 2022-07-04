import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
from pathlib import Path
import pandas as pd
import sys
sys.path.insert(1, '/home/clilje/summer_project')
from illustris_python import groupcat,lhalotree,snapshot,sublink,util
import time

tt = time.time()
size = 50
res = 1
dark = False
#z = 0
snapnum = 99

if dark == True:
    basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
else:
   basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass']

print(round(time.time()-tt,2))
tt = time.time()
subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift
#print(a)
#print(z)
print(round(time.time()-tt,2))
tt = time.time()
#print(groupcat.loadHeader(basePath, snapnum))

x = 33981
pdheader = ['ID','Type','x','y','z','mass','vx','vy','vz']
while x <= (len(num_halo)):
    if dark == False: 
        filename = 'HaloParticles50-1-pd/snap_99_halo_'+str(x)
    else: 
        filename = 'HaloParticles50-1-pd/snap_99_halo_'+str(x)+'-dark'
    
    print(round(time.time()-tt,2))
    tt = time.time()
    
    if dark == False:
        gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
        starparts = snapshot.loadHalo(basePath, snapnum, x, 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
        bhparts = snapshot.loadHalo(basePath, snapnum, x, 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
    dmparts = snapshot.loadHalo(basePath, snapnum, x, 'dm', fields=['Coordinates','ParticleIDs','Velocities'])
    
    print(round(time.time()-tt,2))
    tt = time.time()
    
    with h5py.File(snapshot.snapPath(basePath,snapnum),'r') as f:
        header = dict( f['Header'].attrs.items() )
        print(header['MassTable'][1]) # 10^10 msun/h
        dmmass = [(header['MassTable'][1]).astype('float')]*len(dmparts['ParticleIDs'])
    
    lowerbound = 0
    indexjump = 900000
    if gasparts['count'] >= starparts['count'] and gasparts['count'] >= bhparts['count'] and gasparts['count'] >= dmparts['count']:
        max_num_part = gasparts['count']
        print("gas dominates")
    if dmparts['count'] >= gasparts['count'] and dmparts['count'] >= bhparts['count'] and dmparts['count'] >= starparts['count']:
        max_num_part = dmparts['count']
        print("dm dominates")
    if starparts['count'] >= gasparts['count'] and starparts['count'] >= bhparts['count'] and starparts['count'] >= dmparts['count']:
        max_num_part = starparts['count']
        print("star dominates")
    else: 
        max_num_part = bhparts['count']
        print("bh dominates")
        
        
        
        
    derek = pd.DataFrame(columns=pdheader)
    
    while lowerbound < max_num_part:
        
        #print(lowerbound)
        
        miniderek = pd.DataFrame(columns=pdheader)
        if (lowerbound+indexjump)>len(gasparts['ParticleIDs']):
            upperbound = max_num_part
        else:
            upperbound = lowerbound+indexjump
        
        if dark == False:
            if gasparts['count'] != 0:
                #Gas Particles
                miniderek['ID']=gasparts['ParticleIDs'][lowerbound:upperbound]
                miniderek['Type']=['gas']*len(gasparts['ParticleIDs'][lowerbound:upperbound])
                miniderek['x']=gasparts['Coordinates'][:, 0][lowerbound:upperbound]
                miniderek['y']=gasparts['Coordinates'][:, 1][lowerbound:upperbound]
                miniderek['z']=gasparts['Coordinates'][:, 2][lowerbound:upperbound]
                miniderek['mass']=gasparts['Masses'][lowerbound:upperbound]
                miniderek['vx']=gasparts['Velocities'][:, 0][lowerbound:upperbound]
                miniderek['vy']=gasparts['Velocities'][:, 1][lowerbound:upperbound]
                miniderek['vz']=gasparts['Velocities'][:, 2][lowerbound:upperbound]
                
                
                #print(miniderek)
                
                derek = pd.concat([derek,miniderek])
                miniderek = miniderek[0:0]
            
            if starparts['count'] != 0:
                #Star Particles
                miniderek['ID']=starparts['ParticleIDs'][lowerbound:upperbound]
                miniderek['Type']=['star']*len(starparts['ParticleIDs'][lowerbound:upperbound])
                miniderek['x']=starparts['Coordinates'][:, 0][lowerbound:upperbound]
                miniderek['y']=starparts['Coordinates'][:, 1][lowerbound:upperbound]
                miniderek['z']=starparts['Coordinates'][:, 2][lowerbound:upperbound]
                miniderek['mass']=starparts['Masses'][lowerbound:upperbound]
                miniderek['vx']=starparts['Velocities'][:, 0][lowerbound:upperbound]
                miniderek['vy']=starparts['Velocities'][:, 1][lowerbound:upperbound]
                miniderek['vz']=starparts['Velocities'][:, 2][lowerbound:upperbound]
                
                derek = pd.concat([derek,miniderek])
                miniderek = miniderek[0:0]
            
            #print(bhparts)
            #print(bhparts.keys())
            if bhparts['count'] != 0:
                #Black Holes
                miniderek['ID']=bhparts['ParticleIDs'][lowerbound:upperbound]
                miniderek['Type']=['bh']*len(bhparts['ParticleIDs'][lowerbound:upperbound])
                miniderek['x']=bhparts['Coordinates'][:, 0][lowerbound:upperbound]
                miniderek['y']=bhparts['Coordinates'][:, 1][lowerbound:upperbound]
                miniderek['z']=bhparts['Coordinates'][:, 2][lowerbound:upperbound]
                miniderek['mass']=bhparts['Masses'][lowerbound:upperbound]
                miniderek['vx']=bhparts['Velocities'][:, 0][lowerbound:upperbound]
                miniderek['vy']=bhparts['Velocities'][:, 1][lowerbound:upperbound]
                miniderek['vz']=bhparts['Velocities'][:, 2][lowerbound:upperbound]
                derek = pd.concat([derek,miniderek])
                miniderek =miniderek[0:0]

        #DM Particles
        if dmparts['count'] != 0:
            miniderek['ID']=dmparts['ParticleIDs'][lowerbound:upperbound]
            miniderek['Type']=['dm']*len(dmparts['ParticleIDs'][lowerbound:upperbound])
            miniderek['x']=dmparts['Coordinates'][:, 0][lowerbound:upperbound]
            miniderek['y']=dmparts['Coordinates'][:, 1][lowerbound:upperbound]
            miniderek['z']=dmparts['Coordinates'][:, 2][lowerbound:upperbound]
            miniderek['mass']=dmmass[lowerbound:upperbound]
            miniderek['vx']=dmparts['Velocities'][:, 0][lowerbound:upperbound]
            miniderek['vy']=dmparts['Velocities'][:, 1][lowerbound:upperbound]
            miniderek['vz']=dmparts['Velocities'][:, 2][lowerbound:upperbound]
            derek = pd.concat([derek,miniderek])
            miniderek = miniderek[0:0]

        lowerbound = upperbound
        #derek = derek.reset_index(drop=True)
        derek.to_csv(filename+'.csv', mode='a')
        #derek.to_hdf(filename+'.hdf', mode='a')
        derek = derek[0:0]
    print(x)
    x +=1