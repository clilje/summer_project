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
"""
with open('50-1-subhalo-info.csv', 'w', encoding='UTF8', newline='') as subfile:
    header = ['SubhaloIndex','SubhaloPosX','SubhaloPosY','SubhaloPosZ','SubhaloHalfmassRad','SubhaloMass']
    fwriter = csv.writer(subfile, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], subhalos['SubhaloMass']]).transpose()
    fwriter.writerows(data)
"""
#gasparts = snapshot.loadHalo(basePath, snapnum, 0, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
#print(gasparts['Coordinates'])
x = 0
pdheader = ['ID','Type','x','y','z','mass','vx','vy','vz']
while x <= (len(num_halo)):
    if dark == False: 
        filename = 'HaloParticles50-1/snap_99_halo_'+str(x)
    else: 
        filename = 'HaloParticles50-1/snap_99_halo_'+str(x)+'-dark'
    #with open(filename, 'w', encoding='UTF8', newline='') as f:
    
    print(round(time.time()-tt,2))
    tt = time.time()
    
    if dark == False:
        gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
        #gasparts = pd.DataFrame.from_dict(gasparts)
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
    indexjump = int(len(gasparts['ParticleIDs'])/300)
    
    derek = pd.DataFrame(columns=pdheader)
    
    while lowerbound < len(gasparts['ParticleIDs']):
        
        print(lowerbound)
        
        miniderek = pd.DataFrame(columns=pdheader)
        if (lowerbound+indexjump)>len(gasparts['ParticleIDs']):
            upperbound = len(gasparts['ParticleIDs'])
        else:
            upperbound = lowerbound+indexjump
        
        if dark == False:
            
            miniderek['ID']=gasparts['ParticleIDs'][lowerbound:upperbound]
            miniderek['Type']=['gas']*len(gasparts['ParticleIDs'][lowerbound:upperbound])
            miniderek['x']=gasparts['Coordinates'][:, 0][lowerbound:upperbound]
            miniderek['y']=gasparts['Coordinates'][:, 1][lowerbound:upperbound]
            miniderek['z']=gasparts['Coordinates'][:, 2][lowerbound:upperbound]
            miniderek['mass']=gasparts['Masses'][lowerbound:upperbound]
            miniderek['vx']=gasparts['Velocities'][:, 0][lowerbound:upperbound]
            miniderek['vy']=gasparts['Velocities'][:, 1][lowerbound:upperbound]
            miniderek['vz']=gasparts['Velocities'][:, 2][lowerbound:upperbound]
            
            
            print(miniderek)
            
            """
            gasdata = np.vstack([gasparts['ParticleIDs'][lowerbound:upperbound],
                                 ['gas']*len(gasparts['ParticleIDs'][lowerbound:upperbound]), 
                                 gasparts['Coordinates'][:, 0][lowerbound:upperbound],
                                 gasparts['Coordinates'][:, 1][lowerbound:upperbound],
                                 gasparts['Coordinates'][:, 2][lowerbound:upperbound], 
                                 gasparts['Masses'][lowerbound:upperbound], 
                                 gasparts['Velocities'][:, 0][lowerbound:upperbound],
                                 gasparts['Velocities'][:, 1][lowerbound:upperbound],
                                 gasparts['Velocities'][:, 2][lowerbound:upperbound]]).transpose()
            """
            derek = pd.concat([derek,miniderek])
            miniderek = pd.DataFrame(columns=pdheader)
            
            miniderek['ID']=starparts['ParticleIDs'][lowerbound:upperbound]
            miniderek['Type']=['star']*len(starparts['ParticleIDs'][lowerbound:upperbound])
            miniderek['x']=starparts['Coordinates'][:, 0][lowerbound:upperbound]
            miniderek['y']=starparts['Coordinates'][:, 1][lowerbound:upperbound]
            miniderek['z']=starparts['Coordinates'][:, 2][lowerbound:upperbound]
            miniderek['mass']=starparts['Masses'][lowerbound:upperbound]
            miniderek['vx']=starparts['Velocities'][:, 0][lowerbound:upperbound]
            miniderek['vy']=starparts['Velocities'][:, 1][lowerbound:upperbound]
            miniderek['vz']=starparts['Velocities'][:, 2][lowerbound:upperbound]#
            
            derek = pd.concat([derek,miniderek])
            miniderek = pd.DataFrame(columns=pdheader)
            """
            #fwriter.writerows(gasdata)
            stardata = np.vstack([starparts['ParticleIDs'][lowerbound:upperbound],
                                  ['stars']*len(starparts['ParticleIDs'][lowerbound:upperbound]), 
                                  starparts['Coordinates'][:, 0][lowerbound:upperbound],
                                  starparts['Coordinates'][:, 1][lowerbound:upperbound],
                                  starparts['Coordinates'][:, 2][lowerbound:upperbound], 
                                  starparts['Masses'][lowerbound:upperbound], 
                                  starparts['Velocities'][:, 0][lowerbound:upperbound],
                                  starparts['Velocities'][:, 1][lowerbound:upperbound],
                                  starparts['Velocities'][:, 2][lowerbound:upperbound]]).transpose()
            fwriter.writerows(stardata)
            """
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
            miniderek = pd.DataFrame(columns=pdheader)
            '''
            bhdata = np.vstack([bhparts['ParticleIDs'][lowerbound:upperbound],
                                ['bh']*len(bhparts['ParticleIDs'][lowerbound:upperbound]), 
                                bhparts['Coordinates'][:, 0][lowerbound:upperbound],
                                bhparts['Coordinates'][:, 1][lowerbound:upperbound],
                                bhparts['Coordinates'][:, 2][lowerbound:upperbound], 
                                bhparts['Masses'][lowerbound:upperbound], 
                                bhparts['Velocities'][:, 0][lowerbound:upperbound],
                                bhparts['Velocities'][:, 1][lowerbound:upperbound],
                                bhparts['Velocities'][:, 2][lowerbound:upperbound]]).transpose()
            fwriter.writerows(bhdata)
            '''
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
        miniderek = pd.DataFrame(columns=pdheader)
        """
        dmdata = np.vstack([dmparts['ParticleIDs'][lowerbound:upperbound],
                            ['dm']*len(dmparts['ParticleIDs'][lowerbound:upperbound]), 
                            dmparts['Coordinates'][:, 0][lowerbound:upperbound],
                            dmparts['Coordinates'][:, 1][lowerbound:upperbound],
                            dmparts['Coordinates'][:, 2][lowerbound:upperbound], 
                            dmmass[lowerbound:upperbound], 
                            dmparts['Velocities'][:, 0][lowerbound:upperbound],
                            dmparts['Velocities'][:, 1][lowerbound:upperbound],
                            dmparts['Velocities'][:, 2][lowerbound:upperbound]]).transpose()
        fwriter.writerows(dmdata)
        """
        lowerbound = upperbound
        derek.to_csv(filename+'.csv', mode='a')
        derek.to_hdf(filename+'.hdf', mode='a')
        derek = pd.DataFrame(columns=pdheader)
    print(x)
    x +=1