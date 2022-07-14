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
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass']


subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift
x = 0
while x <= (len(num_halo)):
    if dark == False: 
        filename = 'HaloParticles50-1-sub/snap_99_halo_'+str(x)+'.csv'
    else: 
        filename = 'HaloParticles50-1-sub/snap_99_halo_'+str(x)+'-dark.csv'
    with open(filename, 'w', encoding='UTF8', newline='') as f:
        print(x)
        header = ['ID','Type','x','y','z','mass','vx','vy','vz']
        # Create a writer object
        fwriter = csv.writer(f, delimiter=',')
        # Write the header
        fwriter.writerow(header)
        if dark == False:
            gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
            starparts = snapshot.loadHalo(basePath, snapnum, x, 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
            bhparts = snapshot.loadHalo(basePath, snapnum, x, 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
        dmparts = snapshot.loadHalo(basePath, snapnum, x, 'dm', fields=['Coordinates','ParticleIDs','Velocities'])
        
        with h5py.File(snapshot.snapPath(basePath,snapnum),'r') as f:
            header = dict( f['Header'].attrs.items() )
            #print(header['MassTable'][1]) # 10^10 msun/h
            dmmass = [(header['MassTable'][1]).astype('float')]*len(dmparts['ParticleIDs'])
        
        lowerbound = 0
        indexjump = 999999
        while lowerbound < len(gasparts['ParticleIDs']):
            print(lowerbound)
            if (lowerbound+indexjump)>len(gasparts['ParticleIDs']):
                upperbound = len(gasparts['ParticleIDs'])
            else:
                upperbound = lowerbound+indexjump
            if dark == False:
                gasdata = np.vstack([gasparts['ParticleIDs'][lowerbound:upperbound],
                                     ['gas']*len(gasparts['ParticleIDs'][lowerbound:upperbound]), 
                                     gasparts['Coordinates'][:, 0][lowerbound:upperbound],
                                     gasparts['Coordinates'][:, 1][lowerbound:upperbound],
                                     gasparts['Coordinates'][:, 2][lowerbound:upperbound], 
                                     gasparts['Masses'][lowerbound:upperbound], 
                                     gasparts['Velocities'][:, 0][lowerbound:upperbound],
                                     gasparts['Velocities'][:, 1][lowerbound:upperbound],
                                     gasparts['Velocities'][:, 2][lowerbound:upperbound]]).transpose()
                fwriter.writerows(gasdata)
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
            
            lowerbound = upperbound
    print(x)
    x +=1