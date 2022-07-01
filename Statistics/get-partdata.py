"""import numpy as np
a = np.array([1, 2, 3, 4])

b = np.array([5, 6, 7, 8])

c = np.array([2, 3, 4, 5])

d = np.array([6, 8, 2, 4])

e = np.vstack([a, b, c, d]).transpose()
print(e)
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
z = 0
snapnum = 99

if dark == True:
    basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
else:
   basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass']
subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

with open('50-1-subhalo-info.csv', 'w', encoding='UTF8', newline='') as subfile:
    header = ['SubhaloIndex','SubhaloPosX','SubhaloPosY','SubhaloPosZ','SubhaloHalfmassRad','SubhaloMass']
    fwriter = csv.writer(subfile, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], subhalos['SubhaloMass']]).transpose()
    fwriter.writerows(data)
    

gasparts = snapshot.loadHalo(basePath, snapnum, 0, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
print(gasparts['Coordinates'])
x = 0
while x <= (len(num_halo)):
    with open('HaloParticles50-1/snap_99_halo_'+str(x)+'.csv', 'w', encoding='UTF8', newline='') as f:
        header = ['ID','Type','x','y','z','mass','vx','vy','vz']
        # Create a writer object
        fwriter = csv.writer(f, delimiter=',')
        # Write the header
        fwriter.writerow(header)
        if dark == False:
            gasparts = snapshot.loadHalo(basePath, snapnum, x, 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
            starparts = snapshot.loadHalo(basePath, snapnum, x, 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
            bhparts = snapshot.loadHalo(basePath, snapnum, x, 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
        dmparts = snapshot.loadHalo(basePath, snapnum, x, 'dm', fields=['Coordinates','ParticleIDs','Masses'])
        dmmass = [(groupcat.loadHeader(basePath, snapnum)['MassTable'][1]).astype('float')]*len(dmparts['ParticleIDs'])
        

        gasdata = np.vstack([gasparts['ParticleIDs'],['gas']*len(gasparts['ParticleIDs']), gasparts['Coordinates'][:, 0],gasparts['Coordinates'][:, 1],gasparts['Coordinates'][:, 2], gasparts['Masses'], gasparts['Velocities'][:, 0],gasparts['Velocities'][:, 1],gasparts['Velocities'][:, 2]]).transpose()
        fwriter.writerows(gasdata)
        stardata = np.vstack([starparts['ParticleIDs'],['stars']*len(starparts['ParticleIDs']), starparts['Coordinates'][:, 0],starparts['Coordinates'][:, 1],starparts['Coordinates'][:, 2], starparts['Masses'], starparts['Velocities'][:, 0],starparts['Velocities'][:, 1],starparts['Velocities'][:, 2]]).transpose()
        fwriter.writerows(stardata)
        bhdata = np.vstack([bhparts['ParticleIDs'],['bh']*len(bhparts['ParticleIDs']), bhparts['Coordinates'][:, 0],bhparts['Coordinates'][:, 1],bhparts['Coordinates'][:, 2], bhparts['Masses'], bhparts['Velocities'][:, 0],bhparts['Velocities'][:, 1],bhparts['Velocities'][:, 2]]).transpose()
        fwriter.writerows(bhdata)
        dmdata = np.vstack([dmparts['ParticleIDs'],['dm']*len(dmparts['ParticleIDs']), dmparts['Coordinates'][:, 0],dmparts['Coordinates'][:, 1],dmparts['Coordinates'][:, 2], dmparts['Masses'], dmparts['Velocities'][:, 0],dmparts['Velocities'][:, 1],dmparts['Velocities'][:, 2]]).transpose()
        fwriter.writerows(dmdata)
    print(x)
    x +=1