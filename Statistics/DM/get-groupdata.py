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
import pandas as pd
import os
import sys
sys.path.insert(1, '/home/clilje/summer_project')
from illustris_python import groupcat,lhalotree,snapshot,sublink,util
import time

"""
Reading out info from matching catalogue
"""
def get_matching(sim_size, sim_res):
    with h5py.File("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/subhalo_matching_to_dark.hdf5") as file:
        #print(file.keys())
        matchingarr = np.array(file['Snapshot_99/SubhaloIndexDark_LHaloTree'])
        #print(matchingarr)
    return(matchingarr)





size = 50
res = 1
dark = True
#z = 0
snapnum = 99


matchingarr = get_matching(50, 1)


basePathDark = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 
fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass','SubhaloLen']


subhalos = groupcat.loadSubhalos(basePathDark,99,fields=fields)
subhalosBar= groupcat.loadSubhalos(basePath,99,fields=fields)
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

#a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
#z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift


print(matchingarr)


with open('../50-1-subhalo-info-dark-og.csv', 'w', encoding='UTF8', newline='') as subfile:
    header = ['SubhaloIndex','DMIndex','SubhaloPosX','SubhaloPosY','SubhaloPosZ','SubhaloHalfmassRad','SubhaloMass','SubhaloLen']
    fwriter = csv.writer(subfile, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    g = 0
    while g < len(matchingarr):
        
        if matchingarr[g] != -1:
            print(g,matchingarr[g])
            print(subhalos['SubhaloMass'],subhalosBar['SubhaloMass'])
            data = [g, matchingarr[g], subhalos['SubhaloPos'][:, 0][matchingarr[g]],subhalos['SubhaloPos'][:, 1][matchingarr[g]],subhalos['SubhaloPos'][:, 2][matchingarr[g]], subhalos['SubhaloHalfmassRad'][matchingarr[g]], subhalos['SubhaloMass'][matchingarr[g]], subhalos['SubhaloLen'][matchingarr[g]]]
            fwriter.writerow(data)
        g +=1
            
        
    """    
    for k in range(len(matchingarr)):
        x = matchingarr[k]
        data = [int(indices[k]), subhalos['SubhaloPos'][:, 0][x],subhalos['SubhaloPos'][:, 1][x],subhalos['SubhaloPos'][:, 2][x], subhalos['SubhaloHalfmassRad'][x], subhalos['SubhaloMass'][x], subhalos['SubhaloLen'][x]]
        fwriter.writerow(data)
"""

#header = ['SubhaloIndex','SubhaloPosX','SubhaloPosY','SubhaloPosZ','SubhaloHalfmassRad','SubhaloMass','SubhaloLen']

#derek = pd.DataFrame(columns=header)
#x = 0
#while x <= len(np.array(subhalosBar['SubhaloMass'])):
#if matchingarr[x] != -1:

#data = np.vstack([num_halo, subhalos['SubhaloPos'][:, 0],subhalos['SubhaloPos'][:, 1],subhalos['SubhaloPos'][:, 2], subhalos['SubhaloHalfmassRad'], subhalos['SubhaloMass'], subhalos['SubhaloLen']]).transpose()
#fwriter.writerows(data)
