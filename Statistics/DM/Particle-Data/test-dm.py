# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:41:28 2022

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




#Set up and definition of which simulation is used.

size = 50
res = 1
dark = True
#z = 0
snapnum = 99


matchingarr = get_matching(50, 1)

#specify the base path to the files

darkbasePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 


fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass']


#lead the group catalogues data

subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)

subhalosdark = groupcat.loadSubhalos(darkbasePath,99,fields=fields)

#get the indices of each halo, will have to be adapted for DMO 
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift


x = [0,15,23,200,45200,100000,1200000]
while x <= 10:
    if matchingarr[x] != -1:
        print(subhalos['SubhaloPos'][x])
        print(subhalosdark['SubhaloPos'][matchingarr][x])
    else:
        print('no match')