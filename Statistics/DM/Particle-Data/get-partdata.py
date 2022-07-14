"""
Created on 

This is the code used to extract particle data for each halo and write it to file.


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
if dark == True:
    basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res)+'-dark'
else:
   basePath = '/disk01/rmcg/downloaded/tng/tng'+str(size)+'-'+str(res) 


fields = ['SubhaloPos','SubhaloHalfmassRad','SubhaloMass']


#lead the group catalogues data

subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)

#get the indices of each halo, will have to be adapted for DMO 
num_halo = np.arange(len(np.array(subhalos['SubhaloMass'])))

a = (groupcat.loadHeader(basePath, snapnum)['Time']).astype('int')  #scalefactor
z = (groupcat.loadHeader(basePath, snapnum)['Redshift']).astype('float')  #redshift


#starting halo
x = 0
c = []
pdheader = ['ID','Type','x','y','z','mass','vx','vy','vz']


while x <= len(np.array(subhalos['SubhaloMass'])):
    if matchingarr[x] != -1:
        #create the filename to be written to 
        if dark == False: 
            filename = 'FullRun/snap_99_halo_'+str(x)
        else: 
            filename = '../../FullRunDark/snap_99_halo_'+str(x)+'-dark'
            
            
        '''
        #Optional to check whether faulty file is present if nessecary.
        if(os.path.isfile(filename+'.csv')):
        
            #os.remove() function to remove the file
            os.remove(filename+'.csv')
        
            #Printing the confirmation message of deletion
            print("File Deleted successfully")
        else:
            print("File does not exist")
        '''
        
        
        #load all the particle data for each subhalo
        if dark == False:
            gasparts = snapshot.loadSubhalo(basePath, snapnum, matchingarr[x], 'gas', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
            starparts = snapshot.loadSubhalo(basePath, snapnum, matchingarr[x], 'stars', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
            bhparts = snapshot.loadSubhalo(basePath, snapnum, matchingarr[x], 'bh', fields=['Coordinates','ParticleIDs','Velocities','Masses'])
        dmparts = snapshot.loadSubhalo(basePath, snapnum, matchingarr[x], 'dm', fields=['Coordinates','ParticleIDs','Velocities'])
        
        
        print(dmparts['count'])
        
        #get the mass of the DM particles
        with h5py.File(snapshot.snapPath(basePath,snapnum),'r') as f:
            header = dict( f['Header'].attrs.items() )
            #print(header['MassTable'][1]) # 10^10 msun/h
            dmmass = [(header['MassTable'][1]).astype('float')]*dmparts['count']
        
        #getting a chunk to write to file, since some files are too big for cuillins little memory
        lowerbound = 0
        indexjump = 900000
        
        if dark == False: 
            if gasparts['count'] >= starparts['count'] and gasparts['count'] >= bhparts['count'] and gasparts['count'] >= dmparts['count']:
                max_num_part = gasparts['count']
                print("gas dominates")
            if dmparts['count'] >= gasparts['count'] and dmparts['count'] >= bhparts['count'] and dmparts['count'] >= starparts['count']:
                max_num_part = dmparts['count']
                print("dm dominates")
                c.append(x)
            if starparts['count'] >= gasparts['count'] and starparts['count'] >= bhparts['count'] and starparts['count'] >= dmparts['count']:
                max_num_part = starparts['count']
                print("star dominates")
                c.append(x)
            if bhparts['count'] >= gasparts['count'] and bhparts['count'] >= starparts['count'] and bhparts['count'] >= dmparts['count']:
                max_num_part = bhparts['count']
                print("bh dominates")
                c.append(x)
        else:
            print("dm dominates")
            max_num_part = dmparts['count']
    
     
            
            
            
        #create a data frame that will be finally written to file
        derek = pd.DataFrame(columns=pdheader)
        
        while lowerbound < max_num_part:
            
            #create sub-dataframe for each chunk
            miniderek = pd.DataFrame(columns=pdheader)
            
            #define new bounds for chunk
            if (lowerbound+indexjump)>max_num_part:
                upperbound = max_num_part
            else:
                upperbound = lowerbound+indexjump
            
            #read out particle properties and add to data frame
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
                    
                    derek = pd.DataFrame(np.concatenate([derek.values,miniderek.values]),columns=derek.columns)
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
                    
                    derek = pd.DataFrame(np.concatenate([derek.values,miniderek.values]),columns=derek.columns)
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
                    derek = pd.DataFrame(np.concatenate([derek.values,miniderek.values]),columns=derek.columns)
                    miniderek = miniderek[0:0]
    
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
                derek = pd.DataFrame(np.concatenate([derek.values,miniderek.values]),columns=derek.columns)
                print(miniderek)
                miniderek = miniderek[0:0]
            
            #define new bounds
            lowerbound = upperbound
            #write data frame to csv file
            derek.to_csv(filename+'.csv', mode='a')
            derek = derek[0:0]
            
        #check which halo we are on 
        print(x)
    x +=1
print(c)