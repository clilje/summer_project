# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

@author: clara
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd

def get_filenames(sim_size, sim_res, num_files):
    filename = []
    i = 0
    # Making a list of all possible filenames
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
        i += 1
    return(filename)


def get_pos(filename):
    g = 0
    pos = np.array([])
    while g < len(filename):
        with h5py.File(str(filename[g])) as file:
            if 'Subhalo/SubhaloPos'in file:
                #print(file['Subhalo'])
                subpos = np.array(file['Subhalo/SubhaloPos'])
                pos = np.append(pos, subpos)
            g +=1
        #print(pos)
    pos = np.reshape(pos, [int(len(pos)/3),3])
    return(pos)


"""
Function that gets Subhalo Half Mass Radius from given file list.
"""
        
def get_rad(filename):
    g = 0
    rad = np.array([])
    while g < len(filename):
        with h5py.File(str(filename[g])) as file:
            if 'Subhalo/SubhaloHalfmassRad'in file:
                #print(file['Subhalo'])
                subrad = np.array(file['Subhalo/SubhaloHalfmassRad'])
                rad = np.append(rad, subrad)
            g +=1
        #pos = np.reshape(pos, [int(len(pos)/3),3])
    return(rad)

def distancefromcentre(cx, cy, cz, x, y, z, ):
    
    return (np.sqrt((np.power(np.subtract(x,cx), 2)+ np.power(np.subtract(y,cy), 2) + np.power(np.subtract(z,cz), 2)))) # distance between the centre and given point


def radial_density(partx, party, partz, mass, interval, virrad, halox, haloy, haloz):
    density = []
    rad_lowerbound = []
    lowerbound = interval
    i = 0
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    virV = (4/3)*math.pi*(np.power((virrad-10),3)-np.power((virrad+10),3))
    virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = mass.astype(float)
    virM = np.sum(mass[virindex])
    virdensity = virM/virV
    
    while i < (len(interval)-1):
        #print(lowerbound[i])
        dr = (lowerbound[i+1]-lowerbound[i])
        dV = (4/3)*math.pi*(np.power(lowerbound[i+1],3)-np.power(lowerbound[i],3))
        
        nindex = np.where(np.logical_and(dis.astype(float)>float(lowerbound[i]), dis.astype(float)<float(lowerbound[(i+1)])))[0]
        dn = len(nindex)
        #mass = mass.astype(float)
        M = np.sum(mass[nindex])
        density = np.append(density, (M/(dV))/virdensity)
        rad_lowerbound = np.append(rad_lowerbound, lowerbound[i]/virrad)
        i += 1
    return(density, rad_lowerbound)
    


interval = np.logspace(0.1, 2, 50)
files = get_filenames(50, 4, 11)
positions = get_pos(files)
radius = get_rad(files)
g = 0
numhalos = 2
rad_den = np.zeros(numhalos)
while g < numhalos:
    with open('HaloParticles/50-1_snap_99_halo_'+str(g)+'_rad_mass_100kpc.csv', 'r') as datafile:
        csvFile = csv.reader(datafile)
        data_csv = list(csvFile)
        data_transposed = np.array(data_csv).T
        print((data_transposed[0][1:].astype(float), data_transposed[1][1:].astype(float), data_transposed[2][1:].astype(float),data_transposed[3][1:].astype(float), interval, radius[g], positions[g][0], positions[g][1], positions[g][2]))
        rad_den[g] = radial_density(data_transposed[0][1:].astype(float), data_transposed[1][1:].astype(float), data_transposed[2][1:].astype(float),data_transposed[3][1:].astype(float), interval, radius[g], positions[g][0], positions[g][1], positions[g][2])
    g += 1 
    
with open('snap_99_halo_1_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den1 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
    
"""
with open('snap_99_halo_0_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den0 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
#interval = np.logspace(0.1, 1, 50)
    
with open('snap_99_halo_1_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den1 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)

with open('snap_99_halo_2_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den2 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
    
with open('snap_99_halo_3_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den3 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
    
with open('snap_99_halo_4_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den4 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
    
with open('snap_99_halo_5_rad_mass_100kpc.csv', 'r') as datafile:
    csvFile = csv.reader(datafile)  
    data_csv = list(csvFile)
    data_transposed = np.array(data_csv).T
    rad_den5 = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
"""
plt.loglog(rad_den[0][1], rad_den[0][0], "+", color="black", label="Halo_0_099")
plt.loglog(rad_den[1][1], rad_den[1][0], "+", color="blue", label="Halo_1_099")
#plt.loglog(rad_den2[1], rad_den2[0], "+", color="red", label="Halo_2_099")
#plt.loglog(rad_den3[1], rad_den3[0], "+", color="green", label="Halo_3_099")
#plt.loglog(rad_den4[1], rad_den4[0], "+", color="cyan", label="Halo_4_099")
#plt.loglog(rad_den5[1], rad_den5[0], "+", color="pink", label="Halo_5_099")
plt.xlabel(r'Radius ($ckpc/h}$)')
plt.ylabel(r'$\rho$(r) ($10^{10} M_{\odot} h^{-1} ckpc^{-3}$)')
plt.legend()
plt.savefig('rad-den-all5')
plt.show()
