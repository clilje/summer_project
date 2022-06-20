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



def radial_density(dis,mass,interval):
    density = []
    rad_lowerbound = []
    #lowerbound = np.add(min_mass, mass_interval)
    lowerbound = interval
    i = 0
    while i < (len(interval)-1):
        #print(lowerbound[i])
        dr = (lowerbound[i+1]-lowerbound[i])
        dV = (4/3)*math.pi*(np.power(lowerbound[i+1],3)-np.power(lowerbound[i],3))
        #dm = lowerbound[i]
        #print(dm)
        #print(type(lowerbound[i]))
        #print(type(lowerbound[i+1]))
        #print(dis)
        #dis = [int(x, base=32) for x in dis]
        #nindex = np.where(np.logical_and((dis.astype(float)>5), (dis.astype(float)<10)))[0]
        nindex = np.where(np.logical_and(dis.astype(float)>float(lowerbound[i]), dis.astype(float)<float(lowerbound[(i+1)])))[0]
        dn = len(nindex)
        mass = mass.astype(float)
        M = np.sum(mass[nindex])
        density = np.append(density, (M/(dV)))
        #diff_number_per_mass = np.append(diff_number_per_mass, (len(np.where(np.logical_and(mass_list>=lowerbound[i], mass_list<=(lowerbound[(i+1)])))[0]))/(lowerbound[i+1]-lowerbound[i]))
        rad_lowerbound = np.append(rad_lowerbound, lowerbound[i])
        #lowerbound += mass_interval[i]
        i += 1
    return(density, rad_lowerbound)
    


interval = np.logspace(0.1, 2, 50)

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
#data_transposed = np.array(data_csv).T
#print(data_transposed)
#rad_den = radial_density(data_transposed[0][1:], data_transposed[1][1:], interval)
#print(data_csv[3])

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

plt.loglog(rad_den0[1], rad_den0[0], "+", color="black", label="Halo_0_099")
plt.loglog(rad_den1[1], rad_den1[0], "+", color="blue", label="Halo_1_099")
plt.loglog(rad_den2[1], rad_den2[0], "+", color="red", label="Halo_2_099")
plt.loglog(rad_den3[1], rad_den3[0], "+", color="green", label="Halo_3_099")
plt.loglog(rad_den4[1], rad_den4[0], "+", color="cyan", label="Halo_4_099")
plt.loglog(rad_den5[1], rad_den5[0], "+", color="pink", label="Halo_5_099")
plt.xlabel(r'Radius ($ckpc/h}$)')
plt.ylabel(r'$\rho$(r) ($10^{10} M_{\odot} h^{-1} ckpc^{-3}$)')
plt.legend()
plt.savefig('rad-den-all5')
plt.show()
