# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 22:48:58 2022

@author: clara
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import h5py
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)


subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy().astype(int)


fit_param = pd.read_csv('HaloFitsInfo/50-4_snap_99_fit_param.csv')
nfw_chisquare = fit_param['NFW ChiSquare'].to_numpy()
nfw_scalerad = fit_param['NFW Scale Radius'].to_numpy()
datapoint = fit_param['DataPoints'].to_numpy()
indices = fit_param['Halo Number'].to_numpy().astype(int)
numhalos = len(subhalo_index)
weighted_chisquare = nfw_chisquare/datapoint
#g = 51
concentration = []
for g in indices:
    data_csv = pd.read_csv('HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den.csv')
    
    virrad = data_csv['Virial Radius'][0]
    #den = data_csv['Density']
    #uncer = data_csv['Uncertainty']
    concentration.append(virrad/nfw_scalerad[g-indices[0]])
    #num_datapoints = len(data_csv['Radius'])
    g +=1
    
plt.plot((full_mass[indices]*h),concentration,'.')
plt.xscale('log')
plt.xlabel(r'Total Mass of Halo in $10^{10} M_{\odot}$')
plt.ylabel(r'$c_{200}$')
plt.savefig('cmfunc')
plt.show()

fig = plt.figure()
ax = plt.axes(projection ='3d')
xyz = np.arange(len(positionsX))
index = np.random.choice(xyz,2000)
ax.scatter(positionsX[index], positionsY[index], positionsZ[index],c=weighted_chisquare[index], marker='.',cmap='hot',alpha=0.1)
#ax.scatter(halo_50[0], halo_50[1], halo_50[2], marker='+',color='red')
#ax.scatter(positionsX[index_sub],positionsY[index_sub],positionsZ[index_sub],marker='x', color='black')
#ax.scatter(CM[0], CM[1], CM[2], marker='+',color='pink')

ax.set_xlabel('x [ckpc/h]')

ax.set_ylabel('y [ckpc/h]')
ax.set_zlabel('z [ckpc/h]')
fig.savefig('Heatmap')