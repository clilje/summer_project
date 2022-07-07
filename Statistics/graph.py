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
    
plt.plot(full_mass[indices],concentration,'.')
plt.xscale('log')
plt.xlabel(r'Total Mass of cluster in $10^{10} M_{\odot}$')
plt.ylable('$c_{200}$')
plt.savefig('cmfunc')
plt.show()