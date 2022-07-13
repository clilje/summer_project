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
import os
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)


def nfw(r, density_0, scale_radius):
    return(density_0/((r/scale_radius)*np.power((1+(r/scale_radius)),2)))


def virialRadius(radius, density):
    above_virR = np.where((density).astype(float)>float(p_crit*200))[0]
    virIndex = np.argmax(radius[above_virR])
    virR = radius[virIndex]
    #print(virR)
    #print(radius[-10:-1])
    return(virR,above_virR)

subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy().astype(int)

"""
fit_param = pd.read_csv('HaloFitsInfo/50-4_snap_99_fit_param.csv')
nfw_chisquare = fit_param['NFW ChiSquare'].to_numpy()
nfw_scalerad = fit_param['NFW Scale Radius'].to_numpy()
datapoint = fit_param['DataPoints'].to_numpy()
indices = fit_param['Halo Number'].to_numpy().astype(int)
"""

numhalos = len(subhalo_index)
#weighted_chisquare = nfw_chisquare
#g = 51
concentration = []
chisquare = []
g= 0
while g < 360999:
    filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den'
    if(os.path.isfile(filename+'.csv')):
        print(g)
        
        data_csv = pd.read_csv('HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den.csv')
        
        rad = data_csv['Radius']
        den = data_csv['Density']
        uncer = data_csv['Uncertainty']
        virrad,virial_index = virialRadius(rad, den)
        virial_density = p_crit*200
        if virrad/radius[g] < 15:
            nfwfitp, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[virial_density,virrad], sigma=uncer)
            nfwchi_square_test_statistic =  np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/(nfw(rad, nfwfitp[0], nfwfitp[1])))
            chisquare.append(nfwchi_square_test_statistic)
            #virrad = data_csv['Virial Radius'][0]
            #uncer = data_csv['Uncertainty']
            concentration.append(virrad/nfwfitp[0])
            #num_datapoints = len(data_csv['Radius'])
    g +=1


#lower_end = np.where(weighted_chisquare<np.mean(weighted_chisquare))
#weighted_chisquare = nfw_chisquare/weight

#lower_end = np.where(weighted_chisquare<(10**3))[0]
#print(lower_end)
#print(weighted_chisquare)
plt.plot((full_mass)*h,concentration,'.')
plt.xscale('log')
plt.xlabel(r'Total Mass of Halo in $10^{10} M_{\odot}$')
plt.ylabel(r'$c_{200}$')
plt.savefig('cmfunc-reduced')
plt.show()

fig = plt.figure()
ax = plt.axes(projection ='3d')
#xyz = np.arange(len(positionsX[indices]))
#index = np.random.choice(xyz,2000)

chisquare = np.array(chisquare)/min(chisquare)
p = ax.scatter3D(positionsX, positionsY, positionsZ,c=chisquare, marker='.',cmap='hot',alpha=0.3)
#ax.scatter(halo_50[0], halo_50[1], halo_50[2], marker='+',color='red')
#ax.scatter(positionsX[index_sub],positionsY[index_sub],positionsZ[index_sub],marker='x', color='black')
#ax.scatter(CM[0], CM[1], CM[2], marker='+',color='pink')
fig.colorbar(p, ax=ax)
ax.set_xlabel('x [ckpc/h]')

ax.set_ylabel('y [ckpc/h]')
ax.set_zlabel('z [ckpc/h]')
fig.savefig('Heatmap-reduced')