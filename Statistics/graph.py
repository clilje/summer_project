# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 22:48:58 2022
This code is used to plot the concentration mass function of halos.
Could also be used to plot 3d distribution of halos with their associated chisquare values.

@author: clara
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import csv

import pandas as pd
import os
import statistics
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)


def nfw(r, density_0, scale_radius):
    """
    NFW profile formula
    """
    return(density_0/((r/scale_radius)*np.power((1+(r/scale_radius)),2)))


def virialRadius(radius, density):
    """
    

    Parameters
    ----------
    radius : list of radial distances
    density : list of densities

    Returns
    -------
    radius at which halo is 200*critical density of universe.

    """
    above_virR = np.where((density).astype(float)>float(p_crit*200))[0]
    virIndex = np.argmax(radius[above_virR])
    virR = radius[virIndex]
    return(virR,above_virR)

#Get key info from group catalogue from file
subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex'].to_numpy().astype(int)
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
#length = subhalo_info['SubhaloLen'].to_numpy().astype(int)


#Read in the optimal fit parameters as well as chisquare
fit_param = pd.read_csv('HaloFitsInfo/50-1_snap_99_fit_param.csv')
nfw_chisquare = fit_param['NFW ChiSquare'].to_numpy()
nfw_scalerad = fit_param['NFW Scale Radius'].to_numpy()
datapoint = fit_param['DataPoints'].to_numpy()
virrad = fit_param['Virial Radius'].to_numpy()
true_indices = fit_param['Halo Number'].to_numpy().astype(int)

#lists to store data
numhalos = len(subhalo_index)


#calculate concentration from given arrays
concentration = virrad/nfw_scalerad

#Get list of indices plotted
indices = np.linspace(0,len(true_indices))


#Prepare mass to be binned
bins = np.logspace(-2,5.5,12)

#get lists for binning
mean_mass = []
mean_concentration =[]
stdev = []
lowerbound = 0

#loop over bins
for upperbound in bins:
    #get indices of mass lying inside bin
    massindex = np.where(np.logical_and(full_mass[true_indices]<upperbound,full_mass[true_indices]>lowerbound))[0]
    
    #get indices for which concentration values correspond to this
    conc_index = np.searchsorted(true_indices, subhalo_index[true_indices][massindex])
    
    #ensure no error for stdev or mean
    if conc_index.size ==0:
        break
    
    #append all data to lists
    mean_mass.append(((upperbound-lowerbound)/2)*h)
    mean_concentration.append(statistics.mean(concentration[conc_index]))
    stdev.append(statistics.stdev(concentration[conc_index]))
    lowerbound= upperbound


#plot data obtained in scatterplot
plt.errorbar(np.array(mean_mass),mean_concentration,yerr=stdev,fmt='.')
plt.xscale('log')
plt.yscale('log')
plt.ylim((1,30))
plt.xlabel(r'Total Mass of Halo in $10^{10} M_{\odot}$')
plt.ylabel(r'$c_{200}$')
plt.savefig('cmfunc-6')
plt.show()


"""
#This is optional to plot heat map of chisquare over 3d space

fig = plt.figure()
ax = plt.axes(projection ='3d')
xyz = np.arange(len(positionsX[indices]))
index = np.random.choice(xyz,2000)
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
"""