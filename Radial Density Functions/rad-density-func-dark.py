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
from pathlib import Path

def get_filenames(sim_size, sim_res, num_files):
    """
    

    Parameters
    ----------
    sim_size : Size of Edge of simulation box in kpc.
    sim_res : Resolution of simulation, integer from 1 to 4
    num_files : number of files in given simulation folder

    Returns
    -------
    List of filnames in string formats

    """
    filename = []
    i = 0
    # Making a list of all possible filenames
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"-dark/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
        i += 1
    return(filename)


def get_pos(filename):
    """
    

    Parameters
    ----------
    filename : Array of group file names in string format from which to extract the subhalo position.

    Returns
    -------
    Array of arrays containing x,y,z position of subhalo centres.

    """
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

def get_rad(filename):
    """
    

    Parameters
    ----------
    filename : The names of the group files in which to look for halfmass radii.

    Returns
    -------
    Array containing the halfmass radius for each Subhalo in the order they are saved in the files by TNG.

    """
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

"""
Reading out info from matching catalogue
"""
def get_matching(sim_size, sim_res):
    with h5py.File("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/subhalo_matching_to_dark.hdf5") as file:
        #print(file.keys())
        matchingarr = np.array(file['Snapshot_99/SubhaloIndexDark_LHaloTree'])
        #print(matchingarr)
    return(matchingarr)

def distancefromcentre(cx, cy, cz, x, y, z, ):
    """
    

    Parameters
    ----------
    cx : Halo centre x-coord float
    cy : Halo centre y-coord float
    cz : Halo centre z-coord float
    x : Particle Position x-coord
    y : Particle Position y-coord
    z : Particle Position z-coord

    Returns
    -------
    Distance between particle and halo centre

    """
    return (np.sqrt((np.power(np.subtract(x,cx), 2)+ np.power(np.subtract(y,cy), 2) + np.power(np.subtract(z,cz), 2)))) # distance between the centre and given point


def radial_density(partx, party, partz, mass, interval, virrad, halox, haloy, haloz):
    """
    

    Parameters
    ----------
    partx : Particles in halo x coordinates (array)
    party : Particles in halo y coordinates (array)
    partz : Particles in halo z coordinates (array)
    mass : Particle mass (array)
    interval : Interval of radii over which mass is binned and density calculated
    virrad : Half Mass Radius of halo
    halox : halo centre x coordinate
    haloy : halo centre y coordinate
    haloz : halo centre z coordinate

    Returns
    -------
    list with densities at certain radii and radii at which density is calculated.

    """
    density = []
    rad_lowerbound = []
    #lowerbound = np.logspace(0.1, (3*virrad), 50)
    lowerbound = interval
    i = 0
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    virV = (4/3)*math.pi*(np.power((virrad+10),3)-np.power((virrad-10),3))
    virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = mass.astype(float)
    virM = np.sum(mass[virindex])
    virdensity = virM/virV
    print(virindex)
    print(virV)
    print(virM)
    print(virdensity)
    
    while i < (len(lowerbound)-1):
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

matchingarr = get_matching(50, 4)

interval = np.logspace(0.1, 2.5, 100)
files = get_filenames(50, 4, 4)
positions = get_pos(files)
radius = get_rad(files)
g = 0
numhalos = 5
densities = []
radii = []
halo_number = []

g = 0 
#while g < len(halfmassradii)
while g < (numhalos): 
    path = Path('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass_dark.csv')
    if path.is_file():
        data_csv = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass_dark.csv')
        rad_den = radial_density(data_csv['x'], data_csv['y'], data_csv['z'],data_csv['mass'], interval, radius[matchingarr[g]], positions[matchingarr[g]][0], positions[matchingarr[g]][1], positions[matchingarr[g]][2])
        print(rad_den)
        densities.append(list(rad_den[0]))
        radii.append(list(rad_den[1]))
        halo_number.append(g)
    g += 1
    

hsv = plt.get_cmap('gnuplot')
colors = iter(hsv(np.linspace(0,1,15)))
b = 0
while b < (len(radii)):
    print('loop')
    plt.loglog(radii[b], densities[b], "+", label="Halo_dark"+str(halo_number[b])+"_099", color=next(colors))
    b += 1

plt.xlabel(r'Radius ($ckpc/(h*R_{HalfMass}})}$)')
plt.ylabel(r'$\rho$(r) ($10^{10} M_{\odot} h^{-1} ckpc^{-3} (\rho_{HalfMass})^{-1}$)')
plt.legend()
plt.savefig('rad-den-dark-next5')
plt.show()

