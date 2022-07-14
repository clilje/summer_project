# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

This code is designed to find the difference between the Halfmass Radius
of Gravitationally bound particles and Particles that just happend to be 
within a certain radius of the Halo centre. 

@author: clara
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
from pathlib import Path


h = 0.6774


        
def get_mass(filename):
    """
    Function that gets Subhalo Mass from given file list.
    """
    g = 0
    mass = np.array([])
    while g < len(filename):
        with h5py.File(str(filename[g])) as file:
            if 'Subhalo/SubhaloMass'in file:
                #print(file['Subhalo'])
                submass = np.array(file['Subhalo/SubhaloMass'])
                mass = np.append(mass, submass)
            g +=1
    return(mass)



def get_filenames(sim_size, sim_res, num_files, dark):
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
    if dark == True:
        while i < num_files:
            filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"-dark/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
            i += 1
    else:
        while i < num_files:
            filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
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
                subrad = np.array(file['Subhalo/SubhaloHalfmassRad'])
                rad = np.append(rad, subrad)
            g +=1
        #pos = np.reshape(pos, [int(len(pos)/3),3])
    return(rad)


def get_matching(sim_size, sim_res):
    """
    Reading out info from matching catalogue
    """
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

def radial_density(partx, party, partz, mass, binsize, virrad, halox, haloy, haloz):
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
    
    #create lists in which to store final results
    density = np.array([])
    rad_lowerbound = np.array([])
    uncertainties = np.array([])
    
    
    #calculate radial distance of each particle
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    #index= np.argsort(dis)
    
    #calculate the halfmass radius density, here falsely named virial radius
    virV = (4/3)*math.pi*(np.power((virrad+10),3)-np.power((virrad-10),3))
    virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = np.array(mass)
    virM = np.sum(mass[virindex])
    virdensity = virM/virV
    
    #sort list of radial distances to make binning easier
    bin_index = np.argsort(dis)
    radius_lowerbound = 0
    bin_lowerbound = 0
    
    
    while (bin_lowerbound+binsize) < len(dis):
        """
        This loop iterates through radial shells, each containing the same amount of halos.
        Then determines the density in each bin, storing everything in lists.
        Uncertainties are Poisson errors.
        """
        #radial shell volume
        index_in_bin = bin_index[bin_lowerbound:(bin_lowerbound+binsize)]
        radius_upperbound = dis[index_in_bin][-1]
        dV = (4/3)*math.pi*(np.power(radius_upperbound,3)-np.power(radius_lowerbound,3))
        
        #mass contained in shell
        M = np.sum(mass[index_in_bin])
        
        #density = mass/volume
        subdensity = (M/(dV))
        density = np.append(density, subdensity)
        
        #change indices for next iteration
        rad_lowerbound = np.append(rad_lowerbound, radius_upperbound)
        dn = len(index_in_bin)
        uncertainties = np.append(uncertainties, subdensity/np.sqrt(dn))
        radius_lowerbound = radius_upperbound
        bin_lowerbound = bin_lowerbound+binsize

    return(density, rad_lowerbound, uncertainties, virdensity)


def halfmass(partx, party, partz, mass, binsize, hmrad, halox, haloy, haloz):
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
    mass within the given halfmass radius

    """
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    hmindex = np.where(dis.astype(float)<float(hmrad))[0]
    mass = np.array(mass)
    hmM = np.sum(mass[hmindex])
    
    return(hmM)

#Set up defining simulation
files = get_filenames(50, 4, 9, False)
positions = get_pos(files)
radius = get_rad(files)
masses = get_mass(files)

#Which halos to read out
g = 0
c = 0
numhalos = 30

#lists to store results
halo_number = []
hmclara = []


while c < numhalos: 
    #get data from file containing all particles within specific radius of Halo Centre
    data_csv = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass.csv')
    hmrad = radius[g]
    
    #get half mass 
    hmMass = halfmass(data_csv['x'].to_numpy(), data_csv['y'].to_numpy(), data_csv['z'].to_numpy(), data_csv['mass'].to_numpy(), 10, radius[g], positions[g][0], positions[g][1], positions[g][2])

    halo_number.append(g)
    hmclara.append(hmMass)
    c+=1
    g += 1


hmclara = np.array(hmclara)
mclara = 2*hmclara

#Calculate difference between given full mass in Catalogue and the one found by 
#doubling the mass in given HMrad for all particles within that radius of Halo Centre
mass_difference = masses[:len(hmclara)]-mclara

with open('50-4_mass_difference.csv', 'w', encoding='UTF8', newline='') as f:
    '''
    Write results to file
    '''
    header = ['halo_number','halfmass_clara','mass_clara','mass_tng','mass_difference']
    # Create a writer object
    fwriter = csv.writer(f, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    x = 0
    while x <numhalos:
        data = [halo_number[x],hmclara[x],mclara[x],masses[x],mass_difference[x]]
        fwriter.writerow(data)
        x +=1

print('done')