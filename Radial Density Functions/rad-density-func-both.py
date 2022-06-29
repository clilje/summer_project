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


h = 0.6774

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
    density = []
    rad_lowerbound = []
    uncertainties = []
    #lowerbound = np.logspace(0.1, (3*virrad), 50)
    #lowerbound = interval
    #i = 0
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    
    virV = (4/3)*math.pi*(np.power((virrad+10),3)-np.power((virrad-10),3))
    virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = np.array(mass)
    virM = np.sum(mass[virindex])
    virdensity = virM/virV
    
    
    bin_index = np.argsort(dis)
    radius_lowerbound = 0
    bin_lowerbound = 0
    
    while (bin_lowerbound+binsize) < len(dis):
        index_in_bin = bin_index[bin_lowerbound:(bin_lowerbound+binsize)]
        radius_upperbound = dis[index_in_bin][-1]
        dV = (4/3)*math.pi*(np.power(radius_upperbound,3)-np.power(radius_lowerbound,3))
        
        M = np.sum(mass[index_in_bin])
        subdensity = (M/(dV))
        density = np.append(density, subdensity)
        
        rad_lowerbound = np.append(rad_lowerbound, radius_upperbound)
        dn = len(index_in_bin)
        uncertainties = np.append(uncertainties, subdensity/np.sqrt(dn))
        radius_lowerbound = radius_upperbound
        bin_lowerbound = bin_lowerbound+binsize

    return(density, rad_lowerbound, uncertainties, virdensity)


matchingarr = get_matching(50, 4)

files_dark = get_filenames(50, 4, 4, True)
files = get_filenames(50, 4, 10, False)
positions = get_pos(files)
radius = get_rad(files)
positions_dark = get_pos(files_dark)
radius_dark = get_rad(files_dark)
print(positions, radius)
print(positions_dark, radius_dark)
g = 2
c = 0
numhalos = 10
densities = []
radii = []
uncertainties = []
densities_dark = []
radii_dark = []
uncertainties_dark = []
halo_number = []


while c < 4: 
    path = Path('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass_dark.csv')
    if path.is_file():
        data_csv_dark = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass_dark.csv')
        data_csv = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass.csv')
        rad_den_dark = radial_density(data_csv_dark['x'].to_numpy(), data_csv_dark['y'].to_numpy(), data_csv_dark['z'].to_numpy(), data_csv_dark['mass'].to_numpy(), 20, radius_dark[matchingarr[g]], positions_dark[matchingarr[g]][0], positions_dark[matchingarr[g]][1], positions_dark[matchingarr[g]][2])
        rad_den = radial_density(data_csv['x'].to_numpy(), data_csv['y'].to_numpy(), data_csv['z'].to_numpy(), data_csv['mass'].to_numpy(), 20, radius[g], positions[g][0], positions[g][1], positions[g][2])
        #print(rad_den)
        hmrad = radius[g]
        densities.append(list(rad_den[0]))
        radii.append(list(rad_den[1]))
        uncertainties.append(list(rad_den[2]))
        hmden = rad_den[3]
        
        hmrad_dark = radius_dark[matchingarr[g]]
        densities_dark.append(list(rad_den_dark[0]))
        radii_dark.append(list(rad_den_dark[1]))
        uncertainties_dark.append(list(rad_den_dark[2]))
        hmden_dark = rad_den_dark[3]
        print(hmrad,hmden)
        halo_number.append(g)
        c+=1
    g += 1
    

hsv = plt.get_cmap('hsv')
colors = iter(hsv(np.linspace(0,1,6)))
fig, axs = plt.subplots(2, 2, figsize=(15,15))
radii = np.array(radii, dtype=('float64'))
densities = np.array(densities, dtype=('float64'))
uncertainties = np.array(uncertainties, dtype=('float64'))
radii_dark = np.array(radii_dark, dtype=('float64'))
densities_dark = np.array(densities_dark, dtype=('float64'))
uncertainties_dark = np.array(uncertainties_dark, dtype=('float64'))

axs[0,0].errorbar((radii[0])*(h/hmrad), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(halo_number[0])+"_099", color='black')
axs[0,0].errorbar((radii_dark[0])*(h/hmrad_dark), (densities_dark[0])/(10*(h**2)*hmden_dark), yerr=(uncertainties_dark[0]), fmt='.', label="Halo_"+str(halo_number[0])+"_099_dark", color='red')
axs[0,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,0].legend()
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].set_title("Halo"+str(halo_number[0]))

axs[0,1].errorbar((radii[1])*(h/hmrad), (densities[1])/(10*(h**2)*hmden), yerr=(uncertainties[1]), fmt='.', label="Halo_"+str(halo_number[1])+"_099", color='black')
axs[0,1].errorbar((radii_dark[1])*(h/hmrad_dark), (densities_dark[1])/(10*(h**2)*hmden_dark), yerr=(uncertainties_dark[1]), fmt='.', label="Halo_"+str(halo_number[1])+"_099_dark", color='red')
axs[0,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,1].legend()
axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_title("Halo"+str(halo_number[1]))

axs[1,0].errorbar((radii[2])*(h/hmrad), (densities[2])/(10*(h**2)*hmden), yerr=(uncertainties[2]), fmt='.', label="Halo_"+str(halo_number[2])+"_099", color='black')
axs[1,0].errorbar((radii_dark[2])*(h/hmrad_dark), (densities_dark[2])/(10*(h**2)*hmden_dark), yerr=(uncertainties_dark[2]), fmt='.', label="Halo_"+str(halo_number[2])+"_099_dark", color='red')
axs[1,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[1,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,0].legend()
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title("Halo"+str(halo_number[2]))

axs[1,1].errorbar((radii[3])*(h/hmrad), (densities[3])/(10*(h**2)*hmden), yerr=(uncertainties[3]), fmt='.', label="Halo_"+str(halo_number[3])+"_099", color='black')
axs[1,1].errorbar((radii_dark[3])*(h/hmrad_dark), (densities_dark[3])/(10*(h**2)*hmden_dark), yerr=(uncertainties_dark[3]), fmt='.', label="Halo_"+str(halo_number[3])+"_099_dark", color='red')
axs[1,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[1,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,1].legend()
axs[1,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].set_title("Halo"+str(halo_number[3]))


fig.tight_layout()

fig.savefig('rad-den-dark-comparison')
fig.show()

