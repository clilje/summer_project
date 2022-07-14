# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

This code plots the radial density function, from a file which includes all particles,
their position and mass.
This code has the additional functionality of plotting both the radial density function
for the DMO and the DM + Baryons halo.
Lastly the code calculates the enhancement of each halo and adds this to the plot.

This code contains functions for getting filenames for group catalogs, extracting position
and radius from h5py files, as well as getting the matching indices of DMO halos.
Lastly there are functions for radial density binning, distance from centre of particles and 
calculating enhancement that are easily adaptable. 


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

def enhancement(halopos, halopartpos, halopartmass, cutoffrad, halodmpos, halodmpartpos, halodmpartmass):
    """
    

    Parameters
    ----------
    halopos : contains x,y,z of halo
    halopartpos : contains lists of x,y,z of particles belonging to halo
    halopartmass : list of mass of each particle
    cutoffrad : radius within which enhancement should be calculated.
    halodmpos : DMO halo centre x,y,z
    halodmpartpos : particle positions belonging to DMO halo
    halodmpartmass : List of mass of each particle of DMO halo

    Returns
    -------
    enhancement of DM + Baryonic compared to DMO

    """
    #simulation constants
    omega_m = 0.3089
    omega_b = 0.0486
    
    #extracting data from lists
    halox, haloy, haloz = halopos
    halodmx, halodmy, halodmz = halodmpos
    halopartx, haloparty, halopartz = halopartpos
    halodmpartx, halodmparty, halodmpartz = halodmpartpos
    
    #calculating radial distances
    dis = distancefromcentre(halox, haloy, haloz, halopartx, haloparty, halopartz)
    disdm = distancefromcentre(halodmx, halodmy, halodmz, halodmpartx, halodmparty, halodmpartz)
    
    #cutoff indices
    virindex = np.where(dis.astype(float)<float(cutoffrad))[0]
    virindex_dm = np.where(disdm.astype(float)<float(cutoffrad))[0]
    M = np.sum(halopartmass[virindex])
    M_dm = np.sum(halodmpartmass[virindex_dm])
    
    return((M/M_dm)*(omega_m/(omega_m-omega_b)))
    
#get filenames and group file information
matchingarr = get_matching(50, 4)
files_dark = get_filenames(50, 4, 4, True)
files = get_filenames(50, 4, 9, False)
positions = get_pos(files)
radius = get_rad(files)
positions_dark = get_pos(files_dark)
radius_dark = get_rad(files_dark)

#constants defining which halos are selected and iterated through
g = 1
c = 0
numhalos = 10


#lists to store results
densities = []
radii = []
uncertainties = []
densities_dark = []
radii_dark = []
uncertainties_dark = []
halo_number = []
enhance = []

while c < 4: 
    
    #iterate through each dark halo, because there are more DM+Baryonic matched
    path = Path('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass_dark.csv')
    
    
    if path.is_file():
        #extract data
        data_csv_dark = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass_dark.csv')
        data_csv = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass.csv')
        
        #calc radial densities
        rad_den_dark = radial_density(data_csv_dark['x'].to_numpy(), data_csv_dark['y'].to_numpy(), data_csv_dark['z'].to_numpy(), data_csv_dark['mass'].to_numpy(), 10, radius_dark[matchingarr[g]], positions_dark[matchingarr[g]][0], positions_dark[matchingarr[g]][1], positions_dark[matchingarr[g]][2])
        rad_den = radial_density(data_csv['x'].to_numpy(), data_csv['y'].to_numpy(), data_csv['z'].to_numpy(), data_csv['mass'].to_numpy(), 10, radius[g], positions[g][0], positions[g][1], positions[g][2])
        
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
        
        halo_number.append(g)
        
        #calc enhancement
        enhance.append(enhancement([positions[g][0], positions[g][1], positions[g][2]],[data_csv['x'].to_numpy(), data_csv['y'].to_numpy(), data_csv['z'].to_numpy()],data_csv['mass'].to_numpy(),10,[positions_dark[matchingarr[g]][0], positions_dark[matchingarr[g]][1], positions_dark[matchingarr[g]][2]],[data_csv_dark['x'].to_numpy(), data_csv_dark['y'].to_numpy(), data_csv_dark['z'].to_numpy()], data_csv_dark['mass'].to_numpy()))
        
        c+=1
    g += 1
    
#print(enhance)

#This section plots 4 of the halos both DMO and counterparts on 4 subplots, writing the enhancement on each.
hsv = plt.get_cmap('hsv')
colors = iter(hsv(np.linspace(0,1,6)))
fig, axs = plt.subplots(2, 2, figsize=(16,16))

axs[0,0].errorbar(np.array(radii[0])*(h/hmrad), np.array(densities[0])/(10*(h**2)*hmden), yerr=np.array(uncertainties[0]), fmt='.', label="Halo_"+str(halo_number[0])+"_099", color='black')
axs[0,0].errorbar(np.array(radii_dark[0])*(h/hmrad_dark), np.array(densities_dark[0])/(10*(h**2)*hmden_dark), yerr=np.array(uncertainties_dark[0]), fmt='.', label="Halo_"+str(halo_number[0])+"_099_dark", color='red')
axs[0,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,0].text(0.1,0.2,'$\eta$='+str(enhance[0]), transform=axs[0,0].transAxes)
axs[0,0].legend()
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].set_title("Halo"+str(halo_number[0]))

axs[0,1].errorbar(np.array(radii[1])*(h/hmrad), np.array(densities[1])/(10*(h**2)*hmden), yerr=np.array(uncertainties[1]), fmt='.', label="Halo_"+str(halo_number[1])+"_099", color='black')
axs[0,1].errorbar(np.array(radii_dark[1])*(h/hmrad_dark), np.array(densities_dark[1])/(10*(h**2)*hmden_dark), yerr=np.array(uncertainties_dark[1]), fmt='.', label="Halo_"+str(halo_number[1])+"_099_dark", color='red')
axs[0,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,1].text(0.1,0.2,'$\eta$='+str(enhance[1]), transform=axs[0,1].transAxes)
axs[0,1].legend()
axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_title("Halo"+str(halo_number[1]))

axs[1,0].errorbar(np.array(radii[2])*(h/hmrad), np.array(densities[2])/(10*(h**2)*hmden), yerr=np.array(uncertainties[2]), fmt='.', label="Halo_"+str(halo_number[2])+"_099", color='black')
axs[1,0].errorbar(np.array(radii_dark[2])*(h/hmrad_dark), np.array(densities_dark[2])/(10*(h**2)*hmden_dark), yerr=np.array(uncertainties_dark[2]), fmt='.', label="Halo_"+str(halo_number[2])+"_099_dark", color='red')
axs[1,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[1,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,0].text(0.1,0.2,'$\eta$='+str(enhance[2]), transform=axs[1,0].transAxes)
axs[1,0].legend()
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title("Halo"+str(halo_number[2]))

axs[1,1].errorbar(np.array(radii[3])*(h/hmrad), np.array(densities[3])/(10*(h**2)*hmden), yerr=np.array(uncertainties[3]), fmt='.', label="Halo_"+str(halo_number[3])+"_099", color='black')
axs[1,1].errorbar(np.array(radii_dark[3])*(h/hmrad_dark), np.array(densities_dark[3])/(10*(h**2)*hmden_dark), yerr=np.array(uncertainties_dark[3]), fmt='.', label="Halo_"+str(halo_number[3])+"_099_dark", color='red')
axs[1,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[1,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,1].text(0.1,0.2,'$\eta$='+str(enhance[3]), transform=axs[1,1].transAxes)
axs[1,1].legend()
axs[1,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].set_title("Halo"+str(halo_number[3]))


#fig.tight_layout()

fig.savefig('rad-den-enhancement')
fig.show()

