"""
Created on Thu Jun 16 14:38:14 2022

This code creates binned radial density profiles for each halo that
satisfies more than 500 particles and virrad/hmrad <3.
This is done from files containing all halo particle information.

@author: clara
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import h5py
import os
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)

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


def radial_density(partx, party, partz, mass, binsize, halox, haloy, haloz):
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
    #create lists
    density = []
    rad_lowerbound = []
    uncertainties = []
    
    #find radial distance
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    mass = np.array(mass)
    
    
    #set up for indexing
    bin_index = np.argsort(dis)
    radius_lowerbound = 0
    bin_lowerbound = 0
    upperbound =(bin_lowerbound+binsize)
    j = 0
    
    while j < 1:        
        #ensure no overflow
        if upperbound > len(dis):
            upperbound = len(dis)
            j = 1
        if bin_lowerbound == upperbound:
            break
        
        #Find radial shell volume
        index_in_bin = bin_index[bin_lowerbound:upperbound]
        radius_upperbound = dis[index_in_bin][-1]
        dV = (4/3)*math.pi*(np.power(radius_upperbound,3)-np.power(radius_lowerbound,3))
        
        #find mass and density in shell
        M = np.sum(mass[index_in_bin])
        subdensity = (M/(dV))
        density = np.append(density, subdensity)
        
        #reset bounds for next iteration
        rad_lowerbound = np.append(rad_lowerbound, (radius_upperbound+radius_lowerbound)/2)
        dn = len(index_in_bin)
        uncertainties = np.append(uncertainties, subdensity/np.sqrt(dn))
        radius_lowerbound = radius_upperbound
        bin_lowerbound = upperbound
        upperbound = bin_lowerbound+binsize
    
    return(density, rad_lowerbound, uncertainties)


#Read in basic Subhalo Group Info
subhalo_info = pd.read_csv('../50-1-subhalo-info-dark.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy()

#halo to start with
g = 4000000
numhalos = len(subhalo_index)



pdheader = ['Radius','Density','Uncertainty']

#iterating through every halo
while g < 5000000:
    
    #excluding extremely massive halos for storage reasons
    if g not in [0,63864,96762,117250,184931,198182,143880,208811,229933,220595,167392,253861,242788,264883]:
        
        #read in particle data
        chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['x'],dtype={'x':object})
        chunk = chunk[chunk['x'] != 'x']
        partx = chunk['x'].to_numpy().astype(float)
        #print(partx)
        chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['y'],dtype={'y':object})
        chunk = chunk[chunk['y'] != 'y']
        party = chunk['y'].to_numpy().astype(float)
        chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['z'],dtype={'z':object})
        chunk = chunk[chunk['z'] != 'z']
        partz = chunk['z'].to_numpy().astype(float)
        chunk = pd.read_csv('../FullRunDark/snap_99_halo_'+str(g)+'-dark.csv', usecols=['mass'],dtype={'mass':object})
        chunk = chunk[chunk['mass'] != 'mass']
        mass = chunk['mass'].to_numpy().astype(float)
        
        #create new filename to store rad-den-data
        filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den-dark'
        
        #check if file already exists and delete faulty files
        if(os.path.isfile(filename+'.csv')):
        
            #os.remove() function to remove the file
            os.remove(filename+'.csv')
        
            #Printing the confirmation message of deletion
            print("File Deleted successfully")
        else:
            print("File does not exist")
        
        #condition if halo is large neough
        if len(partx) > 500:
            
            
            #50 bins to bin over
            binsize = int(len(partx)/50)
            
            #get radial density values
            rad_den = radial_density((partx*h), (party*h), (partz*h),(mass*h*(10**10)), binsize, (positionsX[g]*h), (h*positionsY[g]), (h*positionsZ[g]))
            #mass in solar masses
            #distances in kpc
            
            #virden = 200*p_crit
            
            #create data frame which then will be added to file
            miniderek = pd.DataFrame(columns=pdheader)
            miniderek['Radius']=rad_den[1]
            #miniderek['Virial Radius']=[virrad]*len(rad_den[0])
            miniderek['Density']=rad_den[0]
            miniderek['Uncertainty']=rad_den[2]
            print(miniderek)
            miniderek.to_csv(filename+'.csv')
            miniderek = miniderek[0:0]
    g += 1 