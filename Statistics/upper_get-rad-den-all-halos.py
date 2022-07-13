"""
Created on Thu Jun 16 14:38:14 2022

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
    density = []
    rad_lowerbound = []
    uncertainties = []
    dis = distancefromcentre(halox, haloy, haloz, partx, party, partz)
    mass = np.array(mass)
    
    
    bin_index = np.argsort(dis)
    radius_lowerbound = 0
    bin_lowerbound = 0
    upperbound =(bin_lowerbound+binsize)
    j = 0
    
    while j < 1:           
        if upperbound > len(dis):
            upperbound = len(dis)
            j = 1
        if bin_lowerbound == upperbound:
            break
        index_in_bin = bin_index[bin_lowerbound:upperbound]
        radius_upperbound = dis[index_in_bin][-1]
        dV = (4/3)*math.pi*(np.power(radius_upperbound,3)-np.power(radius_lowerbound,3))
        
        M = np.sum(mass[index_in_bin])
        subdensity = (M/(dV))
        density = np.append(density, subdensity)
        
        rad_lowerbound = np.append(rad_lowerbound, (radius_upperbound+radius_lowerbound)/2)
        dn = len(index_in_bin)
        uncertainties = np.append(uncertainties, subdensity/np.sqrt(dn))
        radius_lowerbound = radius_upperbound
        bin_lowerbound = upperbound
        upperbound = bin_lowerbound+binsize
    """
    above_virR = np.where((density).astype(float)>float(p_crit*200))[0]
    virR = np.max(rad_lowerbound[above_virR])
    print(virR)
    print(rad_lowerbound[-10:-1])
    """
    return(density, rad_lowerbound, uncertainties)

subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy()

g = 0
numhalos = len(subhalo_index)



pdheader = ['Radius','Density','Uncertainty']
#derek = pd.DataFrame(columns=pdheader)

while g < numhalos:
    
    if g not in [0,63864,96762,117250,184931,198182,143880,208811,229933,220595,167392,253861,242788,264883]:

        #chunksize = 999999
        #chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', dtype={'':int,'ID':object,'Type':'string','x':float,'y':float,'z':float,'mass':float,'vx':float,'vy':float,'vz':float})
        #y = 0
        #for chunk in data_csv:
        #process(chunk)  
        chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['x'],dtype={'x':object})
        chunk = chunk[chunk['x'] != 'x']
        print('success')
        partx = chunk['x'].to_numpy().astype(float)
        print(partx)
        #exit()
        chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['y'],dtype={'y':object})
        chunk = chunk[chunk['y'] != 'y']
        party = chunk['y'].to_numpy().astype(float)
        chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['z'],dtype={'z':object})
        chunk = chunk[chunk['z'] != 'z']
        partz = chunk['z'].to_numpy().astype(float)
        chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['mass'],dtype={'mass':object})
        chunk = chunk[chunk['mass'] != 'mass']
        mass = chunk['mass'].to_numpy().astype(float)
        filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den'
        if(os.path.isfile(filename+'.csv')):
        
            #os.remove() function to remove the file
            os.remove(filename+'.csv')
        
            #Printing the confirmation message of deletion
            print("File Deleted successfully")
        else:
            print("File does not exist")
        
        
        if len(partx) > 500:
            
            
            
            binsize = int(len(partx)/50)
            
            rad_den = radial_density((partx*h), (party*h), (partz*h),(mass*h*(10**10)), binsize, (positionsX[g]*h), (h*positionsY[g]), (h*positionsZ[g]))
            #mass in solar masses
            #distances in kpc
            #virrad = rad_den[3]
            virden = 200*p_crit
            miniderek = pd.DataFrame(columns=pdheader)
            miniderek['Radius']=rad_den[1]
            #miniderek['Virial Radius']=[virrad]*len(rad_den[0])
            miniderek['Density']=rad_den[0]
            miniderek['Uncertainty']=rad_den[2]
            print(miniderek)
            miniderek.to_csv(filename+'.csv')
            miniderek = miniderek[0:0]
            #y += 1
            #print('chunk'+str(y))
    g += 1 
            
            

    
    
    
    """
    #ax = plt.axes(projection =None)
    plt.errorbar(rad_den[1]/(virrad), rad_den[0]/(virden), yerr=rad_den[2]/virden, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    plt.axhline(1)
    plt.xlabel(r'Radial Distance normalized by $r_{200}$ in kpc')
    plt.ylabel(r'Density normalized by $\rho_{200}$ in $M_{\odot}/kpc^{-3}$')
    plt.xscale('log')
    plt.yscale('log')
    
    plt.savefig('HaloFitsInfo/fit-profiles-halo-'+str(g)+'.jpg')
    plt.clf()
    """
    '''
    fig = plt.figure()
    ax = plt.axes(projection ='3d')
    xyz = np.arange(len(partx))
    index = np.random.choice(xyz,200)
    ax.scatter(partx[index], party[index], partz[index], marker='+',color='blue',alpha=0.1)
    ax.scatter(halo_50[0], halo_50[1], halo_50[2], marker='+',color='red')
    #ax.scatter(positionsX[index_sub],positionsY[index_sub],positionsZ[index_sub],marker='x', color='black')
    #ax.scatter(CM[0], CM[1], CM[2], marker='+',color='pink')
    
    ax.set_xlabel('x [ckpc/h]')
    
    ax.set_ylabel('y [ckpc/h]')
    ax.set_zlabel('z [ckpc/h]')
    fig.savefig('HaloFitsInfo/halocomp-'+str(g))
    '''