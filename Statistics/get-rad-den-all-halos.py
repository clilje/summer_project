"""
Created on Thu Jun 16 14:38:14 2022

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
    #print(np.min(dis))
    #virV = (4/3)*math.pi*(np.power((virrad+10),3)-np.power((virrad-10),3))
    #virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = np.array(mass)
    #virM = np.sum(mass[virindex])
    #virdensity = virM/virV
    
    
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
        
        rad_lowerbound = np.append(rad_lowerbound, (radius_upperbound+radius_lowerbound)/2)
        dn = len(index_in_bin)
        uncertainties = np.append(uncertainties, subdensity/np.sqrt(dn))
        radius_lowerbound = radius_upperbound
        bin_lowerbound = bin_lowerbound+binsize
    
    
    #normalized = 
    below_virR = np.where((density).astype(float)<float(p_crit*200))[0]
    print(density)
    print(below_virR)
    print(density[below_virR[-10:-1]])
    #print(rad_lowerbound)
    virR = np.min(rad_lowerbound[below_virR])
    return(density, rad_lowerbound, uncertainties, virR)

subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy().astype(int)

gg = [59,150,74]
numhalos = len(subhalo_index)



pdheader = ['Radius','Density','Uncertainty','Virial Radius']
#derek = pd.DataFrame(columns=pdheader)

for g in gg:
    intervals = np.arange(0,length,900000)
    print(intervals)
    f = 0
    while f <= len(intervals):
        if f == 0:
            to_exclude = np.arange(intervals[1],length,1)
        else:
            to_exclude = np.arange(0,intervals[f]) + np.arange(intervals[f],length,1)
        data_csv = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', dtype={'':int,'ID':object,'Type':'string','x':float,'y':float,'z':float,'mass':float,'vx':float,'vy':float,'vz':float},skiprows=to_exclude)
        f = f+1
        partx = data_csv['x'].to_numpy().astype(float)
        party = data_csv['y'].to_numpy().astype(float)
        partz = data_csv['z'].to_numpy().astype(float)
        mass = data_csv['mass'].to_numpy().astype(float)
        filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den'
        rad_den = radial_density((data_csv['x'].to_numpy()*h), (data_csv['y'].to_numpy()*h), (data_csv['z'].to_numpy()*h),(data_csv['mass'].to_numpy()*h*(10**10)), 25, (positionsX[g]*h), (h*positionsY[g]), (h*positionsZ[g]))
        #mass in solar masses
        #distances in kpc
        virrad = rad_den[3]
        virden = 200*p_crit
        miniderek = pd.DataFrame(columns=pdheader)
        miniderek['Radius']=rad_den[1]/(virrad)
        miniderek['Virial Radius']=[virrad]*len(rad_den[0])
        miniderek['Density']=rad_den[0]/(virden)
        miniderek['Uncertainty']=rad_den[2]
        print(miniderek)
        miniderek.to_csv(filename+'.csv', mode='a')
        miniderek = miniderek[0:0]
        g += 1 
    #halo_50 = [positionsX[g],positionsY[g],positionsZ[g], radius[g], full_mass[g]]
    #print(full_mass[g])
    #print(np.sum(mass))
    

    
    
    
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