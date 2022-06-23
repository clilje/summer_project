import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
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
                subpos = np.array(file['Subhalo/SubhaloPos'])
                pos = np.append(pos, subpos)
            g +=1
    
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
    return(rad)




def get_filenames_snap(sim_size, sim_res, num_files):
    """
    

    Parameters
    ----------
    sim_size : Size of Edge of simulation box in kpc.
    sim_res : Resolution of simulation, integer from 1 to 4
    num_files : number of files in given simulation folder

    Returns
    -------
    List of filnames in string formats for each snapshot

    """
    filename = []
    i = 0
    # Making a list of all possible filenames
    
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/snapshot_99/snap_099."+str(i)+".hdf5")
        i += 1
    return(filename)


def get_filenames_snap_dark(sim_size, sim_res, num_files):
    """
    

    Parameters
    ----------
    sim_size : Size of Edge of simulation box in kpc.
    sim_res : Resolution of simulation, integer from 1 to 4
    num_files : number of files in given simulation folder

    Returns
    -------
    List of filnames in string formats for each snapshot

    """
    filename = []
    i = 0
    # Making a list of all possible filenames
    
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"-dark/snapshot_99/snap_099."+str(i)+".hdf5")
        i += 1
    return(filename)


"""
Creating function to get filenames for given dark simulation.
"""
def get_filenames_dark(sim_size, sim_res, num_files):
    filename = []
    i = 0
    # Making a list of all possible filenames
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"-dark/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
        i += 1
    return(filename)


"""
Reading out info from matching catalogue
"""
def get_matching(sim_size, sim_res):
    with h5py.File("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/subhalo_matching_to_dark.hdf5") as file:
        #print(file.keys())
        matchingarr = np.array(file['Snapshot_99/SubhaloIndexDark_LHaloTree'])
        #print(matchingarr)
    return(matchingarr)


def distancefromcentre(cx, cy, cz, x, y, z):
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


def particle_from_halo(num_halo, position, halfmassrad, rad_to_check, filenamesnap):
    """
    

    Parameters
    ----------
    num_halo : Number of halo, index style
    position : position of centre of halo, three coordinates
    halfmassrad : halfmass radius of halo
    rad_to_check : radius to find particles within
    filenamesnap : names of files in which particles are stored

    Returns
    -------
    No direct output
    File is written with associated halo positions.

    """
    xhalo, yhalo, zhalo = position
    r = rad_to_check
    g = 0
    
    
    with open('HaloParticles/50-1_snap_99_halo_'+str(num_halo)+'_pos_mass.csv', 'w', encoding='UTF8', newline='') as f:
        
        header = ['x','y','z','mass']
        # Create a writer object
        fwriter = csv.writer(f, delimiter=',')
        # Write the header
        fwriter.writerow(header)
        
        while g < len(filenamesnap):
            
            print(filenamesnap[g])
            path = Path(filenamesnap[g])
            if path.is_file():
                with h5py.File(filenamesnap[g]) as file:
                                    
                    #DMParticles
                    if 'PartType1/Coordinates'in file:
                        partpos = np.array(file['PartType1/Coordinates'])#
                        c = 0
                        header = dict( file['Header'].attrs.items() )
                        dmmass = np.ones(np.size(partpos, axis=0))*header['MassTable'][1]  # 10^10 Msun/h
                        
                        dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                        
                        nindex = (np.where(dis<r))[0]
                        
                        if len(nindex) !=0:
                            data = np.hstack((partpos[nindex], np.atleast_2d(dmmass[nindex]).T))
                            fwriter.writerows(data)
                            print("Number of DM Particles"+str(len(nindex)))
                            c = len(data)
                        else:
                            print("No DM found")
                    
                    
                    #GasParticles
                    if 'PartType0/Coordinates'in file:
                        partpos = np.array(file['PartType0/Coordinates'])
                        mass0 = np.array(file['PartType0/Masses'])
                        
                        dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                        
                        nindex = (np.where(dis<r))[0]
                        
                        if len(nindex) != 0:
                            data = np.hstack((partpos[nindex], np.atleast_2d(mass0[nindex]).T))
                            fwriter.writerows(data)
                            c = len(data)
                            print("Number of Gas Particles"+str(c))
                        else:
                            print("No Gas found")
                
    
                    #Stellar and Wind particles
                    if 'PartType4/Coordinates'in file:
                        partpos = np.array(file['PartType4/Coordinates'])
                        mass4 = np.array(file['PartType4/Masses'])
                        
                        dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                        
                        nindex = (np.where(dis<r))[0]
                        
                        if len(nindex) !=0:
                            data = np.hstack((partpos[nindex], np.atleast_2d(mass4[nindex]).T))
                            fwriter.writerows(data)
                            c = len(data)
                            print("Number of Stellar Particles"+str(c))
                            
                        else:
                            print("No Stellar particles found")
                    
                    
                    #Black Hole Particles
                    if 'PartType5/Coordinates'in file:
                        partpos = file['PartType5/Coordinates']
                        mass5 = np.array(file['PartType5/Masses'])
                        
                        dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                        
                        nindex = (np.where(dis<r))[0]
                        
                        if len(nindex) !=0:
                            data = np.hstack((partpos[nindex], np.atleast_2d(mass5[nindex]).T))
                            fwriter.writerows(data)
                            c = len(data)
                            print("Number of BH Particles"+str(c))
                        else:
                            print("No BH found")
            g+= 1

    
    
    

#filename_group = get_filenames(50, 4, 11)

filename_dark = get_filenames_dark(50, 4, 4)

pos = get_pos(filename_dark)

halfmassradii = get_rad(filename_dark)

filename_snap = get_filenames_snap_dark(50, 4, 4)

matcharr = get_matching(50, 4)
darkindex = np.zeros(20)
i = 0
while i<20:
    if matcharr[i] == -1:
        print("No match for Halo"+str(i))
        darkindex[i]=-1
    else:
        darkindex[i] = matcharr[i]
    i +=1
print(matcharr)
g = 0 
#while g < len(halfmassradii)
while g < (np.max(darkindex)+1): 
    if (g in darkindex):
        particle_from_halo(g, pos[g], halfmassradii[g], (halfmassradii[g]*3), filename_snap)
    g += 1
    