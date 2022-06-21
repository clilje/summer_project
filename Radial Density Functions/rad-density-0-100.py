import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv


"""
Creating function to get filenames for given simulation.
"""
def get_filenames(sim_size, sim_res, num_files):
    filename = []
    i = 0
    # Making a list of all possible filenames
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
        i += 1
    return(filename)

"""
Function that gets Subhalo Position from given file list.
"""
        
def get_pos(filename):
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


"""
Function that gets Subhalo Half Mass Radius from given file list.
"""
        
def get_rad(filename):
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
    filename = []
    i = 0
    # Making a list of all possible filenames
    
    while i < num_files:
        filename.append("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/snapshot_99/snap_099."+str(i)+".hdf5")
        i += 1
    return(filename)


"""
Get associated Particles within 10^5 Parsec for each Halo
"""


def distancefromcentre(cx, cy, cz, x, y, z, ):
    
    return (np.sqrt((np.power(np.subtract(x,cx), 2)+ np.power(np.subtract(y,cy), 2) + np.power(np.subtract(z,cz), 2)))) # distance between the centre and given point


def particle_from_halo(num_halo, position, halfmassrad, rad_to_check, filenamesnap):
    xhalo, yhalo, zhalo = position
    r = rad_to_check
    g = 0
    
    
    with open('HaloParticles/50-1_snap_99_halo_'+str(num_halo)+'_rad_mass_100kpc.csv', 'w', encoding='UTF8', newline='') as f:
        
        header = ['x','y','z','mass']
        # Create a writer object
        fwriter = csv.writer(f, delimiter=',')
        # Write the header
        fwriter.writerow(header)
        
        while g < len(filenamesnap):
            
            print(filenamesnap[g])
            
            
            with h5py.File(filenamesnap[g]) as file:
                            
                #DMParticles
                partpos = np.array(file['PartType1/Coordinates'])#
                c = 0
                header = dict( file['Header'].attrs.items() )
                dmmass = np.ones(np.size(partpos, axis=0))*header['MassTable'][1]  # 10^10 Msun/h
                
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                
                nindex = (np.where(dis<r))[0]
                
                if len(nindex) !=0:
                    data = np.hstack((partpos[nindex], np.atleast_2d(dmmass[nindex]).T))
                    fwriter.writerows(data)
                    print(len(nindex))
                    c = len(data)
                else:
                    print("No DM found")
                
                
                #GasParticles
                partpos = np.array(file['PartType0/Coordinates'])
                mass0 = np.array(file['PartType0/Masses'])
                
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                
                nindex = (np.where(dis<r))[0]
                
                if len(nindex) != 0:
                    data = np.hstack((partpos[nindex], np.atleast_2d(mass0[nindex]).T))
                    fwriter.writerows(data)
                    c = len(data)
                    print(c)
                else:
                    print("No Gas found")
            

                #Stellar and Wind particles
                partpos = np.array(file['PartType4/Coordinates'])
                mass4 = np.array(file['PartType4/Masses'])
                
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                
                nindex = (np.where(dis<r))[0]
                
                if len(nindex) !=0:
                    data = np.hstack((partpos[nindex], np.atleast_2d(mass4[nindex]).T))
                    fwriter.writerows(data)
                    print(len(data))
                    c = len(data)
                    
                else:
                    print("No Stellar part found")
                
                
                #Black Hole Particles
                partpos = file['PartType5/Coordinates']
                mass5 = np.array(file['PartType5/Masses'])
                
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                
                nindex = (np.where(dis<r))[0]
                
                if len(nindex) !=0:
                    data = np.hstack((partpos[nindex], np.atleast_2d(mass5[nindex]).T))
                    fwriter.writerows(data)
                    c = len(data)
                    print(c)
                else:
                    print("No BH found")
            g+= 1

    
    
    

filename_group = get_filenames(50, 4, 11)

pos = get_pos(filename_group)

halfmassradii = get_rad(filename_group)

filename_snap = get_filenames_snap(50, 4, 11)


g = 0 
while g < 20:    
    particle_from_halo(g, pos[g], halfmassradii[g], (halfmassradii[g]*3), filename_snap)
    g += 1
    
'''
particle_from_halo(1, pos[1], halfmassradii[1], (halfmassradii[1]*3), filename_snap)
particle_from_halo(2, pos[2], halfmassradii[2], (halfmassradii[2]*3), filename_snap)
particle_from_halo(3, pos[3], halfmassradii[3], (halfmassradii[3]*3), filename_snap)
particle_from_halo(4, pos[4], halfmassradii[4], (halfmassradii[4]*3), filename_snap)
particle_from_halo(5, pos[5], halfmassradii[5], (halfmassradii[5]*3), filename_snap)

'''