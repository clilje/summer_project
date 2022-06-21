import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv


"""
Reading out the positions from the TNG group files.

"""
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
'''
filename = []
pos = np.array([])
i = 0
# Making a list of all possible filenames
while i < 11:
    filename.append("/disk01/rmcg/downloaded/tng/tng50-4/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
#print(filename)

'''
'''
g = 0

# Getting position for each file
while g < len(filename):
    with h5py.File(str(filename[g])) as file:
        #print(file.keys())
        subpos = np.array((file['Subhalo/SubhaloPos']))
        #print(subpos)
        pos= np.append(pos, subpos)
        g +=1
pos = np.reshape(pos, [int(len(pos)/3),3])
#print(pos)
print(np.shape(np.array(pos)))

'''
"""
Function that gets Subhalo Position from given file list.
"""
        
def get_pos(filename):
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


"""
Function that gets Subhalo Half Mass Radius from given file list.
"""
        
def get_rad(filename):
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
Reading out the particle data from TNG snapshot files

"""
'''
filenamesnap = []
#particleprop = np.zeros(4)
#particlepos = np.array([])
#particlemass = np.array([])
i = 0
# Making a list of all possible filenames
while i < 11:
    filenamesnap.append("/disk01/rmcg/downloaded/tng/tng50-4/snapshot_99/snap_099."+str(i)+".hdf5")
    i += 1
'''

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
     
    #x1 = np.power(np.subtract(x,cx), 2)
    #y1 = np.power(np.subtract(y,cy), 2)
    #z1 = np.power(np.subtract(z,cz), 2)
    return (np.sqrt((np.power(np.subtract(x,cx), 2)+ np.power(np.subtract(y,cy), 2) + np.power(np.subtract(z,cz), 2)))) # distance between the centre and given point


def particle_from_halo(num_halo, position, halfmassrad, rad_to_check, filenamesnap):
    xhalo, yhalo, zhalo = position
    r = rad_to_check
    g = 0
    with open('50-1_snap_99_halo_'+str(num_halo)+'_rad_mass_100kpc.csv', 'w', encoding='UTF8', newline='') as f:
        header = ['x','y','z','mass']
        # Create a writer object
        fwriter = csv.writer(f, delimiter=',')
        #print(fwriter)
        #print(header)
        # Write the header
        fwriter.writerow(header)
        while g < len(filenamesnap):
            print(filenamesnap[g])
            with h5py.File(filenamesnap[g]) as file:
                            
                #DMParticles
                partpos = np.array(file['PartType1/Coordinates'])#
                print(np.shape(partpos))
                #partpos = np.reshape(partpos, [int(len(partpos)/3),3])
                #b = 0
                c = 0
                header = dict( file['Header'].attrs.items() )
                dmmass = np.ones(len(partpos))*header['MassTable'][1]  # 10^10 Msun/h
                #while b < len(partpos):
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                print(dis)
                nindex = np.where(dis<r)
                print(nindex)
                #data = [partpos[nindex][:, 0], partpos[nindex][:, 1], partpos[nindex][:, 2], dmmass]
                print(partpos)
                print(dmmass.T)
                data = np.append(partpos, dmmass.T, axis = 1)
                print(data)
                fwriter.writerows(data)
                print(len(data))
                c = len(data)
                '''
                if dis <(r):
                    data = [partpos[b][0], partpos[b][1], partpos[b][2],dmmass]
                    fwriter.writerow(data)
                    c += 1'''
                    
                #b += 1
                #print(c)
                
                
                if c != 0:
                    #GasParticles
                    partpos = np.array(file['PartType0/Coordinates'])
                    mass0 = np.array(file['PartType0/Masses'])
                    #b = 0
                    #c = 0
                    #while b < len(partpos):
                    dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                    nindex = np.where(dis<r)
                    data = [partpos[nindex][:, 0], partpos[nindex][:, 1], partpos[nindex][:, 2], mass0]
                    '''if dis <(r):
                            data = [partpos[b][0], partpos[b][1], partpos[b][2],mass0[b]]
                            fwriter.writerow(data)
                            print(data)
                            c += 1
                        b += 1'''
                    fwriter.writerows(data)
                    print(len(data))
                    c = len(data)
                    print(c)
                
                

                    #Stellar and Wind particles
                    partpos = np.array(file['PartType4/Coordinates'])
                    #b = 0
                    #c = 0
                    mass4 = np.array(file['PartType4/Masses'])
                    dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                    nindex = np.where(dis<r)
                    data = [partpos[nindex][:, 0], partpos[nindex][:, 1], partpos[nindex][:, 2], mass4]
                    fwriter.writerows(data)
                    print(len(data))
                    c = len(data)
                    '''while b < len(partpos):
                        dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[b][0], partpos[b][1], partpos[b][2])
                        if dis <(r):
                            data = [partpos[b][0], partpos[b][1], partpos[b][2], mass4[b]]
                            fwriter.writerow(data)
                            c += 1
                        b += 1'''
                    
                    print(c)
                    
                    #Black Hole Particles
                    partpos = file['PartType5/Coordinates']
                    #b = 0
                    #c = 0
                    mass5 = np.array(file['PartType5/Masses'])
                    dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[:, 0], partpos[:, 1], partpos[:, 2])
                    nindex = np.where(dis<r)
                    data = [partpos[nindex][:, 0], partpos[nindex][:, 1], partpos[nindex][:, 2], mass5]
                    fwriter.writerows(data)
                    print(len(data))
                    c = len(data)
                    '''
                    while b < len(partpos):
                        dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[b][0], partpos[b][1], partpos[b][2])
                        if dis <(r):
                            data = [partpos[b][0], partpos[b][1], partpos[b][2], mass5[0]]
                            fwriter.writerow(data)
                            c +=1
                        b += 1'''
                    print(c)
            g+= 1

    
    
    

filename_group = get_filenames(50, 4, 11)
print(filename_group)
pos = get_pos(filename_group)
print(pos)
halfmassradii = get_rad(filename_group)
print(halfmassradii)
#First example for one Halo
#xhalo, yhalo, zhalo = pos[0]
#hmrad = halfmassradii[0]
#distance to check in:
filename_snap = get_filenames_snap(50, 4, 11)
print(filename_snap)
particle_from_halo(0, pos[0], halfmassradii[0], 100, filename_snap)
particle_from_halo(1, pos[1], halfmassradii[1], 100, filename_snap)
particle_from_halo(2, pos[2], halfmassradii[2], 100, filename_snap)
particle_from_halo(3, pos[3], halfmassradii[3], 100, filename_snap)
particle_from_halo(4, pos[4], halfmassradii[4], 100, filename_snap)
particle_from_halo(5, pos[5], halfmassradii[5], 100, filename_snap)

"""
r = 100
g = 0

# Getting position for each file
#datafile = open('099_0_halo-rad_mass.txt','w')
with open('snap_99_halo_0_rad_mass_100kpc.csv', 'w', encoding='UTF8', newline='') as f:
    header = ['rad', 'mass']
    # Create a writer object
    fwriter = csv.writer(f, delimiter=',')
    print(fwriter)
    print(header)
    # Write the header
    fwriter.writerow(header)
    #datafile.write(str(header[0])+','+str(header[1]))
    while g < len(filenamesnap):
        print(filenamesnap[g])
        with h5py.File(filenamesnap[g]) as file:
            
            #GasParticles
            partpos = np.array(file['PartType0/Coordinates'])
            mass0 = np.array(file['PartType0/Masses'])
            #print(partpos)
            #print(xhalo,yhalo,zhalo)
            b = 0
            c = 0
            #print(np.shape(partpos))
            #print(len(partpos))
            while b < len(partpos):
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[b][0], partpos[b][1], partpos[b][2])
                #if partpos[b][0] >= (xhalo-(10)) and partpos[b][0] <= (xhalo+(10)) and partpos[b][1] >= (yhalo-(10)) and partpos[b][1] <= (yhalo+(10)) and partpos[b][2] >= (zhalo-(10)) and partpos[b][2] <= (zhalo+(10)):
                if dis <(r):
                    #particlepos = np.append(particlepos, [partpos[b][0],partpos[b][1],partpos[b][2]])
                    #particlemass = np.append(particlemass, np.array(file['PartType0/Masses']))
                    data = [dis,mass0[b]]
                    #datafile.write('\n'+str(dis)+','+str(mass0[b]))
                    fwriter.writerow(data)
                    print(data)
                    #print("Success")
                    #print(partpos[b])
                    #print(b)
                    c += 1
                b += 1
            print(c)
            
            
        
            #DMParticles
            partpos = file['PartType1/Coordinates']
            b = 0
            c = 0
            header = dict( file['Header'].attrs.items() )
            dmmass = header['MassTable'][1]  # 10^10 Msun/h
            while b < len(partpos):
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[b][0], partpos[b][1], partpos[b][2])
                #if partpos[b][0] >= (xhalo-(10)) and partpos[b][0] <= (xhalo+(10)) and partpos[b][1] >= (yhalo-(10)) and partpos[b][1] <= (yhalo+(10)) and partpos[b][2] >= (zhalo-(10)) and partpos[b][2] <= (zhalo+(10)):
                if dis <(r):
                    #particlepos = np.append(particlepos, [partpos[b][0],partpos[b][1],partpos[b][2]])
                    #particlemass = np.append(particlemass, np.array(file['PartType0/Masses']))
                    data = [dis,dmmass]
                    fwriter.writerow(data)
                    #datafile.write('\n'+str(dis)+','+str(dmmass))
                    #print(partpos[b])
                    #print(b)
                    c += 1
                b += 1
            print(c)
            
            #Stellar and Wind particles
            partpos = file['PartType4/Coordinates']
            b = 0
            c = 0
            mass4 = np.array(file['PartType4/Masses'])
            while b < len(partpos):
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[b][0], partpos[b][1], partpos[b][2])
                #if partpos[b][0] >= (xhalo-(10)) and partpos[b][0] <= (xhalo+(10)) and partpos[b][1] >= (yhalo-(10)) and partpos[b][1] <= (yhalo+(10)) and partpos[b][2] >= (zhalo-(10)) and partpos[b][2] <= (zhalo+(10)):
                if dis <(r):
                    #particlepos = np.append(particlepos, [partpos[b][0],partpos[b][1],partpos[b][2]])
                    #particlemass = np.append(particlemass, np.array(file['PartType0/Masses']))
                    data = [dis, mass4[b]]
                    fwriter.writerow(data)
                    #datafile.write('\n'+str(dis)+','+str(mass4[b]))
                    #print(partpos[b])
                    #print(b)
                    c += 1
                b += 1
                
            print(c)
    
            #Black Hole Particles
            partpos = file['PartType5/Coordinates']
            b = 0
            c = 0
            mass5 = np.array(file['PartType5/Masses'])
            
            while b < len(partpos):
                dis = distancefromcentre(xhalo, yhalo, zhalo, partpos[b][0], partpos[b][1], partpos[b][2])
                #if partpos[b][0] >= (xhalo-(10)) and partpos[b][0] <= (xhalo+(10)) and partpos[b][1] >= (yhalo-(10)) and partpos[b][1] <= (yhalo+(10)) and partpos[b][2] >= (zhalo-(10)) and partpos[b][2] <= (zhalo+(10)):
                if dis <(r):
                    #particlepos = np.append(particlepos, [partpos[b][0],partpos[b][1],partpos[b][2]])
                    #particlemass = np.append(particlemass, np.array(file['PartType0/Masses']))
                    data = [dis, mass5[0]]
                    fwriter.writerow(data)
                    #datafile.write('\n'+str(dis)+','+str(mass5[b]))
                    #print(partpos[b])
                    #print(b)
                    c +=1
                b += 1
            print(c)
            g+= 1
"""