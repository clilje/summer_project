import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv


"""
Reading out the positions from the TNG group files.

"""
filename = []
pos = np.array([])
i = 0
# Making a list of all possible filenames
while i < 11:
    filename.append("/disk01/rmcg/downloaded/tng/tng50-4/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
#print(filename)
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


"""
Reading out the particle data from TNG snapshot files

"""
filenamesnap = []
#particleprop = np.zeros(4)
#particlepos = np.array([])
#particlemass = np.array([])
i = 0
# Making a list of all possible filenames
while i < 11:
    filenamesnap.append("/disk01/rmcg/downloaded/tng/tng50-4/snapshot_99/snap_099."+str(i)+".hdf5")
    i += 1



"""
Get associated Particles within 10^5 Parsec for each Halo
"""


def distancefromcentre(cx, cy, cz, x, y, z, ):
     
    x1 = math.pow((x-cx), 2)
    y1 = math.pow((y-cy), 2)
    z1 = math.pow((z-cz), 2)
    return (math.sqrt((x1 + y1 + z1))) # distance between the centre and given point
     

#First example for one Halo
xhalo, yhalo, zhalo = pos[0]
#distance to check in:
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
