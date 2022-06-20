import h5py
import numpy as np
import matplotlib.pyplot as plt

"""
Reading out info from matching catalogue
"""
with h5py.File("/disk01/rmcg/downloaded/tng/tng50-4/subhalo_matching_to_dark.hdf5") as file:
    #print(file.keys())
    matchingarr = np.array(file['Snapshot_99/SubhaloIndexDark_LHaloTree'])
print(matchingarr)
"""
Reading out the masses from the TNG files.

"""
filename = []
mass = np.array([])
i = 0
# Making a list of all possible filenames
while i < 11:
    filename.append("/disk01/rmcg/downloaded/tng/tng50-4/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
#print(filename)
g = 0

# Getting mass for each file
while g < len(filename):
    with h5py.File(str(filename[g])) as file:
        #print(file.keys())
        submass = np.array(file['Subhalo/SubhaloMass'])
        #print(submass)
        mass = np.append(mass, submass)
        g +=1
print(str(mass[0]), str(mass[1]))
print(len(mass))
#print("Done")

"""
Reading out masses from dark files
"""
filenamedark = []
massdark = np.array([])
i = 0
# Making a list of all possible filenames
while i < 4:
    filenamedark.append("/disk01/rmcg/downloaded/tng/tng50-4-dark/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
#print(filename)
g = 0

# Getting mass for each file
while g < len(filenamedark):
    with h5py.File(str(filenamedark[g])) as file:
        #print(file.keys())
        submass = np.array(file['Subhalo/SubhaloMass'])
        #print(submass)
        massdark = np.append(massdark, submass)
        g +=1
print(str(massdark[0]), str(massdark[1]))
print(len(massdark))

"""
Determining matching pairs and noting down their masses as corresponding data points
"""
i = 0
orderdm = np.zeros(len(mass))
while i<len(matchingarr):
    if matchingarr[i] == -1:
        mass[i] = 0
        orderdm[i] = 0
    else:
        orderdm[i] = massdark[matchingarr[i]]
    i +=1
        


plt.loglog(mass, orderdm, "+", color="black")
plt.xlabel(r'Mass of DM + Baryonic in $10^{10}$ $M_{\odot}$')
plt.ylabel(r'Mass of DM only in $10^{10}$ $M_{\odot}$')
plt.savefig('matchedmasses')
plt.show()