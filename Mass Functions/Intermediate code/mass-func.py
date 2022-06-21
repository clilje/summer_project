import h5py
import numpy as np
import matplotlib.pyplot as plt


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
#print(str(mass[0]), str(mass[1]))
#print("Done")

"""
Calculating the mass function.
"""
#Manipulating mass for understanding
max_mass = np.max(mass)
min_mass = np.amin(mass)
#mass.sort()
#print(mass[-2], mass[-1])


def mass_func(mass_list, mass_interval):
    """
    
    Parameters
    ----------
    mass_list : List that includes all possible masses
    mass_interval : Array that gives the intervals at which mass function should be evaluated

    Returns
    -------
    x and y parameters of Mass function
    x: mass at intervals
    y: number of haloes with mass less than this

    """
    number_per_mass = []
    x_mass_lowerbound = []
    #lowerbound = np.add(min_mass, mass_interval)
    lowerbound = mass_interval
    i = 0
    while i < (len(interval)-1):
        #print(lowerbound[i])
        number_per_mass = np.append(number_per_mass, len(np.where(np.logical_and(mass_list>=lowerbound[i], mass_list<=(lowerbound[(i+1)])))[0]))
        x_mass_lowerbound = np.append(x_mass_lowerbound, lowerbound[i])
        #lowerbound += mass_interval[i]
        i += 1
    return(number_per_mass, x_mass_lowerbound)

#Preparing use of mass function
#interval = np.geomspace(min_mass, max_mass, 100)
interval = np.logspace(0.1, 5, 100)
#print(interval)
#print(min_mass)
#print(max_mass)
mass_function = mass_func(mass,interval)
print(mass_function)
print(sum(mass_function[0]))
plt.loglog(mass_function[1], mass_function[0], ".")
plt.xlabel('Mass in 10^10 Solar Units')
plt.ylabel('Number of haloes per mass increment')
plt.savefig('massfunction')
plt.show()