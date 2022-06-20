import h5py
import numpy as np
import matplotlib.pyplot as plt


"""
Reading out the masses from the TNG50-4 files.

"""
filename504 = []
mass504 = np.array([])
i = 0
# Making a list of all possible filenames
while i < 11:
    filename504.append("/disk01/rmcg/downloaded/tng/tng50-4/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
g = 0
print(filename504)

# Getting mass for each file
while g < len(filename504):
    with h5py.File(str(filename504[g])) as file:
        #print(file.keys())
        submass = np.array(file['Subhalo/SubhaloMass'])
        #print(submass)
        mass504 = np.append(mass504, submass)
        g +=1
print(len(mass504))
"""
Reading out the masses from the TNG50-1 files.

"""
"""
filename501 = []
mass501 = np.array([])
i = 0
# Making a list of all possible filenames
while i < 630:
    filename501.append("/disk01/rmcg/downloaded/tng/tng50-1/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
g = 3
print(filename501)
# Getting mass for each file
while g < len(filename501):
    with h5py.File(str(filename501[g])) as file:
        print(filename501[g])
        print(file.keys())
        print(np.array(file['Subhalo']))
        submass = np.array(file['Subhalo/SubhaloMass'])
        #print(submass)
        mass501 = np.append(mass501, submass)
        g +=1
print(len(mass501))

"""

"""
Reading out the masses from the TNG100-3 files.

"""
filename1003 = []
mass1003 = np.array([])
i = 0
# Making a list of all possible filenames
while i < 7:
    filename1003.append("/disk01/rmcg/downloaded/tng/tng100-3/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
g = 0
print(filename1003)

# Getting mass for each file
while g < len(filename1003):
    with h5py.File(str(filename1003[g])) as file:
        #print(file.keys())
        submass = np.array(file['Subhalo/SubhaloMass'])
        #print(submass)
        mass1003 = np.append(mass1003, submass)
        g +=1
print(len(mass1003))



"""
Reading out the masses from the TNG100-2 files.

"""
filename1002 = []
mass1002 = np.array([])
i = 0
# Making a list of all possible filenames
while i < 56:
    filename1002.append("/disk01/rmcg/downloaded/tng/tng100-2/fof_subfind_snapshot_99/fof_subhalo_tab_099."+str(i)+".hdf5")
    i += 1
g = 0
print(filename1002)

# Getting mass for each file
while g < len(filename1002):
    with h5py.File(str(filename1002[g])) as file:
        #print(file.keys())
        submass = np.array(file['Subhalo/SubhaloMass'])
        #print(submass)
        mass1002 = np.append(mass1002, submass)
        g +=1
print(len(mass1002))

"""
Calculating the mass function.
"""
#Manipulating mass for understanding
volume50 = pow(51.7, 3)  #in Mpc
volume100 = pow(110.7,3)


def mass_func(mass_list, mass_interval, volume):
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
    diff_number_per_mass = []
    x_mass_lowerbound = []
    #lowerbound = np.add(min_mass, mass_interval)
    lowerbound = mass_interval
    i = 0
    while i < (len(interval)-1):
        #print(lowerbound[i])
        dm = (lowerbound[i+1]-lowerbound[i])
        #dm = lowerbound[i]
        #print(dm)
        dn = len(np.where(np.logical_and(mass_list>=lowerbound[i], mass_list<=(lowerbound[(i+1)])))[0])
        diff_number_per_mass = np.append(diff_number_per_mass, (dn/(dm*volume)))
        #diff_number_per_mass = np.append(diff_number_per_mass, (len(np.where(np.logical_and(mass_list>=lowerbound[i], mass_list<=(lowerbound[(i+1)])))[0]))/(lowerbound[i+1]-lowerbound[i]))
        x_mass_lowerbound = np.append(x_mass_lowerbound, lowerbound[i])
        #lowerbound += mass_interval[i]
        i += 1
    return(diff_number_per_mass, x_mass_lowerbound)

interval = np.logspace(0.1, 5, 50)
mass_function_504 = mass_func(mass504,interval, volume50)
#mass_function_501 = mass_func(mass501,interval)
mass_function_1003 = mass_func(mass1003, interval, volume100)
mass_function_1002 = mass_func(mass1002, interval, volume100)


#length = len(mass_function_all[0])
plt.loglog(mass_function_504[1], mass_function_504[0], "+", color="black", label="TNG50-4")
#plt.loglog(mass_function_501[1], mass_function_501[0], "x", color="red", label="TNG50-1")
plt.loglog(mass_function_1003[1], mass_function_1003[0], "_", color="blue", label="TNG100-3")
plt.loglog(mass_function_1002[1], mass_function_1002[0], "x", color="green", label="TNG100-2")

#plt.loglog(mass_function_all[1][int(length/2):(length-1)], mass_function_all[0][int(length/2):(length-1)], "v", color="green", label="Subselection 2")
plt.xlabel(r'Mass in $10^{10}$ $M_{\odot}$$h^{-1}$')
plt.ylabel(r'$dn/dlog_{10}(M)$ ($Mpc^{-3}$$h^{-1}$)')
plt.legend()
plt.savefig('resmassfunc')
plt.show()