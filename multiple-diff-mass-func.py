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
#print(str(massdark[0]), str(massdark[1]))
#print(len(massdark))



"""
Calculating the mass function.
"""
#Manipulating mass for understanding
#max_mass = np.max(mass)
#min_mass = np.amin(mass)
volume = pow(51.7, 3)  #in Mpc
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

#Preparing use of mass function
interval = np.logspace(0.1, 4.1, 50)
#interval = np.geomspace(min_mass, max_mass, 50)
#interval = np.linspace(min_mass, max_mass, 100)
#print(interval)
#print(min_mass)
#print(max_mass)
mass_function_all = mass_func(mass,interval)
mass_function_dark = mass_func(massdark,interval)

#print(mass_function)
#print(sum(mass_function[0]))
#normvals = ((mass_function[1]-min_mass)/(max_mass-min_mass))
#print(normvals)
length = len(mass_function_all[0])
plt.loglog(mass_function_all[1], mass_function_all[0], "+", color="black", label="DM + Baryonic")
plt.loglog(mass_function_dark[1], mass_function_dark[0], "x", color="red", label="DM only")
plt.loglog(mass_function_all[1][0:int(length/2)], mass_function_all[0][0:int(length/2)], "_", color="blue", label="Subselection 1")
plt.loglog(mass_function_all[1][int(length/2):(length-1)], mass_function_all[0][int(length/2):(length-1)], "v", color="green", label="Subselection 2")
plt.xlabel(r'Mass in $10^{10}$ $M_{\odot}$')
plt.ylabel(r'$dn/dlog_{10}(M)$ ($Mpc^{-3}$)')
plt.legend()
plt.savefig('selcomparemassfunc')
plt.show()