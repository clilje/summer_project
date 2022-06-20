import h5py
import numpy as np
import matplotlib.pyplot as plt
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
Function that gets Subhalo Mass from given file list.
"""
        
def get_mass(filename):
    g = 0
    mass = np.array([])
    while g < len(filename):
        with h5py.File(str(filename[g])) as file:
            if 'Subhalo'in file:
                print(file['Subhalo'])
                submass = np.array(file['Subhalo/SubhaloMass'])
                mass = np.append(mass, submass)
            g +=1
    return(mass)




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
    while i < (len(mass_interval)-1):
        dm = (lowerbound[i+1]-lowerbound[i])
        dn = len(np.where(np.logical_and(mass_list>=lowerbound[i], mass_list<=(lowerbound[(i+1)])))[0])
        diff_number_per_mass = np.append(diff_number_per_mass, (dn/(dm*volume)))
        x_mass_lowerbound = np.append(x_mass_lowerbound, lowerbound[i])
        i += 1
    return(diff_number_per_mass, x_mass_lowerbound)

interval = np.logspace(0.1, 5, 50)


volume50 = pow(51.7, 3)  #in Mpc
volume100 = pow(110.7,3)

files504 = get_filenames(50, 4, 11)
mass504 = get_mass(files504)
mass_function_504 = mass_func(mass504,interval, volume50)

files501 = get_filenames(50, 1, 680)

mass501 = get_mass(files501)
mass_function_504 = mass_func(mass501,interval, volume50)

files1003 = get_filenames(100, 3, 7)
mass1003 = get_mass(files1003)
mass_function_1003 = mass_func(mass1003, interval, volume100)

files1002 = get_filenames(100, 2, 56)
mass1002 = get_mass(files1002)
mass_function_1002 = mass_func(mass1002, interval, volume100)


plt.loglog(mass_function_504[1], mass_function_504[0], "+", color="black", label="TNG50-4")
#plt.loglog(mass_function_501[1], mass_function_501[0], "x", color="red", label="TNG50-1")
plt.loglog(mass_function_1003[1], mass_function_1003[0], "_", color="blue", label="TNG100-3")
plt.loglog(mass_function_1002[1], mass_function_1002[0], "x", color="green", label="TNG100-2")

#plt.loglog(mass_function_all[1][int(length/2):(length-1)], mass_function_all[0][int(length/2):(length-1)], "v", color="green", label="Subselection 2")
plt.xlabel(r'Mass in $10^{10}$ $M_{\odot}$$h^{-1}$')
plt.ylabel(r'$dn/dlog_{10}(M)$ ($Mpc^{-3}$$h^{-1}$)')
plt.legend()
plt.savefig('resolutionmassfunc')
plt.show()