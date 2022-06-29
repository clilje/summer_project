# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

@author: clara
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774


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
    # Making a list of all possible filenames
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
                #print(file['Subhalo'])
                subpos = np.array(file['Subhalo/SubhaloPos'])
                pos = np.append(pos, subpos)
            g +=1
        #print(pos)
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
                #print(file['Subhalo'])
                subrad = np.array(file['Subhalo/SubhaloHalfmassRad'])
                rad = np.append(rad, subrad)
            g +=1
        #pos = np.reshape(pos, [int(len(pos)/3),3])
    return(rad)

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


def radial_density(partx, party, partz, mass, binsize, virrad, halox, haloy, haloz):
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
    
    virV = (4/3)*math.pi*(np.power((virrad+10),3)-np.power((virrad-10),3))
    virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = np.array(mass)
    virM = np.sum(mass[virindex])
    virdensity = virM/virV
    
    
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

    return(density, rad_lowerbound, uncertainties, virdensity)


interval = np.logspace(0.1, 2.5, 100)
files = get_filenames(50, 4, 11)
positions = get_pos(files)
radius = get_rad(files)
#print(positions[2])
#print(radius[2])
halonumber = []
g = 23
numhalos = 24
densities = []
uncertainties = []
radii = []
while g < numhalos:
    data_csv = pd.read_csv('HaloParticles/50-4_snap_99_halo_'+str(g)+'_pos_mass.csv')
    rad_den = radial_density(data_csv['x'].to_numpy(), data_csv['y'].to_numpy(), data_csv['z'].to_numpy(),data_csv['mass'].to_numpy(), 40, radius[g], positions[g][0], positions[g][1], positions[g][2])
    #print(rad_den)
    hmrad = radius[g]
    densities.append(list(rad_den[0]))
    radii.append(list(rad_den[1]))
    uncertainties.append(list(rad_den[2]))
    hmden = rad_den[3]
    halonumber.append(g)
    #print(hmrad,hmden)
    
    g += 1 

densities = np.array(densities)
radii = np.array(radii)
uncertainties = np.array(uncertainties)
uncertainties[uncertainties == np.nan] = 0



def nfw(r, density_0, scale_radius):
    return(density_0/((r/scale_radius)*np.power((1+(r/scale_radius)),2)))

def einasto(r, density_e, r_e, n):
    d_n = (3*n)-(1/3)+(0.0079/n)
    return(density_e*np.exp((-1*d_n)*(np.power((r/r_e),(1/n))-1)))

def burkert(r, density_0, r_s):
    return((density_0*np.power(r_s,3))/((r+r_s)*(np.power(r,2)+np.power(r_s,2))))

def dehnen_twoparam(r, density_s, r_s):
    return(((2**6)*density_s)/((np.power((r/r_s),(7/9)))*np.power((1+(np.power((r/r_s),(4/9)))),6)))

def dehnen_threeparam(r, density_s, r_s, gamma):
    return(((2**6)*density_s)/((np.power((r/r_s),gamma))*np.power((1+(np.power((r/r_s),((3-gamma)/5)))),6)))

nfwfitp, nfwfitcov = scopt.curve_fit(nfw, radii[0]*h, densities[0]/(10*(h**2)), p0=[0.001,20], sigma=uncertainties[0])
nfwchi_square_test_statistic =  np.sum((np.square(((densities[0])/(10*(h**2)*hmden))-(nfw(radii[0], nfwfitp[0], nfwfitp[1])/hmden)))/(nfw(radii[0], nfwfitp[0], nfwfitp[1])/hmden))
nfwp_value = scipy.stats.distributions.chi2.sf(nfwchi_square_test_statistic,(len(densities[0])-1))
print ('ChiSquare and P values for NFW', nfwchi_square_test_statistic)
print ('Fitted value for NFW', nfwfitp)
print ('Uncertainties for NFW', np.sqrt(np.diag(nfwfitcov)))

einastofitp, einastofitcov = scopt.curve_fit(einasto, radii[0]*h, densities[0]/(10*(h**2)), p0=[0.001,10,4], sigma=uncertainties[0])
einastochi_square_test_statistic =  np.sum((np.square(((densities[0])/(10*(h**2)*hmden))-(einasto(radii[0], einastofitp[0], einastofitp[1],einastofitp[2])/hmden)))/(einasto(radii[0], einastofitp[0], einastofitp[1],einastofitp[2])/hmden))
einastop_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(densities[0])-1))
print ('ChiSquare and P values for Einasto', einastochi_square_test_statistic, einastop_value)
print ('Fitted value for Einasto', einastofitp)
print ('Uncertainties for Einasto', np.sqrt(np.diag(einastofitcov)))

burkertfitp, burkertfitcov = scopt.curve_fit(burkert,radii[0]*h, densities[0]/(10*(h**2)), p0=[0.1,10], sigma=uncertainties[0])
burkertchi_square_test_statistic =  np.sum((np.square(((densities[0])/(10*(h**2)*hmden))-(burkert(radii[0], burkertfitp[0], burkertfitp[1])/hmden)))/(burkert(radii[0], burkertfitp[0], burkertfitp[1])/hmden))
burkertp_value = scipy.stats.distributions.chi2.sf(burkertchi_square_test_statistic,(len(densities[0])-1))
print ('ChiSquare and P values for Burkert', burkertchi_square_test_statistic, burkertp_value)
print ('Fitted value for Burkert', burkertfitp)
print ('Uncertainties for Burkert', np.sqrt(np.diag(burkertfitcov)))


dehnen_twoparamfitp, dehnen_twoparamfitcov = scopt.curve_fit(dehnen_twoparam, radii[0]*h, densities[0]/(10*(h**2)), p0=[0.01,30], sigma=uncertainties[0])
dehnentwochi_square_test_statistic =  np.sum((np.square(((densities[0])/(10*(h**2)*hmden))-(dehnen_twoparam(radii[0], dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/hmden)))/(dehnen_twoparam(radii[0], dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/hmden))
dehnentwop_value = scipy.stats.distributions.chi2.sf(dehnentwochi_square_test_statistic,(len(densities[0])-1))
print ('ChiSquare and P values for dehnentwo', dehnentwochi_square_test_statistic, dehnentwop_value)
print ('Fitted value for Dehnen Two Parameters', dehnen_twoparamfitp)
print ('Uncertainties for Dehnen Two Parameters', np.sqrt(np.diag(dehnen_twoparamfitcov)))

dehnen_threeparamfitp, dehnen_threeparamfitcov = scopt.curve_fit(dehnen_threeparam, radii[0]*h, densities[0]/(10*(h**2)), p0=[0.01,25,0.02], sigma=uncertainties[0])
dehnenthreechi_square_test_statistic =  np.sum((np.square(((densities[0])/(10*(h**2)*hmden))-(dehnen_threeparam(radii[0],dehnen_threeparamfitp[0],dehnen_threeparamfitp[1],dehnen_threeparamfitp[2])/hmden)))/(dehnen_threeparam(radii[0], dehnen_threeparamfitp[0], dehnen_threeparamfitp[1],dehnen_threeparamfitp[2])/hmden))
dehnenthreep_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(densities[0])-1))
print ('ChiSquare and P values for dehnentwo', dehnenthreechi_square_test_statistic, dehnenthreep_value)
print ('Fitted value for Dehnen Three Parameters', dehnen_threeparamfitp)
print ('Uncertainties for Dehnen Three Parameters', np.sqrt(np.diag(dehnen_threeparamfitcov)))


with open('HaloFits/50-4_snap_99_halo_'+str(g-1)+'_fit_param.csv', 'w', encoding='UTF8', newline='') as f:
    
    header = ['Halo Number','NFW Scale Density','NFW Scale Radius','NFW Scale Density Uncertainty',
              'NFW Scale Radius Uncertainty','NFW ChiSquare','NFW P-Value','Burkert Scale Density','Burkert Scale Radius',
              'Burkert Scale Density Uncertainty','Burkert Scale Radius Uncertainty','Burkert ChiSquare','Burkert P-Value', 
              'Dehnen-2 Scale Density','Dehnen-2 Scale Radius','Dehnen-2 Scale Density Uncertainty',
              'Dehnen-2 Scale Radius Uncertainty','Dehnen-2 ChiSquare','Dehnen-2 P-Value','Einasto Scale Density',
              'Einasto Scale Radius','Einasto n', 'Einasto Scale Density Uncertainty',
              'Einasto Scale Radius Uncertainty','Einasto n Uncertainty','Einasto ChiSquare','Einasto P-Value',
              'Dehnen-3 Scale Density','Dehnen-3 Scale Radius','Dehnen-3 gamma', 'Dehnen-3 Scale Density Uncertainty',
              'Dehnen-3 Scale Radius Uncertainty','Dehnen-3 gamma Uncertainty','Dehnen-3 ChiSquare', 'Dehnen-3 P-Value']
    # Create a writer object
    fwriter = csv.writer(f, delimiter=',')
    # Write the header
    fwriter.writerow(header)
    data = [halonumber[0],nfwfitp[0],nfwfitp[1],np.sqrt(np.diag(nfwfitcov))[0],
              np.sqrt(np.diag(nfwfitcov))[1],nfwchi_square_test_statistic,nfwp_value,
              burkertfitp[0],burkertfitp[1],
              np.sqrt(np.diag(burkertfitcov))[0],np.sqrt(np.diag(burkertfitcov))[1],burkertchi_square_test_statistic, burkertp_value,
              dehnen_twoparamfitp[0],dehnen_twoparamfitp[1],np.sqrt(np.diag(dehnen_twoparamfitcov))[0],
              np.sqrt(np.diag(dehnen_twoparamfitcov))[1],dehnentwochi_square_test_statistic,dehnentwop_value,
              einastofitp[0],
              einastofitp[1],einastofitp[2], np.sqrt(np.diag(einastofitcov))[0],
              np.sqrt(np.diag(einastofitcov))[1],np.sqrt(np.diag(einastofitcov))[2],einastochi_square_test_statistic,einastop_value,
              dehnen_threeparamfitp[0],dehnen_threeparamfitp[1],dehnen_threeparamfitp[2], np.sqrt(np.diag(dehnen_threeparamfitcov))[0],
              np.sqrt(np.diag(dehnen_threeparamfitcov))[1],np.sqrt(np.diag(dehnen_threeparamfitcov))[2],dehnenthreechi_square_test_statistic,dehnenthreep_value]
    fwriter.writerow(data)
        



fig, axs = plt.subplots(3, 2, figsize=(15,15))


axs[0,0].errorbar((radii[0])*(h/hmrad), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(1)+"_099", color='green')



axs[0,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,0].legend()
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].set_title("Data from TNG")

axs[0,1].errorbar(radii[0]/hmrad, nfw(radii[0], nfwfitp[0], nfwfitp[1])/hmden, fmt='-', label="NFW fit Halo_"+str(1)+"_099", color='blue')
axs[0,1].errorbar((radii[0])*(h/hmrad), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[0,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,1].set_ylabel(r'($\rho$(r) ($M_{\odot} pc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,1].legend()
axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_title('NFW fit for Data')

axs[1,0].errorbar(radii[0]/hmrad, einasto(radii[0], einastofitp[0], einastofitp[1], einastofitp[2])/hmden, fmt='-', label="Einasto fit Halo_"+str(1)+"_099", color='blue')
axs[1,0].errorbar((radii[0])*(h/hmrad), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[1,0].set_xlabel(r'(Radius ($kpc/(h*R_{HalfMass}})}$))')
axs[1,0].set_ylabel(r'($\rho$(r) ($M_{\odot} pc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,0].legend()
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title('Einasto fit for Data')


axs[1,1].errorbar(radii[0]/hmrad, burkert(radii[0], burkertfitp[0], burkertfitp[1])/hmden, fmt='-', label="Bukert fit Halo_"+str(1)+"_099", color='blue')
axs[1,1].errorbar((radii[0]*(h/hmrad)), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[1,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[1,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,1].legend()
axs[1,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].set_title('Burkert fit for Data')

axs[2,0].errorbar(radii[0]/hmrad, dehnen_twoparam(radii[0], dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/hmden, fmt='-', label="Dehnen-2 fit Halo_"+str(1)+"_099", color='blue')
axs[2,0].errorbar((radii[0])*(h/hmrad), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[2,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[2,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[2,0].legend()
axs[2,0].set_yscale('log')
axs[2,0].set_xscale('log')
axs[2,0].set_title('Denhen-2 fit for Data')


axs[2,1].errorbar(radii[0]/hmrad, dehnen_threeparam(radii[0], dehnen_threeparamfitp[0], dehnen_threeparamfitp[1], dehnen_threeparamfitp[2])/hmden, fmt='-', label="Dehnen-3 fit Halo_"+str(1)+"_099", color='blue')
axs[2,1].errorbar((radii[0])*(h/hmrad), (densities[0])/(10*(h**2)*hmden), yerr=(uncertainties[0]), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[2,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[2,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[2,1].legend()
axs[2,1].set_yscale('log')
axs[2,1].set_xscale('log')
axs[2,1].set_title('Denhen-3 fit for Data')

fig.tight_layout()
fig.savefig('fit-profiles-halo-23')
fig.show()
