# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

@author: clara
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)

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


def radial_density(partx, party, partz, mass, binsize, halox, haloy, haloz):
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
    
    #virV = (4/3)*math.pi*(np.power((virrad+10),3)-np.power((virrad-10),3))
    #virindex = np.where(np.logical_and(dis.astype(float)>float(virrad-10), dis.astype(float)<float(virrad+10)))[0]
    mass = np.array(mass)
    #virM = np.sum(mass[virindex])
    #virdensity = virM/virV
    
    
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
    
    below_virR = np.where((dis[bin_index]).astype(float)<float(p_crit*200))[0]
    virR = dis[bin_index][below_virR][-1]
    return(density, rad_lowerbound, uncertainties, virR)


interval = np.logspace(0.1, 2.5, 100)
#files = get_filenames(50, 4, 11)
subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
#full_mass = subhalo_info['SubhaloMass'].to_numpy()
#print(pd.read_csv('HaloFits/50-4_halodata.csv'))
#print(positions[2])
#print(radius[2])
#halonumber = []
g = 4000
numhalos = len(subhalo_index)
#densities = []
#uncertainties = []
#radii = []

pdheader = ['Radius','Density','Uncertainty','Virial Radius']
#derek = pd.DataFrame(columns=pdheader)

while g < 4003:
    data_csv = pd.read_csv('HaloParticles50-1-pd/snap_99_halo_'+str(g)+'.csv', dtype={'':int,'ID':int,'Type':'string','x':float,'y':float,'z':float,'mass':float,'vx':float,'vy':float,'vz':float})
    print(data_csv)
    filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den'
    rad_den = radial_density((data_csv['x'].to_numpy()*h), (data_csv['y'].to_numpy()*h), (data_csv['z'].to_numpy()*h),(data_csv['mass'].to_numpy()*h*(10**10)), 10, (positionsX[g]*h), (h*positionsY[g]), (h*positionsZ[g]))
    #mass in solar masses
    #distances in kpc
    print(rad_den)
    #hmrad = radius[g]
    virrad = rad_den[3]
    virden = 200*p_crit
    miniderek = pd.DataFrame(columns=pdheader)
    miniderek['Radius']=rad_den[1]/(virrad)
    miniderek['Virial Radius']=[virrad]*len(rad_den[0])
    miniderek['Density']=rad_den[0]/(virden)
    miniderek['Uncertainty']=rad_den[2]
    #densities.append(list(rad_den[0]))
    #radii.append(list(rad_den[1]))
    #uncertainties.append(list(rad_den[2]))
    #hmden = rad_den[3]
    #halonumber.append(g)
    #print(hmrad,hmden)
    #derek = pd.concat([derek,miniderek])
    miniderek.to_csv(filename+'.csv', mode='a')
    miniderek = miniderek[0:0]
    g += 1 
    plt.errorbar(rad_den[1]/(virrad), rad_den[0]/(virden), yerr=rad_den[2], fmt='.', label="Halo_"+str(g)+"_099", color='green')
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('fit-profiles-halo-'+str(g)+'.jpg')

#densities = np.array(densities)
#radii = np.array(radii)
#half_rad_index = int(len(radii[0])/2)
#uncertainties = np.array(uncertainties)
#uncertainties[uncertainties == np.nan] = 0
#indices = np.arange(half_rad_index,int(len(radii[0])))
#print(len(radii[0]))
#print(indices)
#shortrad = np.delete(radii[0],indices)
#shortden = np.delete(densities[0],indices)
#shortuncer = np.delete(uncertainties[0],indices)
#print(len(radii[0]))
"""
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

nfwfitp, nfwfitcov = scopt.curve_fit(nfw, shortrad*h, shortden/(10*(h**2)), p0=[0.001,100], sigma=shortuncer)
nfwchi_square_test_statistic =  np.sum((np.square(((shortden)/(10*(h**2)*hmden))-(nfw(shortrad, nfwfitp[0], nfwfitp[1])/hmden)))/(nfw(shortrad, nfwfitp[0], nfwfitp[1])/hmden))
nfwp_value = scipy.stats.distributions.chi2.sf(nfwchi_square_test_statistic,(len(shortden)-1))
print ('ChiSquare and P values for NFW', nfwchi_square_test_statistic)
print ('Fitted value for NFW', nfwfitp)
print ('Uncertainties for NFW', np.sqrt(np.diag(nfwfitcov)))

einastofitp, einastofitcov = scopt.curve_fit(einasto, shortrad*h, shortden/(10*(h**2)), p0=[0.0001,1000,4], sigma=shortuncer)
einastochi_square_test_statistic =  np.sum((np.square(((shortden)/(10*(h**2)*hmden))-(einasto(shortrad, einastofitp[0], einastofitp[1],einastofitp[2])/hmden)))/(einasto(shortrad, einastofitp[0], einastofitp[1],einastofitp[2])/hmden))
einastop_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(shortden)-1))
print ('ChiSquare and P values for Einasto', einastochi_square_test_statistic, einastop_value)
print ('Fitted value for Einasto', einastofitp)
print ('Uncertainties for Einasto', np.sqrt(np.diag(einastofitcov)))

burkertfitp, burkertfitcov = scopt.curve_fit(burkert,shortrad*h, shortden/(10*(h**2)), p0=[0.1,10], sigma=shortuncer)
burkertchi_square_test_statistic =  np.sum((np.square(((shortden)/(10*(h**2)*hmden))-(burkert(shortrad, burkertfitp[0], burkertfitp[1])/hmden)))/(burkert(shortrad, burkertfitp[0], burkertfitp[1])/hmden))
burkertp_value = scipy.stats.distributions.chi2.sf(burkertchi_square_test_statistic,(len(shortden)-1))
print ('ChiSquare and P values for Burkert', burkertchi_square_test_statistic, burkertp_value)
print ('Fitted value for Burkert', burkertfitp)
print ('Uncertainties for Burkert', np.sqrt(np.diag(burkertfitcov)))


dehnen_twoparamfitp, dehnen_twoparamfitcov = scopt.curve_fit(dehnen_twoparam, shortrad*h, shortden/(10*(h**2)), p0=[0.01,300], sigma=shortuncer)
dehnentwochi_square_test_statistic =  np.sum((np.square(((shortden)/(10*(h**2)*hmden))-(dehnen_twoparam(shortrad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/hmden)))/(dehnen_twoparam(shortrad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/hmden))
dehnentwop_value = scipy.stats.distributions.chi2.sf(dehnentwochi_square_test_statistic,(len(shortden)-1))
print ('ChiSquare and P values for dehnentwo', dehnentwochi_square_test_statistic, dehnentwop_value)
print ('Fitted value for Dehnen Two Parameters', dehnen_twoparamfitp)
print ('Uncertainties for Dehnen Two Parameters', np.sqrt(np.diag(dehnen_twoparamfitcov)))

dehnen_threeparamfitp, dehnen_threeparamfitcov = scopt.curve_fit(dehnen_threeparam, shortrad*h, shortden/(10*(h**2)), p0=[0.01,250,0.02], sigma=shortuncer)
dehnenthreechi_square_test_statistic =  np.sum((np.square(((shortden)/(10*(h**2)*hmden))-(dehnen_threeparam(shortrad,dehnen_threeparamfitp[0],dehnen_threeparamfitp[1],dehnen_threeparamfitp[2])/hmden)))/(dehnen_threeparam(shortrad, dehnen_threeparamfitp[0], dehnen_threeparamfitp[1],dehnen_threeparamfitp[2])/hmden))
dehnenthreep_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(shortden)-1))
print ('ChiSquare and P values for dehnentwo', dehnenthreechi_square_test_statistic, dehnenthreep_value)
print ('Fitted value for Dehnen Three Parameters', dehnen_threeparamfitp)
print ('Uncertainties for Dehnen Three Parameters', np.sqrt(np.diag(dehnen_threeparamfitcov)))


with open('HaloFits/50-4_snap_99_halo_'+str(g)+'_fit_param_shortrad.csv', 'w', encoding='UTF8', newline='') as f:
    
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
   """     


"""
fig, axs = plt.subplots(3, 2, figsize=(15,15))


axs[0,0].errorbar((shortrad)*(h/hmrad), (shortden)/(10*(h**2)*hmden), yerr=(shortuncer), fmt='.', label="Halo_"+str(1)+"_099", color='green')



axs[0,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,0].legend()
axs[0,0].set_yscale('log')
axs[0,0].set_xscale('log')
axs[0,0].set_title("Data from TNG")

axs[0,1].errorbar(shortrad/hmrad, nfw(shortrad, nfwfitp[0], nfwfitp[1])/hmden, fmt='-', label="NFW fit Halo_"+str(1)+"_099", color='blue')
axs[0,1].errorbar((shortrad)*(h/hmrad), (shortden)/(10*(h**2)*hmden), yerr=(shortuncer), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[0,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[0,1].set_ylabel(r'($\rho$(r) ($M_{\odot} pc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[0,1].legend()
axs[0,1].set_yscale('log')
axs[0,1].set_xscale('log')
axs[0,1].set_title('NFW fit for Data')

axs[1,0].errorbar(shortrad/hmrad, einasto(shortrad, einastofitp[0], einastofitp[1], einastofitp[2])/hmden, fmt='-', label="Einasto fit Halo_"+str(1)+"_099", color='blue')
axs[1,0].errorbar((shortrad)*(h/hmrad), (shortden)/(10*(h**2)*hmden), yerr=(shortuncer), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[1,0].set_xlabel(r'(Radius ($kpc/(h*R_{HalfMass}})}$))')
axs[1,0].set_ylabel(r'($\rho$(r) ($M_{\odot} pc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,0].legend()
axs[1,0].set_yscale('log')
axs[1,0].set_xscale('log')
axs[1,0].set_title('Einasto fit for Data')


axs[1,1].errorbar(shortrad/hmrad, burkert(shortrad, burkertfitp[0], burkertfitp[1])/hmden, fmt='-', label="Bukert fit Halo_"+str(1)+"_099", color='blue')
axs[1,1].errorbar((shortrad*(h/hmrad)), (shortden)/(10*(h**2)*hmden), yerr=(shortuncer), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[1,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[1,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[1,1].legend()
axs[1,1].set_yscale('log')
axs[1,1].set_xscale('log')
axs[1,1].set_title('Burkert fit for Data')

axs[2,0].errorbar(shortrad/hmrad, dehnen_twoparam(shortrad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/hmden, fmt='-', label="Dehnen-2 fit Halo_"+str(1)+"_099", color='blue')
axs[2,0].errorbar((shortrad)*(h/hmrad), (shortden)/(10*(h**2)*hmden), yerr=(shortuncer), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[2,0].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[2,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[2,0].legend()
axs[2,0].set_yscale('log')
axs[2,0].set_xscale('log')
axs[2,0].set_title('Denhen-2 fit for Data')


axs[2,1].errorbar(shortrad/hmrad, dehnen_threeparam(shortrad, dehnen_threeparamfitp[0], dehnen_threeparamfitp[1], dehnen_threeparamfitp[2])/hmden, fmt='-', label="Dehnen-3 fit Halo_"+str(1)+"_099", color='blue')
axs[2,1].errorbar((shortrad)*(h/hmrad), (shortden)/(10*(h**2)*hmden), yerr=(shortuncer), fmt='.', label="Halo_"+str(1)+"_099", color='green')
axs[2,1].set_xlabel(r'(Radius ($kpc/(R_{HalfMass}})}$))')
axs[2,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{HalfMass})^{-1}$))')
axs[2,1].legend()
axs[2,1].set_yscale('log')
axs[2,1].set_xscale('log')
axs[2,1].set_title('Denhen-3 fit for Data')

fig.tight_layout()
fig.savefig('shortrad-fit-profiles-halo-3')
print('hello')
fig.show()
"""