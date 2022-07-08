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

def virialRadius(radius, density):
    above_virR = np.where((density).astype(float)>float(p_crit*200))[0]
    virIndex = np.argmax(radius[above_virR])
    virR = radius[virIndex]
    print(virR)
    print(radius[-10:-1])
    return(virR,above_virR)
    

subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
length = subhalo_info['SubhaloLen'].to_numpy().astype(int)

gg = [0,63864,96762,117250,184931,198182,143880,208811,229933,220595,167392,253861,242788,264883]
numhalos = len(subhalo_index)
#densities = []
#uncertainties = []
#radii = []

pdheader = ['Radius','Density','Uncertainty']
#derek = pd.DataFrame(columns=pdheader)

with open('HaloFitsInfo/50-4_snap_99_fit_param.csv', 'x', encoding='UTF8', newline='') as f:
    
    header = ['Halo Number','DataPoints','NFW Scale Density','NFW Scale Radius','NFW Scale Density Uncertainty',
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

for g in gg:
    data_csv = pd.read_csv('HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den.csv')
    
    rad = data_csv['Radius'].to_numpy()
    den = data_csv['Density'].to_numpy()
    uncer = data_csv['Uncertainty'].to_numpy()
    num_datapoints = len(rad)
    
    virial_radius,virial_index = virialRadius(rad, den)
    virial_density = p_crit*200
    filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'param'
    
    rad = rad[virial_index]
    den = den[virial_index]
    uncer = uncer[virial_index]
    nfwfitp, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[1,10], sigma=uncer)
    nfwchi_square_test_statistic =  np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/(nfw(rad, nfwfitp[0], nfwfitp[1])))
    nfwp_value = scipy.stats.distributions.chi2.sf(nfwchi_square_test_statistic,(len(den)-1))
    print ('ChiSquare and P values for NFW', nfwchi_square_test_statistic)
    print ('Fitted value for NFW', nfwfitp)
    print ('uncer/(p_crit*200)tainties for NFW', np.sqrt(np.diag(nfwfitcov)))
    
    einastofitp, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=[1,10,5], sigma=uncer)
    einastochi_square_test_statistic =  np.sum((np.square(((den))-(einasto(rad, einastofitp[0], einastofitp[1],einastofitp[2]))))/(einasto(rad, einastofitp[0], einastofitp[1],einastofitp[2])))
    einastop_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(den)-1))
    print ('ChiSquare and P values for Einasto', einastochi_square_test_statistic, einastop_value)
    print ('Fitted value for Einasto', einastofitp)
    print ('uncertainties for Einasto', np.sqrt(np.diag(einastofitcov)))
    
    burkertfitp, burkertfitcov = scopt.curve_fit(burkert,rad, den, p0=[0.1,100], sigma=uncer)
    burkertchi_square_test_statistic =  np.sum((np.square(((den))-(burkert(rad, burkertfitp[0], burkertfitp[1]))))/(burkert(rad, burkertfitp[0], burkertfitp[1])))
    burkertp_value = scipy.stats.distributions.chi2.sf(burkertchi_square_test_statistic,(len(den)-1))
    print ('ChiSquare and P values for Burkert', burkertchi_square_test_statistic, burkertp_value)
    print ('Fitted value for Burkert', burkertfitp)
    print ('uncertainties for Burkert', np.sqrt(np.diag(burkertfitcov)))
    
    
    dehnen_twoparamfitp, dehnen_twoparamfitcov = scopt.curve_fit(dehnen_twoparam, rad, den, p0=[1,10], sigma=uncer)
    dehnentwochi_square_test_statistic =  np.sum((np.square(((den))-(dehnen_twoparam(rad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1]))))/(dehnen_twoparam(rad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])))
    dehnentwop_value = scipy.stats.distributions.chi2.sf(dehnentwochi_square_test_statistic,(len(den)-1))
    print ('ChiSquare and P values for dehnentwo', dehnentwochi_square_test_statistic, dehnentwop_value)
    print ('Fitted value for Dehnen Two Parameters', dehnen_twoparamfitp)
    print ('uncertainties for Dehnen Two Parameters', np.sqrt(np.diag(dehnen_twoparamfitcov)))
    
    dehnen_threeparamfitp, dehnen_threeparamfitcov = scopt.curve_fit(dehnen_threeparam, rad, den, p0=[0.1,20,0.02], sigma=uncer)
    dehnenthreechi_square_test_statistic =  np.sum((np.square(((den))-(dehnen_threeparam(rad,dehnen_threeparamfitp[0],dehnen_threeparamfitp[1],dehnen_threeparamfitp[2]))))/(dehnen_threeparam(rad, dehnen_threeparamfitp[0], dehnen_threeparamfitp[1],dehnen_threeparamfitp[2])))
    dehnenthreep_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(den)-1))
    print ('ChiSquare and P values for dehnenthree', dehnenthreechi_square_test_statistic, dehnenthreep_value)
    print ('Fitted value for Dehnen Three Parameters', dehnen_threeparamfitp)
    print ('uncertainties for Dehnen Three Parameters', np.sqrt(np.diag(dehnen_threeparamfitcov)))
    
    
    with open('HaloFitsInfo/50-4_snap_99_fit_param.csv', 'a', encoding='UTF8', newline='') as f:
        fwriter = csv.writer(f, delimiter=',')
        data = [subhalo_index[g],num_datapoints,nfwfitp[0],nfwfitp[1],np.sqrt(np.diag(nfwfitcov))[0],
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
    
    
    axs[0,0].errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    
    
    
    axs[0,0].set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs[0,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs[0,0].legend()
    axs[0,0].set_yscale('log')
    axs[0,0].set_xscale('log')
    axs[0,0].set_title("Data from TNG")
    
    axs[0,1].errorbar(rad/virial_radius, nfw(rad, nfwfitp[0], nfwfitp[1])/virial_density, fmt='-', label="NFW fit Halo_"+str(g)+"_099", color='blue')
    axs[0,1].errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    axs[0,1].set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs[0,1].set_ylabel(r'($\rho$(r) ($M_{\odot} pc^{-3} (\rho_{200})^{-1}$))')
    axs[0,1].legend()
    axs[0,1].set_yscale('log')
    axs[0,1].set_xscale('log')
    axs[0,1].set_title('NFW fit for Data')
    
    axs[1,0].errorbar(rad/virial_radius, einasto(rad, einastofitp[0], einastofitp[1], einastofitp[2])/virial_density, fmt='-', label="Einasto fit Halo_"+str(1)+"_099", color='blue')
    axs[1,0].errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    axs[1,0].set_xlabel(r'(Radius ($kpc/(h*R_{200}})}$))')
    axs[1,0].set_ylabel(r'($\rho$(r) ($M_{\odot} pc^{-3} (\rho_{200})^{-1}$))')
    axs[1,0].legend()
    axs[1,0].set_yscale('log')
    axs[1,0].set_xscale('log')
    axs[1,0].set_title('Einasto fit for Data')
    
    
    axs[1,1].errorbar(rad/virial_radius, burkert(rad, burkertfitp[0], burkertfitp[1])/virial_density, fmt='-', label="Bukert fit Halo_"+str(g)+"_099", color='blue')
    axs[1,1].errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    axs[1,1].set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs[1,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs[1,1].legend()
    axs[1,1].set_yscale('log')
    axs[1,1].set_xscale('log')
    axs[1,1].set_title('Burkert fit for Data')
    
    axs[2,0].errorbar(rad/virial_radius, dehnen_twoparam(rad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])/virial_density, fmt='-', label="Dehnen-2 fit Halo_"+str(g)+"_099", color='blue')
    axs[2,0].errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    axs[2,0].set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs[2,0].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs[2,0].legend()
    axs[2,0].set_yscale('log')
    axs[2,0].set_xscale('log')
    axs[2,0].set_title('denhen-2 fit for Data')
    
    
    axs[2,1].errorbar(rad/virial_radius, dehnen_threeparam(rad, dehnen_threeparamfitp[0], dehnen_threeparamfitp[1], dehnen_threeparamfitp[2])/virial_density, fmt='-', label="Dehnen-3 fit Halo_"+str(g)+"_099", color='blue')
    axs[2,1].errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='green')
    axs[2,1].set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs[2,1].set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs[2,1].legend()
    axs[2,1].set_yscale('log')
    axs[2,1].set_xscale('log')
    axs[2,1].set_title('denhen-3 fit for Data')
    
    fig.tight_layout()
    fig.savefig('radius-fit-profiles-halo'+str(g))
    print('hello')
    
    fig.clf()
    fig.show()
    g += 1
    #fig.clf()
    #fig.show()
