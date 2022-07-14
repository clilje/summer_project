# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

This code is to check whether NFW fits depend sensitively on initial parameters

@author: clara
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)

"""
These are all possible profiles that can be tested.
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



def virialRadius(radius, density):
    """
    

    Parameters
    ----------
    radius : radius of the radial shells at which density is calculated
    density : density for each radial shell

    Returns
    -------
    radius at which halo has 200 times the critical density of the univers

    """
    above_virR = np.where((density).astype(float)>float(p_crit*200))[0]
    virIndex = np.argmax(radius[above_virR])
    virR = radius[virIndex]
    return(virR,above_virR)
    
def chiSquareNFW(rad,den, uncertainties, nfwfitp):
    """

    Parameters
    ----------
    radius : radius of the radial shells at which density is calculated
    density : density for each radial shell
    uncertainties : poisson errors
    nfwfitp : optimised parameters
    Returns
    -------
    weighted chisquare for the most optimal fit compared to the given data for nfw profile

    """
    return(np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/np.square(uncertainties)))

def getChiSquareplot(rad,den,uncertainties,nfwfitopt):
    """


    Parameters
    ----------
    radius : radius of the radial shells at which density is calculated
    density : density for each radial shell
    uncertainties : poisson errors
    nfwfitp : optimised parameters

    Returns
    -------
    A NxN grid of the chisquare values of different parameters for the NFW profile

    """
    results = np.zeros((len(nfwfitopt[0]),len(nfwfitopt[1])))
    #print(nfwfitopt)
    for x in np.arange(len(nfwfitopt[0])):
        for y in np.arange(len(nfwfitopt[1])):
            #print(x,y)
            results[x,y] = chiSquareNFW(rad,den,uncertainties,[nfwfitopt[0][x,y],nfwfitopt[1][x,y]])
    return(results)


def plotting(rad, den,uncertainties, virial_radius, virial_density,nfwfitp,burkertfitp,dehnen_threeparamfitp,dehnen_twoparamfitp,einastofitp):
    """
    In case this is required this plots a heat map of Chi Square value with the most optimised value indicated.
    It is possible to add the location of different starting values and the fits they produce.
    This function needs to be hand-tuned to each halo it is operated on to find optimal framing of the plot as well as starting
    parameters.
    
    """
    
    #create axis with sublots
    fig, axs = plt.subplots(1, 2, figsize=(20,20))
    axs2 = axs[1]
    
    #plot the heatmap on the right
    x, y = np.linspace(100000, 1000000000, 1000), np.linspace(0, 5, 1000)
    X, Y = np.meshgrid(x, y)
    Z = getChiSquareplot(rad, den, uncertainties, [X,Y])
    
    normalizeC = np.min(Z)
    Z = Z/normalizeC
    
    pc = axs2.pcolor(X,Y,Z,norm=matplotlib.colors.LogNorm(vmin=1, vmax=15000),cmap='PuBu_r')
    
    #add different initial parameters

    nfwfitp0, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
    axs[1].scatter(nfwfitp0[0],nfwfitp0[1], color='red', label="Optimised Parameters initial parameters: Virrad and Virdensity")
    axs[1].scatter(den[-1],rad[-1], color='red',marker='x', label="Initial parameters: Virrad and Virdensity")
    
    
    init1 = [100000,3]
    nfwfitp1, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=init1, sigma=uncer)
    axs[1].scatter(nfwfitp1[0],nfwfitp1[1], color='pink', label="Optimised Parameters initial parameters: 100000,3")
    axs[1].scatter(100000,3, color='pink',marker='x', label="Initial parameters: 100000,3")
    
    
    init2 = [200000000,0.5]
    nfwfitp2, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[200000000,0.5], sigma=uncer)
    axs[1].scatter(nfwfitp2[0],nfwfitp2[1], color='cyan', label="Optimised Parameters initial parameters:200000000,0.5")
    axs[1].scatter(200000000,0.5, color='cyan',marker='x', label="Initial parameters: 200000000,0.5")
    
    
    init3 = [800000000,0.3]
    nfwfitp3, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[800000000,0.3], sigma=uncer)
    axs[1].scatter(nfwfitp3[0],nfwfitp3[1], color='green', label="Optimised Parameters initial parameters:800000000,0.3")
    axs[1].scatter(800000000,0.3, color='green',marker='x', label="Initial parameters: 800000000,0.3")
    
    
    
    #init4 = [3000000000,1]
    #nfwfitp4, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[3000000000,1], sigma=uncer)
    #axs[1].scatter(nfwfitp4[0],nfwfitp4[1], color='black', label="Optimised Parameters initial parameters: 3000000000,1")
    #axs[1].scatter(3000000000,1, color='black',marker='x', label="Initial parameters: 3000000000,1")
    
    
    
    
    axs[1].set_xlabel("Scale Density")
    axs[1].set_ylabel("Scale Radius")
    axs[1].set_title("Heatmap of ChiSquare values for NFW profile parameters (minimum ChiSquare set to 1)")
    axs[1].legend()
    fig.colorbar(pc, ax=axs2, extend='max')
    
    #Plot fits to original data on the left
    axs1 = axs[0]
    axs1.errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='blue')
    axs1.scatter(rad/virial_radius, nfw(rad, nfwfitp0[0], nfwfitp0[1])/virial_density, marker='x', label="NFW fit optimised param"+str(nfwfitp0), color='red')
    #add fits of different initial parameters
    axs1.scatter(rad/virial_radius, nfw(rad, init1[0], init1[1])/virial_density,marker='x', label="NFW fit initial param"+str(init1), color='pink')
    axs1.scatter(rad/virial_radius, nfw(rad, init2[0], init2[1])/virial_density, marker='x', label="NFW fit initial param"+str(init2), color='cyan')
    axs1.scatter(rad/virial_radius, nfw(rad, init3[0], init3[1])/virial_density, marker='x', label="NFW fit initial param"+str(init3), color='green')
    #axs1.scatter(rad/virial_radius, nfw(rad, init4[0], init4[1])/virial_density, marker='x', label="NFW fit initial param"+str(init4), color='black')
    
    axs1.set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs1.set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs1.legend()
    axs1.set_yscale('log')
    axs1.set_xscale('log')
    axs1.set_title("Fits to Data from TNG")
    
    fig.tight_layout()
    fig.savefig('HaloFitsInfo/zoom-fit-profiles-halo'+str(g))

    
    fig.clf()



subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()

#list of halos to check for
gg =[2000]
numhalos = len(subhalo_index)

pdheader = ['Radius','Density','Uncertainty']

#iterate over each halo
for g in gg:
    
    #Read radial density information in 
    data_csv = pd.read_csv('HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den.csv')
    
    rad = data_csv['Radius'].to_numpy()
    den = data_csv['Density'].to_numpy()
    uncer = data_csv['Uncertainty'].to_numpy()
    num_datapoints = len(rad)
    
    #find virial radius
    virial_radius,virial_index = virialRadius(rad, den)
    virial_density = p_crit*200
    
    
    #criterion for disruption
    if (virial_radius/radius)[g] < 10:
        
        filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'param'
        
        #find optimised fit up to virial radius for nFW
        rad = rad[virial_index]
        den = den[virial_index]
        uncer = uncer[virial_index]
        nfwfitp, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
        nfwchi_square_test_statistic =  np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/(nfw(rad, nfwfitp[0], nfwfitp[1])))
        nfwp_value = scipy.stats.distributions.chi2.sf(nfwchi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for NFW', nfwchi_square_test_statistic)
        print ('Fitted value for NFW', nfwfitp)
        print ('uncertainties for NFW', np.sqrt(np.diag(nfwfitcov)))
        
        #This section can be uncommented if other fits are required
        """
        einastofitp, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=[den[-1],rad[-1],5], sigma=uncer)
        einastochi_square_test_statistic =  np.sum((np.square(((den))-(einasto(rad, einastofitp[0], einastofitp[1],einastofitp[2]))))/(einasto(rad, einastofitp[0], einastofitp[1],einastofitp[2])))
        einastop_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for Einasto', einastochi_square_test_statistic, einastop_value)
        print ('Fitted value for Einasto', einastofitp)
        print ('uncertainties for Einasto', np.sqrt(np.diag(einastofitcov)))
        
        burkertfitp, burkertfitcov = scopt.curve_fit(burkert,rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
        burkertchi_square_test_statistic =  np.sum((np.square(((den))-(burkert(rad, burkertfitp[0], burkertfitp[1]))))/(burkert(rad, burkertfitp[0], burkertfitp[1])))
        burkertp_value = scipy.stats.distributions.chi2.sf(burkertchi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for Burkert', burkertchi_square_test_statistic, burkertp_value)
        print ('Fitted value for Burkert', burkertfitp)
        print ('uncertainties for Burkert', np.sqrt(np.diag(burkertfitcov)))
        
        
        dehnen_twoparamfitp, dehnen_twoparamfitcov = scopt.curve_fit(dehnen_twoparam, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
        dehnentwochi_square_test_statistic =  np.sum((np.square(((den))-(dehnen_twoparam(rad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1]))))/(dehnen_twoparam(rad, dehnen_twoparamfitp[0], dehnen_twoparamfitp[1])))
        dehnentwop_value = scipy.stats.distributions.chi2.sf(dehnentwochi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for dehnentwo', dehnentwochi_square_test_statistic, dehnentwop_value)
        print ('Fitted value for Dehnen Two Parameters', dehnen_twoparamfitp)
        print ('uncertainties for Dehnen Two Parameters', np.sqrt(np.diag(dehnen_twoparamfitcov)))
        
        dehnen_threeparamfitp, dehnen_threeparamfitcov = scopt.curve_fit(dehnen_threeparam, rad, den, p0=[den[-1],rad[-1],0.02], sigma=uncer)
        dehnenthreechi_square_test_statistic =  np.sum((np.square(((den))-(dehnen_threeparam(rad,dehnen_threeparamfitp[0],dehnen_threeparamfitp[1],dehnen_threeparamfitp[2]))))/(dehnen_threeparam(rad, dehnen_threeparamfitp[0], dehnen_threeparamfitp[1],dehnen_threeparamfitp[2])))
        dehnenthreep_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for dehnenthree', dehnenthreechi_square_test_statistic, dehnenthreep_value)
        print ('Fitted value for Dehnen Three Parameters', dehnen_threeparamfitp)
        print ('uncertainties for Dehnen Three Parameters', np.sqrt(np.diag(dehnen_threeparamfitcov)))
        
        
        """
        
        #unnescessary parameters for rad den
        burkertfitp = [0,0]
        dehnen_threeparamfitp = [0,0]
        dehnen_twoparamfitp = [0,0]
        einastofitp = [0,0]
        
        #plot
        plotting(rad, den, uncer, virial_radius, virial_density,nfwfitp,burkertfitp,dehnen_threeparamfitp,dehnen_twoparamfitp,einastofitp)
        
    g += 1
