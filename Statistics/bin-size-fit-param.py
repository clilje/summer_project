# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

@author: clara
"""
import numpy as np
import matplotlib
import math
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats

plt.rc('font', size=20) #controls default text size
plt.rc('axes', titlesize=20) #fontsize of the title
plt.rc('axes', labelsize=20) #fontsize of the x and y labels
plt.rc('xtick', labelsize=20) #fontsize of the x tick labels
plt.rc('ytick', labelsize=20) #fontsize of the y tick labels
plt.rc('legend', fontsize=20) #fontsize of the legend


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
    mass = np.array(mass)
    
    
    bin_index = np.argsort(dis)
    radius_lowerbound = 0
    bin_lowerbound = 0
    upperbound =(bin_lowerbound+binsize)
    j = 0
    
    while j < 1:           
        if upperbound > len(dis):
            upperbound = len(dis)
            j = 1
        if bin_lowerbound == upperbound:
            break
        index_in_bin = bin_index[bin_lowerbound:upperbound]
        radius_upperbound = dis[index_in_bin][-1]
        dV = (4/3)*math.pi*(np.power(radius_upperbound,3)-np.power(radius_lowerbound,3))
        
        M = np.sum(mass[index_in_bin])
        subdensity = (M/(dV))
        density = np.append(density, subdensity)
        
        rad_lowerbound = np.append(rad_lowerbound, (radius_upperbound+radius_lowerbound)/2)
        dn = len(index_in_bin)
        uncertainties = np.append(uncertainties, subdensity/np.sqrt(dn))
        radius_lowerbound = radius_upperbound
        bin_lowerbound = upperbound
        upperbound = bin_lowerbound+binsize
    """
    above_virR = np.where((density).astype(float)>float(p_crit*200))[0]
    virR = np.max(rad_lowerbound[above_virR])
    print(virR)
    print(rad_lowerbound[-10:-1])
    """
    return(density, rad_lowerbound, uncertainties)


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
    #print(virR)
    #print(radius[-10:-1])
    return(virR,above_virR)
    
def chiSquareNFW(rad,den, uncertainties, nfwfitp):
    return(np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/np.square(uncertainties)))

def getChiSquareplot(rad,den,uncertainties,nfwfitopt):
    results = np.zeros((len(nfwfitopt[0]),len(nfwfitopt[1])))
    #print(nfwfitopt)
    for x in np.arange(len(nfwfitopt[0])):
        for y in np.arange(len(nfwfitopt[1])):
            #print(x,y)
            results[x,y] = chiSquareNFW(rad,den,uncertainties,[nfwfitopt[0][x,y],nfwfitopt[1][x,y]])
    return(results)

def plotting(rad, den,uncertainties, virial_radius, virial_density,nfwfitp,burkertfitp,dehnen_threeparamfitp,dehnen_twoparamfitp,einastofitp):
    #fig, axs = plt.subplots(3, 2, figsize=(15,15))
    #fig = plt.figure(figsize=(30,20))
    fig, axs = plt.subplots(1, 2, figsize=(20,20))
    #gs0 = fig.add_gridspec(3, 3)
   # gs0 = fig.add_gridspec(1, 2)
    
    """
    gs = fig.add_gridspec(3, 3)
    ax1 = fig.add_subplot(gs[0, :])
    ax1.set_title('gs[0, :]')
    ax2 = fig.add_subplot(gs[1, :-1])
    ax2.set_title('gs[1, :-1]')
    """
    #axs2 = fig.add_subplot(gs0[0,1])
    axs2 = axs[1]
    print('1')
    
    
    x, y = np.linspace(100000, 1000000000, 1000), np.linspace(0, 5, 1000)
    X, Y = np.meshgrid(x, y)
    Z = getChiSquareplot(rad, den, uncertainties, [X,Y])
    #print(X)
    #print(Y)
    #print(Z)
    normalizeC = np.min(Z)
    #print(normalizeC)
    Z = Z/normalizeC
    #print(Z)
    #Z[np.where(Z==np.inf)[0]] = 10**10
    #print(Z)
    pc = axs2.pcolor(X,Y,Z,norm=matplotlib.colors.LogNorm(vmin=1, vmax=15000),cmap='PuBu_r')
    #pc = axs2.pcolor(Z)
    
    #print('2')
    nfwfitp0, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
    axs[1].scatter(nfwfitp0[0],nfwfitp0[1], color='red', label="Optimised Parameters initial parameters: Virrad and Virdensity")
    axs[1].scatter(den[-1],rad[-1], color='red',marker='x', label="Initial parameters: Virrad and Virdensity")
    #print('3')
    init1 = [100000,3]
    nfwfitp1, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=init1, sigma=uncer)
    axs[1].scatter(nfwfitp1[0],nfwfitp1[1], color='pink', label="Optimised Parameters initial parameters: 100000,3")
    axs[1].scatter(100000,3, color='pink',marker='x', label="Initial parameters: 100000,3")
    #print('4')
    init2 = [200000000,0.5]
    nfwfitp2, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[200000000,0.5], sigma=uncer)
    axs[1].scatter(nfwfitp2[0],nfwfitp2[1], color='cyan', label="Optimised Parameters initial parameters:200000000,0.5")
    axs[1].scatter(200000000,0.5, color='cyan',marker='x', label="Initial parameters: 200000000,0.5")
    #print('5')
    init3 = [800000000,0.3]
    nfwfitp3, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[800000000,0.3], sigma=uncer)
    axs[1].scatter(nfwfitp3[0],nfwfitp3[1], color='green', label="Optimised Parameters initial parameters:800000000,0.3")
    axs[1].scatter(800000000,0.3, color='green',marker='x', label="Initial parameters: 800000000,0.3")
    #print('6')
    #init4 = [3000000000,1]
    #nfwfitp4, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[3000000000,1], sigma=uncer)
    #axs[1].scatter(nfwfitp4[0],nfwfitp4[1], color='black', label="Optimised Parameters initial parameters: 3000000000,1")
    #axs[1].scatter(3000000000,1, color='black',marker='x', label="Initial parameters: 3000000000,1")
    
    
    
    
    axs[1].set_xlabel("Scale Density")
    axs[1].set_ylabel("Scale Radius")
    axs[1].set_title("Heatmap of ChiSquare values for NFW profile parameters (minimum ChiSquare set to 1)")
    axs[1].legend()
    fig.colorbar(pc, ax=axs2, extend='max')
    #print(nfwfitp0, nfwfitp1, nfwfitp2,nfwfitp3,nfwfitp4)
    
    print(rad,den)
    axs1 = axs[0]
    axs1.errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='blue')
    axs1.scatter(rad/virial_radius, nfw(rad, nfwfitp0[0], nfwfitp0[1])/virial_density, marker='x', label="NFW fit optimised param"+str(nfwfitp0), color='red')
    axs1.scatter(rad/virial_radius, nfw(rad, init1[0], init1[1])/virial_density,marker='x', label="NFW fit initial param"+str(init1), color='pink')
    axs1.scatter(rad/virial_radius, nfw(rad, init2[0], init2[1])/virial_density, marker='x', label="NFW fit initial param"+str(init2), color='cyan')
    axs1.scatter(rad/virial_radius, nfw(rad, init3[0], init3[1])/virial_density, marker='x', label="NFW fit initial param"+str(init3), color='green')
    #axs1.scatter(rad/virial_radius, nfw(rad, init4[0], init4[1])/virial_density, marker='x', label="NFW fit initial param"+str(init4), color='black')
    
    #print(chiSquareNFW(rad, den, uncertainties, nfwfitp1))
    #print(chiSquareNFW(rad, den, uncertainties, init1))
    #print(chiSquareNFW(rad, den, uncertainties, init2))
    #print(chiSquareNFW(rad, den, uncertainties, init3))
    #print(chiSquareNFW(rad, den, uncertainties, init4))
    #print(min([chiSquareNFW(rad,uncertainties, den, nfwfitp1),chiSquareNFW(rad, den,uncertainties, init1),chiSquareNFW(rad, den,uncertainties, init2),chiSquareNFW(rad, den,uncertainties, init3),chiSquareNFW(rad, den,uncertainties, init4)]))
    
    axs1.set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs1.set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs1.legend()
    axs1.set_yscale('log')
    axs1.set_xscale('log')
    axs1.set_title("Fits to Data from TNG")
    
    fig.tight_layout()
    fig.savefig('HaloFitsInfo/zoom-fit-profiles-halo'+str(g))
    print('hello')
    
    fig.clf()
    #fig.show()



subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_index = subhalo_info['SubhaloIndex']
positionsX = subhalo_info['SubhaloPosX'].to_numpy()
positionsY = subhalo_info['SubhaloPosY'].to_numpy()
positionsZ = subhalo_info['SubhaloPosZ'].to_numpy()
radius = subhalo_info['SubhaloHalfmassRad'].to_numpy()
full_mass = subhalo_info['SubhaloMass'].to_numpy()
#length = subhalo_info['SubhaloLen'].to_numpy().astype(int)

#gg = [0,4,20,200,2000,198182,19999]
gg =[20,50,200,300,2000,4000,4]
numhalos = len(subhalo_index)
#densities = []
#uncertainties = []
#radii = []

pdheader = ['Radius','Density','Uncertainty']
#derek = pd.DataFrame(columns=pdheader)

for g in gg:
    
    chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['x'],dtype={'x':object})
    chunk = chunk[chunk['x'] != 'x']
    #print('success')
    partx = chunk['x'].to_numpy().astype(float)
    #print(partx)
    #exit()
    chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['y'],dtype={'y':object})
    chunk = chunk[chunk['y'] != 'y']
    party = chunk['y'].to_numpy().astype(float)
    chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['z'],dtype={'z':object})
    chunk = chunk[chunk['z'] != 'z']
    partz = chunk['z'].to_numpy().astype(float)
    chunk = pd.read_csv('FullRun/snap_99_halo_'+str(g)+'.csv', usecols=['mass'],dtype={'mass':object})
    chunk = chunk[chunk['mass'] != 'mass']
    mass = chunk['mass'].to_numpy().astype(float)
    #filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den'
    """
    if len(partx) < 600:
        binsize = 3
    else:
        binsize = int(len(partx)/200)
    """
    if len(mass)/350 > 5:
        binsizes = np.linspace((len(mass)/350),(len(mass)/150),15).astype(int)
        #binsizes = np.linspace((len(mass)/250),(len(mass)/7),15).astype(int)
    else:
        #bin_num = len(mass)/3
        binsizes = np.linspace(3,len(mass)/15,15).astype(int)
        #binsizes = np.linspace(3,len(mass)/7,15).astype(int)
        
    scale_den = []
    scale_rad = []
    chisquare = []
    uncer_den = []
    uncer_rad = []
    for binsize in binsizes:
        rad_den = radial_density((partx*h), (party*h), (partz*h),(mass*h*(10**10)), binsize, (positionsX[g]*h), (h*positionsY[g]), (h*positionsZ[g]))
        #mass in solar masses
        #distances in kpc
        #virrad = rad_den[3]
        virden = 200*p_crit
        #miniderek = pd.DataFrame(columns=pdheader)
        rad=rad_den[1]
        #miniderek['Virial Radius']=[virrad]*len(rad_den[0])
        den=rad_den[0]
        uncer=rad_den[2]
        #print(miniderek)
        #miniderek.to_csv(filename+'.csv')
        #miniderek = miniderek[0:0]
        #y += 1
        #print('chunk'+str(y))
        
        
        #data_csv = pd.read_csv('HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den.csv')
        
        #rad = data_csv['Radius'].to_numpy()
        #den = data_csv['Density'].to_numpy()
        #uncer = data_csv['Uncertainty'].to_numpy()
        num_datapoints = len(rad)
        
        virial_radius,virial_index = virialRadius(rad, den)
        virial_density = p_crit*200
        #print(str((virial_radius/radius)[g]))
        #if (virial_radius/radius)[g] < 10:
        
        #filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'param'
        
        rad = rad[virial_index]
        den = den[virial_index]
        uncer = uncer[virial_index]
        nfwfitp, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[virial_density,virial_radius], sigma=uncer)
        nfwchi_square_test_statistic =  np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/np.square(uncer))
        nfwp_value = scipy.stats.distributions.chi2.sf(nfwchi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for NFW', nfwchi_square_test_statistic)
        print ('Fitted value for NFW', nfwfitp)
        print ('uncertainties for NFW', np.sqrt(np.diag(nfwfitcov)))
        scale_den.append(nfwfitp[0])
        scale_rad.append(nfwfitp[1])
        chisquare.append(nfwchi_square_test_statistic)
        uncer_den.append(np.sqrt(np.diag(nfwfitcov))[0])
        uncer_rad.append(np.sqrt(np.diag(nfwfitcov))[1])
    burkertfitp = [0,0]
    dehnen_threeparamfitp = [0,0]
    dehnen_twoparamfitp = [0,0]
    einastofitp = [0,0]
    print(binsizes)
    print(scale_den)
    print(scale_rad)
    print(chisquare)
    fig, axs = plt.subplots(1, 2, figsize=(30,20))
    cm = plt.cm.get_cmap('RdYlBu')
    #xy = range(20)
    z = chisquare
    sc = axs[0].errorbar(binsizes, scale_den, c=z,yerr=uncer_den, s=300,marker='*')
    
    axs[0].set_xlabel(r'Binsize')
    axs[0].set_ylabel(r'Scale Density')
    #axs1.legend()
    #axs1.set_yscale('log')
    #axs1.set_xscale('log')
    axs[0].set_title("Scale Density at different bin sizes")
    
    sc2 = axs[1].errorbar(binsizes, scale_rad, c=z,yerr=uncer_rad, s=300,marker='*',  cmap=cm)
    
    #axs[1].colorbar(sc)
    axs[1].set_xlabel(r'Binsize')
    axs[1].set_ylabel(r'Scale Radius')
    #axs1.legend()
    #axs1.set_yscale('log')
    #axs1.set_xscale('log')
    axs[1].set_title("Scale Radius at different bin sizes")
    
    fig.tight_layout()
    #fig.savefig('HaloFitsInfo/bins-fit-profiles-halo'+str(g))
    print('hello')
    fig.colorbar(sc, ax=axs[0])
    fig.colorbar(sc2, ax=axs[1])
    fig.savefig('HaloFitsInfo/bins-smallrange-fit-profiles-halo'+str(g))
    #fig.show()
    #plotting(rad, den, uncer, virial_radius, virial_density,nfwfitp,burkertfitp,dehnen_threeparamfitp,dehnen_twoparamfitp,einastofitp)
    
    g += 1
        #fig.clf()
        #fig.show()
