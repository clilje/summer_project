# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 14:38:14 2022

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
n = 1/0.16


def nfw(r, density_0, scale_radius):
    return(density_0/((r/scale_radius)*np.power((1+(r/scale_radius)),2)))

def einasto(r, density_e, r_e):
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
    
def chiSquareeinasto(rad,den, uncertainties, einastofitp):
    return(np.sum((np.square(((den))-(einasto(rad, einastofitp[0], einastofitp[1]))))/np.square(uncertainties)))

def getChiSquareplot(rad,den,uncertainties,einastofitopt):
    results = np.zeros((len(einastofitopt[0]),len(einastofitopt[1])))
    #print(einastofitopt)
    for x in np.arange(len(einastofitopt[0])):
        for y in np.arange(len(einastofitopt[1])):
            #print(x,y)
            results[x,y] = chiSquareeinasto(rad,den,uncertainties,[einastofitopt[0][x,y],einastofitopt[1][x,y]])
    return(results)

def chiSquareeinasto(rad,den, uncertainties, nfwfitp):
    return(np.sum((np.square(((den))-(einasto(rad, nfwfitp[0], nfwfitp[1]))))/np.square(uncertainties)))

def getChiSquareeinasto(rad,den,uncertainties,nfwfitopt):
    results = np.zeros((len(nfwfitopt[0]),len(nfwfitopt[1])))
    #print(nfwfitopt)
    for x in np.arange(len(nfwfitopt[0])):
        for y in np.arange(len(nfwfitopt[1])):
            #print(x,y)
            results[x,y] = chiSquareeinasto(rad,den,uncertainties,[nfwfitopt[0][x,y],nfwfitopt[1][x,y]])
    return(results)



def plotting(rad, den,uncertainties, virial_radius, virial_density,nfwfitp,burkertfitp,dehnen_threeparamfitp,dehnen_twoparamfitp,einastofitp):
    #fig, axs = plt.subplots(3, 2, figsize=(15,15))
    #fig = plt.figure(figsize=(30,20))
    fig, axs = plt.subplots(1, 2, figsize=(20,20))
    #gs0 = fig.add_gridspec(3, 3)
    #   gs0 = fig.add_gridspec(1, 2)
    
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
    
    
    x, y = np.linspace(1000, 700000, 1000), np.linspace(0, 15, 1000)
    X, Y = np.meshgrid(x, y)
    Z = getChiSquareplot(rad, den, uncertainties, [X,Y])
    #print(X)
    #print(Y)
    #print(Z)
    normalizeC = np.min(Z)
    #print(normalizeC)
    Z = Z/normalizeC
    print(np.max(Z))
    #Z[np.where(Z==np.inf)[0]] = 10**10
    #print(Z)
    #pc = axs2.pcolor(X,Y,Z,norm=matplotlib.colors.LogNorm(vmin=1, vmax=15000),cmap='PuBu_r')
    pc = axs2.pcolor(X,Y,Z,norm=matplotlib.colors.LogNorm(vmin=1, vmax=200),cmap='PuBu_r')
    #pc = axs2.pcolor(Z)
    
    #print('2')
    einastofitp0, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
    axs[1].scatter(einastofitp0[0],einastofitp0[1], color='red',s=200, label="Optimised Parameters initial parameters: Virrad and Virdensity")
    axs[1].scatter(den[-1],rad[-1], color='red',marker='x', s=200,label="Initial parameters: Virrad and Virdensity")
    #print('3')
    init1 = [5000,8]
    einastofitp1, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=init1, sigma=uncer)
    axs[1].scatter(einastofitp1[0],einastofitp1[1], color='pink',s=200, label="Optimised Parameters initial parameters: "+str(init1))
    axs[1].scatter(init1[0],init1[1], color='pink',marker='x',s=200, label="Initial parameters: "+str(init1))
    #print('4')
    init2 = [70000,5]
    einastofitp2, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=init2, sigma=uncer)
    axs[1].scatter(einastofitp2[0],einastofitp2[1], color='cyan',s=200, label="Optimised Parameters initial parameters: "+str(init2))
    axs[1].scatter(init2[0],init2[1], color='cyan',marker='x',s=200, label="Initial parameters:  "+str(init2))
    #print('5')
    init3 = [80000,7]
    einastofitp3, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=init3, sigma=uncer)
    axs[1].scatter(einastofitp3[0],einastofitp3[1], color='green', s=200,label="Optimised Parameters initial parameters: "+str(init3))
    axs[1].scatter(init3[0],init3[1], color='green',marker='x',s=200, label="Initial parameters: "+str(init3))
    #print('6')
    init4 = [650000,2.2]
    einastofitp4, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=init4, sigma=uncer)
    axs[1].scatter(einastofitp4[0],einastofitp4[1], color='black',s=200, label="Optimised Parameters initial parameters: "+str(init4))
    axs[1].scatter(init4[0],init4[1], color='black',marker='x',s=200, label="Initial parameters: "+str(init4))
    
    
    
    
    axs[1].set_xlabel("Scale Density")
    axs[1].set_ylabel("Scale Radius")
    axs[1].set_title("Heatmap of ChiSquare values for einasto profile parameters (minimum ChiSquare set to 1)")
    axs[1].legend()
    fig.colorbar(pc, ax=axs2, extend='max')
    #print(einastofitp0, einastofitp1, einastofitp2,einastofitp3,einastofitp4)
    
    print(rad,den)
    axs1 = axs[0]
    axs1.errorbar((rad/virial_radius), (den/virial_density), yerr=uncer/virial_density, fmt='.', label="Halo_"+str(g)+"_099", color='blue')
    axs1.scatter(rad/virial_radius, einasto(rad, einastofitp0[0], einastofitp0[1])/virial_density, marker='x', label="einasto fit optimised param"+str(einastofitp0), color='red')
    axs1.scatter(rad/virial_radius, einasto(rad, init1[0], init1[1])/virial_density,marker='x', label="einasto fit initial param"+str(init1), color='pink')
    axs1.scatter(rad/virial_radius, einasto(rad, init2[0], init2[1])/virial_density, marker='x', label="einasto fit initial param"+str(init2), color='cyan')
    axs1.scatter(rad/virial_radius, einasto(rad, init3[0], init3[1])/virial_density, marker='x', label="einasto fit initial param"+str(init3), color='green')
    #axs1.scatter(rad/virial_radius, einasto(rad, init4[0], init4[1])/virial_density, marker='x', label="einasto fit initial param"+str(init4), color='black')
    
    #print(chiSquareeinasto(rad, den, uncertainties, einastofitp1))
    #print(chiSquareeinasto(rad, den, uncertainties, init1))
    #print(chiSquareeinasto(rad, den, uncertainties, init2))
    #print(chiSquareeinasto(rad, den, uncertainties, init3))
    #print(chiSquareeinasto(rad, den, uncertainties, init4))
    #print(min([chiSquareeinasto(rad,uncertainties, den, einastofitp1),chiSquareeinasto(rad, den,uncertainties, init1),chiSquareeinasto(rad, den,uncertainties, init2),chiSquareeinasto(rad, den,uncertainties, init3),chiSquareeinasto(rad, den,uncertainties, init4)]))
    
    axs1.set_xlabel(r'(Radius ($kpc/(R_{200}})}$))')
    axs1.set_ylabel(r'($\rho$(r) ($M_{\odot} kpc^{-3} (\rho_{200})^{-1}$))')
    axs1.legend()
    axs1.set_yscale('log')
    axs1.set_xscale('log')
    axs1.set_title("Fits to Data from TNG")
    
    fig.tight_layout()
    fig.savefig('HaloFitsInfo/einsto-fit-profiles-halo'+str(g))
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
gg =[6,20,2000,200,19999]
numhalos = len(subhalo_index)
#densities = []
#uncertainties = []
#radii = []

pdheader = ['Radius','Density','Uncertainty']
#derek = pd.DataFrame(columns=pdheader)
"""
with open('HaloFitsInfo/50-1_snap_99_fit_param.csv', 'x', encoding='UTF8', newline='') as f:
    
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
"""
einasto_chiquare = []
nfw_chisquare = []
for g in gg:
    data_csv = pd.read_csv('HaloFitsInfo/snap_99_halo_'+str(g)+'rad-den.csv')
    
    rad = data_csv['Radius'].to_numpy()
    den = data_csv['Density'].to_numpy()
    uncer = data_csv['Uncertainty'].to_numpy()
    num_datapoints = len(rad)
    
    virial_radius,virial_index = virialRadius(rad, den)
    virial_density = p_crit*200
    print(str((virial_radius/radius)[g]))
    print(g)
    if (virial_radius/radius)[g] < 15:
        
        filename = 'HaloFitsInfo/snap_99_halo_'+str(g)+'param'
        
        rad = rad[virial_index]
        den = den[virial_index]
        uncer = uncer[virial_index]
        nfwfitp, nfwfitcov = scopt.curve_fit(nfw, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
        nfwchi_square_test_statistic =  np.sum((np.square(((den))-(nfw(rad, nfwfitp[0], nfwfitp[1]))))/(nfw(rad, nfwfitp[0], nfwfitp[1])))
        nfwp_value = scipy.stats.distributions.chi2.sf(nfwchi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for NFW', nfwchi_square_test_statistic)
        print ('Fitted value for NFW', nfwfitp)
        print ('uncertainties for NFW', np.sqrt(np.diag(nfwfitcov)))

        einastofitp, einastofitcov = scopt.curve_fit(einasto, rad, den, p0=[den[-1],rad[-1]], sigma=uncer)
        einastochi_square_test_statistic =  np.sum(np.square(((den))-(einasto(rad, einastofitp[0], einastofitp[1])))/np.square(uncer))
        einastop_value = scipy.stats.distributions.chi2.sf(einastochi_square_test_statistic,(len(den)-1))
        print ('ChiSquare and P values for Einasto', einastochi_square_test_statistic, einastop_value)
        print ('Fitted value for Einasto', einastofitp)
        print ('uncertainties for Einasto', np.sqrt(np.diag(einastofitcov)))
        
        nfw_chisquare.append(nfwchi_square_test_statistic)
        einasto_chiquare.append(einastochi_square_test_statistic)
        """
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
        '''
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
        
        
        '''
        
        """
        burkertfitp = [0,0]
        dehnen_threeparamfitp = [0,0]
        dehnen_twoparamfitp = [0,0]
        nfwfitp = [0,0]
        #plotting(rad, den, uncer, virial_radius, virial_density,nfwfitp,burkertfitp,dehnen_threeparamfitp,dehnen_twoparamfitp,einastofitp)
        print(g)
    #g += 1
print(nfw_chisquare,einasto_chiquare)
difference = np.array(nfw_chisquare)-np.array(einasto_chiquare)
print(difference/np.array(nfw_chisquare))
fig, ax = plt.subplots(1)
ax.plot(gg,nfw_chisquare,'_', label ="NFW ChiSquare", color="magenta", markersize = 20 )
ax2 = ax.twinx()
ax2.plot(gg,einasto_chiquare,'.',  label ="Einasto ChiSquare", color="cornflowerblue", markersize = 20)
ax.set_xlabel('Halo Number')
ax.set_ylabel('NWF ChiSquare')
ax2.set_ylabel('Einasto ChiSquare')
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
plt.gca().set_xscale('log')
ax.set_yscale('log')
ax2.set_yscale('log')
#ax.legend()
ax.legend(h1+h2,l1+l2)

fig.savefig('chisquarediff.jpg')
fig.show()
        #fig.clf()
        #fig.show()
