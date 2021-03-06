# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:22:05 2022

@author: clara
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
import statistics
#import scikit-learn as sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import sklearn.metrics
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (30, 10),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)


h = 0.6774
p_crit = 127 #m_sun/(kpc^3)

#Get key info from group catalogue from file


#Read in the optimal fit parameters as well as chisquare
fit_param = pd.read_csv('50-1_snap_99_fit_param.csv')
print(fit_param)
fit_param['Halo Number'] = fit_param['Halo Number'].astype(int)
print(fit_param)
fit_param = fit_param.set_index('Halo Number')



subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_info['SubhaloIndex'] = subhalo_info['SubhaloIndex'].astype(int)
subhalo_info = subhalo_info.set_index(['SubhaloIndex'])
#subhalo_info['Df_cat'] = pd.Categorical(subhalo_info['SubhaloIndex'],
#                                             categories = true_indices,
#                                             ordered=True)
#print(subhalo_info_dark.sort_values('Df_cat'))
#sorted_df = subhalo_info.sort_values('Df_cat').dropna()
#print(sorted_df)
sorted_fit, sorted_df = fit_param.align(subhalo_info, join='inner', axis=0)

mass_sorted = sorted_df['SubhaloMass'].to_numpy()


#Similarly for DMO


#Read in the optimal fit parameters as well as chisquare
fit_param_dark = pd.read_csv('50-1_snap_99_fit_param-dark.csv')
fit_param_dark['Halo Number']=fit_param_dark['Halo Number'].astype(int)
fit_param_dark = fit_param_dark.set_index(['Halo Number'])



#Get key info from group catalogue from file
subhalo_info_dark = pd.read_csv('50-1-subhalo-info-dark',index_col=0)
#subhalo_info_dark['SubhaloIndex'] = subhalo_info_dark['SubhaloIndex'].astype(int)

#subhalo_info_dark = subhalo_info_dark.set_index('SubhaloIndex')
#subhalo_info_dark['Df_cat'] = pd.Categorical(subhalo_info_dark['SubhaloIndex'],
#                                             categories = true_indices_dark,
#                                             ordered=True)
#print(subhalo_info_dark.sort_values('Df_cat'))
#sorted_df_dark = subhalo_info_dark.sort_values('Df_cat').dropna()
#print(sorted_df)
sorted_fit_dark, sorted_df_dark = fit_param_dark.align(subhalo_info_dark, join='inner', axis=0)

mass_sorted_dark = sorted_df_dark['SubhaloMass'].to_numpy()


nfw_chisquare = sorted_fit['NFW ChiSquare'].to_numpy()
nfw_scalerad = sorted_fit['NFW Scale Radius'].to_numpy()
datapoint = sorted_fit['DataPoints'].to_numpy()
virrad = sorted_fit['Virial Radius'].to_numpy()
#true_indices = sorted_fit['Halo Number'].to_numpy().astype(int)



nfw_chisquare_dark = sorted_fit_dark['NFW ChiSquare'].to_numpy()
nfw_scalerad_dark = sorted_fit_dark['NFW Scale Radius'].to_numpy()
datapoint_dark = sorted_fit_dark['DataPoints'].to_numpy()
virrad_dark = sorted_fit_dark['Virial Radius'].to_numpy()
#true_indices_dark = sorted_fit_dark['Halo Number'].to_numpy().astype(int)



#lists to store data#
#numhalos = len(subhalo_index)
print(sorted_df)
print(sorted_df_dark)

#calculate concentration from given arrays
concentration = np.array(virrad/nfw_scalerad)
concentration_dark = np.array(virrad_dark/nfw_scalerad_dark)

y = concentration
y_dark = concentration_dark
print(y)
print(y_dark)


#Prepare mass to be binned


mean_mass = []
conc_hist = []
conc_hist_dark = []

conc_hist_ML = []
conc_hist_dark_ML = []


lowerbound = 0.01
bins = [0.1,1,10]


fig = plt.figure(figsize=(30,30))
# add grid specifications
gs = fig.add_gridspec(3, 3)
# open axes/subplots
axs = []
axs.append( fig.add_subplot(gs[0,0]) )
axs.append( fig.add_subplot(gs[0,1]) )   
axs.append( fig.add_subplot(gs[0,2]) )  
axs.append( fig.add_subplot(gs[1,0]) )
axs.append( fig.add_subplot(gs[1,1]) )
axs.append( fig.add_subplot(gs[1,2]) )
axs.append( fig.add_subplot(gs[2,0]) )  
axs.append( fig.add_subplot(gs[2,1]) )
axs.append( fig.add_subplot(gs[2,2]) )
print(len(axs))
i = 0
#loop over bins
for upperbound in bins:
    #get indices of mass lying inside bin
    massindex = np.where(np.logical_and(mass_sorted<upperbound,mass_sorted>lowerbound))[0]
    mass_sorted = np.array(mass_sorted)
    massindex_dark = np.where(np.logical_and(mass_sorted_dark<upperbound,mass_sorted_dark>lowerbound))[0]
    
    print(len(massindex))
    print(len(concentration[massindex]))
    print(len(massindex_dark))
    print(len(concentration_dark[massindex_dark]))
    #append all data to lists
    mean_mass.append(((upperbound-lowerbound)/2)*h)
    
    sorted_df_inter = sorted_df.reset_index()
    sorted_df_inter_dark = sorted_df_dark.reset_index()
    #print(len(sorted_df_inter['SubhaloSpinZ'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])]))
    #print(len(sorted_df_inter_dark['SubhaloSpinZ'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])]))

    conc_hist.append(concentration[massindex])
    conc_hist_dark.append(concentration_dark[massindex_dark])
    
    bin_X_final = pd.DataFrame([sorted_df_inter['SubhaloGasMass'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                                sorted_df_inter['SubhaloStarMass'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloBHMass'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloDMMass'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloSpinX'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloSpinY'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloSpinZ'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloVelDisp'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloVmax'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloBHMdot'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['SubhaloSFR'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['FoFMass'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])],
                             sorted_df_inter['FoFDistanceCenter'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])]]).T
    
    bin_X_dark_final = pd.DataFrame([sorted_df_inter_dark['SubhaloDMMass'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['SubhaloSpinX'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['SubhaloSpinY'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['SubhaloSpinZ'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['SubhaloVelDisp'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['SubhaloVmax'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['FoFMass'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])],
                             sorted_df_inter_dark['FoFDistanceCenter'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])]]).T
        
        
    Xtrain, Xtest, ytrain, ytest = train_test_split(bin_X_final, concentration[massindex],
                                                    random_state=1)
    
    Xtrain_dark, Xtest_dark, ytrain_dark, ytest_dark = train_test_split(bin_X_dark_final, concentration_dark[massindex_dark],
                                                    random_state=1)
    
    #Train Model for DM+Baryons
    model = RandomForestRegressor(n_estimators=1000,n_jobs=50)
    model.fit(Xtrain,ytrain)
    y_pred = model.predict(Xtest)
    importances = model.feature_importances_
    std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)

    print(y)
    print(y_pred)
    #Train Model for DMO
    model_dark = RandomForestRegressor(n_estimators=1000,n_jobs=50)
    model_dark.fit(Xtrain_dark,ytrain_dark)
    y_pred_dark = model_dark.predict(Xtest_dark)
    importances_dark = model_dark.feature_importances_
    std_dark = np.std([tree_dark.feature_importances_ for tree_dark in model_dark.estimators_], axis=0)

    
    print('R_2')
    print(sklearn.metrics.r2_score(ytest, y_pred))
    print(sklearn.metrics.r2_score(ytest_dark, y_pred_dark))
    
    print('Mean squared error')
    print(sklearn.metrics.mean_squared_error(ytest, y_pred))
    print(sklearn.metrics.mean_squared_error(ytest_dark, y_pred_dark))
    
    forest_importances = pd.Series(importances, index=['SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                    'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                    'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter'])

    forest_importances_dark = pd.Series(importances_dark, index=['SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                    'FoFMass','FoFDistanceCenter'])

        
    axs[i].hist(conc_hist[-1][np.where(conc_hist[-1]<30)[0]], alpha = 0.5, color='magenta', label='Full Physics', density = True, bins=100)
    axs[i].hist(y_pred_dark[np.where(y_pred_dark<30)[0]], alpha = 0.5, color='green', label='ML DMO', density = True, bins=100)
    axs[i].hist(y_pred[np.where(y_pred<30)[0]], alpha = 0.5, color='blue', label='ML Full Physics', density = True, bins=100)
    axs[i].hist(conc_hist_dark[-1][np.where(conc_hist_dark[-1]<30)[0]], alpha = 0.5, color='red', label='DMO', density = True, bins=100)
    

    axs[i].set_xlabel(r'$c_{200}$')
    axs[i].set_ylabel(r'Number of Halos')
    axs[i].set_xlim(0,30)
    axs[i].set_title('Mass Bin '+str(round(mean_mass[-1],5))+r' $10^{10} M_{\odot}$')
    axs[i].legend()
    
    forest_importances.plot.bar(yerr=std, ax=axs[i+1])
    axs[i+1].set_xlabel(r'Feature importances using MDI')
    axs[i+1].set_ylabel(r'Mean decrease in impurity')
    axs[i+1].set_title('Feature Importance DM+Baryons')
    
    forest_importances_dark.plot.bar(yerr=std_dark, ax=axs[i+2])
    axs[i+2].set_xlabel(r'Feature importances using MDI')
    axs[i+2].set_ylabel(r'Mean decrease in impurity')
    axs[i+2].set_title('Feature Importance DMO')
        
    i += 3
    
    lowerbound= upperbound
fig.tight_layout()
fig.savefig('ML-per-bin-newindex')


