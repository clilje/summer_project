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

#get the input data from the Group Catalogues
#DM + Baryons
"""
data_csv = pd.read_csv('50-1-subhalo-info.csv',usecols=['SubhaloIndex','SubhaloMass','SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter'])
column_names = ['SubhaloIndex','SubhaloMass','SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter']
#X = data_csv[column_names]

#DMO
data_csv_dark = pd.read_csv('50-1-subhalo-info-dark.csv',usecols=['SubhaloIndex','SubhaloMass','SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 
                     'SubhaloVmax','FoFMass','FoFDistanceCenter', ])
column_names_dark = ['SubhaloIndex','SubhaloMass','SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 
                     'SubhaloVmax','FoFMass','FoFDistanceCenter']
#X_dark = data_csv_dark[column_names_dark]
"""
"""
#Get the nessecary data to calculate the concentration from the fit files

#Prepare the dataframes to be resorted to match Halo Indices
fit_param_dark = pd.read_csv('50-1_snap_99_fit_param-dark.csv')
true_indices_dark = fit_param_dark['Halo Number'].to_numpy().astype(int)

#Import the Fit Parameters for DM+Baryons
#Reorder according to DMO Halo Indices, cut all NaN
fit_param = pd.read_csv('50-1_snap_99_fit_param.csv')
fit_param['Df_cat'] = pd.Categorical(fit_param['Halo Number'],
                                             categories = true_indices_dark,
                                             ordered=True)
sorted_df = fit_param.sort_values('Df_cat').dropna()
nfw_scalerad = sorted_df['NFW Scale Radius'].to_numpy()
virrad = sorted_df['Virial Radius'].to_numpy()
true_indices = sorted_df['Halo Number'].to_numpy().astype(int)

#Reorder DMO Halos according to leftover DM+Baryon indices, cut all NaN
fit_param_dark['Df_cat'] = pd.Categorical(fit_param_dark['Halo Number'],
                                             categories = true_indices,
                                             ordered=True)
sorted_df_dark = fit_param_dark.sort_values('Df_cat').dropna()
nfw_scalerad_dark = sorted_df_dark['NFW Scale Radius'].to_numpy()
virrad_dark = sorted_df_dark['Virial Radius'].to_numpy()


#Reorder the input values according to leftover DM+Baryon indices, cut all NaN
X['Df_cat'] = pd.Categorical(X['SubhaloIndex'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data = X.sort_values('Df_cat').dropna().copy()
sorted_X = pd.DataFrame([sorted_data['SubhaloIndex'],sorted_data['SubhaloMass'],sorted_data['SubhaloGasMass'],sorted_data['SubhaloStarMass'],
                         sorted_data['SubhaloBHMass'],sorted_data['SubhaloDMMass'],
                         sorted_data['SubhaloSpinX'],sorted_data['SubhaloSpinY'],
                         sorted_data['SubhaloSpinZ'],sorted_data['SubhaloVelDisp'],
                         sorted_data['SubhaloVmax'],sorted_data['SubhaloBHMdot'],
                         sorted_data['SubhaloSFR'],sorted_data['FoFMass'],
                         sorted_data['FoFDistanceCenter']]).T

X_dark['Df_cat'] = pd.Categorical(X_dark['SubhaloIndex'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data_dark = X_dark.sort_values('Df_cat').dropna().copy()
sorted_X_dark = pd.DataFrame([sorted_data_dark['SubhaloIndex'],sorted_data_dark['SubhaloMass'],sorted_data_dark['SubhaloDMMass'],
                         sorted_data_dark['SubhaloSpinX'],sorted_data_dark['SubhaloSpinY'],
                         sorted_data_dark['SubhaloSpinZ'],sorted_data_dark['SubhaloVelDisp'],
                         sorted_data_dark['SubhaloVmax'],sorted_data_dark['FoFMass'],
                         sorted_data_dark['FoFDistanceCenter']]).T
"""#

#Calculate concentration from Re-indexed input arrays and set as expected value
#concentration = virrad/nfw_scalerad
#concentration_dark = virrad_dark/nfw_scalerad_dark



#print(sorted_X)
#print(sorted_X_dark)




#Get key info from group catalogue from file
#subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
#subhalo_index = subhalo_info['SubhaloIndex'].to_numpy().astype(int)
#full_mass = subhalo_info['SubhaloMass'].to_numpy()


#Read in the optimal fit parameters as well as chisquare
fit_param = pd.read_csv('50-1_snap_99_fit_param.csv')
nfw_chisquare = fit_param['NFW ChiSquare'].to_numpy()
nfw_scalerad = fit_param['NFW Scale Radius'].to_numpy()
datapoint = fit_param['DataPoints'].to_numpy()
virrad = fit_param['Virial Radius'].to_numpy()
true_indices = fit_param['Halo Number'].to_numpy().astype(int)

subhalo_info = pd.read_csv('50-1-subhalo-info.csv')
subhalo_info['Df_cat'] = pd.Categorical(subhalo_info['SubhaloIndex'],
                                             categories = true_indices,
                                             ordered=True)
#print(subhalo_info_dark.sort_values('Df_cat'))
sorted_df = subhalo_info.sort_values('Df_cat').dropna()
#print(sorted_df)
mass_sorted = sorted_df['SubhaloMass']


#Similarly for DMO


#Read in the optimal fit parameters as well as chisquare
fit_param_dark = pd.read_csv('50-1_snap_99_fit_param-dark.csv')
nfw_chisquare_dark = fit_param_dark['NFW ChiSquare'].to_numpy()
nfw_scalerad_dark = fit_param_dark['NFW Scale Radius'].to_numpy()
datapoint_dark = fit_param_dark['DataPoints'].to_numpy()
virrad_dark = fit_param_dark['Virial Radius'].to_numpy()
true_indices_dark = fit_param_dark['Halo Number'].to_numpy().astype(int)


#Get key info from group catalogue from file
subhalo_info_dark = pd.read_csv('50-1-subhalo-info-dark.csv')
subhalo_info_dark['Df_cat'] = pd.Categorical(subhalo_info_dark['SubhaloIndex'],
                                             categories = true_indices_dark,
                                             ordered=True)
#print(subhalo_info_dark.sort_values('Df_cat'))
sorted_df_dark = subhalo_info_dark.sort_values('Df_cat').dropna()
#print(sorted_df)
mass_sorted_dark = sorted_df_dark['SubhaloMass']
#lists to store data#
#numhalos = len(subhalo_index)
print(sorted_df)
print(sorted_df_dark)

#calculate concentration from given arrays
concentration = virrad/nfw_scalerad
concentration_dark = virrad_dark/nfw_scalerad_dark

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


fig = plt.figure(figsize=(20,50))
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

i = 0
#loop over bins
for upperbound in bins:
    #get indices of mass lying inside bin
    massindex = np.where(np.logical_and(mass_sorted<upperbound,mass_sorted>lowerbound))[0]
    mass_sorted = np.array(mass_sorted)
    massindex_dark = np.where(np.logical_and(mass_sorted_dark<upperbound,mass_sorted_dark>lowerbound))[0]
    
    #get indices for which concentration values correspond to this
    #conc_index = np.searchsorted(true_indices, subhalo_index[true_indices][massindex])
    
    #ensure no error for stdev or mean
    #if conc_index.size ==1 or conc_index.size ==0:
    #    break
    print(len(massindex))
    print(len(concentration[massindex]))
    print(len(massindex_dark))
    print(len(concentration_dark[massindex_dark]))
    #append all data to lists
    mean_mass.append(((upperbound-lowerbound)/2)*h)
    #mean_concentration.append(statistics.mean(concentration[conc_index]))
    #stdev.append(statistics.stdev(concentration[conc_index]))
    
    sorted_df_inter = sorted_df.reset_index()
    sorted_df_inter_dark = sorted_df_dark.reset_index()
    print(len(sorted_df_inter['SubhaloSpinZ'][sorted_df_inter.SubhaloMass.isin(mass_sorted[massindex])]))
    print(len(sorted_df_inter_dark['SubhaloSpinZ'][sorted_df_inter_dark.SubhaloMass.isin(mass_sorted_dark[massindex_dark])]))

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
    
    #set up plotting params
    fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(30, 10))
    
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
fig.savefig('ML-per-bin')
#print(bins_dark)












"""
for upperbound_dark in bins_dark:
    #get indices of mass lying inside bin
    #mass_sorted = []
    #for x in true_indices_dark:
        #mass_sorted.append(full_mass_dark[np.where(x == subhalo_index_dark)[0]][0])
    #print(mass_sorted)
    #print(true_indices_dark)
    mass_sorted = np.array(mass_sorted)
    massindex_dark = np.where(np.logical_and(mass_sorted<upperbound_dark,mass_sorted>lowerbound_dark))[0]
    
    #get indices for which concentration values correspond to this
    #conc_index_dark = np.searchsorted(true_indices_dark, subhalo_index_dark[true_indices_dark][massindex_dark])
    conc_index_dark = massindex_dark
    #ensure no error for stdev or mean
    #print(massindex_dark)
    if conc_index_dark.size == 0 or conc_index_dark.size ==1:
        break
    
    #append all data to lists
    mean_mass_dark.append(((upperbound_dark-lowerbound_dark)/2)*h)
    mean_concentration_dark.append(statistics.mean(concentration_dark[conc_index_dark]))
    stdev_dark.append(statistics.stdev(concentration_dark[conc_index_dark]))
    
    conc_hist_dark.append(concentration_dark[conc_index_dark])
    
    lowerbound_dark= upperbound_dark
"""  
    
"""
#ML

#loop over bins
for upperbound_ML in bins_ML:
    #get indices of mass lying inside bin
    massindex_ML = np.where(np.logical_and(mass_sorted_ML<upperbound_ML,mass_sorted_ML>lowerbound_ML))[0]
    
    #get indices for which concentration values correspond to this
    conc_index_ML = massindex_ML
    #ensure no error for stdev or mean
    if conc_index_ML.size == 0 or conc_index_ML.size == 1:
        break
    
    #append all data to lists
    mean_mass_ML.append(((upperbound_ML-lowerbound_ML)/2)*h)
    mean_concentration_ML.append(statistics.mean(y_pred[conc_index_ML]))
    stdev_ML.append(statistics.stdev(y_pred[conc_index_ML]))
    
    conc_hist_ML.append(y_pred[conc_index_ML])
    
    lowerbound_ML= upperbound_ML
    
#print(bins_dark)

for upperbound_dark_ML in bins_dark_ML:
    #get indices of mass lying inside bin
    massindex_dark_ML = np.where(np.logical_and(mass_sorted_ML_dark<upperbound_dark_ML,mass_sorted_ML_dark>lowerbound_dark_ML))[0]
    conc_index_dark_ML = massindex_dark_ML
    #ensure no error for stdev or mean
    if conc_index_dark_ML.size ==0 or conc_index_dark_ML.size ==1:
        break
    
    #append all data to lists
    mean_mass_dark_ML.append(((upperbound_dark_ML-lowerbound_dark_ML)/2)*h)
    mean_concentration_dark_ML.append(statistics.mean(y_pred_dark[conc_index_dark_ML]))
    stdev_dark_ML.append(statistics.stdev(y_pred_dark[conc_index_dark_ML]))
    
    conc_hist_dark_ML.append(y_pred_dark[conc_index_dark_ML])
    
    lowerbound_dark_ML= upperbound_dark_ML  



#Plotting concentration Histograms

#print(len(conc_hist))
#print(len(conc_hist_dark))
#print(max(conc_hist[1]))
#print(max(conc_hist_dark[1]))
# open figure
fig = plt.figure(figsize=(20,50))
# add grid specifications
gs = fig.add_gridspec(5, 3)
# open axes/subplots
axs = []
axs.append( fig.add_subplot(gs[0,:]) ) 
axs.append( fig.add_subplot(gs[1,0]) )
axs.append( fig.add_subplot(gs[1,1]) )   
axs.append( fig.add_subplot(gs[1,2]) )  
axs.append( fig.add_subplot(gs[2,0]) )
axs.append( fig.add_subplot(gs[2,1]) )
axs.append( fig.add_subplot(gs[2,2]) )
axs.append( fig.add_subplot(gs[3,0]) )  
axs.append( fig.add_subplot(gs[3,1]) )
axs.append( fig.add_subplot(gs[3,2]) )
axs.append( fig.add_subplot(gs[4,0]) )
axs.append( fig.add_subplot(gs[4,1]) )  
axs.append( fig.add_subplot(gs[4,2]) )



#fig = plt.figure(constrained_layout=True, figsize=(10, 50))
#subfigs = fig.subfigures(2, 1, wspace=0.07)
#axsup = subfigs[0]
#axsdown = subfigs[1]
#axsdsub = axsdown.add_subplot(6,2,(1,1))
#plot data obtained in scatterplot
offset = np.logspace(-3,3,10)

axs[0].errorbar(np.array(mean_mass),mean_concentration,yerr=stdev,fmt='.', color='gold', label='Full Physics Run',markersize='10')
axs[0].errorbar(np.array(mean_mass_dark)+np.logspace(-3,3,len(mean_mass_dark)),mean_concentration_dark,yerr=stdev_dark,fmt='x',color='darksalmon',label='DMO',markersize='10')

axs[0].errorbar(np.array(mean_mass_ML)+2*np.logspace(-3,3,len(mean_mass_ML)),mean_concentration_ML,yerr=stdev_ML,fmt='.', color='indigo', label='ML Full Physics Run',markersize='10')
axs[0].errorbar(np.array(mean_mass_dark_ML)+3*np.logspace(-3,3,len(mean_mass_dark_ML)),mean_concentration_dark_ML,yerr=stdev_dark_ML,fmt='x',color='darkslategrey',label='ML DMO',markersize='10')

axs[0].set_xscale('log')
axs[0].set_yscale('log')
#plt.ylim((1,30))
axs[0].set_xlabel(r'Total Mass of Halo in $10^{10} M_{\odot}$')
axs[0].set_ylabel(r'$c_{200}$')
axs[0].set_title('Concentration Mass Relation')
axs[0].legend(loc='lower right')

i = 1
while i < len(conc_hist) and i <12:
    axs[i].hist(conc_hist[i-1][np.where(conc_hist[i-1]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)

j = 1 
while j < len(conc_hist_dark)and j <12:
    axs[j].hist(conc_hist_dark[j-1][np.where(conc_hist_dark[j-1]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)

x = 1
while x < len(conc_hist_ML) and x <12:
    axs[x].hist(conc_hist_ML[x-1][np.where(conc_hist_ML[x-1]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)

y = 1
while y < len(conc_hist_dark_ML) and y <12:
    axs[y].hist(conc_hist_dark_ML[y-1][np.where(conc_hist_dark_ML[y-1]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

#axs[1].hist(conc_hist[0][np.where(conc_hist[0]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[1].hist(conc_hist_dark[0][np.where(conc_hist_dark[0]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[1].hist(conc_hist_ML[0][np.where(conc_hist_ML[0]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[1].hist(conc_hist_dark_ML[0][np.where(conc_hist_dark_ML[0]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)


axs[1].set_xlabel(r'$c_{200}$')
axs[1].set_ylabel(r'Number of Halos')
axs[1].set_xlim(0,30)
axs[1].set_title('Mass Bin '+str(round(mean_mass[0],5))+r' $10^{10} M_{\odot}$')
axs[1].legend()

#axs[2].hist(conc_hist[1][np.where(conc_hist[1]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[2].hist(conc_hist_dark[1][np.where(conc_hist_dark[1]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[2].hist(conc_hist_ML[1][np.where(conc_hist_ML[1]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[2].hist(conc_hist_dark_ML[1][np.where(conc_hist_dark_ML[1]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

axs[2].set_xlabel(r'$c_{200}$')
axs[2].set_ylabel(r'Number of Halos')
axs[2].set_xlim(0,30)
axs[2].set_title('Mass Bin '+str(round(mean_mass[1],4))+r' $10^{10} M_{\odot}$')
axs[2].legend()

#axs[3].hist(conc_hist[2][np.where(conc_hist[2]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[3].hist(conc_hist_dark[2][np.where(conc_hist_dark[2]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[3].hist(conc_hist_ML[2][np.where(conc_hist_ML[2]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[3].hist(conc_hist_dark_ML[2][np.where(conc_hist_dark_ML[2]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

axs[3].set_xlabel(r'$c_{200}$')
axs[3].set_ylabel(r'Number of Halos')
axs[3].set_title('Mass Bin '+str(round(mean_mass[2],4))+r' $10^{10} M_{\odot}$')
axs[3].legend()

#axs[4].hist(conc_hist[3][np.where(conc_hist[3]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[4].hist(conc_hist_dark[3][np.where(conc_hist_dark[3]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[4].hist(conc_hist_ML[3][np.where(conc_hist_ML[3]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[4].hist(conc_hist_dark_ML[3][np.where(conc_hist_dark_ML[3]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

axs[4].set_xlabel(r'$c_{200}$')
axs[4].set_ylabel(r'Number of Halos')
axs[4].set_title('Mass Bin '+str(round(mean_mass[3],4))+r' $10^{10} M_{\odot}$')
axs[4].legend()

#axs[5].hist(conc_hist[4][np.where(conc_hist[4]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[5].hist(conc_hist_dark[4][np.where(conc_hist_dark[4]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[5].hist(conc_hist_ML[4][np.where(conc_hist_ML[4]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[5].hist(conc_hist_dark_ML[4][np.where(conc_hist_dark_ML[4]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

#axs[5].set_xlabel(r'$c_{200}$')
#axs[5].set_ylabel(r'Number of Halos')
#axs[5].set_title('Mass Bin '+str(round(mean_mass[4],4))+r' $10^{10} M_{\odot}$')
#axs[5].legend()

#axs[6].hist(conc_hist[5][np.where(conc_hist[5]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[6].hist(conc_hist_dark[5][np.where(conc_hist_dark[5]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[6].hist(conc_hist_ML[5][np.where(conc_hist_ML[5]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[6].hist(conc_hist_dark_ML[5][np.where(conc_hist_dark_ML[5]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

axs[6].set_xlabel(r'$c_{200}$')
axs[6].set_ylabel(r'Number of Halos')
axs[6].set_title('Mass Bin '+str(round(mean_mass[5],4))+r' $10^{10} M_{\odot}$')
axs[6].legend()

#axs[7].hist(conc_hist[6][np.where(conc_hist[6]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=100)
#axs[7].hist(conc_hist_dark[6][np.where(conc_hist_dark[6]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=100)
#axs[7].hist(conc_hist_ML[6][np.where(conc_hist_ML[6]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=100)
#axs[7].hist(conc_hist_dark_ML[6][np.where(conc_hist_dark_ML[6]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=100)

axs[7].set_xlabel(r'$c_{200}$')
axs[7].set_ylabel(r'Number of Halos')
axs[7].legend()

#axs[8].hist(conc_hist[7][np.where(conc_hist[7]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=20)
#axs[8].hist(conc_hist_dark[7][np.where(conc_hist_dark[7]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=20)
#axs[8].hist(conc_hist_ML[7][np.where(conc_hist_ML[7]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=20)
#axs[8].hist(conc_hist_dark_ML[7][np.where(conc_hist_dark_ML[7]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=20)

axs[8].set_xlabel(r'$c_{200}$')
axs[8].set_ylabel(r'Number of Halos')
axs[8].set_title('Mass Bin '+str(round(mean_mass[7],4))+r' $10^{10} M_{\odot}$')
axs[8].legend()

#axs[9].hist(conc_hist[8][np.where(conc_hist[8]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=15)
#axs[9].hist(conc_hist_dark[8][np.where(conc_hist_dark[8]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=15)
#axs[9].hist(conc_hist_ML[8][np.where(conc_hist_ML[8]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=15)
#axs[9].hist(conc_hist_dark_ML[8][np.where(conc_hist_dark_ML[8]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=15)

axs[9].set_xlabel(r'$c_{200}$')
axs[9].set_ylabel(r'Number of Halos')
axs[9].set_title('Mass Bin '+str(round(mean_mass[8],4))+r' $10^{10} M_{\odot}$')
axs[9].legend()
#
#axs[10].hist(conc_hist[9][np.where(conc_hist[9]<30)[0]], alpha = 0.5, color='gold', label='Full Physics', density = True, bins=10)
#axs[10].hist(conc_hist_dark[9][np.where(conc_hist_dark[9]<30)[0]], alpha = 0.5, color='darksalmon', label='DMO', density = True, bins=10)
#axs[10].hist(conc_hist_ML[9][np.where(conc_hist_ML[9]<30)[0]], alpha = 0.5, color='indigo', label='ML Full Physics', density = True, bins=10)
#axs[10].hist(conc_hist_dark_ML[9][np.where(conc_hist_dark_ML[9]<30)[0]], alpha = 0.5, color='darkslategrey', label='ML DMO', density = True, bins=10)

axs[10].set_xlabel(r'$c_{200}$')
axs[10].set_ylabel(r'Number of Halos')
axs[10].set_title('Mass Bin '+str(round(mean_mass[9],4))+r' $10^{10} M_{\odot}$')
axs[10].legend()


fig.tight_layout()
fig.savefig('cmfunc-bins-ML')
#fig.show()

def fraction(binnum,cut):   
    if binnum <= len(conc_hist):
        frac_bin = len(np.where(conc_hist[binnum]>cut)[0])/ len(conc_hist[binnum])
    else:
        frac_bin = str('No halos in this bin')
    if binnum <= len(conc_hist_dark):
        frac_bin_dark = len(np.where(conc_hist_dark[binnum]>cut)[0])/len(conc_hist_dark[binnum])
    else:
        frac_bin_dark = str('No halos in this bin')
    
    if binnum <= len(conc_hist_ML):
        frac_bin_ML = len(np.where(conc_hist_ML[binnum]>cut)[0])/ len(conc_hist_ML[binnum])
    else:
        frac_bin_ML = str('No halos in this bin')
    if binnum <= len(conc_hist_dark_ML):
        frac_bin_dark_ML = len(np.where(conc_hist_dark_ML[binnum]>cut)[0])/len(conc_hist_dark_ML[binnum])
    else:
        frac_bin_dark_ML = str('No halos in this bin')
    return(frac_bin,frac_bin_dark,frac_bin_ML,frac_bin_dark_ML)


bin0, bind0, bin0_ML, bind0_ML = fraction(0,60)
print('Bin 0: Fraction above 60: '+str(bin0)+' DMO: '+str(bind0)+' ML: '+str(bin0_ML)+' ML DMO: '+str(bind0_ML))
bin1, bind1, bin1_ML, bind1_ML = fraction(1,60)
print('Bin 1: Fraction above 60: '+str(bin1)+' DMO: '+str(bind1)+' ML: '+str(bin1_ML)+' ML DMO: '+str(bind1_ML))
bin2, bind2, bin2_ML, bind2_ML = fraction(2,60)
print('Bin 2: Fraction above 60: '+str(bin2)+' DMO: '+str(bind2)+' ML: '+str(bin2_ML)+' ML DMO: '+str(bind2_ML))
bin3, bind3, bin3_ML, bind3_ML = fraction(3,60)
print('Bin 3: Fraction above 60: '+str(bin3)+' DMO: '+str(bind3)+' ML: '+str(bin3_ML)+' ML DMO: '+str(bind3_ML))
bin4, bind4, bin4_ML, bind4_ML = fraction(4,60)
print('Bin 4: Fraction above 60: '+str(bin4)+' DMO: '+str(bind4)+' ML: '+str(bin4_ML)+' ML DMO: '+str(bind4_ML))
bin5, bind5, bin5_ML, bind5_ML = fraction(5,60)
print('Bin 5: Fraction above 60: '+str(bin5)+' DMO: '+str(bind5)+' ML: '+str(bin5_ML)+' ML DMO: '+str(bind5_ML))
bin6, bind6, bin6_ML, bind6_ML = fraction(6,60)
print('Bin 6: Fraction above 60: '+str(bin6)+' DMO: '+str(bind6)+' ML: '+str(bin6_ML)+' ML DMO: '+str(bind6_ML))
bin7, bind7, bin7_ML, bind7_ML = fraction(7,60)
print('Bin 7: Fraction above 60: '+str(bin7)+' DMO: '+str(bind7)+' ML: '+str(bin7_ML)+' ML DMO: '+str(bind7_ML))
bin8, bind8, bin8_ML, bind8_ML = fraction(8,60)
print('Bin 8: Fraction above 60: '+str(bin8)+' DMO: '+str(bind8)+' ML: '+str(bin8_ML)+' ML DMO: '+str(bind8_ML))
bin9, bind9, bin9_ML, bind9_ML = fraction(9,60)
print('Bin 9: Fraction above 60: '+str(bin9)+' DMO: '+str(bind9)+' ML: '+str(bin9_ML)+' ML DMO: '+str(bind9_ML))
print('Overall Fraction of halos with concentration above 30: '
      +str(bin0+bin1+bin2+bin3+bin4+bin5+bin6+bin7+bin8+bin9)+' DMO: '
      +str(bind0+bind1+bind2+bind3+bind4+bind5+bind6+bind7+bind8+bind9)+' ML: '
      +str(bin0_ML+bin1_ML+bin2_ML+bin3_ML+bin4_ML+bin5_ML+bin6_ML+bin7_ML+bin8_ML+bin9_ML)+' ML DMO: '
      +str(bind0_ML+bind1_ML+bind2_ML+bind3_ML+bind4_ML+bind5_ML+bind6_ML+bind7_ML+bind8_ML+bind9_ML))



bin0, bind0, bin0_ML, bind0_ML = fraction(0,30)
print('Bin 0: Fraction above 30: '+str(bin0)+' DMO: '+str(bind0)+' ML: '+str(bin0_ML)+' ML DMO: '+str(bind0_ML))
bin1, bind1, bin1_ML, bind1_ML = fraction(1,30)
print('Bin 1: Fraction above 30: '+str(bin1)+' DMO: '+str(bind1)+' ML: '+str(bin1_ML)+' ML DMO: '+str(bind1_ML))
bin2, bind2, bin2_ML, bind2_ML = fraction(2,30)
print('Bin 2: Fraction above 30: '+str(bin2)+' DMO: '+str(bind2)+' ML: '+str(bin2_ML)+' ML DMO: '+str(bind2_ML))
bin3, bind3, bin3_ML, bind3_ML = fraction(3,30)
print('Bin 3: Fraction above 30: '+str(bin3)+' DMO: '+str(bind3)+' ML: '+str(bin3_ML)+' ML DMO: '+str(bind3_ML))
bin4, bind4, bin4_ML, bind4_ML = fraction(4,30)
print('Bin 4: Fraction above 30: '+str(bin4)+' DMO: '+str(bind4)+' ML: '+str(bin4_ML)+' ML DMO: '+str(bind4_ML))
bin5, bind5, bin5_ML, bind5_ML = fraction(5,30)
print('Bin 5: Fraction above 30: '+str(bin5)+' DMO: '+str(bind5)+' ML: '+str(bin5_ML)+' ML DMO: '+str(bind5_ML))
bin6, bind6, bin6_ML, bind6_ML = fraction(6,30)
print('Bin 6: Fraction above 30: '+str(bin6)+' DMO: '+str(bind6)+' ML: '+str(bin6_ML)+' ML DMO: '+str(bind6_ML))
bin7, bind7, bin7_ML, bind7_ML = fraction(7,30)
print('Bin 7: Fraction above 30: '+str(bin7)+' DMO: '+str(bind7)+' ML: '+str(bin7_ML)+' ML DMO: '+str(bind7_ML))
bin8, bind8, bin8_ML, bind8_ML = fraction(8,30)
print('Bin 8: Fraction above 30: '+str(bin8)+' DMO: '+str(bind8)+' ML: '+str(bin8_ML)+' ML DMO: '+str(bind8_ML))
bin9, bind9, bin9_ML, bind9_ML = fraction(9,30)
print('Bin 9: Fraction above 30: '+str(bin9)+' DMO: '+str(bind9)+' ML: '+str(bin9_ML)+' ML DMO: '+str(bind9_ML))
print('Overall Fraction of halos with concentration above 30: '
      +str(bin0+bin1+bin2+bin3+bin4+bin5+bin6+bin7+bin8+bin9)+' DMO: '
      +str(bind0+bind1+bind2+bind3+bind4+bind5+bind6+bind7+bind8+bind9)+' ML: '
      +str(bin0_ML+bin1_ML+bin2_ML+bin3_ML+bin4_ML+bin5_ML+bin6_ML+bin7_ML+bin8_ML+bin9_ML)+' ML DMO: '
      +str(bind0_ML+bind1_ML+bind2_ML+bind3_ML+bind4_ML+bind5_ML+bind6_ML+bind7_ML+bind8_ML+bind9_ML))

"""