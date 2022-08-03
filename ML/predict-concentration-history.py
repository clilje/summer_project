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
import h5py
#import scikit-learn as sklearn
from sklearn.ensemble import RandomForestRegressor
import sklearn.metrics
import matplotlib
from sklearn.model_selection import train_test_split
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (30, 10),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

def get_matching(sim_size, sim_res):
    with h5py.File("/disk01/rmcg/downloaded/tng/tng"+str(sim_size)+"-"+str(sim_res)+"/subhalo_matching_to_dark.hdf5") as file:
        #print(file.keys())
        matchingarr = np.array(file['Snapshot_99/SubhaloIndexDark_LHaloTree'])
        #print(matchingarr)
    return(matchingarr)



h = 0.6774
p_crit = 127 #m_sun/(kpc^3)
matchingarr = get_matching(50, 1)

#get the input data from the Group Catalogues
#DM + Baryons
data_csv = pd.read_csv('50-1-subhalo-history.csv')
data_csv['index'] = data_csv['index'].astype(int)
data_csv =data_csv.set_index('index', drop=False)
#DMO
data_csv_dark = pd.read_csv('50-1-subhalo-history-dark.csv')
data_csv_dark['index'] = data_csv_dark['index'].astype(int)
data_csv_dark = data_csv_dark.set_index('index')

sorter = matchingarr.astype(int)
print(sorter)
print(data_csv_dark)
data_csv_dark = data_csv_dark.reindex(sorter)
print(data_csv_dark)
data_csv_dark.reset_index(inplace=True,drop=True)
print(data_csv_dark)
data_csv_dark.dropna(inplace=True)
data_csv_dark['index'] = data_csv_dark.index
print(data_csv_dark)

#Get the nessecary data to calculate the concentration from the fit files

#Prepare the dataframes to be resorted to match Halo Indices
fit_param_dark = pd.read_csv('50-1_snap_99_fit_param-dark.csv')
fit_param_dark['Halo Number'] = fit_param_dark['Halo Number'].astype(int)
fit_param_dark = fit_param_dark.set_index('Halo Number')
#true_indices_dark = fit_param_dark['Halo Number'].to_numpy().astype(int)

#Import the Fit Parameters for DM+Baryons
#Reorder according to DMO Halo Indices, cut all NaN
fit_param = pd.read_csv('50-1_snap_99_fit_param.csv')
fit_param['Halo Number'] = fit_param['Halo Number'].astype(int)
fit_param = fit_param.set_index('Halo Number')

sorted_data, sorted_Y = data_csv.align(fit_param, join='inner', axis=0)
sorted_data_dark, sorted_Y_dark = data_csv_dark.align(fit_param_dark, join='inner', axis=0)
sorted_data, sorted_data_dark = sorted_data.align(sorted_data_dark, join='inner', axis=0)
sorted_Y, sorted_Y_dark = sorted_Y.align(sorted_Y_dark, join='inner', axis=0)

print(sorted_data)
print(sorted_Y)
print(sorted_data_dark)
print(sorted_Y_dark)


all_snap = np.arange(2,100,1)
to_keep = np.arange(9,100,10)
#to_keep = np.array([99])
to_drop = np.setdiff1d(all_snap, to_keep)


column_drop = []
column_drop_dark = column_drop.copy()
column_keep = []
column_keep_dark = column_keep.copy()
for i in to_drop:
    column_drop.extend([str(i)+'gas_mass',str(i)+'dm_mass',str(i)+'stellar_mass',
                        str(i)+'bh_mass',str(i)+'spinX',str(i)+'spinY',
                        str(i)+'spinZ',str(i)+'vel_dispersion',str(i)+'v_max',
                        str(i)+'bh_dot',str(i)+'sfr',str(i)+'fof_mass',
                        str(i)+'fof_distance'])
    column_drop_dark.extend([str(i)+'dm_mass',str(i)+'spinX',str(i)+'spinY',
                        str(i)+'spinZ',str(i)+'vel_dispersion',str(i)+'v_max'])
for j in np.flipud(to_keep):
    column_keep.extend(['index',str(j)+'gas_mass',str(j)+'dm_mass',str(j)+'stellar_mass',
                        str(j)+'bh_mass',str(j)+'spinX',str(j)+'spinY',
                        str(j)+'spinZ',str(j)+'vel_dispersion',str(j)+'v_max',
                        str(j)+'bh_dot',str(j)+'sfr',str(j)+'fof_mass',
                        str(j)+'fof_distance'])
    column_keep_dark.extend([str(j)+'dm_mass_DMO',str(j)+'spinX_DMO',str(j)+'spinY_DMO',
                        str(j)+'spinZ_DMO',str(j)+'vel_dispersion_DMO',str(j)+'v_max_DMO','index_DMO'])


for x in all_snap:
    column_drop.extend([str(x)+'positionX',str(x)+'positionY',str(x)+'positionZ',
                        str(x)+'halfmass_rad',str(x)+'particle_number'])
    column_drop_dark.extend([str(x)+'positionX',str(x)+'positionY',str(x)+'positionZ',
                             str(x)+'halfmass_rad',str(x)+'particle_number'])
column_keep_ratio = column_keep.copy()
column_keep_ratio.extend(column_keep_dark)

#ToDo: deleta all unnecessary columns to retain required data


nfw_scalerad = sorted_Y['NFW Scale Radius'].to_numpy()
virrad = sorted_Y['Virial Radius'].to_numpy()

nfw_scalerad_dark = sorted_Y_dark['NFW Scale Radius'].to_numpy()
virrad_dark = sorted_Y_dark['Virial Radius'].to_numpy()

sorted_X = sorted_data.drop(column_drop, axis=1)
sorted_X_dark = sorted_data_dark.drop(column_drop_dark, axis=1)
sorted_X_dark = sorted_X_dark.add_suffix('_DMO')
#Calculate concentration from Re-indexed input arrays and set as expected value
concentration = virrad/nfw_scalerad
concentration_dark = virrad_dark/nfw_scalerad_dark

y = concentration
y_dark = concentration_dark
print(y)
print(y_dark)
y = np.log10(y)
y_dark  = np.log10(y_dark)

print(sorted_X)
print(sorted_X_dark)

Xtrain, Xtest, ytrain, ytest = train_test_split(sorted_X, y,
                                                random_state=42)

Xtrain_dark, Xtest_dark, ytrain_dark, ytest_dark = train_test_split(sorted_X_dark, y_dark,
                                                random_state=42)


y_ratio = y/y_dark
#y_dark_ratio = concentration_dark_ratio
#y_conc_ratio = y_ratio/y_dark_ratio

#sorted_X_ratio = sorted_data_ratio.drop(column_drop, axis=1)
#sorted_X_dark_ratio = sorted_data_dark_ratio.drop(column_drop_dark, axis=1)
#sorted_X_dark_ratio = sorted_X_dark_ratio.add_suffix('_DMO')
#Predict the ratio using ML
X_ratio = pd.concat([sorted_X,sorted_X_dark],axis=1)
print(X_ratio)

Xtrain_ratio, Xtest_ratio, ytrain_ratio, ytest_ratio = train_test_split(X_ratio, y_ratio,
                                                random_state=42)


#set up plotting params
fig, axs = plt.subplots(1,4,constrained_layout=True, figsize=(40, 10))

#Train Model for DM+Baryons
model = RandomForestRegressor(n_estimators=1000,n_jobs=50)
model.fit(Xtrain,ytrain)
y_pred = model.predict(Xtest)
importances = model.feature_importances_
quantiles_lower = np.quantile([tree.feature_importances_ for tree in model.estimators_],0.25, axis=0)
quantiles_upper = np.quantile([tree.feature_importances_ for tree in model.estimators_],0.75, axis=0)



print(y)
print(y_pred)
#Plot Predicted vs actual values
im = axs[0].hexbin(ytest,y_pred, gridsize = 70,norm=matplotlib.colors.LogNorm())
axs[0].set_xlabel(r'Log Concentration of Halos')
axs[0].set_ylabel(r'Predicted Log Concentration of Halos')
axs[0].set_xlim(0, 3)
axs[0].set_ylim(0, 3)
axs[0].set_title('Predicted Halo Concentration from Mass Contents, Vmax, VelDisp, Spin, FoF Properties')
cb = fig.colorbar(im)

#Train Model for DMO
model_dark = RandomForestRegressor(n_estimators=1000,n_jobs=50)
model_dark.fit(Xtrain_dark,ytrain_dark)
y_pred_dark = model_dark.predict(Xtest_dark)
importances_dark = model_dark.feature_importances_
#std_dark = np.std([tree_dark.feature_importances_ for tree_dark in model_dark.estimators_], axis=0)
quantiles_lower_dark = np.quantile([tree_dark.feature_importances_ for tree_dark in model_dark.estimators_],0.25, axis=0)
quantiles_upper_dark = np.quantile([tree_dark.feature_importances_ for tree_dark in model_dark.estimators_],0.75, axis=0)

print(y_dark)
print(y_pred_dark)
#Plot predicted vs actual
axs[1].hexbin(ytest_dark,y_pred_dark, gridsize = 70, norm=matplotlib.colors.LogNorm())
axs[1].set_xlabel(r'Log Concentration of DMO Halos')
axs[1].set_ylabel(r'Predicted Log Concentration of DMO Halos')
axs[1].set_xlim(0, 2)
axs[1].set_ylim(0, 2)
axs[1].set_title('Predicted Halo Concentration from Mass Contents, Vmax, VelDisp, Spin, FoF Properties')



model_ratio = RandomForestRegressor(n_estimators=1000,n_jobs=50)
model_ratio.fit(Xtrain_ratio,ytrain_ratio)
y_pred_ratio = model_ratio.predict(Xtest_ratio)
importances_ratio = model_ratio.feature_importances_
#std_ratio = np.std([tree_ratio.feature_importances_ for tree_ratio in model_ratio.estimators_], axis=0)
quantiles_lower_ratio = np.quantile([tree_ratio.feature_importances_ for tree_ratio in model_ratio.estimators_],0.25, axis=0)
quantiles_upper_ratio = np.quantile([tree_ratio.feature_importances_ for tree_ratio in model_ratio.estimators_],0.75, axis=0)


#Plot predicted vs actual
axs[2].hexbin(ytest_ratio,y_pred_ratio, gridsize = 70,norm=matplotlib.colors.LogNorm())
axs[2].set_xlabel(r'Log of Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[2].set_ylabel(r'Predicted Log of Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[2].set_xlim(0, 3)
axs[2].set_ylim(0, 3)
axs[2].set_title('Predicted Halo Concentration ratio from Mass Contents, Vmax, VelDisp, Spin, FoF Properties')

ytest_ratio_calc = ytest/ytest_dark
y_pred_ratio_calc = y_pred/y_pred_dark
axs[3].hexbin(ytest_ratio_calc,y_pred_ratio_calc, gridsize = 70,norm=matplotlib.colors.LogNorm())
axs[3].set_xlabel(r'Log of Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[3].set_ylabel(r'Predicted Log of Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[3].set_xlim(0, 3)
axs[3].set_ylim(0, 3)
axs[3].set_title('Predicted Halo Concentration ratio calculated from two left panels')



print('R_2')
print('FP: '+str(sklearn.metrics.r2_score(ytest, y_pred)))
print('DMO: '+str(sklearn.metrics.r2_score(ytest_dark, y_pred_dark)))
print('Ratio ML: '+str(sklearn.metrics.r2_score(ytest_ratio, y_pred_ratio)))
print('Ratio Calculated: '+str(sklearn.metrics.r2_score(ytest_ratio_calc, y_pred_ratio_calc)))

print('Mean squared error')
print('FP: '+str(sklearn.metrics.mean_squared_error(ytest, y_pred)))
print('DMO: '+str(sklearn.metrics.mean_squared_error(ytest_dark, y_pred_dark)))
print('Ratio ML: '+str(sklearn.metrics.mean_squared_error(ytest_ratio, y_pred_ratio)))
print('Ratio Calculated: '+str(sklearn.metrics.mean_squared_error(ytest_ratio_calc, y_pred_ratio_calc)))
fig.savefig('concentration_ratio_history0801.jpg')


fig.clf()


forest_importances = pd.Series(importances, index=column_keep)
forest_quantiles_lower = pd.Series(quantiles_lower, index=column_keep)
forest_quantiles_upper = pd.Series(quantiles_upper, index=column_keep)


forest_importances_dark = pd.Series(importances_dark, index=column_keep_dark)
forest_quantiles_lower_dark = pd.Series(quantiles_lower_dark, index=column_keep_dark)
forest_quantiles_upper_dark = pd.Series(quantiles_upper_dark, index=column_keep_dark)

forest_importances_ratio = pd.Series(importances_ratio, index=column_keep_ratio)
forest_quantiles_lower_ratio = pd.Series(quantiles_lower_ratio, index=column_keep_ratio)
forest_quantiles_upper_ratio = pd.Series(quantiles_upper_ratio, index=column_keep_ratio)



forest_importances.to_csv('forest-importances.csv')
forest_quantiles_lower.to_csv('forest_quantiles_lower.csv')
forest_quantiles_upper.to_csv('forest_quantiles_upper.csv')
forest_importances_dark.to_csv('forest-importances_dark.csv')
forest_quantiles_lower_dark.to_csv('forest_quantiles_lower_dark.csv')
forest_quantiles_upper_dark.to_csv('forest_quantiles_upper_dark.csv')
forest_importances_ratio.to_csv('forest-importances_ratio.csv')
forest_quantiles_lower_ratio.to_csv('forest_quantiles_lower_ratio.csv')
forest_quantiles_upper_ratio.to_csv('forest_quantiles_upper_ratio.csv')


feature_names=['index','gas_mass','dm_mass','stellar_mass',
                    'bh_mass','spinX','spinY',
                    'spinZ','vel_dispersion','v_max',
                    'bh_dot','sfr','fof_mass',
                    'fof_distance']
feature_names_dark= ['dm_mass_DMO','spinX_DMO','spinY_DMO',
                    'spinZ_DMO','vel_dispersion_DMO','v_max_DMO','index']
feature_names_ratio= ['index','gas_mass','dm_mass','stellar_mass',
                    'bh_mass','spinX','spinY',
                    'spinZ','vel_dispersion','v_max',
                    'bh_dot','sfr','fof_mass',
                    'fof_distance',
                    'dm_mass_DMO','spinX_DMO','spinY_DMO',
                    'spinZ_DMO','vel_dispersion_DMO','v_max_DMO','index_dmo']
fig, axs = plt.subplots(3,3,constrained_layout=True, figsize=(30, 30))
#Plot Predicted vs actual values

cmap = plt.get_cmap('plasma',20)

column_keep = np.array(column_keep)
column_keep_dark = np.array(column_keep_dark)
column_keep_ratio = np.array(column_keep_ratio)
for i in np.arange(0,4,1):
    index = (np.arange(i,len(column_keep),len(feature_names)))
    axs[0,0].plot(to_keep,forest_importances[index], 
                    label=feature_names[i],
                    color=cmap(i))
    axs[0,0].fill_between(to_keep, forest_quantiles_lower[index], forest_quantiles_upper[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

for i in np.arange(4,8,1):
    index = (np.arange(i,len(column_keep),len(feature_names)))
    axs[1,0].plot(to_keep,forest_importances[index], 
                    label=feature_names[i],
                    color=cmap(i))
    axs[1,0].fill_between(to_keep, forest_quantiles_lower[index], forest_quantiles_upper[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

for i in np.arange(8,12,1):
    index = (np.arange(i,len(column_keep),len(feature_names)))
    axs[2,0].plot(to_keep,forest_importances[index], 
                    label=feature_names[i],
                    color=cmap(i))
    axs[2,0].fill_between(to_keep, forest_quantiles_lower[index], forest_quantiles_upper[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

#forest_importances.plot.bar(yerr=std, ax=axs[0])
axs[0,0].set_xlabel(r'Snap Number')
axs[0,0].set_ylabel(r'Mean decrease in impurity')
axs[0,0].set_title('Feature Importance using MDI DM+Baryons')
axs[0,0].legend()
axs[1,0].set_xlabel(r'Snap Number')
axs[1,0].set_ylabel(r'Mean decrease in impurity')
axs[1,0].set_title('Feature Importance using MDI DM+Baryons')
axs[1,0].legend()
axs[2,0].set_xlabel(r'Snap Number')
axs[2,0].set_ylabel(r'Mean decrease in impurity')
axs[2,0].set_title('Feature Importance using MDI DM+Baryons')
axs[2,0].legend()



#Plot predicted vs actual
for i in np.arange(0,2,1):
    index = (np.arange(i,len(column_keep_dark),len(feature_names_dark)))
    axs[0,1].plot(to_keep,forest_importances_dark[index], 
                    label=feature_names_dark[i],
                    color=cmap(i))
    axs[0,1].fill_between(to_keep, forest_quantiles_lower_dark[index], forest_quantiles_upper_dark[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

for i in np.arange(2,4,1):
    index = (np.arange(i,len(column_keep_dark),len(feature_names_dark)))
    axs[1,1].plot(to_keep,forest_importances_dark[index], 
                    label=feature_names_dark[i],
                    color=cmap(i))
    axs[1,1].fill_between(to_keep, forest_quantiles_lower_dark[index], forest_quantiles_upper_dark[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

for i in np.arange(4,5,1):
    index = (np.arange(i,len(column_keep_dark),len(feature_names_dark)))
    axs[2,1].plot(to_keep,forest_importances_dark[index], 
                    label=feature_names_dark[i],
                    color=cmap(i))
    axs[2,1].fill_between(to_keep, forest_quantiles_lower_dark[index], forest_quantiles_upper_dark[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

#forest_importances_dark.plot.bar(yerr=std_dark, ax=axs[1])
axs[0,1].set_xlabel(r'Feature importances using MDI')
axs[0,1].set_ylabel(r'Mean decrease in impurity')
axs[0,1].set_title('Feature Importance DMO')
axs[0,1].legend()
axs[1,1].set_xlabel(r'Feature importances using MDI')
axs[1,1].set_ylabel(r'Mean decrease in impurity')
axs[1,1].set_title('Feature Importance DMO')
axs[1,1].legend()
axs[2,1].set_xlabel(r'Feature importances using MDI')
axs[2,1].set_ylabel(r'Mean decrease in impurity')
axs[2,1].set_title('Feature Importance DMO')
axs[2,1].legend()


#Plot predicted vs actual
for i in np.arange(0,6,1):
    index = (np.arange(i,len(column_keep_ratio),len(feature_names_ratio)))
    axs[0,2].plot(to_keep,forest_importances_ratio[index], 
                    label=feature_names_ratio[i],
                    color=cmap(i))
    axs[0,2].fill_between(to_keep, forest_quantiles_lower_ratio[index], forest_quantiles_upper_ratio[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

for i in np.arange(6,12,1):
    index = (np.arange(i,len(column_keep_ratio),len(feature_names_ratio)))
    axs[1,2].plot(to_keep,forest_importances_ratio[index], 
                    label=feature_names_ratio[i],
                    color=cmap(i))
    axs[1,2].fill_between(to_keep, forest_quantiles_lower_ratio[index], forest_quantiles_upper_ratio[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill

for i in np.arange(12,18,1):
    index = (np.arange(i,len(column_keep_ratio),len(feature_names_ratio)))
    axs[2,2].plot(to_keep,forest_importances_ratio[index], 
                    label=feature_names_ratio[i],
                    color=cmap(i))
    axs[2,2].fill_between(to_keep, forest_quantiles_lower_ratio[index], forest_quantiles_upper_ratio[index],
                 facecolor=cmap(i), # The fill color
                 color=cmap(i),       # The outline color
                 alpha=0.2)          # Transparency of the fill
#forest_importances_ratio.plot.bar(yerr=std_ratio, ax=axs[2])
axs[0,2].set_xlabel(r'Feature importances using MDI')
axs[0,2].set_ylabel(r'Mean decrease in impurity')
axs[0,2].set_title('Feature Importance Ratio')
axs[0,2].legend()
axs[1,2].set_xlabel(r'Feature importances using MDI')
axs[1,2].set_ylabel(r'Mean decrease in impurity')
axs[1,2].set_title('Feature Importance Ratio')
axs[1,2].legend()
axs[2,2].set_xlabel(r'Feature importances using MDI')
axs[2,2].set_ylabel(r'Mean decrease in impurity')
axs[2,2].set_title('Feature Importance Ratio')
axs[2,2].legend()

fig.savefig('feature-importance_history-0801.jpg')