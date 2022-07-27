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
#import scikit-learn as sklearn
from sklearn.ensemble import RandomForestRegressor
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


h = 0.6774
p_crit = 127 #m_sun/(kpc^3)

#get the input data from the Group Catalogues
#DM + Baryons
data_csv = pd.read_csv('50-1-subhalo-history.csv')

#DMO
data_csv_dark = pd.read_csv('50-1-subhalo-history-dark.csv')


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
data_csv['Df_cat'] = pd.Categorical(data_csv['index'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data = data_csv.sort_values('Df_cat').dropna().copy()

all_snap = np.arange(2,100,1)
to_keep = np.arange(9,100,10)

to_drop = np.setdiff1d(all_snap, to_keep)

column_drop = ['index','Df_cat']
column_drop_dark = column_drop.copy()
for i in to_drop:
    column_drop.extend([str(i)+'positionX',str(i)+'positionY',str(i)+'positionZ',
                        str(i)+'gas_mass',str(i)+'dm_mass',str(i)+'stellar_mass',
                        str(i)+'bh_mass',str(i)+'spinX',str(i)+'spinY',
                        str(i)+'spinZ',str(i)+'vel_dispersion',str(i)+'v_max',
                        str(i)+'bh_dot',str(i)+'sfr',str(i)+'fof_mass',
                        str(i)+'fof_distance'])
    column_drop_dark.extend([str(i)+'positionX',str(i)+'positionY',str(i)+'positionZ',
                        str(i)+'dm_mass',str(i)+'spinX',str(i)+'spinY',
                        str(i)+'spinZ',str(i)+'vel_dispersion',str(i)+'v_max'])

for x in all_snap:
    column_drop.extend([str(x)+'halfmass_rad',str(x)+'particle_number'])
    column_drop_dark.extend([str(x)+'halfmass_rad',str(x)+'particle_number'])
#ToDo: deleta all unnecessary columns to retain required data

sorted_X = sorted_data.drop(column_drop, axis=1)

data_csv_dark['Df_cat'] = pd.Categorical(data_csv_dark['index'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data_dark = data_csv_dark.sort_values('Df_cat').dropna().copy()

sorted_X_dark = sorted_data_dark.drop(column_drop_dark, axis=1)
sorted_X_dark = sorted_X_dark.add_suffix('_DMO')
#Calculate concentration from Re-indexed input arrays and set as expected value
concentration = virrad/nfw_scalerad
concentration_dark = virrad_dark/nfw_scalerad_dark

y = concentration
y_dark = concentration_dark
print(y)
print(y_dark)


print(sorted_X)
print(sorted_X_dark)

Xtrain, Xtest, ytrain, ytest = train_test_split(sorted_X, y,
                                                random_state=1)

Xtrain_dark, Xtest_dark, ytrain_dark, ytest_dark = train_test_split(sorted_X_dark, y_dark,
                                                random_state=1)

#Calculate the C_Bar/C_DMO ratio
y_conc_ratio = y/y_dark

#Predict the ratio using ML
X_ratio = pd.concat([sorted_X,sorted_X_dark])

Xtrain_ratio, Xtest_ratio, ytrain_ratio, ytest_ratio = train_test_split(X_ratio, y_conc_ratio,
                                                random_state=1)



#set up plotting params
fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(30, 10))

#Train Model for DM+Baryons
model = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model.fit(Xtrain,ytrain)
y_pred = model.predict(Xtest)
importances = model.feature_importances_
std = np.std([tree.feature_importances_ for tree in model.estimators_], axis=0)

print(y)
print(y_pred)
#Plot Predicted vs actual values
im = axs[0].hexbin(ytest,y_pred, gridsize = 70,xscale ='log',yscale='log',norm=matplotlib.colors.LogNorm())
axs[0].set_xlabel(r'Concentration of Halos')
axs[0].set_ylabel(r'Predicted Concentration of Halos')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlim(3*10**0, 3*10)
axs[0].set_ylim(3*10**0, 3*10)
axs[0].set_title('Predicted Halo Concentration from Mass Contents, Vmax, VelDisp, Spin, FoF Properties')
cb = fig.colorbar(im)

#Train Model for DMO
model_dark = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model_dark.fit(Xtrain_dark,ytrain_dark)
y_pred_dark = model_dark.predict(Xtest_dark)
importances_dark = model_dark.feature_importances_
std_dark = np.std([tree_dark.feature_importances_ for tree_dark in model_dark.estimators_], axis=0)

print(y_dark)
print(y_pred_dark)
#Plot predicted vs actual
axs[1].hexbin(ytest_dark,y_pred_dark, gridsize = 70,xscale ='log',yscale='log',norm=matplotlib.colors.LogNorm())
axs[1].set_xlabel(r'Concentration of DMO Halos')
axs[1].set_ylabel(r'Predicted Concentration of DMO Halos')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_xlim(4*10**0, 2*10)
axs[1].set_ylim(4*10**0, 2*10)
axs[1].set_title('Predicted Halo Concentration from Mass Contents, Vmax, VelDisp, Spin, FoF Properties')



model_ratio = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model_ratio.fit(Xtrain_ratio,ytrain_ratio)
y_pred_ratio = model_ratio.predict(Xtest_ratio)
importances_ratio = model_ratio.feature_importances_
std_ratio = np.std([tree_ratio.feature_importances_ for tree_ratio in model_ratio.estimators_], axis=0)


#Plot predicted vs actual
plt.hexbin(ytest_ratio,y_pred_ratio, gridsize = 70,xscale ='log',yscale='log',norm=matplotlib.colors.LogNorm())
axs[2].set_xlabel(r'Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[2].set_ylabel(r'Predicted Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].set_xlim(3*10**(-1), 3*10**0)
axs[2].set_ylim(3*10**(-1), 3*10**0)
axs[2].set_title('Predicted Halo Concentration ratio from Mass Contents, Vmax, VelDisp, Spin, FoF Properties')
fig.savefig('concentration_ratio_fof-final.jpg')


fig.clf()

"""
forest_importances = pd.Series(importances, index=['SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter'])

forest_importances_dark = pd.Series(importances_dark, index=['SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'FoFMass','FoFDistanceCenter'])

forest_importances_ratio = pd.Series(importances_ratio, index=['SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter',
                'SubhaloDMMass - DMO','SubhaloSpinX- DMO','SubhaloSpinY- DMO','SubhaloSpinZ- DMO','SubhaloVelDisp- DMO', 
                                     'SubhaloVmax- DMO','FoFMass- DMO','FoFDistanceCenter- DMO'])


fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(30, 10))
#Plot Predicted vs actual values
forest_importances.plot.bar(yerr=std, ax=axs[0])
axs[0].set_xlabel(r'Feature importances using MDI')
axs[0].set_ylabel(r'Mean decrease in impurity')
axs[0].set_title('Feature Importance DM+Baryons')



#Plot predicted vs actual
forest_importances_dark.plot.bar(yerr=std_dark, ax=axs[1])
axs[1].set_xlabel(r'Feature importances using MDI')
axs[1].set_ylabel(r'Mean decrease in impurity')
axs[1].set_title('Feature Importance DMO')

#Plot predicted vs actual
forest_importances_ratio.plot.bar(yerr=std_ratio, ax=axs[2])
axs[2].set_xlabel(r'Feature importances using MDI')
axs[2].set_ylabel(r'Mean decrease in impurity')
axs[2].set_title('Feature Importance Ratio')
fig.savefig('feature-importance_fof-final.jpg')

"""