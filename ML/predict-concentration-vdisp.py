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
data_csv = pd.read_csv('50-1-subhalo-info.csv')
column_names = ['SubhaloIndex','SubhaloVelDisp']
X = data_csv[column_names]

#DMO
data_csv_dark = pd.read_csv('50-1-subhalo-info-dark.csv')
column_names_dark = ['SubhaloIndex','SubhaloVelDisp']
X_dark = data_csv_dark[column_names_dark]


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
sorted_X = sorted_data['SubhaloVelDisp'].to_numpy()

print(X_dark)

X_dark['Df_cat'] = pd.Categorical(X_dark['SubhaloIndex'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data_dark = X_dark.sort_values('Df_cat').dropna().copy()
print(sorted_data_dark)
sorted_X_dark = sorted_data_dark['SubhaloVelDisp'].to_numpy()


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

#set up plotting params
fig, axs = plt.subplots(1,3,constrained_layout=True, figsize=(30, 10))

#Train Model for DM+Baryons
model = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model.fit(Xtrain.reshape(-1,1),ytrain)
y_pred = model.predict(Xtest.reshape(-1,1))
print(y)
print(y_pred)
#Plot Predicted vs actual values
axs[0].scatter(ytest,y_pred, marker="x",color="black")
axs[0].set_xlabel(r'Concentration of Halos')
axs[0].set_ylabel(r'Predicted Concentration of Halos')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_title('Predicted Halo Concentration from VelDisp')

#Train Model for DMO
model_dark = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model.fit(Xtrain_dark.reshape(-1,1),ytrain_dark)
y_pred_dark = model.predict(Xtest_dark.reshape(-1,1))
print(y_dark)
print(y_pred_dark)
#Plot predicted vs actual
axs[1].scatter(ytest_dark,y_pred_dark, marker="x",color="black")
axs[1].set_xlabel(r'Concentration of DMO Halos')
axs[1].set_ylabel(r'Predicted Concentration of DMO Halos')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_title('Prediced Halo Concentration from VelDisp')


#Calculate the C_Bar/C_DMO ratio
conc_ratio = ytest/ytest_dark
conc_ratio_pred = y_pred/y_pred_dark
#Plot predicted vs actual
plt.scatter(conc_ratio,conc_ratio_pred, marker="x",color="black")
axs[2].set_xlabel(r'Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[2].set_ylabel(r'Predicted Ratio of $\frac{C_{B}}{C_{DMO}}$')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].set_title('Predicted Halo Concentration ratio from VelDisp')
fig.savefig('concentration_ratio_vdisp.jpg')

