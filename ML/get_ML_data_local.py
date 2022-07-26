# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:51:11 2022

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
data_csv = pd.read_csv('50-1-subhalo-info.csv',usecols=['SubhaloIndex','SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter'])
column_names = ['SubhaloIndex','SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass',
                'SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 'SubhaloVmax',
                'SubhaloBHMdot','SubhaloSFR','FoFMass','FoFDistanceCenter']
X = data_csv[column_names]

#DMO
data_csv_dark = pd.read_csv('50-1-subhalo-info-dark.csv',usecols=['SubhaloIndex','SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 
                     'SubhaloVmax','FoFMass','FoFDistanceCenter'])
column_names_dark = ['SubhaloIndex','SubhaloDMMass','SubhaloSpinX','SubhaloSpinY','SubhaloSpinZ','SubhaloVelDisp', 
                     'SubhaloVmax','FoFMass','FoFDistanceCenter']
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
sorted_X = pd.DataFrame([sorted_data['SubhaloGasMass'],sorted_data['SubhaloStarMass'],
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
sorted_X_dark = pd.DataFrame([sorted_data_dark['SubhaloDMMass'],
                         sorted_data_dark['SubhaloSpinX'],sorted_data_dark['SubhaloSpinY'],
                         sorted_data_dark['SubhaloSpinZ'],sorted_data_dark['SubhaloVelDisp'],
                         sorted_data_dark['SubhaloVmax'],sorted_data_dark['FoFMass'],
                         sorted_data_dark['FoFDistanceCenter']]).T


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
model = RandomForestRegressor(n_estimators=1000,n_jobs=50)
model.fit(Xtrain,ytrain)
y_pred = model.predict(Xtest)
print(y)
print(y_pred)
#Train Model for DMO
model_dark = RandomForestRegressor(n_estimators=1000,n_jobs=50)
model_dark.fit(Xtrain_dark,ytrain_dark)
y_pred_dark = model_dark.predict(Xtest_dark)

data = pd.DataFrame([Xtest,y_pred,Xtest_dark,y_pred_dark], columns=['Xtest','y_predicted','Xtest_dark','y_predicted_dark'])
data.to_csv('ML_data')

