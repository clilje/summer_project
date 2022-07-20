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

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)


data_csv = pd.read_csv('50-1-subhalo-info.csv')
column_names = ['SubhaloIndex','SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass','SubhaloDMMass']
X = data_csv[column_names]

data_csv_dark = pd.read_csv('50-1-subhalo-info-dark.csv')
column_names_dark = ['SubhaloIndex','SubhaloDMMass']
X_dark = data_csv_dark[column_names_dark]



fit_param_dark = pd.read_csv('50-1_snap_99_fit_param-dark.csv')
true_indices_dark = fit_param_dark['Halo Number'].to_numpy().astype(int)


fit_param = pd.read_csv('50-1_snap_99_fit_param.csv')
fit_param['Df_cat'] = pd.Categorical(fit_param['Halo Number'],
                                             categories = true_indices_dark,
                                             ordered=True)
sorted_df = fit_param.sort_values('Df_cat').dropna()
nfw_scalerad = sorted_df['NFW Scale Radius'].to_numpy()
virrad = sorted_df['Virial Radius'].to_numpy()
true_indices = sorted_df['Halo Number'].to_numpy().astype(int)

fit_param_dark['Df_cat'] = pd.Categorical(fit_param_dark['Halo Number'],
                                             categories = true_indices,
                                             ordered=True)
sorted_df_dark = fit_param_dark.sort_values('Df_cat').dropna()
nfw_scalerad_dark = sorted_df_dark['NFW Scale Radius'].to_numpy()
virrad_dark = sorted_df_dark['Virial Radius'].to_numpy()


print(X)
X['Df_cat'] = pd.Categorical(X['SubhaloIndex'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data = X.sort_values('Df_cat').dropna().copy()
print(sorted_data)
print(sorted_data['Df_cat'])
print(sorted_data['SubhaloGasMass'])
sorted_X = sorted_data.loc['SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass','SubhaloDMMass']
print(sorted_X)

X_dark['Df_cat'] = pd.Categorical(X_dark['SubhaloIndex'],
                                             categories = true_indices,
                                             ordered=True)
sorted_data_dark = X_dark.sort_values('Df_cat').dropna().copy()
sorted_X_dark = sorted_data_dark.loc['SubhaloDMMass']


concentration = virrad/nfw_scalerad
concentration_dark = virrad_dark/nfw_scalerad_dark

y = concentration
y_dark = concentration_dark


fig, axs = plt.subplots(3,constrained_layout=True, figsize=(10, 30))
model = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model.fit(X,y)
data_csv['predicted'] = model.predict(X)
print(y)
print(data_csv['predicted'])
axs[0].set_xlabel(r'Concentration of Halos')
axs[0].set_ylabel(r'Predicted Concentration of Halos')
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_title('Prediced Halo Concentration from Stellar, Gas, BH and DM Mass')

model_dark = RandomForestRegressor(n_estimators=1000,n_jobs=10)
model_dark.fit(X_dark,y_dark)
data_csv_dark['predicted'] = model_dark.predict(X)
print(y_dark)
print(data_csv_dark['predicted'])
axs[1].scatter(y_dark,data_csv_dark['predicted'], marker="x",color="black")
axs[1].set_xlabel(r'Concentration of DMO Halos')
axs[1].set_ylabel(r'Predicted Concentration of DMO Halos')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
axs[1].set_title('Prediced Halo Concentration from Stellar, Gas, BH and DM Mass')

conc_ratio = y/y_dark
conc_ratio_pred = data_csv['predicted']/data_csv_dark['predicted']
plt.scatter(conc_ratio,conc_ratio_pred, marker="x",color="black")
axs[2].set_xlabel(r'Concentration of DMO Halos')
axs[2].set_ylabel(r'Predicted Concentration of DMO Halos')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
fig.savefig('concentation_ratio_forest.jpg')

