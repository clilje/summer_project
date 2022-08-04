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
data_csv =data_csv.set_index('index', drop=True)
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
#data_csv_dark['index'] = data_csv_dark.index
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
    column_keep.extend([str(j)+'gas_mass',str(j)+'dm_mass',str(j)+'stellar_mass',
                        str(j)+'bh_mass',str(j)+'spinX',str(j)+'spinY',
                        str(j)+'spinZ',str(j)+'vel_dispersion',str(j)+'v_max',
                        str(j)+'bh_dot',str(j)+'sfr',str(j)+'fof_mass',
                        str(j)+'fof_distance'])
    column_keep_dark.extend([str(j)+'dm_mass_DMO',str(j)+'spinX_DMO',str(j)+'spinY_DMO',
                        str(j)+'spinZ_DMO',str(j)+'vel_dispersion_DMO',str(j)+'v_max_DMO'])


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

FP_number = []
DMO_number = []
number_halos_FP = len(sorted_data.index.to_numpy())
number_halos_DMO = len(sorted_data_dark.index.to_numpy())
for x in to_keep:
    FP_number.append(len(np.where(sorted_X[str(x)+'dm_mass'].to_numpy() != 0)[0]))
    DMO_number.append(len(len(np.where(sorted_X_dark[str(x)+'dm_mass_DMO'].to_numpy() != 0)[0])))
                      
FP_number = np.array(FP_number)/number_halos_FP
DMO_number = np.array(DMO_number)/number_halos_DMO
plt.loglog(to_keep, FP_number)
plt.loglog(to_keep, DMO_number)
plt.xlabel('Snapshot number')
plt.ylabel('Percentage of halos traced back')
plt.savefig('halo-percentage')