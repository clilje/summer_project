 # -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 23:17:55 2022

@author: clara
"""

#import h5py
import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import pandas as pd
from pathlib import Path


h = 0.6774

g = 0
data = []
rowlabels= []
while g < 29:
    filename = 'HaloFits/50-4_snap_99_halo_'+str(g)+'_fit_param.csv'
    path = Path('HaloFits/50-4_snap_99_halo_'+str(g)+'_fit_param.csv')
    if path.is_file():
        data_csv = pd.read_csv('HaloFits/50-4_snap_99_halo_'+str(g)+'_fit_param.csv')
        #print(data_csv)
        #print(data_csv.loc[[0]])
        #print(data_csv.to_numpy())
        data.append(data_csv.to_numpy()[0][1:])
        #print(data)
        rowlabels.append(data_csv.to_numpy()[0][0])
        columnlabels = data_csv.columns.values[1:]
    g +=1       

data = np.array(data)
print(data)
print(rowlabels)
print(columnlabels)


#rowlabels = data_csv.columns.values        
#data_csv = pd.read_csv('50-4_mass_difference.csv')
#rowlabels = np.loadtxt(filename)[1][0]
#column_labels = np.loadtxt(filename)[0]
fig, ax = plt.subplots(5,1,figsize=(30,12))
'''data=[[1,2,3],
      [5,6,7],
      [8,9,10]]
'''
#column_labels=["Column 1", "Column 2", "Column 3"]
df=pd.DataFrame(data,columns=columnlabels)
#ax.axis('tight')
ax[0].axis('off')
ax[0].table(cellText=df.values[0:,0:6],
        colLabels=df.columns[0:6],
        rowLabels=rowlabels,
        rowColours =["pink"] * 6,  
        colColours =["pink"] * 6,
        loc="center",
        fontsize=50)

ax[1].axis('off')
ax[1].table(cellText=df.values[0:,6:12],
        colLabels=df.columns[6:12],
        rowLabels=rowlabels,
        rowColours =["pink"] * 6,  
        colColours =["pink"] * 6,
        loc="center",
        fontsize=50)

ax[2].axis('off')
ax[2].table(cellText=df.values[0:,12:18],
        colLabels=df.columns[12:18],
        rowLabels=rowlabels,
        rowColours =["pink"] * 6,  
        colColours =["pink"] * 6,
        loc="center",
        fontsize=50)

print(df.values[0:,18:26])
ax[3].axis('off')
ax[3].table(cellText=df.values[0:,18:26],
        colLabels=df.columns[18:26],
        rowLabels=rowlabels,
        rowColours =["pink"] * 8,  
        colColours =["pink"] * 8,
        loc="center",
        fontsize=200)

ax[4].axis('off')
ax[4].table(cellText=df.values[0:,26:34],
        colLabels=df.columns[26:34],
        rowLabels=rowlabels,
        rowColours =["pink"] * 8,  
        colColours =["pink"] * 8,
        loc="center")
plt.show()
fig.savefig('tablefit')