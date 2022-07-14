# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 12:20:24 2022
plot halo positions
@author: clara
"""

import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits import mplot3d
#sys.path.insert(1, '/home/clilje/summer_project')
#from illustris_python import groupcat,lhalotree,snapshot,sublink,util
import math
import csv
import pandas as pd
import scipy.optimize as scopt
import scipy.linalg
import scipy.stats


data_csv = pd.read_csv('50-1-subhalo-info.csv')
#print(data_csv['SubhaloPosX'])
fig = plt.figure()
 
# syntax for 3-D projection
ax = plt.axes(projection ='3d')

ax.scatter(data_csv['SubhaloPosX'].to_numpy()[0:1500], data_csv['SubhaloPosY'].to_numpy()[0:1500], data_csv['SubhaloPosZ'].to_numpy()[0:1500], marker='+')

ax.set_xlabel('x [ckpc/h]')

ax.set_ylabel('y [ckpc/h]')
ax.set_zlabel('z [ckpc/h]')
fig.savefig('halopos')
