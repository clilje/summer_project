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
from sklearn.linear_model import LinearRegression

h = 0.6774
p_crit = 127 #m_sun/(kpc^3)

data_csv = pd.read_csv('50-1-subhalo-info.csv')
print(data_csv.columns.values)
column_names = ['SubhaloGasMass', 'SubhaloStarMass','SubhaloBHMass']
y = data_csv['SubhaloDMMass']
X = data_csv[column_names]
model = LinearRegression(fit_intercept=False)
model.fit(X,y)
data_csv['predicted'] = model.predict(X)
plt.scatter(y*(10**10)*h,data_csv['predicted']*(10**10)*h)
plt.xlabel(r'DM Mass of Halos in$10^{10} M_{\odot}$')
plt.ylabel(r'predicted DM Mass of Halos in$10^{10} M_{\odot}$')
plt.savefig('massprediction.jpg')
plt.show()
