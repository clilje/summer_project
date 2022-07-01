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
data_csv = pd.read_csv('50-4_mass_difference.csv')
plt.plot(data_csv['halo_number'], data_csv['mass_difference']/data_csv['mass_tng'], "+", color="black")
plt.xlabel(r'Halo Index$')
plt.ylabel(r'Difference in Mass/ Mass given from TNG')
plt.savefig('Halfmassdifference')
plt.show()