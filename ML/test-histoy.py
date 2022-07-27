# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 11:44:08 2022

@author: clara
"""

import os

import h5py
import numpy as np
import pandas as pd

data_csv = pd.read_csv('50-1-subhalo-history.csv')
print(data_csv.loc[data_csv['index']==0])