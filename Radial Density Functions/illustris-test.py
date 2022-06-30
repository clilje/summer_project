# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:23:03 2022

@author: clara
"""
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
sys.path.insert(1, '/home/clilje/summer_project')
from illustris_python import groupcat,lhalotree,snapshot,sublink,util
#from localhost.Ubuntu.home.clara.summer_project import illustris_python as il
#sys.path.insert(1, '/')

basePath = '/disk01/rmcg/downloaded/tng/tng50-1'
fields = ['SubhaloMass','SubhaloSFRinRad']
subhalos = groupcat.loadSubhalos(basePath,99,fields=fields)
#tree = sublink.loadtre()
print(subhalos.keys())
print(subhalos['SubhaloMass'].shape)

dm_pos = snapshot.loadSubset(basePath,99,'dm',['Coordinates']);

plt.hist2d(dm_pos[:,0], dm_pos[:,1], norm=mpl.colors.LogNorm(), bins=64);

plt.xlim([0,75000])

plt.ylim([0,75000])

plt.xlabel('x [ckpc/h]')

plt.ylabel('y [ckpc/h]')

plt.savefig('dmhistogram')