# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 11:23:03 2022

@author: clara
"""
import sys
sys.path.insert(1, '/home/clilje/summer_project')
from illustris_python import groupcat,lhalotree,snapshot,sublink,util
#from localhost.Ubuntu.home.clara.summer_project import illustris_python as il


basePath = './disk01/rmcg/downloaded/tng/tng50-1'
fields = ['SubhaloMass','SubhaloSFRinRad']
subhalos = groupcat.loadSubhalos(basePath,135,fields=fields)
#tree = sublink.loadtre()
subhalos.keys()
subhalos['SubhaloMass'].shape
