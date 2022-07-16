#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:05:04 2022

@author: bayat2
"""
import numpy as np
#import math
import os
import sys
from scipy.interpolate import lagrange
Current_path=os.path.abspath(os.getcwd())
Collocation_Path=Current_path+'/Collocation_Librray';
sys.path.insert(1, Collocation_Path)
#import arrayprint as ap
#import occ
from occ import OrthColloc

#Generic = 0
#Gauss = 1
#Lobatto = 2
#Chebyshev = 3
#RadauR = 4
#RadauL = 5
#Cheb1 = 6

LGR_Nodes=OrthColloc(5,5,Shift=False,Geometry=0).Xpoints()[0:-1]
LG_Nodes=OrthColloc(5,1,Shift=False,Geometry=0).Xpoints()[0:-1]