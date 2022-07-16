#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:28:45 2022

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

class LGR_Class:
    def __init__(self,n,last_section):
        #n : number of collotaion points >=1
        #last section: True or False
        num_intermediate_LGR_nodes=n-1;
        OC=OrthColloc(num_intermediate_LGR_nodes,5,Shift=False,Geometry=0)
        self.LGR_Nodes=OC.Xpoints()[0:-1]
        self.LGR_Weights=OC.WeightQ()[0:-1]
        self.LGR_Diff_Matrix=self.Compute_D(last_section)
        
    def Compute_D(self,last_section):
        if last_section==False:
            Regression_points=self.LGR_Nodes
        elif last_section==True:
            Regression_points=np.append(self.LGR_Nodes,1.0)
        D=np.zeros((len(self.LGR_Nodes),len(Regression_points)))
        for i in np.arange(len(Regression_points)):
            index_vec=np.zeros(len(Regression_points))
            index_vec[i]=1.0
            Li=lagrange(Regression_points,index_vec)
            #Li_dot=[Li[j] * j for j in range(1, len(Li)+1)]
            Li_dot=[Li[j] * j for j in range(len(Li),0,-1)]
            for k in np.arange(len(self.LGR_Nodes)):
                D[k,i]=np.polyval(Li_dot,self.LGR_Nodes[k])
        return D

                
if __name__ == '__main__':  
    cl=LGR_Class(5,False)      
    print(cl.LGR_Nodes)
    print(cl.LGR_Weights)
    print(cl.LGR_Diff_Matrix)