#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:28:45 2022

@author: bayat2
"""

import numpy as np
import math
import os
import sys
from LGR_Data import LGR_Class
from cyipopt import minimize_ipopt
from scipy.optimize import rosen, rosen_der
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib
import pickle as pkl
from matplotlib.font_manager import FontProperties
fontP=FontProperties()
fontP.set_size('x-large')
matplotlib.rcParams['interactive'] == True

class MPC_Class:
    def __init__ (self,Model_data,MPC_params,Mesh,Functions,optionsdict):
        self.Model_data=Model_data
        self.MPC_params=MPC_params
        self.Mesh=Mesh
        self.Functions=Functions
        self.optionsdict=optionsdict
        
        self.I={}
        self.Data_I_Itertations_init()
        #self.Bound(0)
        #self.Guess(0)
        #self.Run_MPC_iter(iter_num=0)
        #self.Run_MPC_Loop()
        k=1
        
        if self.Model_data["plot_iter_flag"]==True:
            self.fig_state, self.ax_state = plt.subplots()
            self.fig_control, self.ax_control = plt.subplots()
            
    
    
    def Data_I_Itertations_init(self):
        self.I["num_ietartion"]=round(self.Model_data["tf"]/self.MPC_params["Ts"])
        self.I["opt_success_flag"]=np.zeros(self.I["num_ietartion"])
        self.I["num_s"]=len(self.Model_data["A_scaled"])
        self.I["num_u"]=len(self.Model_data["C_scaled"])
        self.I["time_h"]=np.arange(self.Model_data["t0"],
                                   self.Model_data["tf"]+0.01*self.MPC_params["Ts"]/self.Mesh["n_h_Sol"],
                                   self.MPC_params["Ts"]/self.Mesh["n_h_Sol"])
        self.I["num_time_h"]=len(self.I["time_h"])
        self.I["s_sc_matrix"]=np.zeros((self.I["num_time_h"],self.I["num_s"]))
        self.I["s_unsc_matrix"]=np.zeros((self.I["num_time_h"],self.I["num_s"]))
        self.I["u_sc_matrix"]=np.zeros((self.I["num_time_h"],self.I["num_u"]))
        self.I["u_unsc_matrix"]=np.zeros((self.I["num_time_h"],self.I["num_u"]))
        self.I["tm_iters"]=np.zeros((self.I["num_ietartion"],2))
        self.I["tp_iters"]=np.zeros((self.I["num_ietartion"],2))
        self.I["num_Ts_in_p_iters"]=np.zeros(self.I["num_ietartion"])
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["tm_iters"][enum_iter,0]=enum_iter*self.MPC_params["Ts"]
            self.I["tp_iters"][enum_iter,0]=enum_iter*self.MPC_params["Ts"]
            self.I["tm_iters"][enum_iter,1]=min(self.Model_data["tf"],
                                                self.I["tm_iters"][enum_iter,0]+
                                                self.MPC_params["m"]*self.MPC_params["Ts"])
            self.I["tp_iters"][enum_iter,1]=min(self.Model_data["tf"],
                                                self.I["tp_iters"][enum_iter,0]+
                                                self.MPC_params["p"]*self.MPC_params["Ts"])
            self.I["num_Ts_in_p_iters"][enum_iter]=round((self.I["tp_iters"][enum_iter,1]-
                                                          self.I["tp_iters"][enum_iter,0])/
                                                         self.MPC_params["Ts"])
        
        
        self.I["time_m_iters"]=[None] * self.I["num_ietartion"]
        self.I["time_p_iters"]=[None] * self.I["num_ietartion"]
        self.I["num_time_m_iters"]=np.zeros(self.I["num_ietartion"])
        self.I["num_time_p_iters"]=np.zeros(self.I["num_ietartion"])
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["time_m_iters"][enum_iter]=np.arange(self.I["tm_iters"][enum_iter,0], 
                                                        self.I["tm_iters"][enum_iter,1]+
                                                        0.01*self.MPC_params["Ts"],
                                                        self.MPC_params["Ts"])
            self.I["time_p_iters"][enum_iter]=np.arange(self.I["tp_iters"][enum_iter,0], 
                                                        self.I["tp_iters"][enum_iter,1]+
                                                        0.01*self.MPC_params["Ts"],
                                                        self.MPC_params["Ts"])
            
            self.I["num_time_m_iters"][enum_iter]=len(self.I["time_m_iters"][enum_iter])
            self.I["num_time_p_iters"][enum_iter]=len(self.I["time_p_iters"][enum_iter])
        
        self.I["num_seg_iters"]=np.zeros(self.I["num_ietartion"],dtype=int)
        self.I["num_nodes_each_seg_iters"]=np.zeros(self.I["num_ietartion"],dtype=int)
        self.I["num_All_nodes_each_seg_iters"]=np.zeros(self.I["num_ietartion"],dtype=int)
        for enum_iter in np.arange(self.I["num_ietartion"]):
            #self.I["num_seg_iters"][enum_iter]=int(self.I["num_time_p_iters"][enum_iter]/20+1)
            #self.I["num_nodes_each_seg_iters"][enum_iter]=min(self.I["num_time_p_iters"][enum_iter]/self.I["num_seg_iters"][enum_iter],7)
            #self.I["num_nodes_each_seg_iters"][enum_iter]=min(self.I["num_time_p_iters"][enum_iter]/self.I["num_seg_iters"][enum_iter]/2,3)
            self.I["num_seg_iters"][enum_iter]=max(1,round((self.I["tp_iters"][enum_iter,1]-self.I["tp_iters"][enum_iter,0])/self.Mesh["delta_t_seg"]))#4
            self.I["num_nodes_each_seg_iters"][enum_iter]=self.Mesh["num_nodes_seg"]#10
            self.I["num_All_nodes_each_seg_iters"][enum_iter]=self.I["num_nodes_each_seg_iters"][enum_iter]+1
            
        self.I["coll_indices_at_each_segment_iters"]=[None] * self.I["num_ietartion"]
        self.I["coll_indices_at_each_segment_iters_vec"]=[None] * self.I["num_ietartion"]
        self.I["All_indices_at_each_segment_iters"]=[None] * self.I["num_ietartion"]
        self.I["Endpoint_indices_at_each_segment_iters"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            num_nodes=self.I["num_nodes_each_seg_iters"][enum_iter]
            num_All_nodes=self.I["num_All_nodes_each_seg_iters"][enum_iter]
            self.I["coll_indices_at_each_segment_iters"][enum_iter]=dict()
            self.I["coll_indices_at_each_segment_iters_vec"][enum_iter]=list()
            self.I["Endpoint_indices_at_each_segment_iters"][enum_iter]=dict()
            self.I["All_indices_at_each_segment_iters"][enum_iter]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][enum_iter]):
                start_indx=enum_segment*num_All_nodes
                end_indx=start_indx+num_nodes#(enum_segment+1)*num_nodes
                end_indx_All_nodes=start_indx+num_All_nodes#(enum_segment+1)*num_nodes
                self.I["coll_indices_at_each_segment_iters"][enum_iter][enum_segment]=np.arange(start_indx,end_indx,1,dtype=int)
                self.I["All_indices_at_each_segment_iters"][enum_iter][enum_segment]=np.arange(start_indx,end_indx_All_nodes,1,dtype=int)
                self.I["Endpoint_indices_at_each_segment_iters"][enum_iter][enum_segment]=end_indx
                self.I["coll_indices_at_each_segment_iters_vec"][enum_iter].append(self.I["coll_indices_at_each_segment_iters"][enum_iter][enum_segment])
            
        
        self.I["time_unscaled_iters_array"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["time_unscaled_iters_array"][enum_iter]=np.linspace(self.I["tp_iters"][enum_iter,0],
                                                                       self.I["tp_iters"][enum_iter,1],
                                                                       self.I["num_seg_iters"][enum_iter]+1)
        
        self.I["Includes_tf_iter_flag"]=np.zeros(self.I["num_ietartion"])
        self.I["tau_iters"]=[None] * self.I["num_ietartion"]
        self.I["tau_iters_All_nodes"]=[None] * self.I["num_ietartion"]
        self.I["tau_iters_unscaled"]=[None] * self.I["num_ietartion"]
        self.I["tau_iters_All_nodes_unscaled"]=[None] * self.I["num_ietartion"]
        self.I["W_iters"]=[None] * self.I["num_ietartion"]
        self.I["D_iters"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            last_section=False
            self.I["tau_iters"][enum_iter]=dict()
            self.I["tau_iters_All_nodes"][enum_iter]=dict()
            self.I["W_iters"][enum_iter]=dict()
            self.I["D_iters"][enum_iter]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][enum_iter]):
                if self.I["tp_iters"][enum_iter,1]==self.Model_data["tf"]: #includes tf
                    if enum_segment==self.I["num_seg_iters"][enum_iter]-1: #last segment
                        #last_section=True
                        self.I["Includes_tf_iter_flag"][enum_iter]=1
                collocation_data=LGR_Class(self.I["num_nodes_each_seg_iters"][enum_iter],True)
                self.I["tau_iters"][enum_iter][enum_segment]=collocation_data.LGR_Nodes
                self.I["W_iters"][enum_iter][enum_segment]=collocation_data.LGR_Weights
                self.I["D_iters"][enum_iter][enum_segment]=collocation_data.LGR_Diff_Matrix
                self.I["tau_iters_All_nodes"][enum_iter][enum_segment]=np.append(collocation_data.LGR_Nodes,1.0)
                
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["tau_iters_unscaled"][enum_iter]=dict()
            self.I["tau_iters_All_nodes_unscaled"][enum_iter]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][enum_iter]):
                self.I["tau_iters_unscaled"][enum_iter][enum_segment]=((1-self.I["tau_iters"][enum_iter][enum_segment])*\
                                                                        self.I["time_unscaled_iters_array"][enum_iter][enum_segment]+\
                                                                      (1+self.I["tau_iters"][enum_iter][enum_segment])*\
                                                                        self.I["time_unscaled_iters_array"][enum_iter][enum_segment+1])/2     
                self.I["tau_iters_All_nodes_unscaled"][enum_iter][enum_segment]=((1-self.I["tau_iters_All_nodes"][enum_iter][enum_segment])*\
                                                                        self.I["time_unscaled_iters_array"][enum_iter][enum_segment]+\
                                                                      (1+self.I["tau_iters_All_nodes"][enum_iter][enum_segment])*\
                                                                        self.I["time_unscaled_iters_array"][enum_iter][enum_segment+1])/2   
        
        self.I["num_coll_nodes_iters"]=np.zeros(self.I["num_ietartion"],dtype=int)
        self.I["num_All_nodes_iters"]=np.zeros(self.I["num_ietartion"],dtype=int)
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["num_coll_nodes_iters"][enum_iter]=self.I["num_seg_iters"][enum_iter]*self.I["num_nodes_each_seg_iters"][enum_iter]
            self.I["num_All_nodes_iters"][enum_iter]=self.I["num_seg_iters"][enum_iter]*self.I["num_All_nodes_each_seg_iters"][enum_iter]
            
        self.I["tau_iters_unscaled_vec"]=[None] * self.I["num_ietartion"]   
        self.I["tau_iters_All_nodes_unscaled_vec"]=[None] * self.I["num_ietartion"]   
        self.I["tau_iters_vec"]=[None] * self.I["num_ietartion"]  
        self.I["tau_iters_All_nodes_vec"]=[None] * self.I["num_ietartion"]  
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["tau_iters_unscaled_vec"][enum_iter]=np.zeros(self.I["num_coll_nodes_iters"][enum_iter])
            self.I["tau_iters_All_nodes_unscaled_vec"][enum_iter]=np.zeros(self.I["num_All_nodes_iters"][enum_iter])
            self.I["tau_iters_vec"][enum_iter]=np.zeros(self.I["num_coll_nodes_iters"][enum_iter])
            self.I["tau_iters_All_nodes_vec"][enum_iter]=np.zeros(self.I["num_All_nodes_iters"][enum_iter])
            for enum_segment in np.arange(self.I["num_seg_iters"][enum_iter]):
                start=(enum_segment)*self.I["num_nodes_each_seg_iters"][enum_iter]
                end=(enum_segment+1)*self.I["num_nodes_each_seg_iters"][enum_iter]      
                start_All_nodes=(enum_segment)*self.I["num_All_nodes_each_seg_iters"][enum_iter]
                end_All_nodes=(enum_segment+1)*self.I["num_All_nodes_each_seg_iters"][enum_iter]   
                self.I["tau_iters_unscaled_vec"][enum_iter][start:end]=self.I["tau_iters_unscaled"][enum_iter][enum_segment]
                self.I["tau_iters_All_nodes_unscaled_vec"][enum_iter][start_All_nodes:end_All_nodes]=self.I["tau_iters_All_nodes_unscaled"][enum_iter][enum_segment]
                self.I["tau_iters_vec"][enum_iter][start:end]=self.I["tau_iters"][enum_iter][enum_segment]
                self.I["tau_iters_All_nodes_vec"][enum_iter][start_All_nodes:end_All_nodes]=self.I["tau_iters_All_nodes"][enum_iter][enum_segment]
            
        self.I["fixed_u_indices_iter"]=[None] * self.I["num_ietartion"]
        self.I["nonfixed_u_indices_iter"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["fixed_u_indices_iter"][enum_iter]=np.where(self.I["tau_iters_unscaled_vec"][enum_iter]>self.I["tm_iters"][enum_iter,1])[0]
            self.I["nonfixed_u_indices_iter"][enum_iter]=np.where(self.I["tau_iters_unscaled_vec"][enum_iter]<=self.I["tm_iters"][enum_iter,1])[0]
        
        self.I["nonfixed_Ts_satisfied_u_indices_iter"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            last_index_not_included=1
            self.I["nonfixed_Ts_satisfied_u_indices_iter"][enum_iter]=[] #define as a list
            time_unscaled_nonfixed_u=self.I["tau_iters_unscaled_vec"][enum_iter][self.I["nonfixed_u_indices_iter"][enum_iter]]
            for enum_non_fixed_u_index in np.arange(len(self.I["nonfixed_u_indices_iter"][enum_iter])):
                if enum_non_fixed_u_index==0:
                    self.I["nonfixed_Ts_satisfied_u_indices_iter"][enum_iter].append(self.I["nonfixed_u_indices_iter"][enum_iter][enum_non_fixed_u_index])
                elif (time_unscaled_nonfixed_u[enum_non_fixed_u_index]-time_unscaled_nonfixed_u[enum_non_fixed_u_index-last_index_not_included])>=self.MPC_params["Ts"]:
                    last_index_not_included=1
                    self.I["nonfixed_Ts_satisfied_u_indices_iter"][enum_iter].append(self.I["nonfixed_u_indices_iter"][enum_iter][enum_non_fixed_u_index])
                else:
                    last_index_not_included+=1 
        
        #--------------------------Modification------just for last iter--------------------
        #self.I["nonfixed_Ts_satisfied_u_indices_iter"][-1]=self.I["nonfixed_u_indices_iter"][-1].copy()
        self.I["nonfixed_Ts_satisfied_u_indices_iter"][-1]=self.I["nonfixed_u_indices_iter"][-1][np.array([0,-1])]
                
        self.I["tau_iters_unscaled_vec_non_fixed_Ts_satisfied"]=[None] * self.I["num_ietartion"]
        self.I["tau_iters_vec_non_fixed_Ts_satisfied"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["tau_iters_unscaled_vec_non_fixed_Ts_satisfied"][enum_iter]=self.I["tau_iters_unscaled_vec"][enum_iter][self.I["nonfixed_Ts_satisfied_u_indices_iter"][enum_iter]]
            self.I["tau_iters_vec_non_fixed_Ts_satisfied"][enum_iter]=self.I["tau_iters_vec"][enum_iter][self.I["nonfixed_Ts_satisfied_u_indices_iter"][enum_iter]]
        
        self.I["num_nodes_each_u_nonfixed_Ts_satisfied"]=np.zeros(self.I["num_ietartion"],dtype=int)
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][enum_iter]=len(self.I["nonfixed_Ts_satisfied_u_indices_iter"][enum_iter])
        
        self.I["num_totall_variables_iters"]=np.zeros(self.I["num_ietartion"],dtype=int)
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["num_totall_variables_iters"][enum_iter]=(self.I["num_All_nodes_iters"][enum_iter])*self.I["num_s"]+\
                                                             (self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][enum_iter])*\
                                                             self.I["num_u"]
        
        self.I["opt_var_state_index"]=[None] * self.I["num_ietartion"]
        self.I["opt_var_state_colloc_index"]=[None] * self.I["num_ietartion"]
        self.I["opt_var_control_index"]=[None] * self.I["num_ietartion"]
        for enum_iter in np.arange(self.I["num_ietartion"]):
            self.I["opt_var_state_index"][enum_iter]=dict()
            self.I["opt_var_state_colloc_index"][enum_iter]=dict()
            self.I["opt_var_control_index"][enum_iter]=dict()
            for enum_state in np.arange(self.I["num_s"]):
                nsv=self.I["num_All_nodes_iters"][enum_iter] #num_each_state_var
                nsv_colloc=self.I["num_coll_nodes_iters"][enum_iter]
                start_indx_state=enum_state*nsv
                end_indx_state=(enum_state+1)*nsv
                end_indx_state_colloc=start_indx_state+nsv_colloc
                self.I["opt_var_state_index"][enum_iter][enum_state]=np.arange(start_indx_state,end_indx_state,1,dtype=int)
                colloc_indices=np.reshape(np.array(self.I["coll_indices_at_each_segment_iters_vec"][enum_iter]),-1)
                self.I["opt_var_state_colloc_index"][enum_iter][enum_state]=self.I["opt_var_state_index"][enum_iter][enum_state][colloc_indices]
            for enum_control in np.arange(self.I["num_u"]):
                nuv=self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][enum_iter] #num_each_control_var
                start_indx_control=self.I["num_s"]*nsv+enum_control*nuv
                end_indx_control=self.I["num_s"]*nsv+(enum_control+1)*nuv
                self.I["opt_var_control_index"][enum_iter][enum_control]=np.arange(start_indx_control,end_indx_control,1,dtype=int)
        k=1
        
    def Bound(self,iter_num):
        opt_var_lb=-1.0*np.ones(self.I["num_totall_variables_iters"][iter_num])
        opt_var_ub=1.0*np.ones(self.I["num_totall_variables_iters"][iter_num])
        Bound_matrix=np.array([opt_var_lb,opt_var_ub]).T
        Bound_Matrix_Tuple=tuple(map(tuple, Bound_matrix))
        return Bound_Matrix_Tuple
    
    def Guess(self,iter_num):
        
        #You can also get access or set state and control states using:
           # self.I["opt_var_state_index"][iter_num][enum_state]
           # self.I["opt_var_control_index"][iter_num][enum_control]
        
        x0=np.zeros(self.I["num_totall_variables_iters"][iter_num])
        if iter_num==0:
            for enum_state in np.arange(self.I["num_s"]):
                nsv=self.I["num_All_nodes_iters"][iter_num] #num_each_state_var
                start_indx_state=int(enum_state*nsv)
                end_indx_state=int((enum_state+1)*nsv)
                x0[start_indx_state:end_indx_state]=np.linspace(self.Model_data["state_lb_guess"][enum_state],\
                                                                self.Model_data["state_ub_guess"][enum_state],\
                                                                    int(nsv))
            for enum_control in np.arange(self.I["num_u"]):
                nuv=self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][iter_num] #num_each_control_var
                start_indx_control=int(self.I["num_s"]*nsv+enum_control*nuv)
                end_indx_control=int(self.I["num_s"]*nsv+(enum_control+1)*nuv)
                x0[start_indx_control:end_indx_control]=np.linspace(self.Model_data["control_lb_guess"][enum_control],\
                                                                self.Model_data["control_ub_guess"][enum_control],\
                                                                    int(nuv))
        elif self.last_iter_success==True:
            x0=self.x_guess_next
            #To DO
        return x0
    
    def Objective(self,iter_num):
        fun = lambda x: self.Obj_seg_summation(x,iter_num)
        return fun
    
    def Constraint(self,iter_num):
         
        cons_list=[]    
        cons_defect = {'type': 'eq', 'fun': lambda x:  self.Defect_Const_over_segments(x,iter_num)}
        cons_list.append(cons_defect)    
            
        boundary_const={'type': 'eq', 'fun': lambda x:  self.Boundary_over_iters(x,iter_num)}
        cons_list.append(boundary_const) 
        
        Path_const={'type': 'ineq', 'fun': lambda x:  self.Path_over_iter(x,iter_num)}
        cons_list.append(Path_const)
        
        continuity_const={'type': 'eq', 'fun': lambda x:  self.Continuity_over_segments(x,iter_num)}
        cons_list.append(continuity_const)
        
        continuity_initial={'type': 'eq', 'fun': lambda x:  self.Continuity_initial(x,iter_num)}
        cons_list.append(continuity_initial)
        
        cons=tuple(cons_list)    
        return cons
    
    def Path_over_iter(self,x,iter_num):
        
        colloc_indices_vec=np.reshape(np.array(self.I["coll_indices_at_each_segment_iters_vec"][iter_num]),-1)
        
        u_extended_matrix_iter=self.u_extended_matrix(x,iter_num)
        u_extended_matrix_iter_colloc=u_extended_matrix_iter[colloc_indices_vec,:]
        
        tau_iters_unscaled_vec=self.I["tau_iters_unscaled_vec"][iter_num]
        tau_iters_All_nodes_unscaled_vec=self.I["tau_iters_All_nodes_unscaled_vec"][iter_num]
        
        num_colloc_nodes=self.I["num_coll_nodes_iters"][iter_num] #[iter,state]
                
        states_colloc_matrix=np.zeros((num_colloc_nodes,self.I["num_s"])) #all_nodes *ns
        for enum_state in np.arange(self.I["num_s"]):
            states_colloc_matrix[:,enum_state]=x[self.I["opt_var_state_colloc_index"][iter_num][enum_state]]
            #states_colloc_matrix[:,enum_state]=x[self.I["opt_var_state_index"][iter_num][enum_state]][colloc_indices_vec]
            
        
        states_colloc_matrix_unscaled=np.zeros((num_colloc_nodes,self.I["num_s"])) #all_nodes *ns
        u_extended_colloc_nodes_matrix_iter_unscaled=np.zeros((num_colloc_nodes,self.I["num_u"])) #all_nodes *ns
        for enum_state in np.arange(self.I["num_s"]):
            states_colloc_matrix_unscaled[:,enum_state]=self.Scaling(Signal=states_colloc_matrix[:,enum_state],\
                                                                      Number=enum_state,\
                                                                      Type=0,
                                                                      S_US=1)
        for enum_control in np.arange(self.I["num_u"]):
            u_extended_colloc_nodes_matrix_iter_unscaled[:,enum_control]=self.Scaling(Signal=u_extended_matrix_iter_colloc[:,enum_control],\
                                                                            Number=enum_control,\
                                                                            Type=1,\
                                                                            S_US=1)
        
        #---------user func---------
        path_vec=self.Functions["Path_func"](tau_iters_unscaled_vec,\
                                             states_colloc_matrix_unscaled,\
                                             u_extended_colloc_nodes_matrix_iter_unscaled,\
                                             self.Model_data)
        
        return path_vec
    
    def Boundary_over_iters(self,x,iter_num):
        
        init_flag=self.Model_data["Init_state_boudnary_flag"]
        End_flag=self.Model_data["End_state_boudnary_flag"]
        
        num_all_nodes=len(self.I["opt_var_state_index"][iter_num][0]) #[iter,state]
        states_all_matrix=np.zeros((num_all_nodes,self.I["num_s"])) #all_nodes *ns
        for enum_state in np.arange(self.I["num_s"]):
            states_all_matrix[:,enum_state]=x[self.I["opt_var_state_index"][iter_num][enum_state]]
        
        Boundary_const_list=[]
        if (init_flag==1):
            if self.I["tp_iters"][iter_num,0]==self.Model_data["t0"]: #initnal
                init_con=states_all_matrix[0,:]-self.Model_data["Initial_State"]
                Boundary_const_list.append(init_con)
        
        if (End_flag==1):
            if self.I["tp_iters"][iter_num,1]==self.Model_data["tf"]: #terminal
                end_con=states_all_matrix[-1,:]-self.Model_data["End_State"]
                Boundary_const_list.append(end_con)
            
        num_const=len(Boundary_const_list) # 1 or 2
        Boundary_const=np.zeros((num_const,self.I["num_s"]))
        for enum_boundary_const in np.arange(num_const):
            Boundary_const[enum_boundary_const,:]=Boundary_const_list[enum_boundary_const]
                
        Boundary_const_vec=np.reshape(Boundary_const,-1)       
        return Boundary_const_vec
    
    def Continuity_initial(self,x,iter_num):
        
        if iter_num>=1:
        
            inital_state=self.I["Cont_initla_val_const_state"] #obtiand frol last iteration
            continuity_initial_state_vec=np.zeros(self.I["num_s"])
            for enum_state in np.arange(self.I["num_s"]):
                continuity_initial_state_vec[enum_state]=inital_state[enum_state]-x[self.I["opt_var_state_index"][iter_num][enum_state][0]] #0:first element
        else:
            continuity_initial_state_vec=[]
        
        return continuity_initial_state_vec
        
    
    def Continuity_over_segments(self,x,iter_num):
        if self.I["num_seg_iters"][iter_num]==1:
            continuity_const=[]
        else:
            continuity_const=np.zeros((self.I["num_seg_iters"][iter_num]-1)*self.I["num_s"])
            colloc_next_init=np.zeros((self.I["num_seg_iters"][iter_num]-1,self.I["num_s"]))
            non_colloc_now_end=np.zeros((self.I["num_seg_iters"][iter_num]-1,self.I["num_s"]))
            for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]-1):
                
                segment_index_colloc_next_init=self.I["coll_indices_at_each_segment_iters"][iter_num][enum_segment+1][0]
                
                segment_indices_regressin_now_end=self.I["Endpoint_indices_at_each_segment_iters"][iter_num][enum_segment]
                
                opt_var_state_segment_colloc_indices_next_init=[]
                for enum_state in np.arange(self.I["num_s"]):
                    index=self.I["opt_var_state_index"][iter_num][enum_state][segment_index_colloc_next_init]
                    opt_var_state_segment_colloc_indices_next_init.append(index) #list: ns
                
                opt_var_state_segment_regression_indices_now_end=[]
                for enum_state in np.arange(self.I["num_s"]):
                    index=self.I["opt_var_state_index"][iter_num][enum_state][segment_indices_regressin_now_end]
                    opt_var_state_segment_regression_indices_now_end.append(index) #list: ns
                    
                for enum_state in np.arange(self.I["num_s"]):
                    colloc_next_init[enum_segment,enum_state]=x[opt_var_state_segment_colloc_indices_next_init[enum_state]]
                    
                for enum_state in np.arange(self.I["num_s"]):
                    non_colloc_now_end[enum_segment,enum_state]=x[opt_var_state_segment_regression_indices_now_end[enum_state]]
                    
                continuity_con=colloc_next_init[enum_segment,:]-non_colloc_now_end[enum_segment,:]
                
                continuity_const[enum_segment*self.I["num_s"]:(enum_segment+1)*self.I["num_s"]]=continuity_con
                
        return continuity_const
                
    def Defect_Const_over_segments(self,x,iter_num):
        defect_const=np.zeros(self.I["num_coll_nodes_iters"][iter_num]*self.I["num_s"])
        for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
            t_start_segment=self.I["time_unscaled_iters_array"][iter_num][enum_segment]
            t_end_segment=self.I["time_unscaled_iters_array"][iter_num][enum_segment+1]
            h_segment=t_end_segment-t_start_segment
            D_segment=self.I["D_iters"][iter_num][enum_segment]
            
            tau_unscaled_segment=self.I["tau_iters_unscaled"][iter_num][enum_segment]
            segment_indices=self.I["coll_indices_at_each_segment_iters"][iter_num][enum_segment]
            segment_indices_All_nodes=self.I["All_indices_at_each_segment_iters"][iter_num][enum_segment]
            
            opt_var_state_segment_indices_colloc=[]
            opt_var_state_segment_indices_all=[]
            for enum_state in np.arange(self.I["num_s"]):
                index_colloc=self.I["opt_var_state_index"][iter_num][enum_state][segment_indices]
                index_all=self.I["opt_var_state_index"][iter_num][enum_state][segment_indices_All_nodes]
                opt_var_state_segment_indices_colloc.append(index_colloc) #list: n
                opt_var_state_segment_indices_all.append(index_all) #list: n
                
            u_extended_matrix_iter=self.u_extended_matrix(x,iter_num) #scaled
            u_extended_matrix_segment=u_extended_matrix_iter[segment_indices,:] #seg*nu
            
            states_segment_matrix_colloc=np.zeros((self.I["num_nodes_each_seg_iters"][iter_num],self.I["num_s"]))
            for enum_state in np.arange(self.I["num_s"]):
                states_segment_matrix_colloc[:,enum_state]=x[opt_var_state_segment_indices_colloc[enum_state]]
                
            states_segment_matrix_all=np.zeros((self.I["num_All_nodes_each_seg_iters"][iter_num],self.I["num_s"]))
            for enum_state in np.arange(self.I["num_s"]):
                states_segment_matrix_all[:,enum_state]=x[opt_var_state_segment_indices_all[enum_state]]
                
            states_segment_matrix_colloc_unscaled=np.zeros((self.I["num_nodes_each_seg_iters"][iter_num],self.I["num_s"]))
            for enum_state in np.arange(self.I["num_s"]):
                states_segment_matrix_colloc_unscaled[:,enum_state]=self.Scaling(Signal=states_segment_matrix_colloc[:,enum_state],\
                                                                          Number=enum_state,\
                                                                          Type=0,
                                                                          S_US=1)

            u_extended_matrix_segment_unscaled=np.zeros(u_extended_matrix_segment.shape)
            for enum_control in np.arange(self.I["num_u"]):
                u_extended_matrix_segment_unscaled[:,enum_control]=self.Scaling(Signal=u_extended_matrix_segment[:,enum_control],\
                                                                                Number=enum_control,\
                                                                                Type=1,\
                                                                                S_US=1)
            
            F=self.F_dynamics(tau_unscaled_segment,\
                              states_segment_matrix_colloc_unscaled,\
                              u_extended_matrix_segment_unscaled) #num_all_nodes*n_s
                
            segment_Defect_const_matrix=(D_segment@states_segment_matrix_all-h_segment/2*F)    
            segment_Defect_const_vec=np.reshape(segment_Defect_const_matrix,-1)
            start_index=self.I["num_s"]*enum_segment*self.I["num_nodes_each_seg_iters"][iter_num]
            end_index=self.I["num_s"]*(enum_segment+1)*self.I["num_nodes_each_seg_iters"][iter_num]
            defect_const[start_index:end_index]=segment_Defect_const_vec
        
        return defect_const
    
    def F_dynamics(self,tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled):
        
        F_unscaled=self.Functions["F_dynamics_user"](tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,self.Model_data)
        
        F=np.zeros(F_unscaled.shape)
        for enum_state in np.arange(self.I["num_s"]):
            F[:,enum_state]=self.Model_data["B_scaled"][enum_state]**(-1)*F_unscaled[:,enum_state]
        
        return F
    
    def Lagrange(self,x,iter_num,seg_num):
        
        tau_unscaled_segment=self.I["tau_iters_unscaled"][iter_num][seg_num]
        segment_indices=self.I["coll_indices_at_each_segment_iters"][iter_num][seg_num]
        
        opt_var_state_segment_indices=[]
        for enum_state in np.arange(self.I["num_s"]):
            index=self.I["opt_var_state_index"][iter_num][enum_state][segment_indices]
            opt_var_state_segment_indices.append(index) #list: ns
        
        u_extended_matrix_iter=self.u_extended_matrix(x,iter_num)
        u_extended_matrix_segment=u_extended_matrix_iter[segment_indices,:] #seg*nu
        #------------------------user-------------------
        
        states_segment_matrix=np.zeros((self.I["num_nodes_each_seg_iters"][iter_num],self.I["num_s"]))
        for enum_state in np.arange(self.I["num_s"]):
            states_segment_matrix[:,enum_state]=x[opt_var_state_segment_indices[enum_state]]
        
        states_segment_matrix_unscaled=np.zeros((self.I["num_nodes_each_seg_iters"][iter_num],self.I["num_s"]))
        for enum_state in np.arange(self.I["num_s"]):
            states_segment_matrix_unscaled[:,enum_state]=self.Scaling(Signal=states_segment_matrix[:,enum_state],\
                                                                      Number=enum_state,\
                                                                      Type=0,
                                                                      S_US=1)

        u_extended_matrix_segment_unscaled=np.zeros(u_extended_matrix_segment.shape)
        for enum_control in np.arange(self.I["num_u"]):
            u_extended_matrix_segment_unscaled[:,enum_control]=self.Scaling(Signal=u_extended_matrix_segment[:,enum_control],\
                                                                            Number=enum_control,\
                                                                            Type=1,\
                                                                            S_US=1)
            
        L=self.Functions["Lagrange_User_Func"](tau_unscaled_segment,\
                                               states_segment_matrix_unscaled,\
                                               u_extended_matrix_segment_unscaled,\
                                               self.Model_data)
        
        return L
        
    def Scaling(self,Signal,Number,Type,S_US):
        #Number: Number of state or control
        #Type: 0 for state, 1 for control
        #S_US: 0 for scaling, 1 for unscaling
        if Type==0: #state
            if S_US==0: #scaling
                Transformed_Signal=(Signal-self.Model_data["A_scaled"][Number])/self.Model_data["B_scaled"][Number]
            elif S_US==1: #unscaling
                Transformed_Signal=self.Model_data["A_scaled"][Number]+self.Model_data["B_scaled"][Number]*Signal
        elif Type==1: #control
            if S_US==0: #scaling
                Transformed_Signal=(Signal-self.Model_data["C_scaled"][Number])/self.Model_data["D_scaled"][Number]
            elif S_US==1: #unscaling
                Transformed_Signal=self.Model_data["C_scaled"][Number]+self.Model_data["D_scaled"][Number]*Signal
        return Transformed_Signal    
    
    def Mayer_Cost(self,x,iter_num):
        Meyr_Cost_val=0.
        if self.I["tp_iters"][iter_num,0]==self.Model_data["t0"]:
            state_t0_scaled=np.zeros(self.I["num_s"])
            state_t0_unscaled=np.zeros(self.I["num_s"])
            for enum_state in np.arange(self.I["num_s"]):
                state_t0_scaled[enum_state]=x[self.I["opt_var_state_index"][iter_num][enum_state]][0]
            for enum_state in np.arange(self.I["num_s"]):
                state_t0_unscaled[enum_state]=self.Scaling(Signal=state_t0_scaled[enum_state],\
                                                                             Number=enum_state,\
                                                                             Type=0,
                                                                             S_US=1)
            Meyr_Cost_val+=self.Functions["User_Mayer_Cost"](self.Model_data["t0"],state_t0_unscaled,self.Model_data)
        if self.I["tp_iters"][iter_num,1]==self.Model_data["tf"]:
            state_tf_scaled=np.zeros(self.I["num_s"])
            state_tf_unscaled=np.zeros(self.I["num_s"])
            for enum_state in np.arange(self.I["num_s"]):
                state_tf_scaled[enum_state]=x[self.I["opt_var_state_index"][iter_num][enum_state]][-1]
            for enum_state in np.arange(self.I["num_s"]):
                state_tf_unscaled[enum_state]=self.Scaling(Signal=state_tf_scaled[enum_state],\
                                                                              Number=enum_state,\
                                                                              Type=0,
                                                                              S_US=1)
            Meyr_Cost_val+=self.Functions["User_Mayer_Cost"](self.Model_data["tf"],state_tf_unscaled,self.Model_data)
        return Meyr_Cost_val
            
    
    def Obj_seg_summation(self,x,iter_num):
        obj_iter=0.0
        for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
            t_start_segment=self.I["time_unscaled_iters_array"][iter_num][enum_segment]
            t_end_segment=self.I["time_unscaled_iters_array"][iter_num][enum_segment+1]
            h_segment=t_end_segment-t_start_segment
            w_segment=self.I["W_iters"][iter_num][enum_segment]
            L=self.Lagrange(x,iter_num,enum_segment)
            obj_iter+=h_segment/2.0*sum(w_segment*L)
        mayer_Cost=self.Mayer_Cost(x,iter_num)
        obj_iter+=mayer_Cost
        return obj_iter
            
    def u_extended_matrix(self,x,iter_num): #scaled
        time_coll_iters=self.I["tau_iters_unscaled_vec"][iter_num]
        time_all_nodes=self.I["tau_iters_All_nodes_unscaled_vec"][iter_num]
        tau_unsclaed_cont_opt_var=self.I["tau_iters_unscaled_vec_non_fixed_Ts_satisfied"][iter_num]
        u_opt_var_matrix=np.zeros((self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][iter_num],self.I["num_u"]))
        
        for enum_control in np.arange(self.I["num_u"]):
            u_opt_var_matrix[:,enum_control]=x[self.I["opt_var_control_index"][iter_num][enum_control]]
        
        if self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][iter_num]!=1:
            u_interp_func=interp1d(tau_unsclaed_cont_opt_var,u_opt_var_matrix,axis=0,kind='previous',fill_value="extrapolate")
            u_extended_matrix=u_interp_func(time_all_nodes)
        else:
            u_extended_matrix=np.tile(u_opt_var_matrix[0,:],(self.I["num_All_nodes_iters"][iter_num],1))
        return u_extended_matrix
    
    def Extraxt_Data_at_each_iter(self,\
                          iter_num,\
                          tau_iters_All_nodes_unscaled_vec,\
                          states_all_matrix,\
                          u_extended_matrix_iter,
                          states_all_matrix_unscaled,\
                          u_extended_matrix_iter_unscaled\
                          ):
        
        
        time_all_nodes_segment=dict()
        states_all_matrix_segments=[None] * self.I["num_s"]
        states_all_matrix_segments_unscaled=[None] * self.I["num_s"]
        u_extended_matrix_iter_segments=[None] * self.I["num_u"]
        u_extended_matrix_iter_segments_unscaled=[None] * self.I["num_u"]
        
        for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
            time_all_nodes_segment[enum_segment]=tau_iters_All_nodes_unscaled_vec[self.I["All_indices_at_each_segment_iters"][iter_num][enum_segment]]
        
        for enum_state in np.arange(self.I["num_s"]):
            states_all_matrix_segments[enum_state]=dict()
            states_all_matrix_segments_unscaled[enum_state]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
                seg_indices=self.I["All_indices_at_each_segment_iters"][iter_num][enum_segment]
                states_all_matrix_segments[enum_state][enum_segment]=states_all_matrix[seg_indices,enum_state]
                states_all_matrix_segments_unscaled[enum_state][enum_segment]=states_all_matrix_unscaled[seg_indices,enum_state]
        
        
        for enum_control in np.arange(self.I["num_u"]):
            u_extended_matrix_iter_segments[enum_control]=dict()
            u_extended_matrix_iter_segments_unscaled[enum_control]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
                seg_indices=self.I["All_indices_at_each_segment_iters"][iter_num][enum_segment]
                u_extended_matrix_iter_segments[enum_control][enum_segment]=u_extended_matrix_iter[seg_indices,enum_control]
                u_extended_matrix_iter_segments_unscaled[enum_control][enum_segment]=u_extended_matrix_iter_unscaled[seg_indices,enum_control]
                
        
        states_all_matrix_segments_polinomial=[None] * self.I["num_s"]
        states_all_matrix_segments_unscaled_polinomial=[None] * self.I["num_s"]
        u_extended_matrix_iter_segments_interp1=[None] * self.I["num_u"]
        u_extended_matrix_iter_segments_unscaled_interp1=[None] * self.I["num_u"]
        
        
        for enum_state in np.arange(self.I["num_s"]):
            states_all_matrix_segments_polinomial[enum_state]=dict()
            states_all_matrix_segments_unscaled_polinomial[enum_state]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
                states_all_matrix_segments_polinomial[enum_state][enum_segment]=lagrange(time_all_nodes_segment[enum_segment],\
                                                                                         states_all_matrix_segments[enum_state][enum_segment])
                states_all_matrix_segments_unscaled_polinomial[enum_state][enum_segment]=lagrange(time_all_nodes_segment[enum_segment],\
                                                                                         states_all_matrix_segments_unscaled[enum_state][enum_segment])
                    
        for enum_control in np.arange(self.I["num_u"]):
            u_extended_matrix_iter_segments_interp1[enum_control]=dict()
            u_extended_matrix_iter_segments_unscaled_interp1[enum_control]=dict()
            for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
                u_extended_matrix_iter_segments_interp1[enum_control][enum_segment]=interp1d(time_all_nodes_segment[enum_segment],\
                                                                                         u_extended_matrix_iter_segments[enum_control][enum_segment],kind='previous',fill_value="extrapolate")
                u_extended_matrix_iter_segments_unscaled_interp1[enum_control][enum_segment]=interp1d(time_all_nodes_segment[enum_segment],\
                                                                                         u_extended_matrix_iter_segments_unscaled[enum_control][enum_segment],kind='previous',fill_value="extrapolate")
        
        
        t_fine_mesh=np.linspace(tau_iters_All_nodes_unscaled_vec[0],tau_iters_All_nodes_unscaled_vec[-1],20*len(tau_iters_All_nodes_unscaled_vec))
        
        states_all_matrix_fine_mesh=np.zeros((len(t_fine_mesh),self.I["num_s"]))
        states_all_matrix_unscaled_fine_mesh=np.zeros((len(t_fine_mesh),self.I["num_s"]))
        controls_all_matrix_fine_mesh=np.zeros((len(t_fine_mesh),self.I["num_u"]))
        controls_all_matrix_unscaled_fine_mesh=np.zeros((len(t_fine_mesh),self.I["num_u"]))
        
        for numel_t_fine in np.arange(len(t_fine_mesh)):
            t_fine=t_fine_mesh[numel_t_fine]
            seg_index=self.index_finder(iter_num,t_fine,time_all_nodes_segment)
            for enum_state in np.arange(self.I["num_s"]):
                state_poly=states_all_matrix_segments_polinomial[enum_state][seg_index]
                state_unsc_poly=states_all_matrix_segments_unscaled_polinomial[enum_state][seg_index]
                states_all_matrix_fine_mesh[numel_t_fine,enum_state]=np.polyval(state_poly,t_fine)
                states_all_matrix_unscaled_fine_mesh[numel_t_fine,enum_state]=np.polyval(state_unsc_poly,t_fine)
            for enum_control in np.arange(self.I["num_u"]):
                controls_all_matrix_fine_mesh[numel_t_fine,enum_control]=u_extended_matrix_iter_segments_interp1[enum_control][seg_index](t_fine)
                controls_all_matrix_unscaled_fine_mesh[numel_t_fine,enum_control]=u_extended_matrix_iter_segments_unscaled_interp1[enum_control][seg_index](t_fine)
         
            
        if self.Model_data["plot_iter_flag"]==True: 
            
            self.ax_state.cla()
            self.ax_control.cla()
            
            #fig_state, ax_state = plt.subplots()        
            for enum_state in np.arange(self.I["num_s"]):
                self.ax_state.plot(tau_iters_All_nodes_unscaled_vec,states_all_matrix[:,enum_state],'o',markersize=4,\
                                   label=fr'$\xi_{enum_state+1}$'+'$\mathrm{\,\,Sol}$')
                self.ax_state.plot(t_fine_mesh,states_all_matrix_fine_mesh[:,enum_state],\
                                   linewidth=2,label=fr'$\xi_{enum_state+1}$'+'$\mathrm{\,\,Sim}$')  
            self.ax_state.set_xlabel('$t\,\mathrm{[s]}$')
            self.ax_state.set_ylabel('$\mathrm{[Scaled\,\,States]}$')
            self.ax_state.legend(prop=fontP)
            self.ax_state.set_ylim(bottom=-1-0.1, top=1+0.1)
            self.ax_state.set_xlim(left=self.Model_data["t0"]-self.Model_data["tf"]/20,right=self.Model_data["tf"]*1.05)
            
            
            #fig_state_unsc, ax_state_unsc = plt.subplots()        
            """for enum_state in np.arange(self.I["num_s"]):
                self.ax_state .plot(tau_iters_All_nodes_unscaled_vec,states_all_matrix_unscaled[:,enum_state],'o-')
                self.ax_state .plot(t_fine_mesh,states_all_matrix_unscaled_fine_mesh[:,enum_state])"""
            
        
            #fig_control, ax_control = plt.subplots()        
            for enum_control in np.arange(self.I["num_u"]):
                self.ax_control.plot(tau_iters_All_nodes_unscaled_vec,u_extended_matrix_iter[:,enum_control],'o-',markersize=4,\
                                   label=fr'$u_{enum_control+1}$'+'$\mathrm{\,\,Sol}$')
                self.ax_control.plot(t_fine_mesh,controls_all_matrix_fine_mesh[:,enum_control],\
                                   linewidth=2,label=fr'$u_{enum_control+1}$'+'$\mathrm{\,\,Sim}$')    
            self.ax_control.set_xlabel('$t\,\mathrm{[s]}$')
            self.ax_control.set_ylabel('$\mathrm{[Scaled\,\,Controls]}$')
            self.ax_control.legend(prop=fontP)
            self.ax_control.set_ylim(bottom=-1-0.1, top=1+0.1)
            self.ax_control.set_xlim(left=self.Model_data["t0"]-self.Model_data["tf"]/20,right=self.Model_data["tf"]*1.05)
            
            #fig_control_unsc, ax_control_unsc = plt.subplots()        
            """for enum_control in np.arange(self.I["num_u"]):
                self.ax_control.plot(tau_iters_All_nodes_unscaled_vec,u_extended_matrix_iter_unscaled[:,enum_control],'o-')
                self.ax_control.plot(t_fine_mesh,controls_all_matrix_unscaled_fine_mesh[:,enum_control]) """  
            
            plt.show()
            plt.pause(0.1)
        
            
        return states_all_matrix_segments_polinomial,\
               states_all_matrix_segments_unscaled_polinomial,\
               u_extended_matrix_iter_segments_interp1,\
               u_extended_matrix_iter_segments_unscaled_interp1,\
               time_all_nodes_segment
               #[state][seg](time_unscaled)
               #[control][seg](time_unscaled)
                
            
    def index_finder(self,iter_num,t_node,time_all_nodes_segment):
        for enum_segment in np.arange(self.I["num_seg_iters"][iter_num]):
            if t_node<=time_all_nodes_segment[enum_segment][-1]:
                seg_index=enum_segment
                break
        return seg_index
        
            
            
    def Run_MPC_iter(self,iter_num=0):
        x0=self.Guess(iter_num)
        optionsdict=self.optionsdict
        Bound_Matrix_Tuple=self.Bound(iter_num)
        cons=self.Constraint(iter_num)
        Objective_function=self.Objective(iter_num)
        
        """Opt_result = minimize_ipopt(fun=Objective_function,\
                                    x0=x0,\
                                    options=optionsdict,\
                                    bounds=Bound_Matrix_Tuple,\
                                    constraints=cons)"""
    
        Opt_result = minimize(fun=Objective_function,\
                                    x0=x0,\
                                    options=optionsdict,\
                                    bounds=Bound_Matrix_Tuple,\
                                    constraints=cons)
            
        #print(Opt_result["x"])
        #print(Opt_result["success"])
        #print(Opt_result["message"])
        
        
        x=Opt_result["x"]
        tau_iters_All_nodes_unscaled_vec=self.I["tau_iters_All_nodes_unscaled_vec"][iter_num]
        u_extended_matrix_iter=self.u_extended_matrix(x,iter_num)
        num_All_nodes_each_state=len(self.I["opt_var_state_index"][iter_num][0]) #[iter,state]=
        states_all_matrix=np.zeros((num_All_nodes_each_state,self.I["num_s"])) #all_nodes *ns
        for enum_state in np.arange(self.I["num_s"]):
            states_all_matrix[:,enum_state]=x[self.I["opt_var_state_index"][iter_num][enum_state]]
        u_extended_matrix_iter_unscaled=np.zeros((num_All_nodes_each_state,self.I["num_u"])) #all_nodes *ns
        states_all_matrix_unscaled=np.zeros((num_All_nodes_each_state,self.I["num_s"])) #all_nodes *ns
        for enum_state in np.arange(self.I["num_s"]):
            states_all_matrix_unscaled[:,enum_state]=self.Scaling(Signal=states_all_matrix[:,enum_state],\
                                                                      Number=enum_state,\
                                                                      Type=0,
                                                                      S_US=1)
        for enum_control in np.arange(self.I["num_u"]):
            u_extended_matrix_iter_unscaled[:,enum_control]=self.Scaling(Signal=u_extended_matrix_iter[:,enum_control],\
                                                                            Number=enum_control,\
                                                                            Type=1,\
                                                                            S_US=1)
                
        self.last_iter_success=Opt_result["success"]  
        self.I["opt_success_flag"][iter_num]=Opt_result["success"]
        #print(self.last_iter_success)
        #Provide nintal guess for the nest iteration       
        time_coll_iters=self.I["tau_iters_unscaled_vec"][iter_num]
        time_all_nodes=self.I["tau_iters_All_nodes_unscaled_vec"][iter_num]
        tau_unsclaed_cont_opt_var=self.I["tau_iters_unscaled_vec_non_fixed_Ts_satisfied"][iter_num]
        u_opt_var_matrix=np.zeros((self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][iter_num],self.I["num_u"]))
        state_opt_var_matrix=np.zeros((self.I["num_All_nodes_iters"][iter_num],self.I["num_s"]))
        
        for enum_control in np.arange(self.I["num_u"]):
            u_opt_var_matrix[:,enum_control]=x[self.I["opt_var_control_index"][iter_num][enum_control]]
            
        for enum_state in np.arange(self.I["num_s"]):
            state_opt_var_matrix[:,enum_state]=x[self.I["opt_var_state_index"][iter_num][enum_state]]
        
            
        if iter_num!=self.I["num_ietartion"]-1:
            num_nodes_current=self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][iter_num]
            num_nodes_next=self.I["num_nodes_each_u_nonfixed_Ts_satisfied"][iter_num+1]
            u_opt_var_matrix_next_iter=np.zeros((num_nodes_next,self.I["num_u"]))
            if num_nodes_next==num_nodes_current:
                u_opt_var_matrix_next_iter=u_opt_var_matrix.copy()
            elif num_nodes_next<num_nodes_current:
                u_opt_var_matrix_next_iter=u_opt_var_matrix[0:num_nodes_next,:]
            elif num_nodes_next>num_nodes_current:
                n_diff=num_nodes_next-num_nodes_current
                u_opt_var_matrix_next_iter=np.append(u_opt_var_matrix,np.tile(np.reshape(u_opt_var_matrix[-1,:],(1,-1)),(n_diff,1)),axis=0)
                
        if iter_num!=self.I["num_ietartion"]-1:
            num_nodes_current=self.I["num_All_nodes_iters"][iter_num]
            num_nodes_next=self.I["num_All_nodes_iters"][iter_num+1]
            state_opt_var_matrix_next_iter=np.zeros((num_nodes_next,self.I["num_s"]))
            if num_nodes_next==num_nodes_current:
                state_opt_var_matrix_next_iter=state_opt_var_matrix.copy()
            elif num_nodes_next<num_nodes_current:
                state_opt_var_matrix_next_iter=state_opt_var_matrix[0:num_nodes_next,:]
            elif num_nodes_next>num_nodes_current:
                n_diff=num_nodes_next-num_nodes_current
                state_opt_var_matrix_next_iter=np.append(state_opt_var_matrix,np.tile(np.reshape(state_opt_var_matrix[-1,:],(1,-1)),(n_diff,1)),axis=0)
        
        if iter_num!=self.I["num_ietartion"]-1:
            x_guess_next=np.zeros(self.I["num_totall_variables_iters"][iter_num+1])
            for enum_state in np.arange(self.I["num_s"]):
                x_guess_next[self.I["opt_var_state_index"][iter_num+1][enum_state]]=state_opt_var_matrix_next_iter[:,enum_state]
            
            
            for enum_control in np.arange(self.I["num_u"]):
                x_guess_next[self.I["opt_var_control_index"][iter_num+1][enum_control]]=u_opt_var_matrix_next_iter[:,enum_control]
                
            
            self.x_guess_next=x_guess_next
        
        return tau_iters_All_nodes_unscaled_vec,\
               states_all_matrix,\
               u_extended_matrix_iter,\
               states_all_matrix_unscaled,\
               u_extended_matrix_iter_unscaled 
               
    def plot_Result(self):
        
        fig_state, ax_state = plt.subplots(nrows=self.I["num_s"], ncols=1,sharex=True)
        if self.I["num_s"]==1:
            ax_state=[ax_state]
        for enum_state in np.arange(self.I["num_s"]):
            ax_state[enum_state].plot(self.I["time_h"],self.I["s_sc_matrix"][:,enum_state],'-',linewidth=2)
            ax_state[enum_state].set_ylabel(fr'$\xi_{enum_state+1}$')
        ax_state[enum_state].set_xlabel('$t\,\mathrm{[s]}$')
        fig_state.suptitle('$Scaled\,\,States$')
        
        fig_uns_state, ax_uns_state = plt.subplots(nrows=self.I["num_s"], ncols=1,sharex=True)
        if self.I["num_s"]==1:
            ax_uns_state=[ax_uns_state]
        for enum_state in np.arange(self.I["num_s"]):
            ax_uns_state[enum_state].plot(self.I["time_h"],self.I["s_unsc_matrix"][:,enum_state],'-',linewidth=2)
            ax_uns_state[enum_state].set_ylabel(fr'$\xi_{enum_state+1}$')
        ax_uns_state[enum_state].set_xlabel('$t\,\mathrm{[s]}$')
        fig_uns_state.suptitle('$Uncaled\,\,States$')
        
        fig_control,ax_control = plt.subplots(nrows=self.I["num_u"], ncols=1,sharex=True)
        if self.I["num_u"]==1:
            ax_control=[ax_control]
        for enum_control in np.arange(self.I["num_u"]):
            ax_control[enum_control].plot(self.I["time_h"],self.I["u_sc_matrix"][:,enum_control],'-',linewidth=2)
            ax_control[enum_control].set_ylabel(fr'$u_{enum_control+1}$')
        ax_control[enum_control].set_xlabel('$t\,\mathrm{[s]}$')
        fig_control.suptitle('$Scaled\,\,Control$')
        
        fig_uns_control,ax_uns_control = plt.subplots(nrows=self.I["num_u"], ncols=1,sharex=True)
        if self.I["num_u"]==1:
            ax_uns_control=[ax_uns_control]
        for enum_control in np.arange(self.I["num_u"]):
            ax_uns_control[enum_control].plot(self.I["time_h"],self.I["u_unsc_matrix"][:,enum_control],'-',linewidth=2)
            ax_uns_control[enum_control].set_ylabel(fr'$u_{enum_control+1}$')
        ax_uns_control[enum_control].set_xlabel('$t\,\mathrm{[s]}$')
        fig_uns_control.suptitle('$Unscaled\,\,Control$')
        
        fig_iters,ax_iters = plt.subplots(nrows=1, ncols=1)
        ax_iters.plot(np.arange(1,self.I["num_ietartion"]+1),self.I["opt_success_flag"],'o',markersize=5)
        ax_iters.set_xlabel('$\# \,\mathrm{Iter}$')
        ax_iters.set_ylabel('$\mathrm{Opt\,\,Success\,\, Flag}$')
        
        
        
               
    def ODE_Simulation(self,iter_num,tau_iters_All_nodes_unscaled_vec,states_all_matrix,u_extended_matrix_iter,states_all_matrix_unscaled,u_extended_matrix_iter_unscaled):
        start_index=iter_num*self.Mesh["n_h_Sol"]
        end_index=(iter_num+1)*self.Mesh["n_h_Sol"]+1
        index_vec=np.arange(start_index,end_index)
        
        time_h_iter=self.I["time_h"][start_index:end_index]
        t_span=(time_h_iter[0],time_h_iter[-1])
        
        U_interp_func=[None]*self.I["num_u"]
        for enum_control in np.arange(self.I["num_u"]):
            U_interp_func[enum_control]=interp1d(tau_iters_All_nodes_unscaled_vec,u_extended_matrix_iter_unscaled[:,enum_control],kind='previous',fill_value="extrapolate")
            
        inital_unscaled_statter_iter=states_all_matrix_unscaled[0,:]   
        
        dyn_fun=lambda t, unc_state: self.dynamic_ode(t,unc_state,U_interp_func)
        
        sol_ode=solve_ivp(dyn_fun, t_span, y0=inital_unscaled_statter_iter, method='RK45', t_eval=time_h_iter)
        
        time_sol_ode=sol_ode.t
        state_sol_ode=sol_ode.y.T
        
        
        for enum_state in np.arange(self.I["num_s"]):
            self.I["s_unsc_matrix"][start_index:end_index,enum_state]=state_sol_ode[:,enum_state]
        
        for enum_control in np.arange(self.I["num_u"]):
            self.I["u_unsc_matrix"][start_index:end_index,enum_control]=U_interp_func[enum_control](time_h_iter)
            
        
        for enum_state in np.arange(self.I["num_s"]):
            self.I["s_sc_matrix"][start_index:end_index,enum_state]=self.Scaling(Signal=self.I["s_unsc_matrix"][start_index:end_index,enum_state],\
                                                                      Number=enum_state,\
                                                                      Type=0,
                                                                      S_US=0)
        for enum_control in np.arange(self.I["num_u"]):
            self.I["u_sc_matrix"][start_index:end_index,enum_control]=self.Scaling(Signal=self.I["u_unsc_matrix"][start_index:end_index,enum_control],\
                                                                            Number=enum_control,\
                                                                            Type=1,\
                                                                            S_US=0)
        
        
        return self.I["s_sc_matrix"][end_index-1,:]
            
        
    def  dynamic_ode(self,t,unc_state,U_interp_func):
        
        nu=len(U_interp_func)
        
        u_extended_matrix_segment_unscaled=np.zeros((1,nu))
        for enum_control in np.arange(nu):
            u_extended_matrix_segment_unscaled[:,enum_control]=U_interp_func[enum_control](t)
        
        F_unscaled=self.Functions["F_dynamics_user"](np.array([t]),np.array([unc_state]),u_extended_matrix_segment_unscaled,self.Model_data)
        
        return F_unscaled.flatten()

    def Save_Result(self):
        Result_data=dict()
        Result_data["I"]=self.I
        #Result_data["Model_data"]=self.Model_data
        Result_data["MPC_params"]=self.MPC_params
        Result_data["Mesh"]=self.Mesh
        Result_data["optionsdict"]=self.optionsdict
        
        file_name_addrees=self.Model_data["File_Name"]
        with open(file_name_addrees,'wb') as file:
            pkl.dump(Result_data,file)
        
               
    def Loads_Data_Iter(self,\
                        iter_num,\
                        states_all_matrix_segments_polinomial,\
                        states_all_matrix_segments_unscaled_polinomial,\
                        u_extended_matrix_iter_segments_interp1,
                        u_extended_matrix_iter_segments_unscaled_interp1,\
                        time_all_nodes_segment):
        start_index=iter_num*self.Mesh["n_h_Sol"]
        end_index=(iter_num+1)*self.Mesh["n_h_Sol"]+1
        index_vec=np.arange(start_index,end_index)
        
        time_h_iter=self.I["time_h"][start_index:end_index]
        
        for time_index in np.arange(len(time_h_iter)):
            time_node=time_h_iter[time_index]
            seg_index=self.index_finder(iter_num,time_node,time_all_nodes_segment)
            for enum_state in np.arange(self.I["num_s"]):
                self.I["s_sc_matrix"][index_vec[time_index],enum_state]=np.polyval(states_all_matrix_segments_polinomial[enum_state][seg_index],time_node)
                self.I["s_unsc_matrix"][index_vec[time_index],enum_state]=np.polyval(states_all_matrix_segments_unscaled_polinomial[enum_state][seg_index],time_node)
            for enum_control in np.arange(self.I["num_u"]):
                self.I["u_sc_matrix"][index_vec[time_index],enum_control]=u_extended_matrix_iter_segments_interp1[enum_control][seg_index](time_node)
                self.I["u_unsc_matrix"][index_vec[time_index],enum_control]=u_extended_matrix_iter_segments_unscaled_interp1[enum_control][seg_index](time_node)
            
        return self.I["s_sc_matrix"][end_index-1,:]
            
    def Run_MPC_Loop(self):
        for iter_num in np.arange(self.I["num_ietartion"]):
            tau_iters_All_nodes_unscaled_vec,\
            states_all_matrix,\
            u_extended_matrix_iter,\
            states_all_matrix_unscaled,\
            u_extended_matrix_iter_unscaled =self.Run_MPC_iter(iter_num)
            
            Cont_initla_val_const_state=self.ODE_Simulation(iter_num,\
                                                            tau_iters_All_nodes_unscaled_vec,\
                                                            states_all_matrix,u_extended_matrix_iter,\
                                                            states_all_matrix_unscaled,\
                                                            u_extended_matrix_iter_unscaled)
            states_all_matrix_segments_polinomial,\
            states_all_matrix_segments_unscaled_polinomial,\
            u_extended_matrix_iter_segments_interp1,\
            u_extended_matrix_iter_segments_unscaled_interp1,\
            time_all_nodes_segment=self.Extraxt_Data_at_each_iter(
                                                                 iter_num,\
                                                                 tau_iters_All_nodes_unscaled_vec,\
                                                                 states_all_matrix,\
                                                                 u_extended_matrix_iter,
                                                                 states_all_matrix_unscaled,\
                                                                 u_extended_matrix_iter_unscaled,\
                                                                 ) 
            #[state][seg](time_unscaled)
            #[control][seg](time_unscaled)
            
            """Cont_initla_val_const_state=self.Loads_Data_Iter(
                                                    iter_num,\
                                                    states_all_matrix_segments_polinomial,\
                                                    states_all_matrix_segments_unscaled_polinomial,\
                                                    u_extended_matrix_iter_segments_interp1,\
                                                    u_extended_matrix_iter_segments_unscaled_interp1,\
                                                    time_all_nodes_segment
                                                    )"""
            self.I["Cont_initla_val_const_state"]=Cont_initla_val_const_state
            
        self.Save_Result()
        self.plot_Result()
        
def Read_Result(file_name):
    with open(file_name,'rb') as file:
        Result_data=pkl.load(file)
    return Result_data
 
def Lagrange_User_Func(tau_unscaled_segment,states_segment_matrix_unscaled,u_extended_matrix_segment_unscaled,Model_data):
    state_1_segment_unscaled=states_segment_matrix_unscaled[:,0]
    state_2_segment_unscaled=states_segment_matrix_unscaled[:,1]
    u_1_segment_unscaled=u_extended_matrix_segment_unscaled[:,0]
    
    L=1/2.*u_1_segment_unscaled**2.
    
    return L    

def Path_func(tau_iters_unscaled_vec,states_colloc_matrix_unscaled,u_extended_colloc_nodes_matrix_iter_unscaled,Model_data):
    path_vec_1=(states_colloc_matrix_unscaled[:,0]+20.)/20.
    path_vec_2=(states_colloc_matrix_unscaled[:,1]+20.)/20.
    path_vec=np.append(path_vec_1,path_vec_2)
    return path_vec      

def F_dynamics_user(tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data):
    F=np.zeros((len(tau_unscaled_segment),states_segment_matrix_colloc_unscaled.shape[1]))
    
    F[:,0]=states_segment_matrix_colloc_unscaled[:,1]
    F[:,1]=-states_segment_matrix_colloc_unscaled[:,0]+u_extended_matrix_segment_unscaled[:,0]
    return F   

def User_Mayer_Cost(time, state_unscaled,Model_data):
    mayer_cost=0.0
    if time==Model_data["t0"]:
        mayer_cost+=state_unscaled[0]*0.
    elif time==Model_data["tf"]:
        mayer_cost+=state_unscaled[0]*0.
    return  mayer_cost                                       
                
if __name__ == '__main__':  
    
    Model_data={
                "t0": 0,
                "tf":20,
                "x0":-1/2,
                "v0":1,
                "A_scaled":np.array([0.0,0.0]),
                "B_scaled":np.array([10.0,10.0]),
                "C_scaled":np.array([0.0]),
                "D_scaled":np.array([4.0]),
                "Initial_State":np.array([-0.05,0.1]), #scaled
                "End_State":np.array([0.0,0.0]),      #scaled
                "state_lb_guess":np.array([-0.1,-0.1]), #scaled
                "state_ub_guess":np.array([0.1,0.1]), #scaled
                "control_lb_guess":np.array([-0.1]), #scaled
                "control_ub_guess":np.array([0.1]),  #scaled
                "Init_state_boudnary_flag":1,
                "End_state_boudnary_flag":1,
                "plot_iter_flag":True,
                "File_Name":os.getcwd()+'/Results/'+'Ex1_Test1',
                }
    
    MPC_params={
                "Ts": 0.1,
                "m": 200,
                "p":200,
                }
    
    """Mesh={
        "num_seg_init":2,
        "num_nodes_init":4,
        "Dyn_Sol":1,
        "n_h_Sol":100}"""
    
    Mesh={
        "delta_t_seg":5,
        "num_nodes_seg":10,
        "n_h_Sol":100}
    
    Functions={
        "Lagrange_User_Func": lambda t_unsc_seg,state_unsc_seg,u_unsc_seg,Model_data :\
                                Lagrange_User_Func(t_unsc_seg,state_unsc_seg,u_unsc_seg,Model_data),\
        "Path_func": lambda t_unsc_iter,state_unsc_iter,u_unsc_iter,Model_data :\
                                Path_func(t_unsc_iter,state_unsc_iter,u_unsc_iter,Model_data),
        "F_dynamics_user": lambda tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data:\
                                F_dynamics_user(tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data),
        "User_Mayer_Cost": lambda time, state_unscaled,Model_data:\
                                  User_Mayer_Cost(time, state_unscaled,Model_data)}
    
    """optionsdict={'max_iter' : 200,
                 'acceptable_iter':20,
                 'print_level': 5,
                 'tol':1e-5,
                 'acceptable_tol':1e-4,
                 'dual_inf_tol':1e-4,
                 'acceptable_dual_inf_tol': 1e-3,
                 'constr_viol_tol':1e-7,
                 'acceptable_constr_viol_tol': 1e-5,
                 'compl_inf_tol': 1e-5,
                 'acceptable_compl_inf_tol':1e-2,
                 'acceptable_obj_change_tol':0.1,
                 'print_frequency_iter':2,
                 'nlp_scaling_method':'gradient-based',
                 'obj_scaling_factor':1.0}"""
    optionsdict={'max_iter' : 200,
                 'disp': True,}
    
    MPC_Class_Instance=MPC_Class(Model_data,MPC_params,Mesh,Functions,optionsdict)
    MPC_Class_Instance.Run_MPC_Loop()
    Result_data=Read_Result(os.getcwd()+'/Results/'+'Ex1_Test1')