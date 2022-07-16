#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 15:22:06 2022

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
from MPC_Class import MPC_Class
from MPC_Class import Read_Result
from scipy.integrate import solve_ivp
from matplotlib.font_manager import FontProperties
import pickle as pkl
fontP=FontProperties()


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
    F[:,1]=u_extended_matrix_segment_unscaled[:,0]
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
                "t0": 0.0,
                "tf":1.0,
                "l":1/9,
                "A_scaled":np.array([0.0,0.0]),
                "B_scaled":np.array([1/9,2.0]),
                "C_scaled":np.array([0.0]),
                "D_scaled":np.array([20.0]),
                "Initial_State":np.array([0.0,0.5]), #scaled
                "End_State":np.array([0.0,-0.5]),      #scaled
                "state_lb_guess":np.array([-0.1,-0.1]), #scaled
                "state_ub_guess":np.array([0.1,0.1]), #scaled
                "control_lb_guess":np.array([-0.1]), #scaled
                "control_ub_guess":np.array([0.1]),  #scaled
                "Init_state_boudnary_flag":1,
                "End_state_boudnary_flag":1,
                "plot_iter_flag":True,
                "File_Name":os.getcwd()+'/Results/'+'Ex2_Test3',
                }
    
    MPC_params={
                "Ts": 0.01,
                "m": 100,
                "p":100,
                }
    
    Mesh={
        "delta_t_seg":0.2,
        "num_nodes_seg":6,
        "n_h_Sol":100}
    
    Functions={
        "Lagrange_User_Func": lambda t_unsc_seg,state_unsc_seg,u_unsc_seg,Model_data :\
                                Lagrange_User_Func(t_unsc_seg,state_unsc_seg,u_unsc_seg,Model_data),\
        "Path_func": lambda t_unsc_iter,state_unsc_iter,u_unsc_iter,Model_data :\
                                Path_func(t_unsc_iter,state_unsc_iter,u_unsc_iter,Model_data),
        "F_dynamics_user": lambda tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data:\
                                F_dynamics_user(tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data),\
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
    optionsdict={'maxiter' : 200,
                 'disp': True,}
    
    MPC_Class_Instance=MPC_Class(Model_data,MPC_params,Mesh,Functions,optionsdict)
    MPC_Class_Instance.Run_MPC_Loop()
    
    #-------------------Read Result---------------------
    
    def u_analytic_func(Result_data,time_h_scalar):
        tf=Result_data["Model_data"]["tf"]
        l=Result_data["Model_data"]["l"]
        
        if time_h_scalar<3*l:
            u_star=-2/(3*l)*(1-time_h_scalar/(3*l))
        elif time_h_scalar<1-3*l:
            u_star=0
        else:
            u_star=-2/(3*l)*(1-(1-time_h_scalar)/(3*l))
        
        return u_star
    
    def dynamic_ode(t,state,u_interp_analytic_func):
        state_dot=np.zeros(len(state))
        state_dot[0]=state[1]
        state_dot[1]=u_interp_analytic_func(t)
        return state_dot
    
    Result_data=Read_Result(os.getcwd()+'/Results/'+'Ex2_Test3')
    Result_data["Model_data"]=Model_data
    
    Ts=Result_data["MPC_params"]["Ts"]
    m=Result_data["MPC_params"]["m"]
    p=Result_data["MPC_params"]["p"]
    
    # Optimal control Analytical Solution
    time_h=Result_data["I"]["time_h"]
    u_opt_analytic=np.zeros(len(time_h))
    for time_index in np.arange(len(Result_data["I"]["time_h"])):
        time_h_scalar=Result_data["I"]["time_h"][time_index]
        u_opt_analytic[time_index]=u_analytic_func(Result_data,time_h_scalar)
        
    u_interp_analytic_func=interp1d(time_h,u_opt_analytic,axis=0)
    initial_condition=np.array([0,1])
    
    dyn_fun=lambda t, state: dynamic_ode(t,state,u_interp_analytic_func)
    t_span=(time_h[0],time_h[-1])
    
    sol_ode=solve_ivp(dyn_fun, t_span, y0=initial_condition, method='RK45', t_eval=time_h)
    
    time_sol_ode=sol_ode.t
    state_sol_ode=sol_ode.y.T
    
    Exact_Result=dict()
    Exact_Result["time_sol_ode"]=time_sol_ode
    Exact_Result["state_sol_ode"]=state_sol_ode
    Exact_Result["u_opt_analytic"]=u_opt_analytic
    
    file_ExactSol_addrees=os.getcwd()+'/Results/'+'Ex2_Exact'
    with open(file_ExactSol_addrees,'wb') as file:
        pkl.dump(Exact_Result,file)
    
    fig_1, ax_1 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_1.plot(time_h,Result_data["I"]["s_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_1.plot(time_h,state_sol_ode[:,0],'--',linewidth=2, label='$\mathrm{Exact}$')
    ax_1.set_xlabel('$t\,\mathrm{[s]}$')
    ax_1.set_ylabel(fr'$\xi_{1}$')
    ax_1.legend(prop=fontP)
    #fig_1.suptitle('$States$')
    
    fig_2, ax_2 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_2.plot(time_h,Result_data["I"]["s_unsc_matrix"][:,1],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_2.plot(time_h,state_sol_ode[:,1],'--',linewidth=2, label='$\mathrm{Exact}$')
    ax_2.set_xlabel('$t\,\mathrm{[s]}$')
    ax_2.set_ylabel(fr'$\xi_{2}$')
    ax_2.legend(prop=fontP)
    #fig_2.suptitle('$States$')
    
    fig_3, ax_3 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_3.plot(time_h,Result_data["I"]["u_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_3.plot(time_h,u_opt_analytic,'--',linewidth=2, label='$\mathrm{Exact}$')
    ax_3.set_xlabel('$t\,\mathrm{[s]}$')
    ax_3.set_ylabel(fr'$u$')
    ax_3.legend(prop=fontP)
    #fig_2.suptitle('$States$')