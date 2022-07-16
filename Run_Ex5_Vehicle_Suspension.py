#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 18:57:28 2022

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
import pickle as pkl
from scipy.interpolate import CubicSpline
from MPC_Class import Read_Result
from scipy.io import loadmat
from scipy.integrate import solve_ivp
from matplotlib.font_manager import FontProperties
import pickle as pkl
fontP=FontProperties()


def Lagrange_User_Func(tau_unscaled_segment,states_segment_matrix_unscaled,u_extended_matrix_segment_unscaled,Model_data):
    state_1_segment_unscaled=states_segment_matrix_unscaled[:,0]
    state_2_segment_unscaled=states_segment_matrix_unscaled[:,1]
    state_3_segment_unscaled=states_segment_matrix_unscaled[:,2]
    state_4_segment_unscaled=states_segment_matrix_unscaled[:,3]
    
    u_1_segment_unscaled=u_extended_matrix_segment_unscaled[:,0]
    
    kt=Model_data["kt"]
    mU=Model_data["mU"]
    k1=Model_data["k1"]
    b1=Model_data["b1"]
    mS=Model_data["mS"]
    
    w1=Model_data["w1"]
    w2=Model_data["w2"]
    w3=Model_data["w3"]
    
    dq4=b1/mS*state_2_segment_unscaled-k1/mS*state_3_segment_unscaled-b1/mS*state_4_segment_unscaled+1/mS*u_1_segment_unscaled
    
    L=w1*state_1_segment_unscaled**2+w2*dq4**2+w3*u_1_segment_unscaled**2

    return L    

def Path_func(tau_iters_unscaled_vec,states_colloc_matrix_unscaled,u_extended_colloc_nodes_matrix_iter_unscaled,Model_data):
    #path_vec_1=(states_colloc_matrix_unscaled[:,0]+20.)/20.
    #path_vec_2=(states_colloc_matrix_unscaled[:,1]+20.)/20.
    #path_vec=np.append(path_vec_1,path_vec_2)
    path_vec=[]
    return path_vec      

def F_dynamics_user(tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data):
    F=np.zeros((len(tau_unscaled_segment),states_segment_matrix_colloc_unscaled.shape[1]))
    
    kt=Model_data["kt"]
    mU=Model_data["mU"]
    k1=Model_data["k1"]
    b1=Model_data["b1"]
    mS=Model_data["mS"]
    
    s1=states_segment_matrix_colloc_unscaled[:,0]
    s2=states_segment_matrix_colloc_unscaled[:,1]
    s3=states_segment_matrix_colloc_unscaled[:,2]
    s4=states_segment_matrix_colloc_unscaled[:,3]
    
    u=u_extended_matrix_segment_unscaled[:,0]
    
    zdot=Model_data["road_zdot_interp_func"](tau_unscaled_segment)
    
    F[:,0]=s2-zdot
    F[:,1]=-kt/mU*s1-b1/mU*s2+k1/mU*s3+b1/mU*s4-1/mU*u
    F[:,2]=-s2+s4
    F[:,3]=b1/mS*s2-k1/mS*s3-b1/mS*s4+1/mS*u
    
    return F                                  
                
                
def User_Mayer_Cost(time, state_unscaled,Model_data):
    mayer_cost=0.0
    if time==Model_data["t0"]:
        mayer_cost+=state_unscaled[0]*0.
    elif time==Model_data["tf"]:
        mayer_cost+=0
    return  mayer_cost 

if __name__ == '__main__':  
    
    Road_data=os.getcwd()+'/Road_profile/' + 'Road_data'
    with open(Road_data, 'rb') as file:
        road_t = pkl.load(file)
        road_x = pkl.load(file)
        road_z = pkl.load(file)
        road_zdot = pkl.load(file)
        v=10 # Vehicle speed [m/s]
        road_z_interp_func=CubicSpline(road_t,road_z)
        road_zdot_interp_func=CubicSpline(road_t,road_zdot)
        
    optimal_control_file='optimal_control'
    with open(optimal_control_file, 'rb') as file:
        t_opt = pkl.load(file)
        u_opt = pkl.load(file)
        
    u_opt_func=CubicSpline(t_opt,u_opt)
        
    r_max=0.04

    
    Model_data={
                "t0": 0.0,
                "tf":3.0,
                "mU":65.,
                "mS":325.0, 
                "kt":232.5e3, 
                "k1":37.497e3, 
                "b1":367.5926, 
                "w1":1e5,
                "w2":0.5,
                "w3":1e-5,
                "r_max":r_max,
                "road_z_interp_func": road_z_interp_func,
                "road_zdot_interp_func":road_zdot_interp_func,
                "t_opt":t_opt,
                "u_opt":u_opt,
                "u_opt_func":u_opt_func,
                "A_scaled":np.array([0.0,0.0,0.0,0.0]),
                "B_scaled":np.array([1.0,1.0,r_max,1.0]),
                "C_scaled":np.array([0.0]),
                "D_scaled":np.array([1000.0]),
                "Initial_State":np.array([0.0,0.0,0.0,0.0]), #scaled
                "End_State":np.array([0.0,0.0,0.0,0.0]),      #scaled
                "state_lb_guess":np.array([-0.1,-0.1,-0.1,-0.1]), #scaled
                "state_ub_guess":np.array([0.1,0.1,0.1,0.1]), #scaled
                "control_lb_guess":np.array([-0.1]), #scaled
                "control_ub_guess":np.array([0.1]),  #scaled
                "Init_state_boudnary_flag":1,
                "End_state_boudnary_flag":0,
                "plot_iter_flag":False,
                "File_Name":os.getcwd()+'/Results/'+'Ex5_Test2',
                }
    
    MPC_params={
                "Ts": 0.1,
                "m": 30,
                "p":30,
                }
    
    Mesh={
        "delta_t_seg":0.1,
        "num_nodes_seg":7,
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
    Result_data=Read_Result(os.getcwd()+'/Results/'+'Ex5_Test2')
    Result_data["Model_data"]=Model_data 
    
    Ts=Result_data["MPC_params"]["Ts"]
    m=Result_data["MPC_params"]["m"]
    p=Result_data["MPC_params"]["p"]
    
    GPOPS_Result=loadmat("optimal_control.mat")
    t_GPOPS=GPOPS_Result["t"].flatten()
    u_GPOPS=GPOPS_Result["u"].flatten()
    
    u_interp_analytic_func=interp1d(t_GPOPS,u_GPOPS,axis=0,fill_value="extrapolate")
    initial_condition=np.array([0.0,0.0,0.0,0.0])
    
    def dynamic_ode(t,state,u_interp_analytic_func,Model_data):
        state_dot=np.zeros(len(state))
        
        kt=Model_data["kt"]
        mU=Model_data["mU"]
        k1=Model_data["k1"]
        b1=Model_data["b1"]
        mS=Model_data["mS"]
        
        s1=state[0]
        s2=state[1]
        s3=state[2]
        s4=state[3]
        
        u=u_interp_analytic_func(t)
        
        zdot=Model_data["road_zdot_interp_func"](t)
        
        state_dot[0]=s2-zdot
        state_dot[1]=-kt/mU*s1-b1/mU*s2+k1/mU*s3+b1/mU*s4-1/mU*u
        state_dot[2]=-s2+s4
        state_dot[3]=b1/mS*s2-k1/mS*s3-b1/mS*s4+1/mS*u
        return state_dot
    
    time_h=Result_data["I"]["time_h"]
    
    dyn_fun=lambda t, state: dynamic_ode(t,state,u_interp_analytic_func,Model_data)
    t_span=(time_h[0],time_h[-1])
    
    
    sol_ode=solve_ivp(dyn_fun, t_span, y0=initial_condition, method='RK45', t_eval=time_h)
    
    time_sol_ode=sol_ode.t
    state_sol_ode=sol_ode.y.T
    
    Exact_Result=dict()
    Exact_Result["time_sol_ode"]=time_sol_ode
    Exact_Result["state_sol_ode"]=state_sol_ode
    Exact_Result["u_opt_analytic"]=u_interp_analytic_func(time_sol_ode)
    
    file_ExactSol_addrees=os.getcwd()+'/Results/'+'Ex5_Exact'
    with open(file_ExactSol_addrees,'wb') as file:
        pkl.dump(Exact_Result,file)
    
    fig_1, ax_1 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_1.plot(time_h,Result_data["I"]["s_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_1.plot(time_h,state_sol_ode[:,0],'-',linewidth=1, label='$\mathrm{Exact}$')
    ax_1.set_xlabel('$t\,\mathrm{[s]}$')
    ax_1.set_ylabel(fr'$\xi_1$')
    ax_1.legend(prop=fontP)
    #fig_1.suptitle('$States$')
    
    fig_2, ax_2 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_2.plot(time_h,Result_data["I"]["s_unsc_matrix"][:,1],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_2.plot(time_h,state_sol_ode[:,1],'-',linewidth=1, label='$\mathrm{Exact}$')
    ax_2.set_xlabel('$t\,\mathrm{[s]}$')
    ax_2.set_ylabel(fr'$\xi_2$')
    ax_2.legend(prop=fontP)
    
    fig_3, ax_3 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_3.plot(time_h,Result_data["I"]["s_unsc_matrix"][:,2],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_3.plot(time_h,state_sol_ode[:,2],'-',linewidth=1, label='$\mathrm{Exact}$')
    ax_3.set_xlabel('$t\,\mathrm{[s]}$')
    ax_3.set_ylabel(fr'$\xi_3$')
    ax_3.legend(prop=fontP)
    
    fig_4, ax_4 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_4.plot(time_h,Result_data["I"]["s_unsc_matrix"][:,3],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_4.plot(time_h,state_sol_ode[:,3],'-',linewidth=1, label='$\mathrm{Exact}$')
    ax_4.set_xlabel('$t\,\mathrm{[s]}$')
    ax_4.set_ylabel(fr'$\xi_4$')
    ax_4.legend(prop=fontP)
    
    
    fig_5, ax_5 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_5.plot(time_h,Result_data["I"]["u_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_5.plot(time_h,u_interp_analytic_func(time_sol_ode),'-',linewidth=1, label='$\mathrm{Exact}$')
    ax_5.set_xlabel('$t\,\mathrm{[s]}$')
    ax_5.set_ylabel(fr'$u$')
    ax_5.legend(prop=fontP)
    #fig_2.suptitle('$States$')
    