#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 22:50:40 2022

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
from scipy.interpolate import CubicSpline
from scipy.interpolate import RectBivariateSpline
from MPC_Class import Read_Result
from scipy.io import loadmat
from scipy.integrate import solve_ivp
from matplotlib.font_manager import FontProperties
import pickle as pkl
fontP=FontProperties()

def Lagrange_User_Func(tau_unscaled_segment,states_segment_matrix_unscaled,u_extended_matrix_segment_unscaled,Model_data):
    state_1_segment_unscaled=states_segment_matrix_unscaled[:,0]
    u_1_segment_unscaled=u_extended_matrix_segment_unscaled[:,0]
    u_2_segment_unscaled=u_extended_matrix_segment_unscaled[:,1]
    
    
    R_rotor=Model_data["R_rotor"]
    rho_air=Model_data["rho_air"]
    Irx=Model_data["Irx"]
    P_max=Model_data["P_max"]
    etha=Model_data["etha"]
    v=Model_data["v"]
    Cpfunc=Model_data["Cpfunc"]
    Ctfunc=Model_data["Ctfunc"]

    wind_speed=v(tau_unscaled_segment)
    lambda_w=R_rotor*state_1_segment_unscaled/wind_speed
    
    Cp_vec=np.zeros(len(tau_unscaled_segment))
    for time_index in np.arange(len(tau_unscaled_segment)):
        Cp_vec[time_index]=Cpfunc(lambda_w[time_index],u_2_segment_unscaled[time_index])
    
    tau_a=1/2*rho_air*np.pi*R_rotor**3*Cp_vec/lambda_w*wind_speed**2  #%aerodynamic torque
    P_a=Cp_vec*0.5*rho_air*np.pi*R_rotor**2*wind_speed**3
    
    L=-((P_a)-1.*(1.0*1e-9*u_1_segment_unscaled**2+0*1e8*u_2_segment_unscaled**2))*etha/P_max/Model_data["tf"] #[W/W]
    
    return L    

def Path_func(tau_iters_unscaled_vec,states_colloc_matrix_unscaled,u_extended_colloc_nodes_matrix_iter_unscaled,Model_data):
    #path_vec_1=(states_colloc_matrix_unscaled[:,0]+20.)/20.
    #path_vec_2=(states_colloc_matrix_unscaled[:,1]+20.)/20.
    #path_vec=np.append(path_vec_1,path_vec_2)
    path_vec=[]
    return path_vec      

def F_dynamics_user(tau_unscaled_segment,states_segment_matrix_colloc_unscaled,u_extended_matrix_segment_unscaled,Model_data):
    
    F=np.zeros((len(tau_unscaled_segment),states_segment_matrix_colloc_unscaled.shape[1]))
    state_1_segment_unscaled=states_segment_matrix_colloc_unscaled[:,0]
    u_1_segment_unscaled=u_extended_matrix_segment_unscaled[:,0]
    u_2_segment_unscaled=u_extended_matrix_segment_unscaled[:,1]
    
    R_rotor=Model_data["R_rotor"]
    rho_air=Model_data["rho_air"]
    Irx=Model_data["Irx"]
    P_max=Model_data["P_max"]
    etha=Model_data["etha"]
    v=Model_data["v"]
    Cpfunc=Model_data["Cpfunc"]
    Ctfunc=Model_data["Ctfunc"]
    
    wind_speed=v(tau_unscaled_segment)
    lambda_w=R_rotor*state_1_segment_unscaled/wind_speed
    
    Cp_vec=np.zeros(len(tau_unscaled_segment))
    for time_index in np.arange(len(tau_unscaled_segment)):
        Cp_vec[time_index]=Cpfunc(lambda_w[time_index],u_2_segment_unscaled[time_index])
    
    tau_a=1/2*rho_air*np.pi*R_rotor**3*Cp_vec/lambda_w*wind_speed**2 #aerodynamic torque
    P_a=Cp_vec*0.5*rho_air*np.pi*R_rotor**2*wind_speed**3

    F[:,0]=(tau_a-u_1_segment_unscaled)/Irx
    
    return F   

def User_Mayer_Cost(time, state_unscaled,Model_data):
    mayer_cost=0.0
    if time==Model_data["t0"]:
        mayer_cost+=state_unscaled[0]*0.
    elif time==Model_data["tf"]:
        mayer_cost+=state_unscaled[0]*0.
    return  mayer_cost                                                
                
if __name__ == '__main__':  
    
    wind_profile_30 = loadmat('wind_profile_30.mat')
    t_interp=wind_profile_30['t_interp'][:,0]
    v_interp=wind_profile_30['v_interp'][:,0]
    
    wind_avg=6.
    v=CubicSpline(t_interp/4,wind_avg*1/14.4607*v_interp)
    
    CP_CT_MAT = loadmat('CP_CT_MAT.mat')
    CP_Mat=CP_CT_MAT['CP_Mat']
    CT_Mat=CP_CT_MAT['CT_Mat']
    Lambda=CP_CT_MAT['Lambda']
    Theta_p=CP_CT_MAT['Theta_p']
    THETA_P,LAMBDA = np.meshgrid(Theta_p,Lambda)
    Cpfunc=RectBivariateSpline(Lambda, Theta_p, CP_Mat, kx=3, ky=3, s=0)
    Ctfunc=RectBivariateSpline(Lambda, Theta_p, CT_Mat, kx=3, ky=3, s=0)
    
    Model_data={
                "t0": 0.0,
                "tf":100.0,
                "v":v,
                "Cpfunc":Cpfunc,
                "Ctfunc":Ctfunc,
                "Irx":38759228,
                "rho_air":1.225,
                "omega_max":1.2671,
                "omega_min":0.3,
                "etha":0.944,
                "P_max":5.0e6,
                "R_rotor":63,
                "A_scaled":np.array([0.7835]),
                "B_scaled":np.array([0.4835]),
                "C_scaled":np.array([4.18e6/2,0.6807/2]),
                "D_scaled":np.array([4.18e6/2,0.6807/2]),
                "Initial_State":np.array([-0.2]), #scaled -0.2 0.4759
                "End_State":np.array([0.0]),      #scaled
                "state_lb_guess":np.array([0.5]), #scaled
                "state_ub_guess":np.array([0.8]), #scaled
                "control_lb_guess":np.array([0.5,-1.0]), #scaled
                "control_ub_guess":np.array([0.8,-0.8]),  #scaled
                "Init_state_boudnary_flag":1,
                "End_state_boudnary_flag":0,
                "plot_iter_flag":False,
                "File_Name":os.getcwd()+'/Results/'+'Ex4_Test2',
                }
    
    MPC_params={
                "Ts": 1,
                "m": 5,
                "p":5,
                }
    
    Mesh={
        "delta_t_seg":5,
        "num_nodes_seg":8,
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
    Result_data=Read_Result(os.getcwd()+'/Results/'+'Ex4_Test2')
    Result_data["Model_data"]=Model_data 
    
    s1=Result_data["I"]["s_unsc_matrix"][:,0]
    
    u1=Result_data["I"]["u_unsc_matrix"][:,0]
    u2=Result_data["I"]["u_unsc_matrix"][:,1]
    
    
    Ts=Result_data["MPC_params"]["Ts"]
    m=Result_data["MPC_params"]["m"]
    p=Result_data["MPC_params"]["p"]
    
    time_h=Result_data["I"]["time_h"]
    
    R_rotor=Model_data["R_rotor"]
    rho_air=Model_data["rho_air"]
    Irx=Model_data["Irx"]
    P_max=Model_data["P_max"]
    etha=Model_data["etha"]
    v=Model_data["v"]
    Cpfunc=Model_data["Cpfunc"]
    Ctfunc=Model_data["Ctfunc"]

    wind_speed=v(time_h)
    lambda_w=R_rotor*s1/wind_speed
    
    Cp_vec=np.zeros(len(time_h))
    for time_index in np.arange(len(time_h)):
        Cp_vec[time_index]=Cpfunc(lambda_w[time_index],u2[time_index])
    
    tau_a=1/2*rho_air*np.pi*R_rotor**3*Cp_vec/lambda_w*wind_speed**2  #%aerodynamic torque
    P_a=Cp_vec*0.5*rho_air*np.pi*R_rotor**2*wind_speed**3
    
    L=-((P_a)-1.*(1e-9*u1**2+0*1e8*u2**2))*etha/P_max/Model_data["tf"] #[W/W]
    
    Objective=np.trapz(L,time_h)
    
    fig_1, ax_1 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_1.plot(time_h,s1,'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_1.set_xlabel('$t\,\mathrm{[s]}$')
    ax_1.set_ylabel(r'$\omega \mathrm{[rad/s]}$')
    ax_1.legend(prop=fontP)
    
    fig_2, ax_2 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_2.plot(time_h,Cp_vec,'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_2.set_xlabel('$t\,\mathrm{[s]}$')
    ax_2.set_ylabel(r'$C_p$')
    ax_2.legend(prop=fontP)
    
    print(Objective)
    
    fig_3, ax_3 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_3.plot(time_h,u1,'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_3.set_xlabel('$t\,\mathrm{[s]}$')
    ax_3.set_ylabel(r'$T_{\mathrm{gen}}\,\mathrm{[N.m]}$')
    ax_3.legend(prop=fontP)
    
    fig_4, ax_4 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_4.plot(time_h,u2,'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts}\, ,m={m}\,\, ,p={p}$')
    ax_4.set_xlabel('$t\,\mathrm{[s]}$')
    ax_4.set_ylabel(r'$\theta_{\mathrm{blade}}\, \mathrm{[rad]}$')
    ax_4.legend(prop=fontP)
    
    fig_5, ax_5 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
    ax_5.plot(time_h,wind_speed,'-',linewidth=2)
    ax_5.set_xlabel('$t\,\mathrm{[s]}$')
    ax_5.set_ylabel(r'$V_{\mathrm{wind}}\,\mathrm{[m/s]}$')
    