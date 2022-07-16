#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 13:06:28 2022

@author: bayat2
"""

import numpy as np
import math
import os
import sys
#from LGR_Data import LGR_Class
from cyipopt import minimize_ipopt
from scipy.optimize import rosen, rosen_der
from scipy.interpolate import interp1d
from scipy.optimize import minimize
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from matplotlib.font_manager import FontProperties
import pickle as pkl
fontP=FontProperties()
#fontP.set_size('x-large')
from MPC_Class import Read_Result

Result_data_Test1=Read_Result(os.getcwd()+'/Results/'+'Ex2_Test1')
Result_data_Test2=Read_Result(os.getcwd()+'/Results/'+'Ex2_Test2')
Result_data_Test3=Read_Result(os.getcwd()+'/Results/'+'Ex2_Test3')
Exact=os.getcwd()+'/Results/'+'Ex2_Exact'
with open(Exact,'rb') as file:
    Exact_sol=pkl.load(file)
    
state_sol_ode=Exact_sol["state_sol_ode"]
u_opt_analytic=Exact_sol["u_opt_analytic"]


    
Ts1=Result_data_Test1["MPC_params"]["Ts"]
m1=Result_data_Test1["MPC_params"]["m"]
p1=Result_data_Test1["MPC_params"]["p"]

Ts2=Result_data_Test2["MPC_params"]["Ts"]
m2=Result_data_Test2["MPC_params"]["m"]
p2=Result_data_Test2["MPC_params"]["p"]

Ts3=Result_data_Test3["MPC_params"]["Ts"]
m3=Result_data_Test3["MPC_params"]["m"]
p3=Result_data_Test3["MPC_params"]["p"]

time_h_1=Result_data_Test1["I"]["time_h"]
time_h_2=Result_data_Test2["I"]["time_h"]
time_h_3=Result_data_Test3["I"]["time_h"]

fig_1, ax_1 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
ax_1.plot(time_h_1,Result_data_Test1["I"]["s_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts1}\, ,m={m1}\,\, ,p={p1}$')
ax_1.plot(time_h_2,Result_data_Test2["I"]["s_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts2}\, ,m={m2}\,\, ,p={p2}$')
ax_1.plot(time_h_3,Result_data_Test3["I"]["s_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts3}\, ,m={m3}\,\, ,p={p3}$')
ax_1.plot(time_h_3,state_sol_ode[:,0],'--',linewidth=2, label='$\mathrm{Exact}$')
ax_1.set_xlabel('$t\,\mathrm{[s]}$')
ax_1.set_ylabel(fr'$\xi_{1}$')
ax_1.legend(prop=fontP)

fig_2, ax_2 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
ax_2.plot(time_h_1,Result_data_Test1["I"]["s_unsc_matrix"][:,1],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts1}\, ,m={m1}\,\, ,p={p1}$')
ax_2.plot(time_h_2,Result_data_Test2["I"]["s_unsc_matrix"][:,1],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts2}\, ,m={m2}\,\, ,p={p2}$')
ax_2.plot(time_h_3,Result_data_Test3["I"]["s_unsc_matrix"][:,1],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts3}\, ,m={m3}\,\, ,p={p3}$')
ax_2.plot(time_h_3,state_sol_ode[:,1],'--',linewidth=2, label='$\mathrm{Exact}$')
ax_2.set_xlabel('$t\,\mathrm{[s]}$')
ax_2.set_ylabel(fr'$\xi_{2}$')
ax_2.legend(prop=fontP)

fig_3, ax_3 = plt.subplots(nrows=1, ncols=1,figsize=[6.4,4.8],dpi=200)
ax_3.plot(time_h_1,Result_data_Test1["I"]["u_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts1}\, ,m={m1}\,\, ,p={p1}$')
ax_3.plot(time_h_2,Result_data_Test2["I"]["u_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts2}\, ,m={m2}\,\, ,p={p2}$')
ax_3.plot(time_h_3,Result_data_Test3["I"]["u_unsc_matrix"][:,0],'-',linewidth=2,label=r'$\mathrm{MPC}\,\, ,T_s$'+f'$:{Ts3}\, ,m={m3}\,\, ,p={p3}$')
ax_3.plot(time_h_3,u_opt_analytic,'--',linewidth=2, label='$\mathrm{Exact}$')
ax_3.set_xlabel('$t\,\mathrm{[s]}$')
ax_3.set_ylabel(fr'$u$')
ax_3.legend(prop=fontP)