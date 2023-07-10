#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 17:53:20 2023

@author: aaa
"""

import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
plt.rcParams.update({'figure.max_open_warning': 0})
from PIL import Image as im
from scipy.optimize import curve_fit as fit

alpha=22.5
w_ps=8.005
a21=1
def fit_w1p(x,th,x0):
    return alpha*((1/(1+a21*np.tan(th*np.pi/4)**2*np.exp(-1j*(w_ps*(x-x0)))))).imag

def exp_w1p(x,x0):
    return alpha*((1/(1+a21*np.exp(-1j*(w_ps*(x-x0)))))).imag

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*(x-D))

rad=np.pi/180
inf_file_name="path1pi4cb_g_13Apr1502"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt" 
beta_fold_clean=cleandata+"/Beta"
plots_fold=sorted_fold_path+"/Plots/"
i=0
for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
    files=np.sort(files)[:-1]
    for name in files:
        if i==0:
            tot_data=np.loadtxt(os.path.join(root, name))
            c_amp=tot_data[:,0]
            i=1
        else:
            data=np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos=tot_data[::len(c_amp),-1]
ps_i=109
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 

matrix=np.zeros((len(ps_pos),len(c_amp)))
max_c_amp=np.zeros(len(ps_pos))
beta=np.zeros(len(ps_pos))
w=np.zeros(len(ps_pos))
fit_res0=[503.8204133,  483.67092569,   2.79781313,  -0.5601]
err_res0=[7.24106172, 10.30876345, 0.01913184, 0.01161413]
err_g=np.zeros(len(ps_pos))
max0=fit_res0[-1]
max_int=np.zeros(len(ps_pos))
min_int=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
for i in range(len(ps_pos)):
    P0=[0, np.amax(matrix[i]), 0.1,-0.5]
    B0=([0,0,0.001,-2],[np.inf,np.inf,np.inf,2])
    p,cov=fit(fit_cos,c_amp,matrix[i], p0=P0, bounds=B0)
    
    err=np.diag(cov)**0.5
    # print(p)
    x_plt = np.linspace(c_amp[0], c_amp[-1],100)
    max_int[i]=p[0]+p[1]
    min_int[i]=p[0]-p[1]
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # fig.suptitle("ps_pos="+str(ps_pos[i]))
    # ax.errorbar(c_amp,matrix[i],yerr=np.sqrt(matrix[i]),fmt="ko",capsize=5)
    # ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed",color="b",label="$\beta$="+str("%.3f" % (max_c_amp[i]),))
    # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # ax.set_ylim([0, P0[1]+P0[1]/10])
    # ax.legend(loc=4)
    
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
P0=[0, np.amax(max_int), 8,109.3]
B0=([0,0,0.001,100],[np.inf,np.inf,np.inf,150])
p,cov=fit(fit_cos,ps_pos,max_int, p0=P0, bounds=B0)
print(p)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
fig.suptitle("Max")
ax.errorbar(ps_pos,max_int,yerr=np.sqrt(max_int),fmt="ko",capsize=5)
ax.plot(x_plt,fit_cos(x_plt, *p), "b")

P0=[0, np.amax(max_int), 8,109.3]
B0=([0,0,0.001,100],[np.inf,np.inf,np.inf,150])
p,cov=fit(fit_cos,ps_pos,min_int, p0=P0, bounds=B0)
print(p)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
fig.suptitle("Min")
ax.errorbar(ps_pos,min_int,yerr=np.sqrt(min_int),fmt="ko",capsize=5)
ax.plot(x_plt,fit_cos(x_plt, *p), "b")

fig = plt.figure(figsize=(5,5))
fig.suptitle("Max/min")
ax = fig.add_subplot(111)
# ax.plot(ps_pos, max_int)
# ax.plot(ps_pos, min_int)
ax.plot(ps_pos, max_int/min_int)
    