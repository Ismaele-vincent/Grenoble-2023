#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 12:15:08 2023

@author: aaa
"""

import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
plt.rcParams.update({'figure.max_open_warning': 0})
from PIL import Image as im
from scipy.optimize import curve_fit as fit
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

w_ps=8.002

a1=1/5**0.5
a2=2*a1

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)

def I_px_co(beta, chi, C, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    return C*((a1*np.cos(d/2))**2+(a2*np.cos(e/2))**2+2*a1*a2*np.cos(d/2)*np.cos(e/2)*np.cos(chi))/4

def I_px_in(beta, chi, eta, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    return eta*(np.cos(d/2)**2+(a2/a1)**2*np.cos(e/2)**2)/4

def I_mx_co(beta, chi, C, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    r=np.arctan(gamma/beta)
    p=np.arctan(gamma/(beta+gamma))
    return C*((a1*np.sin(d/2))**2+(a2*np.sin(e/2))**2+2*a1*a2*np.sin(d/2)*np.sin(e/2)*np.cos(chi+r-p))/4

def I_mx_in(beta, chi, eta, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    return eta*(np.sin(d/2)**2+(a2/a1)**2*np.sin(e/2)**2)/4

inf_file_name="path1pi4cb_g_13Apr1502"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt" 
beta_fold_clean=cleandata+"/Beta"
plots_fold=sorted_fold_path+"/Plots/"
i=0
for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
    files=np.sort(files)
    for name in files[:-1]:
        if i==0:
            tot_data=np.loadtxt(os.path.join(root, name))
            c_pos=tot_data[:,0]
            i=1
        else:
            data=np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos=tot_data[::len(c_pos),-1]
# ps_i=109
# ps_f=ps_pos[-1]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

beta=c_pos.copy()#np.linspace(-3*np.pi,3*np.pi,500)
alpha=np.pi/8
chi=ps_pos.copy()#np.linspace(-3*np.pi,3*np.pi,500)
gamma=0
C=0.8
eta=1-C

def fit_I_px(x,beta0, chi0, C):
    X, Y = np.meshgrid(beta+beta0, chi+chi0)
    fit_I_px=I_px_co(X, Y, C, alpha, gamma)+I_px_in(X, Y, C, alpha, gamma)
    return fit_I_px
# print(I_px_co(beta, chi, C, alpha, gamma))
beta, chi = np.meshgrid(beta, chi)

I_px=I_px_co(beta, chi, C, alpha, gamma)+I_px_in(beta, chi, eta, alpha, gamma)
I_mx=I_mx_co(beta, chi, C, alpha, gamma)+I_mx_in(beta, chi, eta, alpha, gamma)
I_x=(I_px-I_mx)/(I_px+I_mx)
I_x_th=(I_px_co(beta, chi, C, alpha, gamma)-I_mx_co(beta, chi, C, alpha, gamma))/(I_px_co(beta, chi, C, alpha, gamma)+I_mx_co(beta, chi, C, alpha, gamma))
x=range(len(matrix.ravel()))
Z= fit_I_px(x,0, 0, C)

# p,cov= fit(fit_I_px,range(len(matrix.ravel())),matrix.ravel()) 

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
# Z=I_px_co(beta, chi, C, alpha, gamma)+I_px_in(beta, chi, eta, alpha, gamma)
ax.contour3D(beta, chi, Z, 30, cmap='binary')
ax.set_xlabel('$\\beta$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 0)

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.plot(beta, I_x, "b")
# ax.plot(beta, I_x_th, "r")



