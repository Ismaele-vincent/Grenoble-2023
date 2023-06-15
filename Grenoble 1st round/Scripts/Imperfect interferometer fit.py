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
ps_i=109
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

ps_data=np.sum(matrix,axis=1)
P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 8.12,ps_pos[0]*8.11]
B0=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0, bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
# ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_ps=p[-2]
ps_0=p[-1]
print(ps_0)

c_data=np.sum(matrix,axis=0)
P0=[(np.amax(c_data)+np.amin(c_data))/2, np.amax(c_data)-np.amin(c_data), 0.1,c_pos[0]]
B0=([0,0,0,0],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,c_pos,c_data, p0=P0, bounds=B0)
x_plt = np.linspace(c_pos[0], c_pos[-1],100)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.errorbar(c_pos,c_data,yerr=np.sqrt(c_data),fmt="ko",capsize=5)  
# ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_c=p[-2]
print(p[1]/p[0])
c_0=p[-1]
# beta=np.linspace(-3*np.pi,3*np.pi,500)#c_pos.copy()#
# chi=np.linspace(-3*np.pi,3*np.pi,500)#ps_pos.copy()#
beta=w_c*c_pos-c_0
chi=w_ps*ps_pos-ps_0
alpha=0#np.pi/8
gamma=0
C=0.8
eta=1-C

def fit_I_px(x, beta0, chi0, w_c, w_ps, C, eta):
    beta=w_c*c_pos-beta0
    chi=w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px=I_px_co(beta, chi, C, alpha, 0) + I_px_in(beta, chi, eta, alpha, 0)
    # print(fit_I_px)
    return fit_I_px.ravel()

P0=(c_0, ps_0, w_c, w_ps, 1, 1)
B0=([0,700,0,7,0,0],[1000,1000,0.3,10,np.inf,np.inf])
p,cov= fit(fit_I_px,range(len(matrix.ravel())),matrix.ravel()/np.amax(matrix.ravel()), bounds=B0) 
print(p)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(fit_I_px(0,*p), "b")
ax.plot(matrix.ravel()/np.amax(matrix.ravel()), "r--")
ax.set_xlim([0,150])
def I_px(x, beta0, chi0, w_c, w_ps, C, eta):
    beta=w_c*c_pos-beta0
    chi=w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px=I_px_co(beta, chi, C, alpha, 0) + I_px_in(beta, chi, eta, alpha, 0)
    # print(fit_I_px)
    return fit_I_px

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
beta, chi = np.meshgrid(beta, chi)
Z= matrix
Z1= I_px(0,*p)*np.amax(matrix)
# Z=I_px_co(beta, chi, C, alpha, gamma)+I_px_in(beta, chi, eta, alpha, gamma)
ax.contour3D(beta, chi, Z, 30, cmap='binary')
ax.contour3D(beta, chi, Z1, 30, cmap='plasma')#cmap='Blues')
ax.set_xlabel('$\\beta$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 45)
plt.show()