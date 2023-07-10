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
from scipy.stats import chisquare
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

w_ps=8.002

a1=1#1/5**0.5
a2=1#2*a1

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

inf_file_name="beta_alpi8off_gamma_chicoarse_BETA_30Jun1447"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 2nd round/exp_CRG-3033/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt" 
beta_fold_clean=cleandata+"/Beta/AlphaON"
plots_fold=sorted_fold_path+"/Plots/"
correct_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 2nd round/exp_CRG-3033/Corrected data/"+inf_file_name+"/Beta/AlphaON"

if not os.path.exists(correct_fold_path):
    os.makedirs(correct_fold_path)
else:
    shutil.rmtree(correct_fold_path)
    os.makedirs(correct_fold_path)

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
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]/tot_data[:,1][tot_data[:,-1]==ps_pos[i]]
    matrix_err=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**-0.5*matrix[i]

ps_data=np.sum(matrix,axis=1)
P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 3,ps_pos[0]*3]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0, bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_ps=p[-2]
ps_0=p[-1]
print(w_ps)

c_data=np.sum(matrix,axis=0)
P0=[(np.amax(c_data)+np.amin(c_data))/2, np.amax(c_data)-np.amin(c_data), 48,c_pos[0]*48]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,c_pos,c_data, p0=P0, bounds=B0)
x_plt = np.linspace(c_pos[0], c_pos[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(c_pos,c_data,yerr=np.sqrt(c_data),fmt="ko",capsize=5)  
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_c=p[-2]
print(w_c)
c_0=p[-1]
# beta=np.linspace(-3*np.pi,3*np.pi,500)#c_pos.copy()#
# chi=np.linspace(-3*np.pi,3*np.pi,500)#ps_pos.copy()#
beta=w_c*c_pos-c_0
chi=w_ps*ps_pos-ps_0
# beta=c_pos
# chi=ps_pos
alpha=np.pi/8
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
print(w_c, w_ps)
P0=(c_0, ps_0, w_c, w_ps, 1, 1)
B0=([-10,-10,40,1,0,0],[10,10,55,4,np.inf,np.inf])
p,cov= fit(fit_I_px,range(len(matrix.ravel())),matrix.ravel()/np.amax(matrix.ravel()), bounds=B0) 
print(p)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(fit_I_px(0,*p), "b")
ax.plot(matrix.ravel()/np.amax(matrix.ravel()), "r--")
# ax.set_xlim([0,150])

print((np.sum(f_obs)-np.sum(f_exp))/np.sum(f_obs))
print(chisquare(f_obs=f_obs, f_exp=f_exp, ddof=7))

def I_px(x, beta0, chi0, w_c, w_ps, C, eta):
    beta=w_c*c_pos-beta0
    chi=w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px=I_px_co(beta, chi, C, alpha, 0) + I_px_in(beta, chi, eta, alpha, 0)
    # print(fit_I_px)
    return fit_I_px

def I_px_corr(x, beta0, chi0, w_c, w_ps, C, eta):
    beta = w_c*c_pos-beta0
    chi = w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px = I_px_co(beta, chi, C, alpha, 0) 
    # print(fit_I_px)
    return fit_I_px


fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
beta, chi = np.meshgrid(beta, chi)
Z= matrix
Z1= I_px(0,*p)*np.amax(matrix)
Z2 = I_px_corr(0, *p)*np.amax(matrix)
# Z=I_px_co(beta, chi, C, alpha, gamma)+I_px_in(beta, chi, eta, alpha, gamma)
ax.contour3D(beta, chi, Z, 30, cmap='binary')
ax.contour3D(beta, chi, Z1, 30, cmap='plasma')#cmap='Blues')
ax.set_xlabel('$\\beta$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 20)
plt.show()

# corrected_matrix=Z2#I_px_new(*p)*np.amax(matrix)
# corrected_matrix_err=np.sqrt(corrected_matrix)
# for i in range(len(ps_pos)):
#     data_txt=np.array([c_pos, corrected_matrix[i], corrected_matrix_err[i], np.ones(len(c_pos))*ps_pos[i]])
#     with open(correct_fold_path+"/beta_ps_"+str("%02d" % (i,))+".txt", 'w') as f:
#         np.savetxt(f, np.transpose(data_txt),  header= "Coil_pos O-Beam err ps_pos", fmt='%.7f %.7f %.7f %.7f' )

