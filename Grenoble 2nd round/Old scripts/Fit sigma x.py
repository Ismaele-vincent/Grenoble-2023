#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 10:55:48 2023

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
def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)

"""
Reference no indium alpha 0
"""
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
"""
Sigma+
"""
max_pos_i=np.where(matrix==np.amax(matrix))[0][0]
max_pos=ps_pos[max_pos_i]
I_xp=matrix[max_pos_i,:]
I_xp_err=matrix_err[max_pos_i,:]

P0=[0, np.amax(I_xp), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
x_0p=p[-1]
p_a0p=p.copy()
err_a0p=np.diag(cov)**0.5
x_0p_err=err_a0p[-1]
sig_xp_0=(fit_cos(x_plt_0, *p_a0p)-fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0p)+fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.vlines(p[-1]-x_0p,0,fit_cos(p[-1]/p[-2], *p))
ax.errorbar(c_pos*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
# ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r")
ax.set_xlabel("$\\beta$")


"""
Sigma-
"""
min_pos_i=np.where(matrix==np.amin(matrix))[0][0]
min_pos=ps_pos[min_pos_i]
I_xm=matrix[min_pos_i,:]
I_xm_err=matrix_err[min_pos_i,:]
P0=[0, np.amax(I_xm), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
p_a0m=p.copy()
err_a0m=np.diag(cov)**0.5
x_0m=p[-1]
x_0m_err=err_a0m[-1]
sig_xm_0=(fit_cos(x_plt_0, *p_a0m)-fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0m)+fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.vlines(p[-1]-x_0m,0,fit_cos(p[-1]/p[-2], *p))

ax.errorbar(c_pos*p[-2]-x_0m,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p), "b")
ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
# ax.plot((x_plt_0)*p[-2]-x_0m,sig_xp_0, "r")
ax.set_xlabel("$\\beta$")


fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection='3d')
Z=matrix
x=c_pos
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.plot3D(c_pos,max_pos+c_pos*0, I_xp, 'red', lw=4)
ax.plot3D(c_pos,min_pos+c_pos*0, I_xm, 'blue', lw=4)
ax.set_xlabel('Coil')
ax.set_ylabel('PS')
ax.set_zlabel('z')
ax.view_init(45, 0)
plt.show()
#%%
"""
alpha=pi/4 path1
"""
#%%
inf_file_name="path1pi4cb_g_08Apr1819"
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
ps_i=109.4
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

max_pos_i=np.where(matrix==np.amax(matrix))[0][0]+1
max_pos=ps_pos[max_pos_i]
I_xp=matrix[max_pos_i,:]
I_xp_err=matrix_err[max_pos_i,:]
P0=[0, np.amax(I_xp), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
p_a0p=p.copy()
err_a0p=np.diag(cov)**0.5
bp_err=(err_a0p[-1]**2 + x_0p_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xp_pi4=(fit_cos(x_plt_0, *p_a0p)-fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0p)+fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta+$ $\\alpha=pi/4$ path 1")
ax.vlines(p[-1]-x_0p,np.amin(sig_xp_pi4),np.amax(sig_xp_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")

print("beta+_pi4_path1=",x_0p-p[-1],"+-", bp_err)

min_pos_i=np.where(matrix==np.amin(matrix))[0][0]-1
min_pos=ps_pos[min_pos_i]
I_xm=matrix[min_pos_i,:]
I_xm_err=matrix_err[min_pos_i,:]
P0=[0, np.amax(I_xm), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
p_a0m=p.copy()
err_a0m=np.diag(cov)**0.5
bm_err=(err_a0m[-1]**2 + x_0m_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xm_pi4=(fit_cos(x_plt_0, *p_a0m)-fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0m)+fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta-$ $\\alpha=pi/4$ path 1")
ax.vlines(p[-1]-x_0m,np.amin(sig_xm_pi4),np.amax(sig_xm_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")


print("beta-_pi4_path1=",x_0m-p_a0m[-1],"+-", bm_err)

fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection='3d')
Z=matrix
x=c_pos
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.plot3D(c_pos,max_pos+c_pos*0, I_xp, 'red', lw=4)
ax.plot3D(c_pos,min_pos+c_pos*0, I_xm, 'blue', lw=4)
ax.set_xlabel('Coil')
ax.set_ylabel('PS')
ax.set_zlabel('z')
ax.view_init(45, 0)

plt.show()
#%%
"""
alpha=pi/4 path2
"""
#%%
inf_file_name="path2pi4cb_g_10Apr1831"
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
ps_i=109.4
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

max_pos_i=np.where(matrix==np.amax(matrix))[0][0]
max_pos=ps_pos[max_pos_i]
I_xp=matrix[max_pos_i,:]
I_xp_err=matrix_err[max_pos_i,:]
P0=[0, np.amax(I_xp), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
p_a0p=p.copy()
err_a0p=np.diag(cov)**0.5
bp_err=(err_a0p[-1]**2 + x_0p_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xp_pi4=(fit_cos(x_plt_0, *p_a0p)-fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0p)+fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta+$ $\\alpha=pi/4$ path 2")
ax.vlines(p[-1]-x_0p,np.amin(sig_xp_pi4),np.amax(sig_xp_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")

print("beta+_pi4_path2=",x_0p-p[-1],"+-", bp_err)

min_pos_i=np.where(matrix==np.amin(matrix))[0][0]
min_pos=ps_pos[min_pos_i]
I_xm=matrix[min_pos_i,:]
I_xm_err=matrix_err[min_pos_i,:]
P0=[0, np.amax(I_xm), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
p_a0m=p.copy()
err_a0m=np.diag(cov)**0.5
bm_err=(err_a0m[-1]**2 + x_0m_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xm_pi4=(fit_cos(x_plt_0, *p_a0m)-fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0m)+fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta-$ $\\alpha=pi/4$ path 2")
ax.vlines(p[-1]-x_0m,np.amin(sig_xm_pi4),np.amax(sig_xm_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")


print("beta-_pi4_path2=",x_0m-p_a0m[-1],"+-", bm_err)

fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection='3d')
Z=matrix
x=c_pos
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.plot3D(c_pos,max_pos+c_pos*0, I_xp, 'red', lw=4)
ax.plot3D(c_pos,min_pos+c_pos*0, I_xm, 'blue', lw=4)
ax.set_xlabel('Coil')
ax.set_ylabel('PS')
ax.set_zlabel('z')
ax.view_init(45, 0)

plt.show()
#%%
"""
alpha=pi/8 path1
"""
#%%
inf_file_name="path1pi8cb_g_09Apr1441"
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
ps_i=109.4
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

max_pos_i=np.where(matrix==np.amax(matrix))[0][0]
max_pos=ps_pos[max_pos_i]
I_xp=matrix[max_pos_i,:]
I_xp_err=matrix_err[max_pos_i,:]
P0=[0, np.amax(I_xp), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
p_a0p=p.copy()
err_a0p=np.diag(cov)**0.5
bp_err=(err_a0p[-1]**2 + x_0p_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xp_pi4=(fit_cos(x_plt_0, *p_a0p)-fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0p)+fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta+$ $\\alpha=pi/8$ path 1")
ax.vlines(p[-1]-x_0p,np.amin(sig_xp_pi4),np.amax(sig_xp_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")

print("beta+_pi8_path1=",x_0p-p[-1],"+-", bp_err)

min_pos_i=np.where(matrix==np.amin(matrix))[0][0]-1
min_pos=ps_pos[min_pos_i]
I_xm=matrix[min_pos_i,:]
I_xm_err=matrix_err[min_pos_i,:]
P0=[0, np.amax(I_xm), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
p_a0m=p.copy()
err_a0m=np.diag(cov)**0.5
bm_err=(err_a0m[-1]**2 + x_0m_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xm_pi4=(fit_cos(x_plt_0, *p_a0m)-fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0m)+fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta-$ $\\alpha=pi/8$ path 1")
ax.vlines(p[-1]-x_0m,np.amin(sig_xm_pi4),np.amax(sig_xm_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")


print("beta-_pi8_path1=",x_0m-p_a0m[-1],"+-", bm_err)

fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection='3d')
Z=matrix
x=c_pos
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.plot3D(c_pos,max_pos+c_pos*0, I_xp, 'red', lw=4)
ax.plot3D(c_pos,min_pos+c_pos*0, I_xm, 'blue', lw=4)
ax.set_xlabel('Coil')
ax.set_ylabel('PS')
ax.set_zlabel('z')
ax.view_init(45, 0)

plt.show()
#%%
"""
alpha=pi/8 path2
"""
#%%
inf_file_name="path2pi8cb_g_12Apr1724"
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
ps_i=109.4
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

max_pos_i=np.where(matrix==np.amax(matrix))[0][0]
max_pos=ps_pos[max_pos_i]
I_xp=matrix[max_pos_i,:]
I_xp_err=matrix_err[max_pos_i,:]
P0=[0, np.amax(I_xp), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
p_a0p=p.copy()
err_a0p=np.diag(cov)**0.5
bp_err=(err_a0p[-1]**2 + x_0p_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xp_pi4=(fit_cos(x_plt_0, *p_a0p)-fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0p)+fit_cos(x_plt_0, *p_a0p[:-1],p_a0p[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta+$ $\\alpha=pi/8$ path 2")
ax.vlines(p[-1]-x_0p,np.amin(sig_xp_pi4),np.amax(sig_xp_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")

print("beta+_pi8_path2=",x_0p-p[-1],"+-", bp_err)

min_pos_i=np.where(matrix==np.amin(matrix))[0][0]+2
min_pos=ps_pos[min_pos_i]
I_xm=matrix[min_pos_i,:]
I_xm_err=matrix_err[min_pos_i,:]
P0=[0, np.amax(I_xm), 0.1,1]
p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
p_a0m=p.copy()
err_a0m=np.diag(cov)**0.5
bm_err=(err_a0m[-1]**2 + x_0m_err**2)**0.5
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
sig_xm_pi4=(fit_cos(x_plt_0, *p_a0m)-fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))/(fit_cos(x_plt_0, *p_a0m)+fit_cos(x_plt_0, *p_a0m[:-1],p_a0m[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("$\\beta-$ $\\alpha=pi/8$ path 2")
ax.vlines(p[-1]-x_0m,np.amin(sig_xm_pi4),np.amax(sig_xm_pi4))
# ax.errorbar(c_pos*p[-2]-x_0p,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_pi4, "r")
ax.plot((x_plt_0)*p[-2]-x_0m,sig_xm_0, "g")
ax.set_xlim([-1,1])
ax.set_xlabel("$\\beta$")


print("beta-_pi8_path2=",x_0m-p_a0m[-1],"+-", bm_err)

fig = plt.figure(figsize=(6,5))
ax = plt.axes(projection='3d')
Z=matrix
x=c_pos
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.plot3D(c_pos,max_pos+c_pos*0, I_xp, 'red', lw=4)
ax.plot3D(c_pos,min_pos+c_pos*0, I_xm, 'blue', lw=4)
ax.set_xlabel('Coil')
ax.set_ylabel('PS')
ax.set_zlabel('z')
ax.view_init(45, 0)

plt.show()

