#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:57:13 2023

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
    for name in files:
        if i==0:
            tot_data=np.loadtxt(os.path.join(root, name))
            c_pos=tot_data[:,0]
            i=1
        else:
            data=np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos=tot_data[::len(c_pos),-1]
ps_i=108.8
ps_f=ps_pos[-1]
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

x_plt_ifg=np.linspace(ps_pos[0],ps_pos[-1],1000)
ifg=np.sum(matrix,axis=1)
P0=[np.amax(ifg), np.amax(ifg), 8,ps_pos[ifg==np.amax(ifg)][0]*8]
p,cov=fit(fit_cos,ps_pos,ifg, p0=P0)
ps_0=ps_pos[abs(ps_pos-p[-1]/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]))][0]
ps_pi=ps_pos[abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2]))][0]
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.vlines(ps_0, 0, ifg[ps_pos==ps_0])
# ax.vlines(ps_pi, 0, ifg[ps_pos==ps_pi])
# ax.errorbar(ps_pos,ifg,yerr=np.sqrt(ifg),fmt="ko",capsize=5)
# ax.plot(x_plt_ifg,fit_cos(x_plt_ifg, *p), "b")

"""
Sigma+
"""
max_pos_i=np.where(ps_pos==ps_0)[0][0]
max_pos=ps_pos[max_pos_i]
I_xp=matrix[max_pos_i,:]
I_xp_err=matrix_err[max_pos_i,:]
P0=[np.amax(I_xp)/2, np.amax(I_xp), 0.1,0]
p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
p_a0p=p.copy()
err_a0p=np.diag(cov)**0.5
x_0p=p[-1]
print(p[-1]/p[-2])
x_0p_err=err_a0p[-1]
sig_xp_0=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.vlines(0,0,fit_cos(p[-1]/p[-2], *p))
ax.errorbar((c_pos)*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
# ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r")
ax.set_xlabel("$\\beta$")

"""
Sigma-
"""
min_pos_i=np.where(ps_pos==ps_pi)[0][0]
min_pos=ps_pos[min_pos_i]
I_xm=matrix[min_pos_i,:]
I_xm_err=matrix_err[min_pos_i,:]
P0=[np.amax(I_xm)/2, np.amax(I_xm), 0.1,0]
p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
p_a0m=p.copy()
err_a0m=np.diag(cov)**0.5
x_0m=p[-1]
x_0m_err=err_a0m[-1]
sig_xm_0=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.vlines(0,0,fit_cos(p[-1]/p[-2], *p))
# ax.errorbar((c_pos)*p[-2]-x_0p,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
# ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
# # ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r")
# ax.set_xlabel("$\\beta$")

"""
alpha=pi/4 path1
"""

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
# ps_i=108.8
# ps_f=ps_pos[-1]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

x_plt_ifg=np.linspace(ps_pos[0],ps_pos[-1],1000)
ifg=np.sum(matrix,axis=1)
P0=[np.amax(ifg), np.amax(ifg), 8,ps_pos[ifg==np.amax(ifg)][0]]
p,cov=fit(fit_cos,ps_pos,ifg, p0=P0)
ps_pos_rad=ps_pos*p[-2]-p[-1]
ps_0=ps_pos[abs(ps_pos-p[-1]/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]))][0]
ps_pi=ps_pos[abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2]))][0]
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.vlines(ps_0, 0, ifg[ps_pos==ps_0])
# ax.vlines(ps_pi, 0, ifg[ps_pos==ps_pi])
# ax.errorbar(ps_pos,ifg,yerr=np.sqrt(ifg),fmt="ko",capsize=5)
# ax.plot(x_plt_ifg,fit_cos(x_plt_ifg, *p), "b")
bp=np.zeros(len(ps_pos))
bp_err=np.zeros(len(ps_pos))
bm=np.zeros(len(ps_pos))
bm_err=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):

    """
    Sigma+
    """
    max_pos=ps_pos[i]
    I_xp=matrix[i,:]
    I_xp_err=matrix_err[i,:]
    P0=[np.amax(I_xp)/2, np.amax(I_xp), 0.1,0]
    p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
    x_pi4p=p[-1]
    errp=np.diag(cov)**0.5
    x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
    sig_xp_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
    bp[i]= x_0p-x_pi4p
    bp_err[i]=(errp[-1]**2+x_0p_err**2)**0.5
    
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # ax.set_title(str(ps_pos[i]))
    # # ax.vlines(x_pi4p-x_0p,0,fit_cos(p[-1]/p[-2], *p))
    # # ax.errorbar((c_pos)*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
    # # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
    # # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
    # ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r", label=str(p[-1]/p[-2]))
    # ax.legend()
    # ax.set_xlabel("$\\beta$")
    
    """
    Sigma-
    """
    I_xm=matrix[i,:]
    I_xm_err=matrix_err[i,:]
    P0=[np.amax(I_xm)/2, np.amax(I_xm), 0.1,0]
    p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
    x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
    errm=np.diag(cov)**0.5
    x_pi4m=p[-1]
    sig_xm_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
    bm[i]=x_0m-x_pi4m
    bm_err[i]=(errm[-1]**2+x_0m_err**2)**0.5
    
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # ax.vlines(x_pi4m-x_0m,0,fit_cos(p[-1]/p[-2], *p))
    # ax.errorbar((c_pos)*p[-2]-x_0m,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
    # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p), "b")
    # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
    # # ax.plot((x_plt_0)*p[-2]-x_0m,sig_xp_0, "r")
    # ax.set_xlabel("$\\beta$")
    

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.set_title("pi/4 path1")
ax.errorbar(ps_pos_rad-ps_pos_rad[0],bp/np.pi*4,yerr=bp_err/np.pi*4,fmt="ko",capsize=5)
# ax.errorbar(ps_pos_rad-ps_pos_rad[0]+np.pi,bm/np.pi*4,yerr=bm_err/np.pi*4,fmt="ro",capsize=5)
ax.set_xlabel("$\chi$")
#ax.set_xlim([3*np.pi,5*np.pi])


# """
# alpha=pi/4 path2
# """

# inf_file_name="path2pi4cb_g_10Apr1831"
# sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
# cleandata=sorted_fold_path+"/Cleantxt" 
# beta_fold_clean=cleandata+"/Beta"
# plots_fold=sorted_fold_path+"/Plots/"
# i=0
# for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
#     files=np.sort(files)
#     for name in files[:-1]:
#         if i==0:
#             tot_data=np.loadtxt(os.path.join(root, name))
#             c_pos=tot_data[:,0]
#             i=1
#         else:
#             data=np.loadtxt(os.path.join(root, name))
#             tot_data = np.vstack((tot_data, data))
# ps_pos=tot_data[::len(c_pos),-1]
# # ps_i=108.8
# # ps_f=ps_pos[-1]
# # ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
# matrix=np.zeros((len(ps_pos),len(c_pos)))
# matrix_err=np.zeros((len(ps_pos),len(c_pos)))
# w=np.zeros(len(ps_pos))
# err_b=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):
#     matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
#     matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

# x_plt_ifg=np.linspace(ps_pos[0],ps_pos[-1],1000)
# ifg=np.sum(matrix,axis=1)
# P0=[np.amax(ifg), np.amax(ifg), 8,ps_pos[ifg==np.amax(ifg)][0]]
# p,cov=fit(fit_cos,ps_pos,ifg, p0=P0)
# ps_pos_rad=ps_pos*p[-2]-p[-1]
# ps_0=ps_pos[abs(ps_pos-p[-1]/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]))][0]
# ps_pi=ps_pos[abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2]))][0]
# # fig = plt.figure(figsize=(5,5))
# # ax = fig.add_subplot(111)
# # ax.vlines(ps_0, 0, ifg[ps_pos==ps_0])
# # ax.vlines(ps_pi, 0, ifg[ps_pos==ps_pi])
# # ax.errorbar(ps_pos,ifg,yerr=np.sqrt(ifg),fmt="ko",capsize=5)
# # ax.plot(x_plt_ifg,fit_cos(x_plt_ifg, *p), "b")
# bp=np.zeros(len(ps_pos))
# bp_err=np.zeros(len(ps_pos))
# bm=np.zeros(len(ps_pos))
# bm_err=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):

#     """
#     Sigma+
#     """
#     max_pos=ps_pos[i]
#     I_xp=matrix[i,:]
#     I_xp_err=matrix_err[i,:]
#     P0=[np.amax(I_xp)/2, np.amax(I_xp), 0.1,0]
#     p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
#     x_pi4p=p[-1]
#     errp=np.diag(cov)**0.5
#     x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
#     sig_xp_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
#     bp[i]= x_0p-x_pi4p
#     bp_err[i]=(errp[-1]**2+x_0p_err**2)**0.5
    
#     # fig = plt.figure(figsize=(5,5))
#     # ax = fig.add_subplot(111)
#     # ax.vlines(x_pi4p-x_0p,0,fit_cos(p[-1]/p[-2], *p))
#     # ax.errorbar((c_pos)*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
#     # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
#     # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
#     # # ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r")
#     # ax.set_xlabel("$\\beta$")
    
#     """
#     Sigma-
#     """
#     I_xm=matrix[i,:]
#     I_xm_err=matrix_err[i,:]
#     P0=[np.amax(I_xm)/2, np.amax(I_xm), 0.1,0]
#     p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
#     x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
#     errm=np.diag(cov)**0.5
#     x_pi4m=p[-1]
#     sig_xm_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
#     bm[i]=x_0m-x_pi4m
#     bm_err[i]=(errm[-1]**2+x_0m_err**2)**0.5
    
#     # fig = plt.figure(figsize=(5,5))
#     # ax = fig.add_subplot(111)
#     # ax.vlines(x_pi4m-x_0m,0,fit_cos(p[-1]/p[-2], *p))
#     # ax.errorbar((c_pos)*p[-2]-x_0m,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
#     # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p), "b")
#     # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
#     # # ax.plot((x_plt_0)*p[-2]-x_0m,sig_xp_0, "r")
#     # ax.set_xlabel("$\\beta$")
    

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.set_title("pi/4 path2")
# ax.errorbar(ps_pos_rad-ps_pos_rad[0],bp/np.pi*4,yerr=bp_err/np.pi*4,fmt="ko",capsize=5)
# ax.errorbar(ps_pos_rad-ps_pos_rad[0]+np.pi,bm/np.pi*4,yerr=bm_err/np.pi*4,fmt="ro",capsize=5)
# ax.set_xlabel("$\chi$")
# #ax.set_xlim([3*np.pi,5*np.pi])
# """
# alpha=pi/8 path1
# """

# inf_file_name="path1pi8cb_g_09Apr1441"
# sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
# cleandata=sorted_fold_path+"/Cleantxt" 
# beta_fold_clean=cleandata+"/Beta"
# plots_fold=sorted_fold_path+"/Plots/"
# i=0
# for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
#     files=np.sort(files)
#     for name in files[:-1]:
#         if i==0:
#             tot_data=np.loadtxt(os.path.join(root, name))
#             c_pos=tot_data[:,0]
#             i=1
#         else:
#             data=np.loadtxt(os.path.join(root, name))
#             tot_data = np.vstack((tot_data, data))
# ps_pos=tot_data[::len(c_pos),-1]
# # ps_i=108.8
# # ps_f=ps_pos[-1]
# # ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
# matrix=np.zeros((len(ps_pos),len(c_pos)))
# matrix_err=np.zeros((len(ps_pos),len(c_pos)))
# w=np.zeros(len(ps_pos))
# err_b=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):
#     matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
#     matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

# x_plt_ifg=np.linspace(ps_pos[0],ps_pos[-1],1000)
# ifg=np.sum(matrix,axis=1)
# P0=[np.amax(ifg), np.amax(ifg), 8,ps_pos[ifg==np.amax(ifg)][0]]
# p,cov=fit(fit_cos,ps_pos,ifg, p0=P0)
# ps_pos_rad=ps_pos*p[-2]-p[-1]
# ps_0=ps_pos[abs(ps_pos-p[-1]/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]))][0]
# ps_pi=ps_pos[abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2]))][0]
# # fig = plt.figure(figsize=(5,5))
# # ax = fig.add_subplot(111)
# # ax.vlines(ps_0, 0, ifg[ps_pos==ps_0])
# # ax.vlines(ps_pi, 0, ifg[ps_pos==ps_pi])
# # ax.errorbar(ps_pos,ifg,yerr=np.sqrt(ifg),fmt="ko",capsize=5)
# # ax.plot(x_plt_ifg,fit_cos(x_plt_ifg, *p), "b")
# bp=np.zeros(len(ps_pos))
# bp_err=np.zeros(len(ps_pos))
# bm=np.zeros(len(ps_pos))
# bm_err=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):

#     """
#     Sigma+
#     """
#     max_pos=ps_pos[i]
#     I_xp=matrix[i,:]
#     I_xp_err=matrix_err[i,:]
#     P0=[np.amax(I_xp)/2, np.amax(I_xp), 0.1,0]
#     p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
#     x_pi4p=p[-1]
#     errp=np.diag(cov)**0.5
#     x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
#     sig_xp_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
#     bp[i]= x_0p-x_pi4p
#     bp_err[i]=(errp[-1]**2+x_0p_err**2)**0.5
    
#     # fig = plt.figure(figsize=(5,5))
#     # ax = fig.add_subplot(111)
#     # ax.vlines(x_pi4p-x_0p,0,fit_cos(p[-1]/p[-2], *p))
#     # ax.errorbar((c_pos)*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
#     # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
#     # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
#     # # ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r")
#     # ax.set_xlabel("$\\beta$")
    
#     """
#     Sigma-
#     """
#     I_xm=matrix[i,:]
#     I_xm_err=matrix_err[i,:]
#     P0=[np.amax(I_xm)/2, np.amax(I_xm), 0.1,0]
#     p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
#     x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
#     errm=np.diag(cov)**0.5
#     x_pi4m=p[-1]
#     sig_xm_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
#     bm[i]=x_0m-x_pi4m
#     bm_err[i]=(errm[-1]**2+x_0m_err**2)**0.5
    
#     # fig = plt.figure(figsize=(5,5))
#     # ax = fig.add_subplot(111)
#     # ax.vlines(x_pi4m-x_0m,0,fit_cos(p[-1]/p[-2], *p))
#     # ax.errorbar((c_pos)*p[-2]-x_0m,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
#     # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p), "b")
#     # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
#     # # ax.plot((x_plt_0)*p[-2]-x_0m,sig_xp_0, "r")
#     # ax.set_xlabel("$\\beta$")
    

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.set_title("pi/8 path1")
# ax.errorbar(ps_pos_rad-ps_pos_rad[0],bp/np.pi*8,yerr=bp_err/np.pi*8,fmt="ko",capsize=5)
# ax.errorbar(ps_pos_rad-ps_pos_rad[0]+np.pi,bm/np.pi*8,yerr=bm_err/np.pi*8,fmt="ro",capsize=5)
# ax.set_xlabel("$\chi$")
# #ax.set_xlim([3*np.pi,5*np.pi])
# """
# alpha=pi/8 path2
# """

# inf_file_name="path2pi8cb_g_12Apr1724"
# sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
# cleandata=sorted_fold_path+"/Cleantxt" 
# beta_fold_clean=cleandata+"/Beta"
# plots_fold=sorted_fold_path+"/Plots/"
# i=0
# for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
#     files=np.sort(files)
#     for name in files[:-1]:
#         if i==0:
#             tot_data=np.loadtxt(os.path.join(root, name))
#             c_pos=tot_data[:,0]
#             i=1
#         else:
#             data=np.loadtxt(os.path.join(root, name))
#             tot_data = np.vstack((tot_data, data))
# ps_pos=tot_data[::len(c_pos),-1]
# # ps_i=108.8
# # ps_f=ps_pos[-1]
# # ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
# matrix=np.zeros((len(ps_pos),len(c_pos)))
# matrix_err=np.zeros((len(ps_pos),len(c_pos)))
# w=np.zeros(len(ps_pos))
# err_b=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):
#     matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
#     matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5

# x_plt_ifg=np.linspace(ps_pos[0],ps_pos[-1],1000)
# ifg=np.sum(matrix,axis=1)
# P0=[np.amax(ifg), np.amax(ifg), 8,ps_pos[ifg==np.amax(ifg)][0]]
# p,cov=fit(fit_cos,ps_pos,ifg, p0=P0)
# ps_pos_rad=ps_pos*p[-2]-p[-1]
# ps_0=ps_pos[abs(ps_pos-p[-1]/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]))][0]
# ps_pi=ps_pos[abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2])==np.amin(abs(ps_pos-p[-1]/p[-2]-np.pi/p[-2]))][0]
# # fig = plt.figure(figsize=(5,5))
# # ax = fig.add_subplot(111)
# # ax.vlines(ps_0, 0, ifg[ps_pos==ps_0])
# # ax.vlines(ps_pi, 0, ifg[ps_pos==ps_pi])
# # ax.errorbar(ps_pos,ifg,yerr=np.sqrt(ifg),fmt="ko",capsize=5)
# # ax.plot(x_plt_ifg,fit_cos(x_plt_ifg, *p), "b")
# bp=np.zeros(len(ps_pos))
# bp_err=np.zeros(len(ps_pos))
# bm=np.zeros(len(ps_pos))
# bm_err=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):

#     """
#     Sigma+
#     """
#     max_pos=ps_pos[i]
#     I_xp=matrix[i,:]
#     I_xp_err=matrix_err[i,:]
#     P0=[np.amax(I_xp)/2, np.amax(I_xp), 0.1,0]
#     p,cov=fit(fit_cos,c_pos,I_xp, p0=P0)
#     x_pi4p=p[-1]
#     errp=np.diag(cov)**0.5
#     x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
#     sig_xp_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
#     bp[i]= x_0p-x_pi4p
#     bp_err[i]=(errp[-1]**2+x_0p_err**2)**0.5
    
#     # fig = plt.figure(figsize=(5,5))
#     # ax = fig.add_subplot(111)
#     # ax.vlines(x_pi4p-x_0p,0,fit_cos(p[-1]/p[-2], *p))
#     # ax.errorbar((c_pos)*p[-2]-x_0p,I_xp,yerr=np.sqrt(I_xp),fmt="ko",capsize=5)
#     # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p), "b")
#     # ax.plot((x_plt_0)*p[-2]-x_0p,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
#     # # ax.plot((x_plt_0)*p[-2]-x_0p,sig_xp_0, "r")
#     # ax.set_xlabel("$\\beta$")
    
#     """
#     Sigma-
#     """
#     I_xm=matrix[i,:]
#     I_xm_err=matrix_err[i,:]
#     P0=[np.amax(I_xm)/2, np.amax(I_xm), 0.1,0]
#     p,cov=fit(fit_cos,c_pos,I_xm, p0=P0)
#     x_plt_0= np.linspace(c_pos[0], c_pos[-1],100)
#     errm=np.diag(cov)**0.5
#     x_pi4m=p[-1]
#     sig_xm_pi4=(fit_cos(x_plt_0, *p)-fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))/(fit_cos(x_plt_0, *p)+fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi))
    
#     bm[i]=x_0m-x_pi4m
#     bm_err[i]=(errm[-1]**2+x_0m_err**2)**0.5
    
#     # fig = plt.figure(figsize=(5,5))
#     # ax = fig.add_subplot(111)
#     # ax.vlines(x_pi4m-x_0m,0,fit_cos(p[-1]/p[-2], *p))
#     # ax.errorbar((c_pos)*p[-2]-x_0m,I_xm,yerr=np.sqrt(I_xm),fmt="ko",capsize=5)
#     # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p), "b")
#     # ax.plot((x_plt_0)*p[-2]-x_0m,fit_cos(x_plt_0, *p[:-1],p[-1]+np.pi), "r")
#     # # ax.plot((x_plt_0)*p[-2]-x_0m,sig_xp_0, "r")
#     # ax.set_xlabel("$\\beta$")
    

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.set_title("pi/8 path2")
# ax.errorbar(ps_pos_rad-ps_pos_rad[0],bp/np.pi*8,yerr=bp_err/np.pi*8,fmt="ko",capsize=5)
# ax.errorbar(ps_pos_rad-ps_pos_rad[0]+np.pi,bm/np.pi*8,yerr=bm_err/np.pi*8,fmt="ro",capsize=5)
# ax.set_xlabel("$\chi$")
#ax.set_xlim([3*np.pi,5*np.pi])