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
plt.rcParams.update({'figure.max_open_warning': 0})
from PIL import Image as im
from scipy.optimize import curve_fit as fit

def fit_func(x,A,B,C,D):
    return A+B*np.cos(C*x-D)
rad=np.pi/180
inf_file_name="moveDC2_08Apr1342"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt" 
beta0_fold_clean=cleandata+"/Beta0"
i=0
for root, dirs, files in os.walk(beta0_fold_clean, topdown=False):
    max_c_pos=np.zeros((len(files)))
    w=np.zeros((len(files)))
    p_tot=np.zeros((len(files),4))
    cov_tot=np.zeros((len(files),4))
    for name in files:
        data=np.loadtxt(os.path.join(root, name))
        c_pos=data[:,0]
        counts=data[:,2]    
        P0=[0, np.amax(counts), 0.1,0]
        # print(P0)
        p,cov=fit(fit_func,c_pos,counts, p0=P0)
        p_tot[i]=p
        cov_tot[i]=np.diag(cov)**0.5
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        fig.suptitle(name)
        x_plt = np.linspace(c_pos[0], c_pos[-1],100)
        x_plt1 = np.linspace(0, 40,100)
        w[i]=p[2]
        max_c_pos[i]=p[-1]/p[-2]#x_plt1[fit_func(x_plt1, *p)==np.amax(fit_func(x_plt1, *p))]
        ax.errorbar(c_pos,counts,yerr=np.sqrt(counts),fmt="ko",capsize=5)
        ax.vlines(max_c_pos[i],0,fit_func(max_c_pos[i], *p),ls="dashed",label="Max coil pos=\n"+str("%.3f" % (max_c_pos[i]),))
        ax.plot(x_plt,fit_func(x_plt, *p), "b")
        ax.set_ylim([0, P0[1]+P0[1]/10])
        ax.legend(loc=1)
        i+=1
        
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(max_c_pos)
print(np.average(p_tot,axis=0))
print(np.average(cov_tot,axis=0))
print(2*np.pi/np.average(w))
plt.show()

counts*=0

for root, dirs, files in os.walk(beta0_fold_clean, topdown=False):
    for name in files:
        data=np.loadtxt(os.path.join(root, name))
        counts+=data[:,2]    
P0=[0, np.amax(counts), 0.1,0]
# print(P0)
p,cov=fit(fit_func,c_pos,counts, p0=P0)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
fig.suptitle("average")
x_plt = np.linspace(c_pos[0], c_pos[-1],100)
x_plt1 = np.linspace(0, 40,100)
#x_plt1[fit_func(x_plt1, *p)==np.amax(fit_func(x_plt1, *p))]
ax.errorbar(c_pos,counts,yerr=np.sqrt(counts),fmt="ko",capsize=5)
ax.vlines(p[-1]/p[-2],0,fit_func(p[-1]/p[-2], *p),ls="dashed",label="Max coil pos=\n"+str("%.3f" % (p[-1]/p[-2]),))
ax.plot(x_plt,fit_func(x_plt, *p), "b")
ax.set_ylim([0, P0[1]+P0[1]/10])
ax.legend(loc=1)
        
print(p)
print(np.diag(cov)**0.5)
print(2*np.pi/np.average(w))
plt.show()