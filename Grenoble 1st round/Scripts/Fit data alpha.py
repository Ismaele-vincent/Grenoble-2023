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
    return A+B*np.cos(C*(x-D))
rad=np.pi/180
inf_file_name="S2Z_10Apr1536"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt" 
beta_fold_clean=cleandata+"/Alpha0"
i_neg=0
i_pos=0
for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
    # max_c_pos=np.zeros((len(files)))
    # w=np.zeros((len(files)))
    for name in files:
        if "neg" in name:
            if i_neg==0:
                data_neg=np.loadtxt(os.path.join(root, name))            
                c_pos_neg = -data_neg[:,0].copy()
                i_neg=1
            else:
                data_neg+=np.loadtxt(os.path.join(root, name))
            # fig = plt.figure(figsize=(5,5))
            # ax = fig.add_subplot(111)
            # fig.suptitle(name)
            # ax.errorbar(c_pos_neg,data_neg[:,2] ,yerr=np.sqrt(data_neg[:,2]),fmt="ko",capsize=5)
        else:
            if i_pos==0:
                data_pos=np.loadtxt(os.path.join(root, name))
                c_pos_pos=data_pos[:,0].copy()
                i_pos=1
            else:
                data_pos+=np.loadtxt(os.path.join(root, name))
            # fig = plt.figure(figsize=(5,5))
            # ax = fig.add_subplot(111)
            # fig.suptitle(name)
            # ax.errorbar(c_pos_pos,data_pos[:,2] ,yerr=np.sqrt(data_pos[:,2]),fmt="ko",capsize=5)
data_neg=data_neg[np.argsort(c_pos_neg)]
c_pos_neg=c_pos_neg[np.argsort(c_pos_neg)]
counts_neg=data_neg[:,2]   
counts_pos=data_pos[:,2]     
counts1=np.append(counts_neg,counts_pos)
c_pos1=np.append(c_pos_neg,c_pos_pos)
c_pos=c_pos1[abs(c_pos1)<6]
counts=counts1[abs(c_pos1)<6]
P0=[0, np.amax(counts), 1,0]
# print(P0)Pos0
p,cov=fit(fit_func,c_pos,counts, p0=P0)
err=np.diag(cov)**0.5
x_plt = np.linspace(c_pos[0], c_pos[-1],100)
x_plt1 = np.linspace(0, 40,100)
w=p[2]
max_c_pos=p[3]#x_plt1[fit_func(x_plt1, *p)==np.amax(fit_func(x_plt1, *p))]
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
fig.suptitle(name)
ax.errorbar(c_pos,counts,yerr=np.sqrt(counts),fmt="ko",capsize=5)
ax.vlines(max_c_pos,0,fit_func(max_c_pos, *p),ls="dashed",label="Max coil pos=\n"+str("%.3f" % (max_c_pos),))
ax.vlines(np.pi/(4*w),0,fit_func(np.pi/(4*w), *p),ls="dashed",label="$\pi/4$ pos=\n"+str("%.3f" % (np.pi/(4*w)),))
ax.plot(x_plt,fit_func(x_plt, *p), "b")
ax.set_ylim([0, P0[1]+P0[1]/10])
ax.legend(loc=1)
# print(max_c_pos/(2*np.pi/w)*360)
# print(2*np.pi/w)
# print(p[-1])
print(p)
print("pi/4=",2*np.pi/w/8,"err=", 2*np.pi/w**2/8*err[2])
print("pi/8=",2*np.pi/w/16,"err=", 2*np.pi/w**2/16*err[2])
print("Flip eff. = ", (p[0]+p[1])/(p[0]-p[1]))
plt.show()