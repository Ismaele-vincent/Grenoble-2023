#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 17 15:07:43 2023

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
plt.rcParams.update(plt.rcParamsDefault)
alpha=45
w_ps=8.002
a21=2 #1.375

def exp_w1p(x,x0):
    return alpha*(1-1/(1+a21*np.exp(-1j*(x-x0)))).real

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)


rad=np.pi/180
inf_file_name="path2pi4cb_g_10Apr1831"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
correct_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Corrected data/"+inf_file_name
beta_fold_clean=correct_fold_path+"/Beta"

i=0
for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
    files=np.sort(files)
    for name in files:
        if i==0:
            tot_data=np.loadtxt(os.path.join(root, name))
            coil=tot_data[:,0]
            i=1
        else:
            data=np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos=tot_data[::len(coil),-1]
# ps_i=109
# ps_f=ps_pos[-10]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
# print(tot_data)
matrix=np.zeros((len(ps_pos),len(coil)))
matrix_err=np.zeros((len(ps_pos),len(coil)))
b=np.zeros(len(ps_pos))
beta=np.zeros(len(ps_pos))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
fit_res0=[1.32214076e+03, 1.17144255e+03, 1.31186510e-01, 1.8086273e+00]
err_res0=[8.10516385e+00, 1.16683697e+01, 3.20540294e-04, 1.75421950e-02]
b0=fit_res0[-1]
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,1][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
for i in range(len(ps_pos)):
    P0=[5, np.amax(matrix[i]), 0.1,0]
    B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
    p,cov=fit(fit_cos,coil,matrix[i], p0=P0, bounds=B0)#, sigma=matrix_err[i])
    err=np.diag(cov)**0.5
    # print(p[3], err[3])
    x_plt = np.linspace(coil[0], coil[-1],100)
    w[i]=p[2]
    #x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
    b[i]=p[3]
    err_b[i]=(err[3]**2+err_res0[3]**2)**0.5/rad
    # err_b[i]=((w[i]*err_eps)**2 +(eps*err_res0[2])**2)**0.5/rad
    beta[i]=(b0-b[i])/rad
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # fig.suptitle("ps_pos="+str(ps_pos[i]))
    # ax.errorbar(coil,matrix[i],yerr=matrix_err[i],fmt="ko",capsize=5)
    # ax.vlines(b[i]/w[i],0,fit_cos(b[i]/w[i], *p),ls="dashed",color="b",label="$\\beta$="+str("%.3f" % (b[i]/w[i]),))
    # # ax.vlines(g0[i]/w0[i],0,fit_cos(g0[i]/w0[i], p[0],p[1],*p0[2:]),ls="dashed",color="r",label="$\\beta_0$="+str("%.3f" % (g0[i]/w0[i]),))
    # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # ax.set_ylim([0, P0[1]+P0[1]/10])
    # # ax.plot(x_plt,fit_cos(x_plt, p[0],p[1],*p0[2:]), "r")
    # ax.legend(loc=4)

ps_data=np.sum(matrix,axis=1)
P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 8,ps_pos[0]*8]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0,bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
# ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
print(p)
w_ps=p[-2]
ps_0=p[-1]
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.plot(ps_pos,max_coil)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.plot(ps_pos,b)


chi=(ps_pos*w_ps-ps_0)
x_plt = np.linspace(chi[0], chi[-1],100)
fig = plt.figure(figsize=(5,5))
gs_t = GridSpec(4,1, figure=fig,hspace=0, bottom=0.1,top=0.98)
gs_b =GridSpec(4,1, figure=fig, wspace=0, top=0.5)
ax = fig.add_subplot(111)
ax.set_title(inf_file_name)
ax.set_xlabel("$\chi$ ($\pi$)")
# ax.set_ylim([-1,1])
ax.errorbar(chi/np.pi, beta, yerr=err_b,fmt="ko",capsize=5)
ax.plot(x_plt/np.pi,exp_w1p(x_plt, 0),"g", label="Exp Re{"+"$\omega_{1+}$}")
ax.legend()
plt.show()

datatxt= np.array([chi,beta,err_b])

with open(correct_fold_path+"/"+inf_file_name[:8]+"_Beta_corrected.txt","w") as f:
    np.savetxt(f,np.transpose(datatxt), header="chi w+ err", fmt='%.7f %.7f %.7f')

# data=np.loadtxt("/home/aaa/Desktop/path1pi8cb_g_09Apr1441.txt")

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.errorbar(data[:,0],data[:,1],yerr=data[:,2], fmt="o")