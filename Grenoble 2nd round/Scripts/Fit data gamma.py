#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 12:29:03 2023

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

alpha=22.5
w_ps=8.002
a21=1/2
def fit_w1p(x,A,th,x0):
    return alpha*(A+(1/(1+a21*np.tan(th*np.pi/4)**2*np.exp(-1j*(w_ps*(x-x0)))))).imag

def exp_w1p(x,x0):
    return alpha*((1/(1+a21*np.exp(-1j*(w_ps*(x-x0)))))).imag

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)

rad=np.pi/180
inf_file_name="beta_alpi8off_gamma_chicoarse_GAMMA_02Jul1638"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 2nd round/exp_CRG-3033/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt"
gamma_fold_clean=cleandata+"/Gamma"
plots_fold=sorted_fold_path+"/Plots/"
i=0
for root, dirs, files in os.walk(gamma_fold_clean, topdown=False):
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
# print(tot_data)
matrix=np.zeros((len(ps_pos),len(coil)))
gamma=np.zeros(len(ps_pos))
w=np.zeros(len(ps_pos))
g=np.zeros(len(ps_pos))
err_g=np.zeros(len(ps_pos))
g0=np.zeros(len(ps_pos))
err_g0=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]/tot_data[:,1][tot_data[:,-1]==ps_pos[i]]
for i in range(len(ps_pos)):
    P0=[(np.amax(matrix[i])+np.amin(matrix[i]))/2, np.amax(matrix[i])-np.amin(matrix[i]), 0.1,0]
    B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
    # print(P0)
    p,cov=fit(fit_cos,coil,matrix[i], p0=P0, bounds=B0)
    err=np.diag(cov)**0.5
    # print(p[3], err[3])
    x_plt = np.linspace(coil[0], coil[-1],100)
    w[i]=p[2]
    #x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
    g[i]=p[3]
    err_g[i]=(err[3]**2+err_g0[i]**2)**0.5/rad
    # err_g[i]=((w[i]*err_eps)**2 +(eps*err_res0[2])**2)**0.5/rad
    gamma[i]=(g0[i]-g[i])/rad
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # fig.suptitle("ps_pos="+str(ps_pos[i]))
    # ax.errorbar(coil,matrix[i],yerr=np.sqrt(matrix[i])/10,fmt="ko",capsize=5)
    # ax.vlines(g[i]/w[i],0,fit_cos(g[i]/w[i], *p),ls="dashed",color="b",label="$\\gamma$="+str("%.3f" % (g[i]/w[i]),))
    # # ax.vlines(g0[i]/w0[i],0,fit_cos(g0[i]/w0[i], p[0],p[1],*p0[2:]),ls="dashed",color="r",label="$\\gamma_0$="+str("%.3f" % (g0[i]/w0[i]),))
    # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # ax.set_ylim([0, P0[1]+P0[1]/10])
    # # ax.plot(x_plt,fit_cos(x_plt, p[0],p[1],*p0[2:]), "r")
    # ax.legend(loc=4)

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.plot(ps_pos,g)

ps_data=np.sum(matrix,axis=1)
P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 3,-3]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0,bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
ax.vlines(p[-1]/p[-2],0,fit_cos(p[-1]/p[-2], *p),ls="dashed")
# print(p)
w_ps=p[-2]

fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
Z=matrix
x=coil
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.set_xlabel('Coil')
ax.set_ylabel('Phase Shifter')
ax.set_zlabel('z')
ax.view_init(45, 20)


P0=[0,0.01,1.5]
B0=([-1,0.01,0],[1,1.5,2*np.pi])
p1,cov1=fit(fit_w1p,ps_pos,-gamma + np.average(gamma), p0=P0, bounds=B0, sigma=err_g/alpha)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
fig = plt.figure(figsize=(5,5))
gs_t = GridSpec(4,1, figure=fig,hspace=0, bottom=0.1,top=0.98)
gs_b =GridSpec(4,1, figure=fig, wspace=0, top=0.5)
ax = fig.add_subplot(111)
ax.set_title(inf_file_name)
ax.set_xlabel("$\chi$ ($\pi$)")
# ax.set_ylim([-1,1])
# ax.plot(ps_pos,gamma/alpha, "ko")
ax.errorbar((ps_pos-ps_pos[0]-p1[-1])*w_ps/np.pi,-gamma + np.average(gamma), yerr=err_g,fmt="ko",capsize=5)
# ax.plot((x_plt-x_plt[0]-p1[-1])*w_ps/np.pi,fit_w1p(x_plt, *p1),"b", label="Fit Re{"+"$\omega_{1+}$}")
ax.plot((x_plt-x_plt[0]-p1[-1])*w_ps/np.pi,exp_w1p(x_plt, p1[-1]),"g", label="Exp Re{"+"$\omega_{1+}$}")
ax.legend()