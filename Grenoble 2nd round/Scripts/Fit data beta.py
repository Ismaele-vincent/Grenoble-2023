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

alpha=-22.5
w_ps=3
a21=1/2
def fit_w1p(x,A,th,x0):
    return alpha*(A+(1/(1+a21*np.tan(th*np.pi/4)**2*np.exp(-1j*(w_ps*(x-x0)))))).real

def exp_w1p(x,x0):
    return alpha*((1/(1+a21*np.exp(-1j*(w_ps*(x-x0)))))).real

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)

rad=np.pi/180
inf_file_name="beta_alpi8off_gamma_chicoarse_BETA_01Jul1712"
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 2nd round/exp_CRG-3033/Sorted data/"+inf_file_name
cleandata=sorted_fold_path+"/Cleantxt"
beta_fold_clean_OFF=cleandata+"/Beta/AlphaOFF"
beta_fold_clean_ON=cleandata+"/Beta/AlphaON"
plots_fold=sorted_fold_path+"/Plots/"
i=0
for root, dirs, files in os.walk(beta_fold_clean_OFF, topdown=False):
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
print(coil)
matrix=np.zeros((len(ps_pos),len(coil)))
matrix_err=np.zeros((len(ps_pos),len(coil)))
w0=np.zeros(len(ps_pos))
err_b0=np.zeros(len(ps_pos))
b0=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]/tot_data[:,1][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**-0.5*matrix[i]
for i in range(len(ps_pos)):
    P0=[(np.amax(matrix[i])+np.amin(matrix[i]))/2, np.amax(matrix[i])-np.amin(matrix[i]), 48,0]
    B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
    # print(P0)
    p0,cov0=fit(fit_cos,coil,matrix[i], p0=P0, sigma=matrix_err[i])
    err=np.diag(cov0)**0.5
    # print(p[3], err[3])
    x_plt = np.linspace(coil[0], coil[-1],100)
    w0[i]=p0[2]
    b0[i]=p0[3]
    err_b0[i]=err[3]
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # fig.suptitle("ps_pos="+str(ps_pos[i]))
    # ax.errorbar(coil,matrix[i],yerr=matrix_err[i],fmt="ko",capsize=5)
    # ax.vlines(b0[i]/w0[i],0,fit_cos(b0[i]/w0[i], *p0),ls="dashed",color="b",label="$\\beta$="+str("%.3f" % (b0[i]/w0[i]),))
    # ax.plot(x_plt,fit_cos(x_plt, *p0), "b")
    # ax.set_ylim([0, P0[1]+P0[1]/10])
    # ax.legend(loc=4)

ps_data=np.sum(matrix,axis=1)
# ps_data=matrix[:,15]
P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 3,-3]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0,bounds=B0)
print(p[1]/p[0])
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(ps_pos,ps_data,yerr=np.sum(matrix_err,axis=1),fmt="ko",capsize=5)  
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
ax.vlines(p[-1]/p[-2],0,fit_cos(p[-1]/p[-2], *p),ls="dashed")
# print(p)
w_ps=p[-2]
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(ps_pos,b0)
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


i=0
for root, dirs, files in os.walk(beta_fold_clean_ON, topdown=False):
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
matrix_err=np.zeros((len(ps_pos),len(coil)))
b=np.zeros(len(ps_pos))
beta=np.zeros(len(ps_pos))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]/tot_data[:,1][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**-0.5*matrix[i]
for i in range(len(ps_pos)):
    P0=[(np.amax(matrix[i])+np.amin(matrix[i]))/2, np.amax(matrix[i])-np.amin(matrix[i]), 48,0]
    B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
    # print(P0)
    p,cov=fit(fit_cos,coil,matrix[i], p0=P0, sigma=matrix_err[i])
    err=np.diag(cov)**0.5
    # print(p[3], err[3])
    x_plt = np.linspace(coil[0], coil[-1],100)
    w[i]=p[2]
    #x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
    b[i]=p[3]
    err_b[i]=(err[3]**2+err_b0[i]**2)**0.5/rad
    # err_b[i]=((w[i]*err_eps)**2 +(eps*err_res0[2])**2)**0.5/rad
    beta[i]=(b0[i]-b[i])/rad
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # fig.suptitle("ps_pos="+str(ps_pos[i]))
    # ax.errorbar(coil,matrix[i],yerr=matrix_err[i],fmt="ko",capsize=5)
    # ax.vlines(b[i]/w[i],0,fit_cos(b[i]/w[i], *p),ls="dashed",color="b",label="$\\beta$="+str("%.3f" % (b[i]/w[i]),))
    # ax.vlines(b0[i]/w0[i],0,fit_cos(b0[i]/w0[i], p[0],p[1],*p0[2:]),ls="dashed",color="r",label="$\\beta_0$="+str("%.3f" % (b0[i]/w0[i]),))
    # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # # ax.set_ylim([0, P0[1]+P0[1]/10])
    # ax.plot(x_plt,fit_cos(x_plt, p[0],p[1],*p0[2:]), "r")
    # ax.legend(loc=4)

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.plot(ps_pos,b)

P0=[0,1,0.2]
B0=([-5,0.5,0],[1,1.5,2*np.pi])
p1,cov1=fit(fit_w1p,ps_pos, beta, p0=P0, bounds=B0, sigma=err_b/alpha)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot()

ax.set_title(inf_file_name)
ax.set_xlabel("$\chi$ ($\pi$)")
# ax.set_ylim([-1,1])
# ax.plot(ps_pos,beta/alpha, "ko")
ax.errorbar((ps_pos-ps_pos[0]-p1[-1])*w_ps/np.pi, beta, yerr=err_b,fmt="ko",capsize=5)
# ax.plot((x_plt-x_plt[0]-p1[-1])*w_ps/np.pi,fit_w1p(x_plt, *p1),"b", label="Fit Re{"+"$\omega_{1+}$}")
ax.plot((x_plt-x_plt[0]-p1[-1])*w_ps/np.pi,exp_w1p(x_plt, p1[-1]),"g", label="Exp Re{"+"$\omega_{1+}$}")