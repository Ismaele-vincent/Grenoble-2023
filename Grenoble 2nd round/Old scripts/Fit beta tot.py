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
from matplotlib.gridspec import GridSpec
plt.rcParams.update({'figure.max_open_warning': 0})
from PIL import Image as im
from scipy.optimize import curve_fit as fit

rad=np.pi/180
inf_file_names_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/Inf file names"
inf_file_names = np.genfromtxt(inf_file_names_path, dtype=str,usecols=(0),comments="#")

# print(inf_file_names)
theta_fit=np.zeros((len(inf_file_names)))
k=0
for inf_file_name in inf_file_names:#["path1pi4_noIn_cb_g_14Apr1720"]:#["path2pi4_In2_cb_g_11Apr1932"]:#
    max0=13.7265133
    if "pi8" in inf_file_name:
        alpha=22.5
    if "pi4" in inf_file_name:
        alpha=45

    if "In2" in inf_file_name:
        a21=0.5
    elif "noIn" in inf_file_name:
        a21=1
        max0=14.0003456
    elif "In08" in inf_file_name:
        a21=1.375
        max0=14.15
    else:
        a21=2
        
    if "path1" in inf_file_name:
        path=1
        def fit_w1p(x,A,th,x0):
            return alpha*(A+(1/(1+a21*np.tan(th*np.pi/2)**2*np.exp(-1j*(w_ps*x-x0))))).real

        def exp_w1p(x,x0):
            return alpha*((1/(1+a21*np.exp(-1j*(w_ps*x-x0))))).real
        
    if "path2" in inf_file_name:
        path=2
        def fit_w1p(x,A,th,x0):
            return alpha*(A+1-(1/(1+a21*np.tan(th*np.pi/2)**2*np.exp(-1j*(w_ps*x-x0))))).real

        def exp_w1p(x,x0):
            return alpha*(1-(1/(1+a21*np.exp(-1j*(w_ps*x-x0))))).real
   
    w_ps=8.002

    def fit_cos(x,A,B,C,D):
        return A+B*np.cos(C*(x-D))
    
    
    
    sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
    cleandata=sorted_fold_path+"/Cleantxt" 
    beta_fold_clean=cleandata+"/Beta"
    plots_fold="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Plots/"
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
    # print(tot_data)
    matrix=np.zeros((len(ps_pos),len(c_pos)))
    max_c_pos=np.zeros(len(ps_pos))
    beta=np.zeros(len(ps_pos))
    w=np.zeros(len(ps_pos))
    err_b=np.zeros(len(ps_pos))
    
    fit_res0=[2.31111864e+03, 2.07213380e+03, 1.30490949e-01, max0]#[2.20354126e+02,  1.95310536e+02,  1.31182401e-01, 1.37265133e+01]#
    err_res0=[2.77104601e+01, 3.96757459e+01, 6.26339602e-04, 2.07305597e-01]
    for i in range(len(ps_pos)):
        matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
    for i in range(len(ps_pos)):
        P0=[5, np.amax(matrix[i]), 0.1,0]
        # print(P0)
        p,cov=fit(fit_cos,c_pos,matrix[i], p0=P0, sigma=np.sqrt(matrix[i]))
        err=np.diag(cov)**0.5
        # print(p[3], err[3])
        x_plt = np.linspace(c_pos[0], c_pos[-1],100)
        w[i]=p[2]
        max_c_pos[i]=p[3]#x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
        eps=max0-max_c_pos[i]
        err_eps=err[3]+err_res0[3] 
        err_b[i]=(err[3]**2+err_res0[3]**2)**0.5*w[i]/rad
        # err_b[i]=((w[i]*err_eps)**2 +(eps*err_res0[2])**2)**0.5/rad
        beta[i]=(max0-max_c_pos[i])*w[i]/rad
        # fig = plt.figure(figsize=(5,5))
        # ax = fig.add_subplot(111)
        # fig.suptitle("ps_pos="+str(ps_pos[i]))
        # ax.errorbar(c_pos,matrix[i],yerr=np.sqrt(matrix[i]),fmt="ko",capsize=5)
        # ax.vlines(max_c_pos[i],0,fit_cos(max_c_pos[i], *p),ls="dashed",color="b",label="$\\beta$="+str("%.3f" % (max_c_pos[i]),))
        # ax.vlines(max0,0,fit_cos(max0, p[0],p[1],*fit_res0[2:]),ls="dashed",color="r",label="$\\beta_0$="+str("%.3f" % (max0),))
        # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
        # ax.set_ylim([0, P0[1]+P0[1]/10])
        # ax.plot(x_plt,fit_cos(x_plt, p[0],p[1],*fit_res0[2:]), "r")
        # ax.legend(loc=4)
        # plt.savefig(plots_fold+str(i)+".png",dpi=200)
        
    ps_data=np.sum(matrix,axis=1)/len(matrix[:,0])
    P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), w_ps,ps_pos[0]]
    p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0)
    x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
    # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
    w_ps=p[-2]
    ps_pos=ps_pos
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # ax.plot(ps_pos,max_c_pos)
    P01=[0,0.5,0.2]
    B0=([-1,0,0],[1,1,2*np.pi])
    p1,cov1=fit(fit_w1p,ps_pos,beta, p0=P01, bounds=B0, sigma=err_b/alpha)
    x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
    fig = plt.figure(figsize=(5,5),dpi=150)
    gs_t = GridSpec(4,1, figure=fig,hspace=0, bottom=0.1,top=0.98)
    gs_b =GridSpec(4,1, figure=fig, wspace=0, top=0.5)
    ax = [fig.add_subplot(gs_t[:-1]), 
          fig.add_subplot(gs_b[-1])]
    ax[-1].tick_params(axis="x", labelbottom=False, bottom = False)
    ax[-1].tick_params(axis="y", labelleft=False, left = False)
    ax[0].set_title(inf_file_name)
    ax[0].set_xlabel("$\chi$ ($\pi$)")
    # ax[0].set_ylim([-1,1])
    # ax[0].plot(ps_pos,beta/alpha, "ko")
    ps_pos_pi=(ps_pos*w_ps-p1[-1])/np.pi
    x_plt_pi=(x_plt*w_ps-p1[-1])/np.pi
    ax[0].errorbar(ps_pos_pi-ps_pos_pi[0],beta/alpha, yerr=err_b/alpha,fmt="ko",capsize=5)
    # ax[0].plot(x_plt_pi-x_plt_pi[0],fit_w1p(x_plt, *p1)/alpha,"b", label="Fit Re{"+"$\omega_{1+}$}")
    ax[0].plot(x_plt_pi-x_plt_pi[0],exp_w1p(x_plt, p1[-1])/alpha,"g", label="Exp Re{"+"$\omega_{1+}$}")
    ax[0].legend()
    theta_fit[k]=p1[1]
    k+=1
    # print("max-min=",np.amax(beta)-np.amin(beta))
    # print("max=",np.amax(beta),"\nmin=",np.amin(beta))
    # print("period=",2*np.pi/np.average(w))
    fit_param_names= ["$off$=","$\\theta$=","$\omega x_0$="]
    text = "Fit results:\t"
    for i in range(len(p1)):
        text+= fit_param_names[i]+str("%.4f"%(p1[i],))+"\t"
    text=text[:-2]
    ax[-1].text(0.5,0.5,text,va="center", ha="center")
    print(inf_file_name,"alpha=",alpha,"a21=",a21,"path=",path)
    # print("avg err=",np.average(err_b/alpha))
    # print("max err=",np.amax(err_b/alpha))
    # print("min err=",np.amin(err_b/alpha))
    # plt.savefig(plots_fold+inf_file_name+"_beta.png",bbox_inches='tight')
    plt.show()
    plt.close(fig)
print(theta_fit)