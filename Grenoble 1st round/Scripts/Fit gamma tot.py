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
plt.rcParams.update({'figure.autolayout': True})
from PIL import Image as im
from scipy.optimize import curve_fit as fit

rad=np.pi/180
inf_file_names_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/Inf file names"
inf_file_names = np.genfromtxt(inf_file_names_path, dtype=str,usecols=(0),comments="#")
data_txt_path ="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Im weak value/"
# print(inf_file_names)
theta_fit=np.zeros((len(inf_file_names)))
k=0
for inf_file_name in inf_file_names:
    if "pi8" in inf_file_name:
        alpha=22.5
    if "pi4" in inf_file_name:
        alpha=45

    if "In2" in inf_file_name:
        a21=0.5
    elif "noIn" in inf_file_name:
        a21=1
    elif "In08" in inf_file_name:
        a21=1.375
    else:
        a21=2
        
    if "path1" in inf_file_name:
        path=1
        def fit_w1p(x,th,x0):
            return alpha*((1/(1+a21*np.tan(th*np.pi/2)**2*np.exp(-1j*(w_ps*(x-x0)))))).imag

        def exp_w1p(x,x0):
            return alpha*((1/(1+a21*np.exp(-1j*(w_ps*(x-x0)))))).imag
        
    if "path2" in inf_file_name:
        path=2
        def fit_w1p(x,th,x0):
            return alpha*(1-(1/(1+a21*np.tan(th*np.pi/2)**2*np.exp(-1j*(w_ps*(x-x0)))))).imag

        def exp_w1p(x,x0):
            return alpha*(1-(1/(1+a21*np.exp(-1j*(w_ps*(x-x0)))))).imag
   
    w_ps=8.002

    def fit_cos(x,A,B,C,D):
        return A+B*np.cos(C*(x-D))
    
    
    
    sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
    cleandata=sorted_fold_path+"/Cleantxt" 
    gamma_fold_clean=cleandata+"/Gamma"
    plots_fold="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Plots/"
    i=0
    for root, dirs, files in os.walk(gamma_fold_clean, topdown=False):
        files=np.sort(files)[:-1]
        for name in files:
            if i==0:
                tot_data=np.loadtxt(os.path.join(root, name))
                c_amp=tot_data[:,0]
                i=1
            else:
                data=np.loadtxt(os.path.join(root, name))
                tot_data = np.vstack((tot_data, data))
    ps_pos=tot_data[::len(c_amp),-2]
    matrix=np.zeros((len(ps_pos),len(c_amp)))
    max_c_amp=np.zeros(len(ps_pos))
    gamma=np.zeros(len(ps_pos))
    w=np.zeros(len(ps_pos))
    fit_res0=[503.8204133,  483.67092569,   2.79781313,  -0.5601]
    err_res0=[7.24106172, 10.30876345, 0.01913184, 0.01161413]
    err_g=np.zeros(len(ps_pos))
    max0=fit_res0[-1]
    for i in range(len(ps_pos)):
        matrix[i]=tot_data[:,2][tot_data[:,-2]==ps_pos[i]]
    for i in range(len(ps_pos)):
        P0=[0, np.amax(matrix[i]), 3,-0.6]
        p,cov=fit(fit_cos,c_amp,matrix[i], p0=P0)
        err=np.diag(cov)**0.5
        # print(p)
        x_plt = np.linspace(c_amp[0], c_amp[-1],100)
        x_plt1 = np.linspace(-2, 0,100)
        w[i]=p[2]
        max_c_amp[i]=p[3]#x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
        gamma[i]=(max_c_amp[i]-max0)*360/(2*np.pi/w[i])
        eps=max0-max_c_amp[i]
        err_eps=err[3]+err_res0[3] 
        err_g[i]=((w[i]*err_eps)**2 +(eps*err_res0[2])**2)**0.5/rad
        # fig = plt.figure(figsize=(5,5))
        # ax = fig.add_subplot(111)
        # fig.suptitle("ps_pos="+str(ps_pos[i]))
        # ax.errorbar(c_amp,matrix[i],yerr=np.sqrt(matrix[i]),fmt="ko",capsize=5)
        # ax.vlines(max_c_amp[i],0,fit_cos(max_c_amp[i], *p),ls="dashed",color="b",label="$\gamma$="+str("%.3f" % (max_c_amp[i]),))
        # ax.vlines(max0,0,fit_cos(max0, p[0],p[1],*fit_res0[2:]),ls="dashed",color="r",label="$\gamma_0$="+str("%.3f" % (max0),))
        # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
        # ax.set_ylim([0, P0[1]+P0[1]/10])
        # ax.plot(x_plt,fit_cos(x_plt, p[0],p[1],*fit_res0[2:]), "r")
        # ax.legend(loc=4)
        # plt.savefig(plots_fold+str(i)+".png",dpi=200)

    ps_data=np.sum(matrix,axis=1)
    P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 8,ps_pos[0]]
    p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0)
    # x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
    # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
    w_ps=p[-2]
    print(p)
    # print(w_ps)
    # fig = plt.figure(figsize=(5,5))
    # ax = fig.add_subplot(111)
    # ax.plot(ps_pos,max_c_pos)
    P01=[0.6,0.5]
    B0=([-1,0],[1,2*np.pi])
    p1,cov1=fit(fit_w1p,ps_pos,gamma, p0=P01, bounds=B0)
    # alpha=p[0]
    x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
    # fig = plt.figure(figsize=(5,5),dpi=150)
    # gs_t = GridSpec(4,1, figure=fig,hspace=0, bottom=0.1,top=0.98)
    # gs_b =GridSpec(4,1, figure=fig, wspace=0, top=0.5)
    # ax = [fig.add_subplot(gs_t[:-1]), 
    #       fig.add_subplot(gs_b[-1])]
    # ax[-1].tick_params(axis="x", labelbottom=False, bottom = False)
    # ax[-1].tick_params(axis="y", labelleft=False, left = False)
    # ax[0].set_title(inf_file_name)
    # ax[0].set_xlabel("$\chi$ ($\pi$)")
    # ax[0].set_ylim([-1,1])
    # x_plt_pi=(x_plt-x_plt[0]+p1[-1])*w_ps/np.pi
    # # ax[0].plot(ps_pos,gamma/alpha, "ko")
    # ax[0].errorbar((ps_pos-ps_pos[0]+p1[-1])*w_ps/np.pi,gamma/alpha, yerr=err_g/alpha,fmt="ko",capsize=5)
    # ax[0].plot(x_plt_pi,exp_w1p(x_plt, p1[-1])/alpha,"g", label="Exp Im{"+"$\omega_{+}$}")
    # ax[0].legend()
    # fit_param_names= ["$\\theta$=","$\omega x_0$="]
    # text = "Fit results:\t"
    # for i in range(len(p1)):
    #     text+= fit_param_names[i]+str("%.4f"%(p1[i],))+"\t"
    # text=text[:-2]
    # ax[-1].text(0.5,0.5,text,va="center", ha="center")
    # ax[0].plot(x_plt_pi,fit_w1p(x_plt, *p1)/alpha,"b", label="Fit Im{"+"$\omega_{+}$}")
    # plt.show()
    # plt.close(fig)
    theta_fit[k]=p1[0]
    k+=1
    # ax[0].plot(x_plt_pi,np.average(exp_w1p(x_plt[x_plt_pi<6], p1[-1])/alpha)+x_plt*0,"r", label="Avg")

    # print("max-min=",np.amax(gamma)-np.amin(gamma))
    # print("max=",np.amax(gamma),"\nmin=",np.amin(gamma))
    # print("period=",2*np.pi/np.average(w))
    # print(inf_file_name,"alpha=",alpha,"a21=",a21,"path=",path)
    # print("avg err=",np.average(err_b/alpha))
    # print("max err=",np.amax(err_b/alpha))
    # print("min err=",np.amin(err_b/alpha))
    # plt.savefig(plots_fold+inf_file_name+"_gamma.png",bbox_inches='tight')
    data_txt=np.array([ps_pos,gamma/alpha,err_g/alpha])
    with open(data_txt_path+inf_file_name+".txt", 'w') as f:
        np.savetxt(f, np.transpose(data_txt), header= "chi(pi) gamma err")


print(theta_fit)