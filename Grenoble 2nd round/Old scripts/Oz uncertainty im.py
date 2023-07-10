#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 16:53:52 2023

@author: aaa
"""
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from PIL import Image as im
from scipy.optimize import curve_fit as fit
from scipy.interpolate import UnivariateSpline as spline

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)

rad=np.pi/180
C=0.7914937937064894 
C_err=0.015566869131704356

fit_res_plus=np.loadtxt("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Scripts/p_plus_fit_param.txt")
data_txt_path ="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Im weak value/"
inf_file_name="path2pi4cb_g_10Apr1831"
p_p_fit=fit_res_plus[0]
p_p_fit_err=fit_res_plus[1]
p_p_fit_rel_err= p_p_fit_err/p_p_fit

data_p= np.loadtxt("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Scripts/p_plus.txt")
ps_pos = data_p[:,0]
p_p=data_p[:,1]
p_p_err=data_p[:,2]
x0=p_p_fit[-1]+2*np.pi
chi=(ps_pos*p_p_fit[-2]-x0)/np.pi

p_m=0.5+C/2*np.cos(chi*np.pi+np.pi)
p_m_err=np.array([*p_p_err[5:],*p_p_err[:5]])
# plt.plot(fit_cos(ps_pos,*p_p_fit))
datatxt=np.array([chi,p_p,p_p_err,p_m,p_m_err])
with open("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Plots data/P_p_p_m.txt","w") as f:
    np.savetxt(f,np.transpose(datatxt), header="chi p_p p_p_err p_m p_m_err")
x_plt=np.linspace(chi[0], chi[-1],1000)
x_plt_pi=x_plt*np.pi
fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(111)
ax.errorbar(chi,p_p, p_p_err, fmt="ko", capsize=5, label="Scaled data $p_+$")
# ax.errorbar(chi,p_m, p_m_err, fmt="go", capsize=5, label="Scaled data $p_-$\n(extracted from fit)")
ax.plot(x_plt, 0.5+C/2*np.cos(x_plt*np.pi), "r", label="$p_+$")
# ax.plot(x_plt, 0.5+C/2*np.cos(x_plt*np.pi+np.pi), "b", label="$p_-$")
ax.set_xlim([-2,2])
ax.set_ylim([-0.1,1])
ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
ax.set_xlabel("$\chi$ ($\pi$)")
ax.legend(loc=4, ncol=2)
# plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Pppm DB.pdf", format="pdf", bbox_inches="tight")
plt.show()

Im_weak_data=np.loadtxt(data_txt_path+inf_file_name+".txt")
ps_pos1=Im_weak_data[:,0]
Im_weak=Im_weak_data[:,1]
Im_weak_err=Im_weak_data[:,2]

p_p_w=np.zeros((len(Im_weak)))
p_p_w_err=np.zeros((len(Im_weak)))

for i in range(len(Im_weak)):
    p_p_w[i]=p_p[abs(ps_pos-ps_pos[i])==np.amin(abs(ps_pos-ps_pos[i]))]
    p_p_w_err[i]=p_p_err[abs(ps_pos-ps_pos[i])==np.amin(abs(ps_pos-ps_pos[i]))]
dpdx=abs(C/2*np.sin((ps_pos1*p_p_fit[-2]-x0)))

fig = plt.figure(figsize=(5,6))
ax = fig.add_subplot(111)
# ax.errorbar(ps_pos1,Im_weak,Im_weak_err, fmt="ko", capsize=5, label="Im weak value")
# ax.errorbar(ps_pos1,p_p_w,p_p_w_err, fmt="ro", capsize=5, label="P+")
ax.plot(ps_pos1,2*p_p_w*abs(Im_weak), "go")
ax.plot(ps_pos1,dpdx, "g")
ax.set_xlabel("$\chi$ ($\pi$)")
ax.legend(loc=4, ncol=2)
# plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Pppm DB.pdf", format="pdf", bbox_inches="tight")
plt.show()












# DB=(1-(p_p-p_m)**2)**0.5
# DB_err = (p_m-p_p)*(p_m_err**2+p_p_err**2)**0.5/DB

# # fig = plt.figure(figsize=(5,4))
# # ax = fig.add_subplot(111)
# # ax.errorbar(chi,DB, DB_err, fmt="ko", capsize=4, label="data")
# # ax.plot(x_plt, (1-(4/5*np.cos(x_plt*np.pi))**2)**0.5,"b", label="theory curve")
# # ax.set_xlim([0,4])
# # ax.set_xlabel("$\chi$")
# # ax.set_xticks([0,0.5,1,1.5,2,2.5,3,3.5,4])
# # ax.set_xticklabels(["0","$\pi/2$","$\pi$","$3/2\pi$","$2\pi$","$5/2\pi$","$3\pi$","$7/2\pi$","$4\pi$"])
# # ax.set_yticks([0.6,0.7,0.8,0.9,1.0])
# # ax.set_ylim([0.55,1.05])
# # ax.legend(ncol=2)
# # datatxt=np.array([chi,DB,DB_err])
# # with open("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Plots data/DB.txt","w") as f:
# #     np.savetxt(f,np.transpose(datatxt), header="chi DB DB_err")
# # plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/DB.pdf", format="pdf", bbox_inches="tight")

# dpdx=abs(C/2*np.sin(chi*np.pi))#abs(C/2*np.sin(x_plt_pi))# abs(p_p_fit[1]*p_p_fit[2]*np.sin(p_p_fit[2]*ps_pos+p_p_fit[3]))#abs(C/2*np.sin(x_plt_pi))#
# dpdx_max=abs((C_err)/2*np.sin(chi*np.pi))
# # dpdx_max=abs((C+C_err)/2*np.sin(chi*np.pi))
# # dpdx_min=abs((C-C_err)/2*np.sin(chi))
# b=C/2
# c=1
# d=np.pi
# b_err=b*p_p_fit_rel_err[1]
# c_err=c*p_p_fit_rel_err[2]
# d_err=d*p_p_fit_rel_err[3]
# dpdx_err=(b**2*c**2*d_err**2*np.cos(d-x_plt*np.pi)**2+c**2*b_err**2*np.sin(d-x_plt_pi)**2+b**2*c_err**2*(-x_plt_pi*c*np.cos(d-x_plt_pi)+np.sin(d-x_plt_pi))**2)**0.5

# fig = plt.figure(figsize=(5,4))
# ax = fig.add_subplot(111)
# # ax.fill_between(
# #         x_plt,
# #         abs(dpdx_min), 
# #         abs(dpdx_max), 
# #         color= "b",
# #         alpha= 0.3)
# ax.errorbar(chi,dpdx,dpdx_max, fmt="ko", capsize=4, label="data") #dpdx_max=err
# ax.plot(x_plt, abs(2/5*np.sin(x_plt_pi)),"b",label="theory curve")
# ax.set_xlim([0,4])
# ax.set_xlabel("$\chi$")
# ax.set_xticks([0,0.5,1,1.5,2,2.5,3,3.5,4])
# ax.set_xticklabels(["0","$\pi/2$","$\pi$","$3/2\pi$","$2\pi$","$5/2\pi$","$3\pi$","$7/2\pi$","$4\pi$"])
# # ax.set_yticks([0.6,0.7,0.8,0.9,1.0])
# ax.set_ylim([0,0.5])
# ax.legend(ncol=2,loc=1)
# datatxt=np.array([chi,dpdx,dpdx_max])
# with open("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Plots data/dpdx.txt","w") as f:
#     np.savetxt(f,np.transpose(datatxt), header="chi dpdx dpdx_err")
# # plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Ozawa_limit.pdf", format="pdf", bbox_inches="tight")


# # f = spline(chi, p_p, k=4, s=0)
# # f_min = spline(chi, p_p-p_p_err, k=4, s=0)
# # f_max = spline(chi, p_p+p_p_err, k=4, s=0)
# # fig = plt.figure(figsize=(5,6))
# # ax = fig.add_subplot(111)
# # ax.errorbar(chi,p_p, p_p_err, fmt="ko", capsize=4)
# # ax.plot(x_plt, f(x_plt))
# # ax.plot(x_plt, f_min(x_plt))
# # ax.plot(x_plt, f_max(x_plt))
# # ax.set_xlim([0,2])
# # ax.set_xlabel("$\chi$")
# # ax.set_xticks([0,0.5,1,1.5,2])
# # ax.set_xticklabels(["0","$\pi/2$","$\pi$","$3/2\pi$","$2\pi$"])
# # # ax.set_ylim([0.7,1])

# # f = spline(chi, p_p, k=4, s=0)
# # f_min = spline(chi, p_p-p_p_err, k=4, s=0)
# # f_max = spline(chi, p_p+p_p_err, k=4, s=0)
# # fig = plt.figure(figsize=(5,6))
# # ax = fig.add_subplot(111)
# # # ax.fill_between(
# # #         x_plt,
# # #         abs(f_max.derivative(1)(x_plt)), 
# # #         abs(f_min.derivative(1)(x_plt)), 
# # #         color= "b",
# # #         alpha= 0.2)
# # ax.plot(x_plt, abs(C/2*np.sin(x_plt*np.pi)))
# # # ax.plot(x_plt, abs(f.derivative()(x_plt)))
# # # ax.plot(x_plt, abs(f_min.derivative()(x_plt)))
# # # ax.plot(x_plt, abs(f_max.derivative()(x_plt)))
# # ax.set_xlim([0,2])
# # ax.set_xlabel("$\chi$")
# # ax.set_xticks([0,0.5,1,1.5,2])
# # ax.set_xticklabels(["0","$\pi/2$","$\pi$","$3/2\pi$","$2\pi$"])
