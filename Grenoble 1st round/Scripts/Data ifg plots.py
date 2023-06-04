#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:11:21 2023

@author: aaa
"""

import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.max_open_warning': 0})
from PIL import Image as im
from scipy.optimize import curve_fit as fit
from matplotlib.gridspec import GridSpec

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)
rad=np.pi/180

"""
Contrast no Indium
"""

# inf_file_name="path1pi4cb_g_13Apr1502"
# sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
# cleandata=sorted_fold_path+"/Cleantxt" 
# beta_fold_clean=cleandata+"/Beta"
# plots_fold=sorted_fold_path+"/Plots/"
# i=0
# for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
#     files=np.sort(files)
#     for name in files:
#         if i==0:
#             tot_data=np.loadtxt(os.path.join(root, name))
#             c_pos=tot_data[:,0]
#             i=1
#         else:
#             data=np.loadtxt(os.path.join(root, name))
#             tot_data = np.vstack((tot_data, data))
# ps_pos=tot_data[::len(c_pos),-1]
# ps_pos1=ps_pos.copy()
# # ps_i=108.8
# # ps_f=ps_pos[-1]
# # ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
# matrix=np.zeros((len(ps_pos),len(c_pos)))
# matrix1=np.zeros((len(ps_pos1),len(c_pos)))
# w=np.zeros(len(ps_pos))
# err_b=np.zeros(len(ps_pos))
# for i in range(len(ps_pos)):
#     matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]
# for i in range(len(ps_pos1)):
#     matrix1[i]=tot_data[:,2][tot_data[:,-1]==ps_pos1[i]]

# max_pos=c_pos[np.where(matrix==np.amax(matrix))[1][0]]
# counts=matrix[:,np.where(matrix==np.amax(matrix))[1][0]]
# P0=[1969.26882216, 1360.72510923, 8.02284756, 10]#[(np.amax(counts)+np.amin(counts))/2, np.amax(counts)-np.amin(counts), 0.1,ps_pos[0]]
# p_no_in,cov=fit(fit_cos,ps_pos,counts, p0=P0)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# fig.suptitle("coil pos="+str(max_pos)+" mm")
# x_plt_0 = np.linspace(ps_pos[0], ps_pos[-1],100)
# ax.errorbar((ps_pos),counts,yerr=np.sqrt(counts),fmt="ko",capsize=5)
# ax.plot((x_plt_0),fit_cos(x_plt_0, *p_no_in), "b")
# ax.set_xlabel("$\\beta$ (mm)")
# plt.show()

# err_no_in=np.diag(cov)**0.5
# print(p_no_in, err_no_in)
# C_no_in_err=(err_no_in[1]**2/p_no_in[0]**2+err_no_in[0]**2*p_no_in[1]**2/p_no_in[0]**4)**0.5
# C_no_in= p_no_in[1]/p_no_in[0]
# fig = plt.figure(figsize=(5,6))
# # fig.suptitle(inf_file_name+".inf")
# ax = plt.axes(projection='3d')
# Z=matrix1
# x=c_pos
# y=ps_pos1
# X, Y = np.meshgrid(x, y)
# ax.contour3D(X, Y, Z, 30, cmap='binary')
# ax.plot3D(15+ps_pos*0,ps_pos, counts, 'r', lw=5,label="Interferogram (coil pos=15)")
# ax.set_xlabel('Coil ($\propto \\beta$)')
# ax.set_ylabel('PS ($\propto \chi$)')
# ax.set_zlabel('counts')
# ax.view_init(50, -20)
# ax.legend(loc=2)

# # plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/IFGs_alpha0.pdf", format="pdf", bbox_inches="tight")
# plt.show()
# print("Contrast no indium=", C_no_in)

"""
Contrast Indium
"""

# inf_file_name="ifg_indium_04Mar1627ff"
# sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
# cleandata=sorted_fold_path+"/Cleantxt" 
# ifg_fold_clean=cleandata+"/IFG"
# i=0
# for root, dirs, files in os.walk(ifg_fold_clean, topdown=False):
#     max_ps_pos=np.zeros((len(files)))
#     w=np.zeros((len(files)))
#     p_tot=np.zeros((len(files),4))
#     cov_tot=np.zeros((len(files),4))
#     fig = plt.figure(figsize=(11,5.5))
#     gs = GridSpec(2,len(files)//2, figure=fig)
#     ax = []
#     for j in range(len(files)//2):
#         ax=np.append(ax, fig.add_subplot(gs[0,j]))
#         ax=np.append(ax, fig.add_subplot(gs[1,j]))
#     for name in files:
#         data=np.loadtxt(os.path.join(root, name))
#         ps_pos=data[:,0]
#         counts=data[:,2]
#         counts_err=np.sqrt(counts)
#         ps_i=108.4
#         ps_f=110.8
#         counts=counts[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
#         counts_err=counts_err[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
#         ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
#         if i==0:
#             x0=ps_pos[5:-5][counts[5:-5]==np.amax(counts[5:-5])][0]
#             # counts_tot=counts.copy()
#         P0=[1922.0919973, 1281.10863081, 8.11851234, 14.18822852]#[np.amax(counts), np.amax(counts), 7,x0]
#         B0=([0,0,0,0],[5000,5000,15,100])
#         # print(P0)
#         p,cov=fit(fit_cos,ps_pos,counts, p0=P0, bounds=B0)
#         p_tot[i]=p
#         cov_tot[i]=np.diag(cov)**0.5
#         x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
#         w[i]=p[2]
#         max_ps_pos[i]=p[3]
#         ax[i].errorbar(ps_pos,counts,yerr=counts_err,fmt="ko",capsize=5)
#         # ax[i].vlines(p[3],0,fit_cos(p[3], *p),ls="dashed",label="Max ps_pos=\n"+str("%.3f" % (p[3]),))
#         ax[i].plot(x_plt,fit_cos(x_plt, *p), "b")
#         # ax[i].set_ylim([0, P0[1]+P0[1]/10])
#         # ax[i].legend(loc=1)
#         i+=1
# ax[0].set_ylabel("Counts")
# ax[1].set_ylabel("Counts")
# ax[3].set_xlabel("Phase shifter position")
# # ax[5].text(ps_pos[0],0,"Phase shifter position", va="top",ha="right")
# # plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/IFGs_indium.pdf", format="pdf", bbox_inches="tight")
# plt.show()

"""
Contrast no Indium
"""

# inf_file_name="ifg_empty_04Mar2142ff"
# sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
# cleandata=sorted_fold_path+"/Cleantxt" 
# ifg_fold_clean=cleandata+"/IFG"
# i=0
# for root, dirs, files in os.walk(ifg_fold_clean, topdown=False):
#     max_ps_pos=np.zeros((len(files)))
#     w=np.zeros((len(files)))
#     p_tot=np.zeros((len(files),4))
#     cov_tot=np.zeros((len(files),4))
#     fig = plt.figure(figsize=(11,5.5))
#     gs = GridSpec(2,len(files)//5, figure=fig)
#     ax = []
#     for j in range(len(files)//5):
#         ax=np.append(ax, fig.add_subplot(gs[0,j]))
#         ax=np.append(ax, fig.add_subplot(gs[1,j]))
#     for name in files:
#         data=np.loadtxt(os.path.join(root, name))
#         ps_pos=data[:,0]
#         counts=data[:,2]
#         counts_err=np.sqrt(counts)
#         ps_i=108.4
#         ps_f=110.8
#         counts=counts[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
#         counts_err=counts_err[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
#         ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
#         if i==0:
#             x0=ps_pos[5:-5][counts[5:-5]==np.amax(counts[5:-5])][0]
#             # counts_tot=counts.copy()
#         P0=[1922.0919973, 1281.10863081, 8.11851234, 30.18822852]#[np.amax(counts), np.amax(counts), 7,x0]
#         B0=([0,0,0,0],[5000,5000,15,100])
#         # print(P0)
#         p,cov=fit(fit_cos,ps_pos,counts, p0=P0, bounds=B0)
#         p_tot[i]=p
#         cov_tot[i]=np.diag(cov)**0.5
#         x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
#         w[i]=p[2]
#         max_ps_pos[i]=p[3]
#         if i<6:
#             ax[i].errorbar(ps_pos,counts,yerr=counts_err,fmt="ko",capsize=5)
#             # ax[i].vlines(p[3],0,fit_cos(p[3], *p),ls="dashed",label="Max ps_pos=\n"+str("%.3f" % (p[3]),))
#             ax[i].plot(x_plt,fit_cos(x_plt, *p), "b")
#             # ax[i].set_ylim([0, P0[1]+P0[1]/10])
#             # ax[i].legend(loc=1)
#         i+=1
# p_in=np.average(p_tot,axis=0)
# print(p_in)
# err_in=np.average(cov_tot,axis=0)
# print(err_in)
# ax[0].set_ylabel("Counts")
# ax[1].set_ylabel("Counts")
# ax[3].set_xlabel("Phase shifter position")
# # ax[2].set_ylabel("Counts")
# # ax[5].text(ps_pos[0],-1500,"Phase shifter position", va="center",ha="right")
# colors=["r","g","b"]
# C_in_err=(err_in[1]**2/p_in[0]**2+err_in[0]**2*p_in[1]**2/p_in[0]**4)**0.5
# C_in=p_in[1]/p_in[0]
# print(C_in,C_in_err)
# # plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/IFGs_empty.pdf", format="pdf", bbox_inches="tight")
# plt.show()

# p_in=np.average(p_tot,axis=0)
# err_in=np.average(cov_tot,axis=0)

# C_in_err=(err_in[1]**2/p_in[0]**2+err_in[0]**2*p_in[1]**2/p_in[0]**4)**0.5

# fig = plt.figure(figsize=(5,6))
# ax = fig.add_subplot(111)
# fig.suptitle("Contrast no Indium vs Indium, $\\alpha=0$")
# C_in=p_in[1]/p_in[0]
# C=C_in/C_no_in
# C_err=(C_in_err**2/C_no_in**2+C_no_in_err**2*C_in**2/C_no_in**4)**0.5
# ax.plot(x_plt,(fit_cos(x_plt, *p_no_in)+p_no_in[1]-p_no_in[0])/(2*p_no_in[1]), "b", label="No Indium")
# ax.plot(x_plt,C*(fit_cos(x_plt, *p_in)+p_in[1]-p_in[0])/(2*p_in[1])+(1-C)/2, "r", label="Indium")
# ax.set_ylim([0,1])
# # ax.legend(loc=1)
# # print(np.average(p_tot,axis=0))
# # print(np.average(cov_tot,axis=0))
# # print("ifg=", np.pi/2/p[2]+p[3])
# ax.legend(loc=3)
# plt.show()

# print("Contrast indium=", C_in)
# print("Contrast ratio=", C,"+/- ",C_err)