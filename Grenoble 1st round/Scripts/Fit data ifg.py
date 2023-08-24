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

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)
rad=np.pi/180
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
#         P0=[1922.0919973, 1281.10863081, 8.11851234, 2]#[np.amax(counts), np.amax(counts), 7,x0]
#         B0=([0,0,0,0],[5000,5000,15,100])
#         # print(P0)
#         p,cov=fit(fit_cos,ps_pos,counts, p0=P0, bounds=B0)
#         p_tot[i]=p
#         cov_tot[i]=np.diag(cov)**0.5
#         x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
#         x_plt1 = np.linspace(0, 40,100)
#         w[i]=p[2]
#         max_ps_pos[i]=p[3]
#         x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
#         fig = plt.figure(figsize=(5,5))
#         ax = fig.add_subplot(111)
#         fig.suptitle(name)
#         ax.errorbar(ps_pos,counts,yerr=counts_err,fmt="ko",capsize=5)
#         ax.vlines(p[-1]/p[-2]+ps_pos[0],0,fit_cos(p[-1]/p[-2]+ps_pos[0], *p),ls="dashed",label="Max ps_pos=\n"+str("%.3f" % (p[3]),))
#         ax.plot(x_plt,fit_cos(x_plt, *p), "b")
#         # ax.set_ylim([0, P0[1]+P0[1]/10])
#         ax.legend(loc=1)
#         i+=1

# p_no_in=np.average(p_tot,axis=0)
# err_no_in=np.average(cov_tot,axis=0)
# C_no_in_err=(err_no_in[1]**2/p_no_in[0]**2+err_no_in[0]**2*p_no_in[1]**2/p_no_in[0]**4)**0.5
# C_no_in= p_no_in[1]/p_no_in[0]

"""
Contrast no Indium alpha=0
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
ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<=(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(c_pos)))
matrix_err=np.zeros((len(ps_pos),len(c_pos)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]#/tot_data[:,5][tot_data[:,-1]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-1]==ps_pos[i]]**0.5#/tot_data[:,5][tot_data[:,-1]==ps_pos[i]]

max_pos=c_pos[np.where(matrix==np.amax(matrix))[1][0]]
counts=matrix[:,np.where(matrix==np.amax(matrix))[1][0]]
counts_err=matrix_err[:,np.where(matrix==np.amax(matrix))[1][0]]
P0=[1922.0919973, 1281.10863081, 8.11851234, 14.18822852]#[(np.amax(counts)+np.amin(counts))/2, np.amax(counts)-np.amin(counts), 0.1,ps_pos[0]]
B0=([0,0,0,0],[5000,5000,15,100])
p_no_in,cov=fit(fit_cos,ps_pos,counts, p0=P0, bounds=B0)
# print(p_no_in)
# print(np.diag(cov)**0.5)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
fig.suptitle("coil pos="+str(max_pos)+" mm")
x_plt_0 = np.linspace(ps_pos[0], ps_pos[-1],100)
ax.errorbar((ps_pos),counts,yerr=counts_err,fmt="ko",capsize=5)
ax.plot((x_plt_0),fit_cos(x_plt_0, *p_no_in), "b")
ax.set_xlabel("$\\beta$ (mm)")
plt.show()

err_no_in=np.diag(cov)**0.5
C_no_in_err=(err_no_in[1]**2/p_no_in[0]**2+err_no_in[0]**2*p_no_in[1]**2/p_no_in[0]**4)**0.5
C_no_in= p_no_in[1]/p_no_in[0]
fig = plt.figure()
fig.suptitle(inf_file_name+".inf")
ax = plt.axes(projection='3d')
Z=matrix
x=c_pos
y=ps_pos
X, Y = np.meshgrid(x, y)
ax.contour3D(X, Y, Z, 30, cmap='binary')
ax.plot3D(15+ps_pos*0,ps_pos, counts, 'red', label="Contrast curve (coil pos=15)")
ax.set_xlabel('Coil')
ax.set_ylabel('PS')
ax.set_zlabel('z')
ax.view_init(45, 45)
ax.legend()
plt.show()

print("Contrast no indium=", C_no_in,"+/- ",C_no_in_err)

# """
# Contrast Indium
# """

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
#     for name in files:
#         data=np.loadtxt(os.path.join(root, name))
#         ps_pos=data[:,0]
#         ps_pos_tot=ps_pos.copy()
#         counts=data[:,2]
#         counts_err_rel=np.sqrt(counts)/counts
#         counts/=data[:,7]
#         counts_err=counts_err_rel*counts
#         if i==0:
#             x0=ps_pos[counts==np.amax(counts)][0]
#             counts_avg=counts.copy()*0
#             counts_err_avg=counts_err.copy()*0
#         counts_avg+=counts.copy()/len(files)
#         counts_err_avg+=counts_err.copy()/len(files)
        
#         ps_i=108.4
#         ps_f=110.8
#         counts=counts[abs(ps_pos-(ps_i+ps_f)/2)<=(ps_f-ps_i)/2]
#         counts_err=counts_err[abs(ps_pos-(ps_i+ps_f)/2)<=(ps_f-ps_i)/2]
#         ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<=(ps_f-ps_i)/2] 
#             # counts_tot=counts.copy()
#         P0=[np.amax(counts), np.amax(counts), 8.11,ps_pos[0]*8.11]#[1922.0919973, 1281.10863081, 8.11851234, 14.18822852]#
#         B0=([0,0,0,0],[5000,5000,15,1000])
#         # print(P0)
#         p,cov=fit(fit_cos,ps_pos,counts, p0=P0, bounds=B0)
#         p_tot[i]=p
#         cov_tot[i]=np.diag(cov)**0.5
#         x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
#         x_plt1 = np.linspace(0, 40,100)
#         w[i]=p[2]
#         max_ps_pos[i]=p[3]
#         #x_plt1[fit_cos(x_plt1, *p)==np.amax(fit_cos(x_plt1, *p))]
#         # fig = plt.figure(figsize=(5,5))
#         # ax = fig.add_subplot(111)
#         # fig.suptitle(name)
#         # ax.errorbar(ps_pos,counts,yerr=counts_err,fmt="ko",capsize=5)
#         # # ax.vlines(p[3],0,fit_cos(p[3], *p),ls="dashed",label="Max ps_pos=\n"+str("%.3f" % (p[3]),))
#         # ax.plot(x_plt,fit_cos(x_plt, *p), "b")
#         # ax.set_ylim([0, P0[1]+P0[1]/10])
#         # ax.legend(loc=1)
#         i+=1
# counts_err_avg_rel=counts_err_avg/counts_avg
# p_in=np.average(p_tot,axis=0)
# err_in=np.average(cov_tot,axis=0)
# with open("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Scripts/p_plus_fit_param.txt", 'w') as f:
#          np.savetxt(f,(p_in, err_in), header="a b c d")
# print(p_in)
# print(err_in)
# C_in_err=(err_in[1]**2/p_in[0]**2+err_in[0]**2*p_in[1]**2/p_in[0]**4)**0.5
# C_in=p_in[1]/p_in[0]
# C=C_in/C_no_in
# C_err=(C_in_err**2/C_no_in**2+C_no_in_err**2*C_in**2/C_no_in**4)**0.5
# x_plt=np.linspace(ps_pos_tot[0],ps_pos_tot[-1],1000)

# scaled=C*(counts_avg-p_in[0]+p_in[1])/(2*p_in[1])+(1-C)/2
# counts_err_avg=counts_err_avg_rel*scaled
# x0=p_in[-1]+2*np.pi
# x_plt_pi=(ps_pos_tot*p_in[-2]-x0)/np.pi
# data_txt=np.transpose(np.array([ps_pos_tot,scaled,counts_err_avg]))

# with open("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Scripts/p_plus.txt", 'w') as f:
#           np.savetxt(f,data_txt, header="chi(pi) scaled_values error")


# fig = plt.figure(figsize=(5,6))
# ax = fig.add_subplot(111)

# # fig.suptitle("Contrast no Indium vs Indium, $\\alpha=0$")
# ax.errorbar((ps_pos_tot*p_in[-2]-x0)/np.pi,scaled,yerr=counts_err_avg,fmt="ko", capsize=5, label="Scaled data $p_+$")
# ax.errorbar((ps_pos_tot*p_in[-2]-x0+np.pi)/np.pi,scaled,yerr=counts_err_avg,fmt="go", capsize=5, label="Scaled data $p_-$")
# # ax.plot(x_plt,(fit_cos(x_plt, *p_no_in)+p_no_in[1]-p_no_in[0])/(2*p_no_in[1]), "b", label="No Indium")
# # ax.plot(x_plt,C/2*fit_cos(x_plt, 0,1,1,0)+1/2, "r", label="$P_+$")
# ax.plot((x_plt*p_in[-2]-x0)/np.pi,C*(fit_cos(x_plt, *p_in)+p_in[1]-p_in[0])/(2*p_in[1])+(1-C)/2, "r", label="$p_+$")
# ax.plot((x_plt*p_in[-2]-x0)/np.pi,C*(fit_cos(x_plt, p_in[0],-p_in[1],*p_in[2:])+p_in[1]-p_in[0])/(2*p_in[1])+(1-C)/2, "b", label="$p_-$")


# print("p+max=",(1-C)/2)
# print("err=",counts_err_avg[abs(scaled-0.1)==np.amin(abs(scaled-0.1))])
# ax.set_ylim([-0.1,1])
# ax.set_xlim([(ps_pos[0]*p_in[-2]-x0)/np.pi,(ps_pos[-1]*p_in[-2]-x0)/np.pi])
# ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
# ax.set_xlabel("$\chi$ ($\pi$)")
# ax.legend(loc=4, ncol=2)
# # plt.savefig("/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/Report/Images/Pppm.pdf", format="pdf", bbox_inches="tight")
# plt.show()

# print("Contrast indium=", C_in,"+/- ",C_in_err)
# print("Contrast ratio=", C,"+/- ",C_err)


