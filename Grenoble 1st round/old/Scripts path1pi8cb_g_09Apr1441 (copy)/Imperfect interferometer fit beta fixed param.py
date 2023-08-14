#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 15:03:44 2023

@author: aaa
"""


import warnings
from scipy.optimize import curve_fit as fit
from PIL import Image as im
import os
import numpy as np
import shutil
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
from scipy.stats import chisquare
plt.rcParams.update({'figure.max_open_warning': 0})
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

w_ps = 8.002

a1 = 1/5**0.5
a2 = 2*a1


def fit_cos(x, A, B, C, D):
    return A+B*np.cos(C*x-D)

def I_px_co(beta, chi, C, alpha, gamma):
    px=((a1*(np.cos(gamma/2)*np.cos((alpha+beta)/2)+1j*np.sin(gamma/2)*np.sin((alpha+beta)/2))+a2*np.exp(-1j*chi)*(np.cos(gamma/2)*np.cos(beta/2)+1j*np.sin(gamma/2)*np.sin(beta/2))))/(2**0.5)
    return C*np.abs(px)**2

# def I_px_in(beta, chi, eta, alpha, gamma):
#     return eta*(np.cos((alpha+beta)/2)**2+(a2/a1)**2*np.cos(beta/2)**2)/4
def I_px_in(beta, chi, eta, alpha, gamma):
    px1=np.cos(gamma/2)*np.cos((alpha+beta)/2)+1j*np.sin(gamma/2)*np.sin((alpha+beta)/2)
    px2=np.cos(gamma/2)*np.cos(beta/2)+1j*np.sin(gamma/2)*np.sin(beta/2)
    return eta*((a1)**2*np.abs(px1)**2+(a2)**2*np.abs(px2)**2)/2




inf_file_name = "path1pi8cb_g_09Apr1441"
sorted_fold_path = "/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
cleandata = sorted_fold_path+"/Cleantxt"
beta_fold_clean = cleandata+"/Beta"
plots_fold = sorted_fold_path+"/Plots/"
correct_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Corrected data/"+inf_file_name+"/Beta"

if not os.path.exists(correct_fold_path):
    os.makedirs(correct_fold_path)
else:
    shutil.rmtree(correct_fold_path)
    os.makedirs(correct_fold_path)


i = 0
for root, dirs, files in os.walk(beta_fold_clean, topdown=False):
    files = np.sort(files)
    for name in files:
        if i == 0:
            tot_data = np.loadtxt(os.path.join(root, name))
            c_pos = tot_data[:, 0]
            i = 1
        else:
            data = np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos = tot_data[::len(c_pos), -1]
# ps_i=109
# ps_f=ps_pos[-1]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
matrix = np.zeros((len(ps_pos), len(c_pos)))
matrix_err = np.zeros((len(ps_pos), len(c_pos)))
w = np.zeros(len(ps_pos))
err_b = np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i] = tot_data[:, 2][tot_data[:, -1] == ps_pos[i]]
    matrix_err[i] = tot_data[:, 2][tot_data[:, -1] == ps_pos[i]]**0.5
ps_data = np.sum(matrix, axis=1)
P0 = [(np.amax(ps_data)+np.amin(ps_data))/2,
      np.amax(ps_data)-np.amin(ps_data), 8.12, ps_pos[0]*8.11]
B0 = ([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf])
p, cov = fit(fit_cos, ps_pos, ps_data, p0=P0, bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1], 100)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)
# ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_ps = p[-2]
print(w_ps)
ps_0 = p[-1]
# print(w_ps)
c_data = np.sum(matrix, axis=0)
P0 = [(np.amax(c_data)+np.amin(c_data))/2, np.amax(c_data)/2, 0.1, c_pos[0]*0.1]
B0 = ([10, 10, 0.001, 0], [np.inf, np.inf, np.inf, 100])
p, cov = fit(fit_cos, c_pos, c_data, p0=P0, bounds=B0)
x_plt = np.linspace(c_pos[0], c_pos[-1], 100)
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.errorbar(c_pos,c_data,yerr=np.sqrt(c_data),fmt="ko",capsize=5)
# ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_c = p[-2]
# print(w_c)
# print(p[1]/p[0])
c_0 = p[-1]
# beta=np.linspace(-3*np.pi,3*np.pi,500)#c_pos.copy()#
# chi=np.linspace(-3*np.pi,3*np.pi,500)#ps_pos.copy()#
beta = w_c*c_pos-c_0
chi = w_ps*ps_pos-ps_0
alpha = np.pi/8
gamma = 0
C = 0.67
eta = 1-C

def fit_I_px(x, beta0, chi0, w_c, A, B):
    beta = (w_c*c_pos-beta0)
    chi = w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px = A*(I_px_co(beta, chi, C, alpha, 0) + I_px_in(beta, chi, eta, alpha, 0)) + B
    # print(fit_I_px)
    return fit_I_px.ravel()

P0 = (c_0, ps_0, w_c, 1, 0)
B0 = ([0, 0, 0, 0, 0], [1000, 1000, 0.3, 10, 10])
p, cov = fit(fit_I_px, range(len(matrix.ravel())), matrix.ravel()/np.amax(matrix.ravel()), bounds=B0)
print(p)
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.plot(fit_I_px(0, *p)*np.amax(matrix.ravel()), "b")
ax.plot(matrix.ravel(), "r--")
# ax.errorbar(np.arange(len(matrix.ravel())),matrix.ravel(), yerr=matrix_err.ravel(), fmt="ko", capsize=3, lw=1)
# ax.set_xlim([0,100])
f_obs=matrix.ravel()
f_exp=fit_I_px(0,*p)*np.amax(matrix.ravel())

# f_obs/=np.sum(f_obs)
# f_exp/=np.sum(f_exp)

# print((np.sum(f_obs)-np.sum(f_exp))/np.sum(f_obs))
# print(chisquare(f_obs=f_obs, f_exp=f_exp, ddof=7))

def I_px(x, beta0, chi0, w_c, A, B):
    beta = w_c*c_pos-beta0
    chi = w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px = A*(I_px_co(beta, chi, C, alpha, 0) + I_px_in(beta, chi, eta, alpha, 0)) + B
    # print(fit_I_px)
    return fit_I_px

def I_px_corr_co(x, beta0, chi0, w_c, A , B):
    beta = w_c*c_pos-beta0
    chi = w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px = I_px_co(beta, chi, C, alpha, 0) 
    # print(fit_I_px)
    return fit_I_px

def I_px_corr_in(x, beta0, chi0, w_c, A, B):
    beta = w_c*c_pos-beta0
    chi = w_ps*ps_pos-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px = I_px_in(beta, chi, eta, alpha, 0) + B
    # print(fit_I_px)
    return fit_I_px


fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
beta, chi = np.meshgrid(beta, chi)
Z = matrix
Z1 = I_px(0, *p)*np.amax(matrix)
Z2 = I_px_corr_co(0, *p)*np.amax(matrix)
Z3 = Z-I_px_corr_in(0, *p)*np.amax(matrix)
# Z=I_px_co(beta, chi, C, alpha, beta)+I_px_in(beta, chi, eta, alpha, beta)
ax.contour3D(beta, chi, Z, 40, cmap='binary')
ax.contour3D(beta, chi, Z1, 40, cmap='plasma')  # cmap='Blues')
ax.set_xlabel('$\\beta$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 45)
plt.show()


corrected_matrix=Z-(Z1-Z2)#
corrected_matrix_err=matrix_err
for i in range(len(ps_pos)):
    data_txt=np.array([c_pos, corrected_matrix[i], corrected_matrix_err[i], np.ones(len(c_pos))*ps_pos[i]])
    with open(correct_fold_path+"/beta_ps_"+str("%02d" % (i,))+".txt", 'w') as f:
        np.savetxt(f, np.transpose(data_txt),  header= "Coil_pos O-Beam err ps_pos", fmt='%.7f %.7f %.7f %.7f' )

