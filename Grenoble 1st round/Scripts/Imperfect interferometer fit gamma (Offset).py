#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 12:15:08 2023

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
plt.rcParams.update({'figure.max_open_warning': 0})
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

rad=np.pi/180
w_ps = 8.002

a2 = 1/5**0.5
a1 = 2*a2


def fit_cos(x, A, B, C, D):
    return A+B*np.cos(C*x-D)


def I_px_co(gamma, chi, C, alpha, beta):
    d = ((alpha+beta)**2+gamma**2)**0.5
    e = (beta**2+gamma**2)**0.5
    return C*((a1*np.cos(d/2))**2+(a2*np.cos(e/2))**2+2*a1*a2*np.cos(d/2)*np.cos(e/2)*np.cos(chi))/4


def I_px_in(gamma, chi, eta, alpha, beta):
    d = ((alpha+beta)**2+gamma**2)**0.5
    e = (beta**2+gamma**2)**0.5
    return eta*(np.cos(d/2)**2+(a2/a1)**2*np.cos(e/2)**2)/4


def I_mx_co(gamma, chi, C, alpha, beta):
    d = ((alpha+beta)**2+gamma**2)**0.5
    e = (beta**2+gamma**2)**0.5
    r = np.arctan(gamma/beta)
    p = np.arctan(gamma/(beta+gamma))
    return C*((a1*np.sin(d/2))**2+(a2*np.sin(e/2))**2+2*a1*a2*np.sin(d/2)*np.sin(e/2)*np.cos(chi+r-p))/4


def I_mx_in(gamma, chi, eta, alpha, beta):
    d = ((alpha+beta)**2+gamma**2)**0.5
    e = (beta**2+gamma**2)**0.5
    return eta*(np.sin(d/2)**2+(a2/a1)**2*np.sin(e/2)**2)/4


inf_file_name = "path1pi8cb_g_09Apr1441"
sorted_fold_path = "/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Sorted data/"+inf_file_name
cleandata = sorted_fold_path+"/Cleantxt"
gamma_fold_clean = cleandata+"/Gamma"
plots_fold = sorted_fold_path+"/Plots/"
correct_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Corrected data/"+inf_file_name+"/Gamma"

if not os.path.exists(correct_fold_path):
    os.makedirs(correct_fold_path)
else:
    shutil.rmtree(correct_fold_path)
    os.makedirs(correct_fold_path)

i = 0
for root, dirs, files in os.walk(gamma_fold_clean, topdown=False):
    files = np.sort(files)
    for name in files:
        if i == 0:
            tot_data = np.loadtxt(os.path.join(root, name))
            c_amp = tot_data[:, 0]
            i = 1
        else:
            data = np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos = tot_data[::len(c_amp), -2]
# ps_i=109
# ps_f=ps_pos[-1]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
matrix = np.zeros((len(ps_pos), len(c_amp)))
matrix_err = np.zeros((len(ps_pos), len(c_amp)))
w = np.zeros(len(ps_pos))
err_b = np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i] = tot_data[:, 2][tot_data[:, -2] == ps_pos[i]]
    matrix_err[i] = tot_data[:, 2][tot_data[:, -2] == ps_pos[i]]**0.5
ps_data = np.sum(matrix, axis=1)
P0 = [(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 8.12, ps_pos[0]*8.11]
B0 = ([0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf])
p, cov = fit(fit_cos, ps_pos, ps_data, p0=P0, bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1], 100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_ps = p[-2]
ps_0 = p[-1]
# print(w_ps)
c_data = np.sum(matrix, axis=0)
P0 = [(np.amax(c_data)+np.amin(c_data))/2, np.amax(c_data)/2, 2.8, c_amp[0]]
B0 = ([100, 100, 0.001, -10], [np.inf, np.inf, np.inf, 10])
p, cov = fit(fit_cos, c_amp, c_data, p0=P0, bounds=B0)
x_plt = np.linspace(c_amp[0], c_amp[-1], 100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(c_amp,c_data,yerr=np.sqrt(c_data),fmt="ko",capsize=5)
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
# ax.vlines(p[-1],0,fit_cos(p[-1], *p),ls="dashed")
w_c = p[-2]
# print(w_c)
# print(p[1]/p[0])
c_0 = p[-1]
# gamma=np.linspace(-3*np.pi,3*np.pi,500)#c_amp.copy()#
# chi=np.linspace(-3*np.pi,3*np.pi,500)#ps_pos.copy()#
alpha = np.pi/8
beta = 0
C = 0.8
eta = 1-C


def fit_I_px(x, gamma0, chi0, w_c, w_ps, C, eta, A):
    gamma = w_c*c_amp-gamma0
    chi = w_ps*ps_pos-chi0
    gamma, chi = np.meshgrid(gamma, chi)
    fit_I_px = (I_px_co(gamma, chi, C, alpha, 0) + eta/4*A*(0.5+0*a1*a2*np.cos(chi)) + I_px_in(gamma, chi, eta, alpha, 0))
    # print(fit_I_px)
    return fit_I_px.ravel()


P0 = (-30, ps_0, w_c, w_ps, 2, 1, 0)
B0 = ([-50, 800, 0, 7, 0, 0, 0], [50, 1050, 5, 10, np.inf, np.inf, 10])
p, cov = fit(fit_I_px, range(len(matrix.ravel())), matrix.ravel()/np.amax(matrix.ravel()), bounds=B0)
print(p)
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
ax.plot(fit_I_px(0, *p), "b")
ax.plot(matrix.ravel()/np.amax(matrix.ravel()), "r--", lw=1)
# ax.plot(fit_I_px(0, *p)-matrix.ravel()/np.amax(matrix.ravel()), "b")
# ax.set_xlim([0,150])


def I_px(x, gamma0, chi0, w_c, w_ps, C, eta, A):
    gamma = w_c*c_amp-gamma0
    chi = w_ps*ps_pos-chi0
    gamma, chi = np.meshgrid(gamma, chi)
    fit_I_px = I_px_co(gamma, chi, C, alpha, 0) +  eta/4*A*(0.5+0*a1*a2*np.cos(chi))+I_px_in(gamma, chi, eta, alpha, 0)
    # print(fit_I_px)
    return fit_I_px

def I_px_corr(x, gamma0, chi0, w_c, w_ps, C, eta, A):
    g = w_c*c_amp-gamma0
    c = w_ps*ps_pos-chi0
    gamma, chi = np.meshgrid(g, c)
    fit_I_px = I_px_co(gamma, chi, C, alpha, 0) 
    # print(fit_I_px)
    return fit_I_px

chi = w_ps*ps_pos-p[1]
gamma = w_c*c_amp-p[0]
fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
x, y= np.meshgrid(gamma, chi)
Z = matrix
Z1 = I_px(0, *p)*np.amax(matrix)
Z2 = I_px_corr(0, *p)*np.amax(matrix)
# Z=I_px_co(beta, chi, C, alpha, beta)+I_px_in(beta, chi, eta, alpha, beta)
ax.contour3D(x, y, Z, 40, cmap='binary')
ax.contour3D(x, y, Z1, 40, cmap='plasma')  # cmap='Blues')
ax.set_xlabel('$\\gamma$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 45)
plt.show()

corrected_matrix=Z1#I_px_new(*p)*np.amax(matrix)
corrected_matrix_err=matrix_err.copy()
for i in range(len(ps_pos)):
    data_txt=np.array([c_amp, corrected_matrix[i], corrected_matrix_err[i], np.ones(len(c_amp))*ps_pos[i]])
    with open(correct_fold_path+"/gamma_ps_"+str("%02d" % (i,))+".txt", 'w') as f:
        np.savetxt(f, np.transpose(data_txt),  header= "Coil_pos O-Beam err ps_pos")