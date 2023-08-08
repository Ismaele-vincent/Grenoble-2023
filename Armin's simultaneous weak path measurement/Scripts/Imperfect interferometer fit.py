#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 17:55:56 2023

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
rad = np.pi/180
alpha1 = np.pi/9
alpha2 = 0
C = 0.7


# def fit_cos(x, A, B, C, D):
#     return A+B*np.cos(C*x+D)

def fit_cos(x, A, B, C, D):
    return A+B*np.cos(C*x+D)


def fit_ph1(x, A, x0):
    return A+np.arctan(C*np.sin(x+x0)/(2**-0.5+C*np.cos(x+x0)))


def fit_ph2(x, A, x0):
    return -A-np.arctan(C*np.sin(x+x0)/(2**-0.5+C*np.cos(x+x0)))


def fit_amp(x, A, B, D, x0):
    return A+B*(C**2*np.sin(D*x+x0)**2 + (np.cos(alpha1/2)+C*np.cos(D*x+x0))**2)**0.5


def fit_mean(x, A, B, D, x0):
    return A+B*np.cos(D*x+x0)*np.cos(alpha1/2)*np.cos(alpha2/2)


# "TOF_PS_RF1_pi2_1h_10Jun1748" #"TOF_PS_RF2_weak_1h.tof_18Jun1231"#"TOF_PS_RF2_pi2_1h_12Jun0717"#
inf_file_name = "TOF_PS_RF1_weak_1h_17Jun1717"
sorted_fold_path = "/home/aaa/Desktop/Fisica/PhD/2023//Armin's simultaneous weak path measurement/Sorted data/"+inf_file_name
cleandata = sorted_fold_path+"/Cleantxt"
# correct_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 1st round/exp_3-16-13/Corrected data/"+inf_file_name+"/Beta"

# if not os.path.exists(correct_fold_path):
#     os.makedirs(correct_fold_path)
# else:
#     shutil.rmtree(correct_fold_path)
#     os.makedirs(correct_fold_path)


i = 0
for root, dirs, files in os.walk(cleandata, topdown=False):
    files = np.sort(files)
    for name in files:
        if i == 0:
            tot_data = np.loadtxt(os.path.join(root, name))
            time = tot_data[:, 0]
            i = 1
        else:
            data = np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos = tot_data[::len(time), -1]
# print(ps_pos)
# ps_i=109
# ps_f=ps_pos[-1]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2]
matrix = np.zeros((len(ps_pos), len(time)))
matrix_err = np.zeros((len(ps_pos), len(time)))
w1 = np.zeros(len(ps_pos))
w2 = np.zeros(len(ps_pos))
w3 = np.zeros(len(ps_pos))
w1_err = np.zeros(len(ps_pos))
w2_err = np.zeros(len(ps_pos))
w3_err = np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i] = tot_data[:, 3][tot_data[:, -1] == ps_pos[i]]
    matrix_err[i] = tot_data[:, 3][tot_data[:, -1] == ps_pos[i]]**0.5
    P0 = [(np.amax(matrix[i])+np.amin(matrix[i]))/2, (np.amax(matrix[i]) - np.amin(matrix[i]))/2, 0.39, fit_ph1(ps_pos[i]*rad, 0.5, 0.2)]
    B0 = ([0, 10, 0, -6], [np.inf, np.inf, np.inf, 6])
    p, cov = fit(fit_cos, time[2:-2], matrix[i][2:-2],
                 p0=P0, bounds=B0, sigma=matrix_err[i][2:-2])
    err = np.diag(cov)**0.5
    w1[i] = p[1]
    w1_err[i] = err[1]
    w2[i] = p[-1]
    w2_err[i] = err[-1]
    w3[i] = p[0]
    w3_err[i] = err[0]
    x_plt = np.linspace(time[0], time[-1], 100)
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    ax.errorbar(time,matrix[i],yerr=matrix_err[i],fmt="ko",capsize=5)
    # ax.plot(time,matrix[i],"k-")
    ax.plot(x_plt,fit_cos(x_plt, *p), "b", lw=5)
    ax.plot(x_plt,fit_cos(x_plt, *p), "b")
    # ax.set_xlim([100,200])
# P0 = [30, 80, 0.7, 1.5, -1]
# B0 = ([0, 0, 0.5, 0.01, -5], [100, 200, 1, 50, 5])
P0 = [30, 80, 1.5, -1]
B0 = ([0, 0, 0.5,  -5], [200, 200, 50, 5])
p, cov = fit(fit_amp, ps_pos*rad, w1, p0=P0, bounds=B0)
print(p)
print(np.diag(cov)**0.5)
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.set_title("Amplitudes")
ax.errorbar(ps_pos, w1, yerr=w1_err, fmt="ko", capsize=5)
x_plt = np.linspace(ps_pos[0]-20, ps_pos[-1]+20, 500)
ax.plot(x_plt, fit_amp(x_plt*rad, *p), "b")

# P0 = [30, 80, 1,  -1]
# B0 = ([0, 0, 0.01, -5], [100, 200, 50, 5])

P0 = [30, 80, 1,  -1]
B0 = ([0, 0, 0.01, -5], [300, 200, 50, 5])
p, cov = fit(fit_mean, ps_pos*rad, w3, p0=P0, bounds=B0)
print(p)
print(np.diag(cov)**0.5)
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.set_title("Mean")
ax.errorbar(ps_pos, w3, yerr=w3_err, fmt="ko", capsize=5)
ax.plot(x_plt, fit_mean(x_plt*rad, *p), "b")

# P0 = [0, 0.7, 0]
# B0 = ([-5, 0, -5], [5, 1, 5])
P0 = [0, 0]
B0 = ([-5, -5], [5, 5])
p, cov = fit(fit_ph1, ps_pos*rad, w2, p0=P0, bounds=B0)
print(p)
print(np.diag(cov)**0.5)
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.set_title("Phase")
ax.errorbar(ps_pos, w2, yerr=w2_err, fmt="ko", capsize=5)
ax.plot(x_plt, fit_ph1(x_plt*rad, *p), "b")

w4=w1/w3
w4_err=(w1_err**2+w3_err**2*w2**2/w3**2)**0.5/w3

# P0 = [30, 80, 1,  -1]
# B0 = ([0, 0, 0.01, -5], [300, 200, 50, 5])
# p, cov = fit(fit_mean, ps_pos*rad, w3, p0=P0, bounds=B0)
# print(p)
# print(np.diag(cov)**0.5)
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)
ax.set_title("Mean")
ax.errorbar(ps_pos, w1/w3, yerr=w4_err, fmt="ko", capsize=5)
# ax.plot(x_plt, fit_mean(x_plt*rad, *p), "b")



# fig = plt.figure(figsize=(10, 10))
# ax = plt.axes(projection='3d')
# time, chi = np.meshgrid(time, ps_pos)
# Z = matrix
# ax.contour3D(time, chi, Z, 40, cmap='binary')
# ax.set_xlabel('$time$')
# ax.set_ylabel('$\chi$')
# ax.set_zlabel('z')
# ax.view_init(0, 0)
# plt.show()
