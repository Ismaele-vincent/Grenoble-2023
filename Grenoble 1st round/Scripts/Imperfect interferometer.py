#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 17 14:19:05 2023

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
import warnings
warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

w_ps=8.002

a1=1/5**0.5
a2=2*a1

def fit_cos(x,A,B,C,D):
    return A+B*np.cos(C*x-D)

def I_px_co(beta, chi, C, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    return C*((a1*np.cos(d/2))**2+(a2*np.cos(e/2))**2+2*a1*a2*np.cos(d/2)*np.cos(e/2)*np.cos(chi))/4

def I_px_in(beta, chi, eta, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    return eta*(np.cos(d/2)**2+(a2/a1)**2*np.cos(e/2)**2)/4

def I_mx_co(beta, chi, C, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    r=np.arctan(gamma/beta)
    p=np.arctan(gamma/(beta+gamma))
    return C*((a1*np.sin(d/2))**2+(a2*np.sin(e/2))**2+2*a1*a2*np.sin(d/2)*np.sin(e/2)*np.cos(chi+r-p))/4

def I_mx_in(beta, chi, eta, alpha, gamma):
    d=((alpha+beta)**2+gamma**2)**0.5
    e=(beta**2+gamma**2)**0.5
    return eta*(np.sin(d/2)**2+(a2/a1)**2*np.sin(e/2)**2)/4

beta=np.linspace(-3*np.pi,3*np.pi,500)
alpha=np.pi/8
chi=np.linspace(-3*np.pi,3*np.pi,500)
gamma=0
C=0.8
eta=1-C
# print(I_px_co(beta, chi, C, alpha, gamma))
beta, chi = np.meshgrid(beta, chi)


I_px=I_px_co(beta, chi, C, alpha, gamma)+I_px_in(beta, chi, eta, alpha, gamma)
I_mx=I_mx_co(beta, chi, C, alpha, gamma)+I_mx_in(beta, chi, eta, alpha, gamma)
I_x=(I_px-I_mx)/(I_px+I_mx)
I_x_th=(I_px_co(beta, chi, C, alpha, gamma)-I_mx_co(beta, chi, C, alpha, gamma))/(I_px_co(beta, chi, C, alpha, gamma)+I_mx_co(beta, chi, C, alpha, gamma))

fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection='3d')
Z=I_px_co(beta, chi, C, alpha, gamma)+I_px_in(beta, chi, eta, alpha, gamma)
ax.contour3D(beta, chi, Z, 30, cmap='binary')
ax.set_xlabel('$\\beta$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 45)

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.plot(beta, I_x, "b")
# ax.plot(beta, I_x_th, "r")



