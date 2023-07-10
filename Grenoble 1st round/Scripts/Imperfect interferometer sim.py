#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 12:15:08 2023

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
a21=2
def exp_w1p(x,x0):
    return alpha*((1/(1+a21*np.exp(-1j*(w_ps*(x-x0)))))).real

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
beta=0
g=np.linspace(-np.pi/2,np.pi/2,100)
c=np.linspace(-2*np.pi,2*np.pi,25)
alpha=np.pi/2

def I_px(x, gamma0, chi0, C, eta):
    gamma=g-gamma0
    chi=c-chi0
    gamma, chi = np.meshgrid(gamma, chi)
    fit_I_px=I_px_co(beta, chi, C, alpha, gamma) #+ I_px_in(0, chi, eta, alpha, gamma)
    return fit_I_px

def I_mx(x, gamma0, chi0, C, eta):
    gamma=g-gamma0
    chi=c-chi0
    gamma, chi = np.meshgrid(gamma, chi)
    fit_I_mx=I_mx_co(beta, chi, C, alpha, gamma) #+ I_mx_in(0, chi, eta, alpha, gamma)
    return fit_I_mx

p=[0,0,1,0]
# beta = exp_w1p(c[i], 0)
Z= I_px(0,*p)
# Z0= I_px(0,0,0,1,0)
# Z2p= I_px(0,2*np.pi,0,1,0.2)
# Z= (I_px(0,*p)-I_mx(0,*p))/(I_px(0,*p)+I_mx(0,*p))

# fig = plt.figure(figsize=(10,10))
# ax = plt.axes(projection='3d')
# x, y = np.meshgrid(g, c)
# ax.contour3D(x, y, Z, 30, cmap='plasma')#cmap='Blues')
# ax.set_xlabel('$\\gamma$')
# ax.set_ylabel('$\chi$')
# ax.set_zlabel('z')
# ax.view_init(40, 45)
wg=np.zeros((len(c)))
for i in range(len(c)):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    ax.set_title(str(c[i]))
    ax.plot(g,Z[0])
    ax.vlines(0,0,np.amax(Z[0]))
    ax.plot(g,Z[i])
    # ax.set_ylim([0,0.05])
    # ax.set_xlim([2.5,5])
    # ax.plot(g,g*0+1)
#     wg[i]=Z2p[i,0]-Z0[i,0]

# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# ax.plot(c,g)
plt.show()