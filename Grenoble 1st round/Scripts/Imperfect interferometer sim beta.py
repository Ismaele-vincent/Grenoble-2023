#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:24:24 2023

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

# def I_px_co(beta, chi, C, alpha, gamma):
#     d=((alpha+beta)**2+gamma**2)**0.5
#     e=(beta**2+gamma**2)**0.5
#     return C*((a1*np.cos(d/2))**2+(a2*np.cos(e/2))**2+2*a1*a2*np.cos(d/2)*np.cos(e/2)*np.cos(chi))/4

def I_px_co(beta, chi, C, alpha, gamma):
    d=(alpha+beta)
    e=beta
    return C*((a1*np.cos(d/2))**2+(a2*np.cos(e/2))**2+2*a1*a2*np.cos(d/2)*np.cos(e/2)*np.cos(chi))/2


def I_px_co1(beta, chi, C, alpha, gamma):
    # beta*=-1
    # alpha*=-1
    return C/8*(2+2*np.cos(alpha/2)*np.cos(chi)*a2+np.cos(gamma)*(np.cos(beta)+np.cos(alpha+beta)-2*a1*np.sin(alpha/2)*np.sin(alpha/2+beta)+2*np.cos(alpha/2+beta)*np.cos(chi)*a2)+2*np.sin(alpha/2)*a2*np.sin(gamma)*np.sin(chi))

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
gamma=0
b=np.linspace(-2*np.pi,2*np.pi,100)
c=np.linspace(0,2*np.pi,25)
alpha=np.pi/4

def I_px(x, beta0, chi0, C, eta):
    beta=b-beta0
    chi=c-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px=I_px_co(beta, chi, C, alpha, 0) #+ I_px_in(beta, chi, eta, alpha, 0)
    # print(fit_I_px)
    return fit_I_px

def I_px1(x, beta0, chi0, C, eta):
    beta=b-beta0
    chi=c-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_px=I_px_co1(beta, chi, C, alpha, 0) #+ I_px_in(beta, chi, eta, alpha, 0)
    # print(fit_I_px)
    return fit_I_px

def I_mx(x, beta0, chi0, C, eta):
    beta=b-beta0
    chi=c-chi0
    beta, chi = np.meshgrid(beta, chi)
    fit_I_mx=I_mx_co(beta, chi, C, alpha, 0) #+ I_mx_in(beta, chi, eta, alpha, 0)
    # print(fit_I_mx)
    return fit_I_mx

p=[0,0,1,0]

# Z= (I_px(0,*p)-I_mx(0,*p))/(I_px(0,*p)+I_mx(0,*p))
Z= I_px(0,*p)
Z1= I_px1(0,*p)
# fig = plt.figure(figsize=(10,10))
# ax = plt.axes(projection='3d')
# x, y = np.meshgrid(g, c)
# ax.contour3D(x, y, Z, 30, cmap='plasma')#cmap='Blues')
# ax.set_xlabel('$\\gamma$')
# ax.set_ylabel('$\chi$')
# ax.set_zlabel('z')
# ax.view_init(40, 45)
w_im=np.zeros((len(c)))
for i in range(len(c)):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    ax.plot(b,Z[0],"r")
    ax.plot(b,Z[i],"g")
    ax.plot(b,Z1[0], "b")
    ax.plot(b,Z1[i], "k")
    # ax.set_ylim([0,0.05])
    # w_im[i]=(np.amax(Z[i])-np.amin(Z[i]))/2
# fig = plt.figure(figsize=(5,5))
# ax = fig.add_subplot(111)
# # ax.plot(c,w_im)
# y=np.arccosh(1/w_im)#2**0.5*(1-w_im)**0.5/alpha
# y[c>np.pi]*=-1
# ax.plot(c,y)
# ax.plot(c,np.arccosh(1/w_im)/alpha)
plt.show()