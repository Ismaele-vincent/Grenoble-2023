#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 17:14:55 2023

@author: aaa
"""
import numpy as np
import shutil
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.gridspec import GridSpec
from scipy.special import jv

mu_N=-9.6623651*1e-27 #J/T
hbar= 6.62607015*1e-34/(2*np.pi)  #J s
B1=10
w1=2*np.pi*1e3
T=19.4e-6

n=5
alpha1=mu_N*B1/(hbar*w1)
alpha=2*alpha1*np.sin(w1*T/2)
eta1=phi1+(w1*T+np.pi)/2
psi=0+1j*0

for i in range(-n,n):
    w_n=w_0+n*w_1
    k_n=(k_0**2-2*m*mu*B_0/hbar+2*m*n*w1/hbar+0*1j)**0.5
    psi+=jv(i, alpha)*np.exp(-1j*n*eta1)*np.exp(1j*k_n*x)*np.exp(-1j*w_n*t)
