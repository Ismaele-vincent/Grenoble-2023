#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 10:54:48 2023

@author: aaa
"""

import numpy as np
import inspect,os,time
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import socket
import shutil
from scipy.optimize import curve_fit as fit
from scipy.stats import chisquare as cs
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
from datetime import datetime
from multiprocessing import Pool
from scipy.stats import norm
from scipy.stats import cosine
from scipy.stats import exponnorm
fit_name="bcr_1_2_phi_no_zeta_0"
p_name=["$(b_c \\rho)_1$","$(b_c \\rho)_2$", "$\mu$", "$\sigma$","$\\tau$", "$x_0$","$\phi$"]
p_units=[" $1/\mu m^2$"," $1/\mu m^2$"," nm", " nm", "", " deg", "$ \pi$"]

pi=np.pi
rad=pi/180

sorted_fold_path="/home/aaa/Desktop/Thesis2/Paper/Sorted data/" #insert folder of sorted meausements files
tiltangles=[0,40,48,61,69,71,"74_52","75_64","76_76","77_88",79,80,81]
foldername=[]
for i in range(len(tiltangles)):
    foldername.append(str(tiltangles[i])+"deg")

tilt=[0,40,48,61,69,71,74.52,75.64,76.76,77.88,79,80,81]
n_pixel = 16384 #number of pixels in one measurement
krange=range(13)
integ_m1=np.zeros((2,len(krange)))
integ_m2=np.zeros((2,len(krange)))
integ_m1[0]=tilt
integ_m2[0]=tilt
for k in krange:
    print(foldername[k])
    data_analysis = sorted_fold_path+foldername[k]+"/Data Analysis/"
    diff_eff =  np.loadtxt(data_analysis+foldername[k]+"_diff_eff_new.mpa",skiprows=1)
    fit_res =  np.loadtxt(data_analysis+foldername[k]+"_fit_results_"+fit_name+".mpa",skiprows=1)
    p=fit_res[0]
    th=np.linspace(diff_eff[0,0],diff_eff[-1,0],100)
    fig = plt.figure(figsize=(5,5),dpi=200)
    for i in range(1,len(diff_eff[0,:])//2):
        spl=UnivariateSpline(diff_eff[:,0], diff_eff[:,2*i],k=5,s=0)
        ax = fig.add_subplot(111)
        ax.set_title(str(tilt[k])+" deg")
        ax.plot(th, spl(th))
        ax.plot(diff_eff[:,0], diff_eff[:,2*i], "k.")
        ax.set_xlim([-0.7,0.7])
        if i==2:
            integ_m1[1,k]=spl.antiderivative()(th)[-1]
        if i==1:
            integ_m2[1,k]=spl.antiderivative()(th)[-1]
        plt.savefig("/home/aaa/Desktop/Thesis2/Paper/Plot int intensities/"+str(tilt[k])+"deg.png")
fig = plt.figure(figsize=(5,5), dpi=200)
ax = fig.add_subplot(111)
ax.plot(integ_m1[0],integ_m1[1], "k^-", label="Int. intensities -1")
ax.plot(integ_m2[0],integ_m2[1], "r^-", label="Int. intensities -2")
ax.legend()
plt.savefig("/home/aaa/Desktop/Thesis2/Paper/Plot int intensities/Pendel.png")
        
        
        
        