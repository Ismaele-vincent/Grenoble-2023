#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 17:39:52 2023

@author: aaa
"""

import numpy as np
import matplotlib.pyplot as plt

def w1p(x,a1,a2):
    return (1/(1+a2/a1*np.exp(-1j*x)))

x=np.linspace(0,4*np.pi, 100)

alpha=1
y1=w1p(x,1,2)*alpha
y2=(1-w1p(x,1,2))*alpha
# plt.plot(x/np.pi,y1.real,label="Re{"+"$\omega_{1+}$}")
# plt.plot(x/np.pi,y1.imag,label="Im{"+"$\omega_{1+}$}")
plt.plot(x/np.pi,y2.real,label="Re{"+"$\omega_{2+}$}")
plt.plot(x/np.pi,y2.imag,label="Im{"+"$\omega_{2+}$}")
plt.legend()
# print(np.amax(y1.real)-np.amin(y1.real))