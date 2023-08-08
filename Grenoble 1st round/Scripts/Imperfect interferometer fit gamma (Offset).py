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

rad=np.pi/180
w_ps=8

a1 = 1/5**0.5
a2 = 2*a1

alpha=np.pi/8
a21=2

def exp_w1p(x):
    return alpha*(1-(1/(1+a21*np.exp(-1j*x)))).real


def fit_cos(x, A, B, C, D):
    return A+B*np.cos(C*x-D)

def I_px_co(beta, chi, C, alpha, gamma):
    px=((a1*(np.cos(gamma/2)*np.cos((alpha+beta)/2)+1j*np.sin(gamma/2)*np.sin((alpha+beta)/2))+a2*np.exp(-1j*chi)*(np.cos(gamma/2)*np.cos(beta/2)+1j*np.sin(gamma/2)*np.sin(beta/2))))/(2**0.5)
    return C*np.abs(px)**2


def I_px_in(beta, chi, eta, alpha, gamma):
    px1=np.cos(gamma/2)*np.cos((alpha+beta)/2)+1j*np.sin(gamma/2)*np.sin((alpha+beta)/2)
    px2=np.cos(gamma/2)*np.cos(beta/2)+1j*np.sin(gamma/2)*np.sin(beta/2)
    return eta*(np.abs(px1)**2+(a2/a1)**2*np.abs(px2)**2)/4


inf_file_name = "path2pi8cb_g_12Apr1724"
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

i=0
for root, dirs, files in os.walk(gamma_fold_clean, topdown=False):
    files=np.sort(files)
    for name in files:
        if i==0:
            tot_data=np.loadtxt(os.path.join(root, name))
            coil=tot_data[:,0]
            i=1
        else:
            data=np.loadtxt(os.path.join(root, name))
            tot_data = np.vstack((tot_data, data))
ps_pos=tot_data[::len(coil),-2]
# ps_i=109
# ps_f=ps_pos[-1]
# ps_pos=ps_pos[abs(ps_pos-(ps_i+ps_f)/2)<(ps_f-ps_i)/2] 
matrix=np.zeros((len(ps_pos),len(coil)))
matrix_err=np.zeros((len(ps_pos),len(coil)))
w=np.zeros(len(ps_pos))
err_b=np.zeros(len(ps_pos))
for i in range(len(ps_pos)):
    matrix[i]=tot_data[:,2][tot_data[:,-2]==ps_pos[i]]
    matrix_err[i]=tot_data[:,2][tot_data[:,-2]==ps_pos[i]]**0.5

ps_data=np.sum(matrix,axis=1)
P0=[(np.amax(ps_data)+np.amin(ps_data))/2, np.amax(ps_data)-np.amin(ps_data), 8,ps_pos[0]*8]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,ps_pos,ps_data, p0=P0, bounds=B0)
x_plt = np.linspace(ps_pos[0], ps_pos[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(ps_pos,ps_data,yerr=np.sqrt(ps_data),fmt="ko",capsize=5)  
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
ax.vlines(p[-1]/p[-2],0,fit_cos(p[-1]/p[-2], *p),ls="dashed")
w_ps=p[-2]
ps_0=p[-1]
print(w_ps)
print(ps_0)

c_data=np.sum(matrix,axis=0)
P0=[(np.amax(c_data)+np.amin(c_data))/2, np.amax(c_data)-np.amin(c_data), 3,0]
B0=([0,0,0,-10],[np.inf,np.inf,np.inf,np.inf])
p,cov=fit(fit_cos,coil,c_data, p0=P0, bounds=B0)
print(p)
print(np.diag(cov)**0.5)
x_plt = np.linspace(coil[0], coil[-1],100)
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
ax.errorbar(coil,c_data,yerr=np.sqrt(c_data),fmt="ko",capsize=5)  
ax.plot(x_plt,fit_cos(x_plt, *p), "b")
ax.vlines(p[-1]/p[-2],0,fit_cos(p[-1]/p[-2], *p),ls="dashed")
w_c=p[-2]
print("w_c=",w_c)
print("contrast gamma=",p[1]/p[0])
c_0=p[-1]
print(c_0)
# gamma=np.linspace(-3*np.pi,3*np.pi,500)#coil.copy()#
# chi=np.linspace(-3*np.pi,3*np.pi,500)#ps_pos.copy()#
gamma=w_c*coil-c_0
chi=w_ps*ps_pos-ps_0
# gamma=coil
# chi=ps_pos
alpha=np.pi/8
C=0.8
eta=1-C

# def fit_I_px(x, gamma0, chi0, w_c, w_ps, C, eta, A):
#     gamma = w_c*coil-gamma0
#     chi = w_ps*ps_pos-chi0
#     gamma, chi = np.meshgrid(gamma, chi)
#     fit_I_px = I_px_co(0, chi, C, alpha, gamma) + eta/4*A + I_px_in(0, chi, eta, alpha, gamma)
#     # print(fit_I_px)
#     return fit_I_px.ravel()

def fit_I_px(x, gamma0, chi0, w_c, w_ps, C, eta, A):
    gamma = w_c*coil-gamma0
    chi = w_ps*ps_pos-chi0
    # gamma, chi = np.meshgrid(gamma, chi)
    beta=exp_w1p(chi)
    for i in range(len(chi)):
        if i==0:
            fit_I_px = I_px_co(beta[i], chi[i], C, alpha, gamma) + eta/4*A + I_px_in(beta[i], chi[i], eta, alpha, gamma)
        else:
            fit_I_px = np.vstack((fit_I_px, I_px_co(beta[i], chi[i], C, alpha, gamma) + eta/4*A+ I_px_in(beta[i], chi[i], eta, alpha, gamma)))
    # print(fit_I_px)
    return fit_I_px.ravel()


P0=(c_0, ps_0, w_c, w_ps, 1, 0.3, 0.5)
B0=([-5,700,0.01,6,0,0,0],[5,1000,4,10,np.inf,np.inf,5])
p, cov = fit(fit_I_px, range(len(matrix.ravel())), matrix.ravel()/np.amax(matrix.ravel()), bounds=B0)
print(p)
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(111)
ax.plot(fit_I_px(0, *p)*np.amax(matrix.ravel()), "b")
# ax.plot(matrix.ravel()/np.amax(matrix.ravel()), "r--", lw=1)
ax.errorbar(np.arange(len(matrix.ravel())),matrix.ravel(), yerr=matrix_err.ravel(), fmt="ko", capsize=3, lw=1)
ax.set_xlim([50,150])
f_obs=matrix.ravel()
f_exp=fit_I_px(0,*p)*np.amax(matrix.ravel())

# f_obs/=np.sum(f_obs)
# f_exp/=np.sum(f_exp)

# print((np.sum(f_obs[:])-np.sum(f_exp[:]))/(np.sum(f_obs[:])))
print(chisquare(f_obs=f_obs, f_exp=f_exp, ddof=7))
# ax.set_xlim([0,150])
# ax.set_ylim([0,0.2])

def I_px(x, gamma0, chi0, w_c, w_ps, C, eta, A):
    gamma = w_c*coil-gamma0
    chi = w_ps*ps_pos-chi0
    # gamma, chi = np.meshgrid(gamma, chi)
    beta=exp_w1p(chi)
    for i in range(len(chi)):
        if i==0:
            fit_I_px = I_px_co(beta[i], chi[i], C, alpha, gamma)
        else:
            fit_I_px = np.vstack((fit_I_px, I_px_co(beta[i], chi[i], C, alpha, gamma) + eta/4*A + I_px_in(beta[i], chi[i], eta, alpha, gamma)))
    # print(fit_I_px)
    return fit_I_px

def I_px_corr(x, gamma0, chi0, w_c, w_ps, C, eta, A):
    gamma = w_c*coil-gamma0
    chi = w_ps*ps_pos-chi0
    # gamma, chi = np.meshgrid(gamma, chi)
    beta=exp_w1p(chi)
    for i in range(len(chi)):
        if i==0:
            fit_I_px = I_px_co(beta[i], chi[i], C, alpha, gamma)
        else:
            fit_I_px = np.vstack((fit_I_px, I_px_co(beta[i], chi[i], C, alpha, gamma)))
    return fit_I_px


fig = plt.figure(figsize=(10, 10))
ax = plt.axes(projection='3d')
gamma, chi = np.meshgrid(gamma, chi)
Z = matrix
Z1 = I_px(0, *p)*np.amax(matrix)
Z2 = I_px_corr(0, *p)*np.amax(matrix)
ax.contour3D(gamma, chi, Z, 40, cmap='binary')
ax.contour3D(gamma, chi, Z1, 40, cmap='plasma')  # cmap='Blues')
ax.set_xlabel('$\\gamma$')
ax.set_ylabel('$\chi$')
ax.set_zlabel('z')
ax.view_init(40, 45)
plt.show()

corrected_matrix=Z-(Z1-Z2)#Z2#
corrected_matrix_err=matrix_err
for i in range(len(ps_pos)):
    data_txt=np.array([coil, corrected_matrix[i], corrected_matrix_err[i], np.ones(len(coil))*ps_pos[i]])
    with open(correct_fold_path+"/gamma_ps_"+str("%02d" % (i,))+".txt", 'w') as f:
        np.savetxt(f, np.transpose(data_txt),  header= "Coil O-Beam err ps_pos", fmt='%.7f %.7f %.7f %.7f' )
