#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 10:03:06 2023

@author: aaa
"""
import os
import numpy as np
import shutil


fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 2nd round/exp_CRG-3033/rawdata/sc/"
inf_file_name="beta_alpi8off_gamma_chicoarse_BETA_30Jun1447"#"path1pi8_noIn_cb_g_17Apr1557"
inf_path = fold_path+inf_file_name+".inf"
# beta0_files=["moveDC2_08Apr1342.dat","moveDC2_08Apr1424.dat","moveDC2_08Apr1514.dat","moveDC2_08Apr1551.dat","moveDC2_08Apr1627.dat","moveDC2_08Apr1703.dat"]
# inf_files_alpha0=["S2Z+_17Jun1652","S2Z-_17Jun1600"]
# inf_path_alpha0 = [fold_path+inf_files_alpha0[0]+".inf",fold_path+inf_files_alpha0[1]+".inf"]
sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Grenoble 2nd round/exp_CRG-3033/Sorted data/"+inf_file_name
if not os.path.exists(sorted_fold_path):
    os.makedirs(sorted_fold_path)

rawdata=sorted_fold_path+"/Rawdata" 
if not os.path.exists(rawdata):
    os.makedirs(rawdata)

cleandata=sorted_fold_path+"/Cleantxt" 
if not os.path.exists(cleandata):
    os.makedirs(cleandata)
    
plots=sorted_fold_path+"/Plots" 
if not os.path.exists(plots):
    os.makedirs(plots)

beta_fold_raw=rawdata+"/Beta"
beta_fold_clean=cleandata+"/Beta"
if (not os.path.exists(beta_fold_raw)) or (not os.path.exists(beta_fold_clean)):
    os.makedirs(beta_fold_raw)
    os.makedirs(beta_fold_clean)
    
alphaON_fold_raw=beta_fold_raw+"/AlphaON"
alphaON_fold_clean=beta_fold_clean+"/AlphaON"
if (not os.path.exists(alphaON_fold_raw)) or (not os.path.exists(alphaON_fold_clean)):
    os.makedirs(alphaON_fold_raw)
    os.makedirs(alphaON_fold_clean)
    
alphaOFF_fold_raw=beta_fold_raw+"/AlphaOFF"
alphaOFF_fold_clean=beta_fold_clean+"/AlphaOFF"
if (not os.path.exists(alphaOFF_fold_raw)) or (not os.path.exists(alphaOFF_fold_clean)):
    os.makedirs(alphaOFF_fold_raw)
    os.makedirs(alphaOFF_fold_clean)

# gamma_fold_raw=rawdata+"/Gamma"
# gamma_fold_clean=cleandata+"/Gamma"
# if (not os.path.exists(gamma_fold_raw)) or (not os.path.exists(gamma_fold_clean)):
#     os.makedirs(gamma_fold_raw)
#     os.makedirs(gamma_fold_clean)
    
# alphaON_fold_raw=gamma_fold_raw+"/AlphaON"
# alphaON_fold_clean=gamma_fold_clean+"/AlphaON"
# if (not os.path.exists(alphaON_fold_raw)) or (not os.path.exists(alphaON_fold_clean)):
#     os.makedirs(alphaON_fold_raw)
#     os.makedirs(alphaON_fold_clean)
    
# alphaOFF_fold_raw=gamma_fold_raw+"/AlphaOFF"
# alphaOFF_fold_clean=gamma_fold_clean+"/AlphaOFF"
# if (not os.path.exists(alphaOFF_fold_raw)) or (not os.path.exists(alphaOFF_fold_clean)):
#     os.makedirs(alphaOFF_fold_raw)
#     os.makedirs(alphaOFF_fold_clean)

# beta0_fold_raw=rawdata+"/Beta0"
# beta0_fold_clean=cleandata+"/Beta0"
# if (not os.path.exists(beta0_fold_raw)) or (not os.path.exists(beta0_fold_clean)):
#     os.makedirs(beta0_fold_raw)
#     os.makedirs(beta0_fold_clean)

# alpha0_fold_raw=rawdata+"/Alpha0"
# alpha0_fold_clean=cleandata+"/Alpha0"
# if (not os.path.exists(alpha0_fold_raw)) or (not os.path.exists(alpha0_fold_clean)):
#     os.makedirs(alpha0_fold_raw)
#     os.makedirs(alpha0_fold_clean)

# gamma0_fold_raw=rawdata+"/Gamma0"
# gamma0_fold_clean=cleandata+"/Gamma0"
# if (not os.path.exists(gamma0_fold_raw)) or (not os.path.exists(gamma0_fold_clean)):
#     os.makedirs(gamma0_fold_raw)
#     os.makedirs(gamma0_fold_clean)

# ifg_fold_raw=rawdata+"/IFG"
# ifg_fold_clean=cleandata+"/IFG"
# if (not os.path.exists(ifg_fold_raw)) or (not os.path.exists(ifg_fold_clean)):
#     os.makedirs(ifg_fold_raw)
#     os.makedirs(ifg_fold_clean)


inf_files = np.genfromtxt(inf_path, dtype=str,usecols=(0,2))
# coil_amp= np.genfromtxt(inf_path, dtype=str,usecols=(0,4), skip_header=1)
ps_num=np.array([range(len(inf_files[::2,0])),inf_files[::2,1].astype(float)])
# ps_num=np.array([range(len(inf_files[:,0])),inf_files[:,1].astype(float)])
print(ps_num)
# inf_files_beta0 = np.genfromtxt(inf_path, dtype=str,usecols=(0,1))
# inf_files_alpha0 = np.genfromtxt(inf_path_alpha0[0], dtype=str,usecols=(0))
# inf_files_alpha0 = np.append(inf_files_alpha0,np.genfromtxt(inf_path_alpha0[1], dtype=str,usecols=(0)))

# inf_files_ifg = np.genfromtxt(inf_path, dtype=str,usecols=(0))

for root, dirs, files in os.walk(fold_path, topdown=False):

    i=0
    for name in files:
        """
        Beta scan sorting
        """
        if (name in inf_files[:,0]):
            print(name)
            if "ON" in name:
                # print(name)
                shutil.copy(os.path.join(root, name), alphaON_fold_raw+"/"+name)
                txt=np.loadtxt(os.path.join(root, name), comments="*", delimiter="\t")[:,1:]
                ps=inf_files[:,1][inf_files[:,0]==name].astype(float) #phashifter position
                cleantxt =np.hstack((txt,np.ones((len(txt[:,0]),1))*ps[0]))
                with open(alphaON_fold_clean+"/beta_ps_"+str("%02d" % (ps_num[0][ps_num[1]==ps][0],)[0])+".txt", 'w') as f:
                    np.savetxt(f, cleantxt,  header= "Position(mm) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f %.7f' )

            if "OFF" in name:
                # print(name)
                shutil.copy(os.path.join(root, name), alphaOFF_fold_raw+"/"+name)
                txt=np.loadtxt(os.path.join(root, name), comments="*", delimiter="\t")[:,1:]
                ps=inf_files[:,1][inf_files[:,0]==name].astype(float) #phashifter position
                cleantxt =np.hstack((txt,np.ones((len(txt[:,0]),1))*ps[0]))
                with open(alphaOFF_fold_clean+"/beta_ps_"+str("%02d" % (ps_num[0][ps_num[1]==ps][0],)[0])+".txt", 'w') as f:
                    np.savetxt(f, cleantxt,  header= "Position(mm) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f %.7f' )
        
        """
        Gamma scan sorting
        """
        
        # if (name in inf_files[:,0]):
        #     print(name)
        #     shutil.copy(os.path.join(root, name), gamma_fold_raw+"/"+name)
        #     txt=np.loadtxt(os.path.join(root, name), skiprows=6, comments="*", delimiter="\t")[:,1:]
        #     ps=inf_files[:,1][inf_files[:,0]==name].astype(float) #phashifter position
        #     print(ps)
        #     cleantxt =np.hstack((txt,np.ones((len(txt[:,0]),1))*ps[0]))
        #     with open(gamma_fold_clean+"/gamma_ps_"+str("%02d" % (ps_num[0][ps_num[1]==ps][0],)[0])+".txt", 'w') as f:
        #         np.savetxt(f, cleantxt,  header= "Position(mm) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f %.7f' )

        
        #     if name[:4] == "DC2Y":
        #         shutil.copy(os.path.join(root, name), gamma_fold_raw+"/"+name)
        #         txt=np.loadtxt(os.path.join(root, name), skiprows=6, comments="*", delimiter="\t")[:,1:]
        #         ps=inf_files[:,1][inf_files[:,0]==name].astype(float)
        #         cleantxt =np.hstack((txt,np.ones((len(txt[:,0]),1))*ps[0]))
        #         c_amp=coil_amp[:,1][coil_amp[:,0]==name].astype(float)
        #         cleantxt =np.hstack((cleantxt,np.ones((len(txt[:,0]),1))*c_amp[0]))
        #         with open(gamma_fold_clean+"/gamma_ps_"+str("%02d" % (ps_num[0][ps_num[1]==ps][0],)[0])+".txt", 'w') as f:
        #                 np.savetxt(f, cleantxt,  header= "Current(A) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos coil_amp(A)", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f %.7f %.7f' )
        
        """
        \beta0 sorting
        """
        # if (name in inf_files_beta0[:,0]):
        #     # print(name)
        #     shutil.copy(os.path.join(root, name), beta0_fold_raw+"/"+name)
        #     # print(np.loadtxt(os.path.join(root, name), skiprows=22, comments="*", delimiter="\t")[0])
        #     cleantxt=np.loadtxt(os.path.join(root, name), skiprows=22, comments="*", delimiter="\t")[:,1:]
        #     with open(beta0_fold_clean+"/"+name, 'w') as f:
        #         np.savetxt(f, cleantxt,  header= "Position(mm) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f' )
        """
        \alpha0 scan sorting
        """  
        # if (name in inf_files_alpha0):
        #     shutil.copy(os.path.join(root, name), alpha0_fold_raw+"/"+name)
        #     cleantxt=np.loadtxt(os.path.join(root, name), skiprows=25, comments="*", delimiter="\t")[:,1:]
        #     with open(alpha0_fold_clean+"/"+name, 'w') as f:
        #         np.savetxt(f, cleantxt,  header= "Current(A) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f' )
        """
        \gamma0 scan sorting
        """  
        # if (name[:len(inf_file_name)] == inf_file_name) and name[-4:]==".dat":
        #     shutil.copy(os.path.join(root, name), gamma0_fold_raw+"/"+name)
        #     cleantxt=np.loadtxt(os.path.join(root, name), skiprows=25, comments="*", delimiter="\t")[:,1:]
        #     with open(gamma0_fold_clean+"/"+name, 'w') as f:
        #         np.savetxt(f, cleantxt,  header= "Current(A) exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f' )
        
        """
        \ifg scan sorting
        """  
        # if name in inf_files_ifg:
        #     shutil.copy(os.path.join(root, name), ifg_fold_raw+"/"+name)
        #     cleantxt=np.loadtxt(os.path.join(root, name), skiprows=20, comments="*", delimiter="\t")[:,1:]
        #     with open(ifg_fold_clean+"/"+name, 'w') as f:
        #         np.savetxt(f, cleantxt,  header= "ps_pos exposure_time(s) O-Beam H-Beam Monitor H2-Beam time(s) O+H encod1 encod2 phaseshif_pos", fmt='%.7f %.1f %i %i %i %i %i %i %.7f %.7f' )