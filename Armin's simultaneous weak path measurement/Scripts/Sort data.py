#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 10:03:06 2023

@author: aaa
"""
import os
import numpy as np
import shutil

TOF_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Armin's simultaneous weak path measurement/sc/TOF/"
for root, dirs, files in os.walk(TOF_fold_path, topdown=False):
    fold_names = dirs.copy()
    

for fold_name in fold_names:
            fold_path=TOF_fold_path+fold_name+"/"
            for root, dirs, files in os.walk(fold_path, topdown=False):
                for name in files:
                    if "TOF_PS_RF2_weak_tof_09Jun1333.inf" in name:
                        # print(name)
                        inf_file_name=name[:-4]
                        inf_path = fold_path+inf_file_name+".inf"
                        sorted_fold_path="/home/aaa/Desktop/Fisica/PhD/2023/Armin's simultaneous weak path measurement/Sorted data/"+inf_file_name
                        
                        if not os.path.exists(sorted_fold_path):
                            os.makedirs(sorted_fold_path)

                        rawdata=sorted_fold_path+"/Rawdata" 
                        if not os.path.exists(rawdata):
                            os.makedirs(rawdata)

                        cleandata=sorted_fold_path+"/Cleantxt"
                        if not os.path.exists(cleandata):
                            os.makedirs(cleandata)
                            

                        inf_files = np.genfromtxt(inf_path, dtype=str,usecols=(0,1))
                        ps_num=np.array([range(len(inf_files[:,0])),inf_files[:,1].astype(float)])
                        # print(inf_files)
                        # print(ps_num)
                        for root, dirs, files in os.walk(fold_path, topdown=False):
                            i=0
                            for name in files:
                                if (name in inf_files[:,0]):
                                    # print(name)
                                    # print(name)
                                    shutil.copy(os.path.join(root, name), rawdata+"/"+name)
                                    # with open(os.path.join(root, name)) as f:
                                    #     lines = f.readlines()
                                    # print(lines[0])
                                    txt=np.loadtxt(os.path.join(root, name), encoding='windows-1252', skiprows=9, comments="*", delimiter="\t")[:,1:]
                                    ps=inf_files[:,1][inf_files[:,0]==name].astype(float) #phashifter position
                                    # if len(txt[0])<4:
                                    #     print(name)
                                    #     cleantxt =np.hstack((np.ones((len(txt[:,0]),1))*txt[:,0],np.ones((len(txt[:,0]),1))*7200.028,txt[:,1:],np.ones((len(txt[:,0]),1))*ps[0]))
                                    # else:
                                    cleantxt =np.hstack((txt,np.ones((len(txt[:,0]),1))*ps[0]))
                                    with open(cleandata+"/beta_ps_"+str("%02d" % (ps_num[0][ps_num[1]==ps][0],)[0])+".txt", 'w') as f:
                                        if len(cleantxt[0])>3:    
                                            np.savetxt(f, cleantxt, header= "BinNR microsec O-Beam AUX chi(deg)", fmt='%.3f %i %i %.3f')