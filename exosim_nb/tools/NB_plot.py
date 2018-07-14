# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 18:08:52 2018

@author: c1341133
"""

import numpy as np
import matplotlib.pyplot as plt

import os

import glob




ch_list= [ 'AIRS CH1', 'NIR Spec', 'NIR Phot', 'FGS Red', 'FGS Prime']
#ch_list= ['FGS Red', 'FGS Prime']
 
 
pl = 'GJ 1214 b'
#pl = 'HD 209458 b'
#pl = 'HD 219134 b'
pl = 'Fake'  


No_real = 2
len_jitterless = 7
len_jitterful = 4

noise_group = 'jitterless'
#noise_group = 'jitterful'
ch = ch_list[0]
               
for filename in glob.iglob('/Volumes/My Passport for Mac/ExoSimOutput/%s_%s--%s*.npy'%(pl, ch, noise_group)):
    f = filename
    
indices = [i for i, s in enumerate(f) if 'x' in s]
idx1 = indices[-1]
indices = [i for i, s in enumerate(f) if '_' in s]
idx0 = indices[1]; idx2 = indices[2]
len_group = int(f[idx0+1:idx1])
No_real = int(f[idx1+1:idx2])

print len_group, No_real

no  = np.load(f)

if ch == 'FGS Red' or ch == 'FGS Prime' or ch == 'NIR Phot':
    wl = no[0]
    sig = no[1]
    aa = no.shape[0]
else:
    wl = no[:,0]
    sig = no[:,1]
    aa = no.shape[1]
 
if noise_group == 'jitterful':
    bb = len_group 
    no_type = ['Total noise', 'Spatial jitter', 'Spectral jitter', 'Combined jitter']

if noise_group == 'jitterless':
    bb = len_group-1 
    no_type = ['All shot noise', 'Source', 'Dark current', 'Zodi', 'Emission', 'Read out']

    




#plt.plot(wl,sig)

mean = np.zeros((len(wl), 2+bb))
std = np.zeros((len(wl), 2+bb))
mean[:,0] = wl
mean[:,1] = sig

for j in range(bb):
    
    no_stack =[]
    for i in range(2+j, aa, bb):
        
        no_stack.append(no[:,i].tolist())
        
    no_stack = np.array(no_stack)
    no_mean = no_stack.mean(axis =0) 
    no_std = no_stack.std(axis =0) 
    
    mean[:,j+2] = no_mean
    std[:,j+2] = no_std

for i in range(1,bb+2):   
    plt.semilogy(wl, mean[:,i], 'o')   
plt.grid()  
 
#plt.plot(wl,mean[:,7])    