# -*- coding: utf-8 -*-
"""
Created on Tue May 15 19:02:44 2018

@author: c1341133

signal

Returns either signal and noise for OOT simulation or
transit depth for IT simulation with noise estimate
"""

import numpy as np
import pytransit
import multiprocessing as mp
from scipy.optimize import fmin
import matplotlib.pyplot as plt


def getSNR(spectra):
    
    if len(spectra.shape) >1:
        signal = np.zeros( spectra.shape[1])
        noise = np.zeros( spectra.shape[1])
        for i in range ( spectra.shape[1]):
            signal[i] =  spectra[:,i].mean()
            noise[i]  = np.std( spectra[:,i])
    else:
        # fix for photometric channel
        signal = spectra.mean()
        noise =  spectra.std()
        
    return signal, noise

def getSpectrum(spectra, info):
    
    z = info['Z']
    multiaccum = info['MACCUM']
    nExp = info['NEXP']
    wl = info['WL']
    signal = np.zeros(spectra.shape[1])
    noise = np.zeros(spectra.shape[1])

#    # if using multiprocessing"
#                            
#    processes=[]
#    work_queue = mp.Queue()
#    done_queue = mp.Queue()
#    wl0 = np.zeros((spectra.shape[1])) 
#    print spectra.shape[1]
#    ct=0
#    for j in range (spectra.shape[1]):
#            print j
#            p = mp.Process(target = curvefit,
#                			 args=(spectra[j], z, multiaccum, nExp, j, wl[j], 
#                                      work_queue))
#            p.start()
#            processes.append(p)
#            p.join()
#    ct = 0        
#    for p in processes:
#        print ct
#        signal[ct], noise[ct], wl0[ct] = work_queue.get()
#        print signal[ct]
#        ct+=1
#    done_queue.put('STOP') 
#    signal = signal[np.argsort(wl0)][::-1]
#    noise = noise[np.argsort(wl0)][::-1]
##    wl0 = wl0[np.argsort(wl0)][::-1]


 
    for j in range (0,spectra.shape[1]):        
        
            signal[j], noise[j], _ = curvefit2(spectra[:,j], z, multiaccum, nExp, j, wl[j])

            print signal[j]
   
                
    return signal, noise 
    
