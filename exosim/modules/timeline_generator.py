#Time line generator module
"""
Created on Wed Mar 11 12:51:06 2015

@author: Subi
"""

import numpy as np
import scipy
import exosim
import matplotlib.pyplot as plt

#set observation parameters

def run():
 
    int_time = 12.0 #enter intergration time in sec for one "exposure" = NDR
    obs_time_hours = 6.0 #enter total observing time in hours
    obs_time = obs_time_hours*60*60.0
    osr = 100
    rms = 4.0e-12
     
    #Call Jitter function
     
    ra_jitter, dec_jitter, time = exosim.lib.exolib.jitter(obs_time,int_time,osr,rms,mode=2)
    
    print "std ra_jitter", np.std(ra_jitter)
    print "std dec_jitter", np.std(dec_jitter)
    
     
    #plot figures
     
    plt.figure(1)
    plt.clf()
    plt.subplot(2,1,1)
    plt.plot(time,ra_jitter,'ro',markersize=0.6)
    plt.subplot(2,1,2)
    plt.plot(time,dec_jitter,'ro',markersize=0.6)





    
    
    
