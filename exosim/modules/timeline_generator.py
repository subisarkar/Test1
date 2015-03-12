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

def run(opt,channel):
    
    
    
    int_time = 12.0 #enter intergration time in sec for one "exposure" = NDR
    obs_time_hours = 6.0 #enter total observing time in hours
    obs_time = obs_time_hours*60*60.0
    osr = 100 # no of osr samples per integration time
    rms = 4.0e-5
     
    #Call Jitter function
     
    ra_jitter, dec_jitter, time = exosim.lib.exolib.jitter(obs_time,int_time,osr,rms,mode=2)
    
    focal_length = opt.channel['SWIR']['wfno'].val
    pixel_size = opt.channel['SWIR']['pixel_size'].val
    osf = channel['SWIR'].osf*channel['SWIR'].ad_osf
    conv_fpa = channel['SWIR'].conv_fpa
    fpa = channel['SWIR'].fpa

    
    plate_scale =  206265/focal_length #arcsec per metre
    plate_scale /= 1e6 #arcsec per micron
    plate_scale /= 60*60 #degrees per micron
                
    ra_jitter /= plate_scale  # convert jitter from degrees to microns
    dec_jitter /= plate_scale

    print "std ra_jitter in microns", np.std(ra_jitter)
    print "std dec_jitter in microns", np.std(dec_jitter)
        
    ra_jitter *= osf/pixel_size  # convert jitter to no of osf 'units'
    dec_jitter *= osf/pixel_size  
    
    pix_count = conv_fpa[int(osf/2)::osf, int(osr/2)::osf]
    
    print pix_count.shape, fpa.shape
    print pix_count.sum(), fpa.sum()

    
    
    plt.figure(999)
    plt.imshow(pix_count)

    
    
    
    
    
    
     
    #plot figures
     
    plt.figure(1)
    plt.clf()
    plt.subplot(2,1,1)
    plt.plot(time,ra_jitter,'ro',markersize=0.6)
    plt.subplot(2,1,2)
    plt.plot(time,dec_jitter,'ro',markersize=0.6)





    
    
    
