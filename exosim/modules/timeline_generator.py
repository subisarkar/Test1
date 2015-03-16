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
    
    pad_y = osf*(np.round(ra_jitter.max()/osr))*2
    pad_x= osf*(np.round(dec_jitter.max()/osr))*2
        
#    conv_fpa0 = conv_fpa
    XX = pad_x*2 + conv_fpa.shape[1]  
    YY = pad_y*2 + conv_fpa.shape[0] 
    
    y_pixels = fpa.shape[0]/channel['SWIR'].osf
    x_pixels = fpa.shape[1]/channel['SWIR'].osf
    
    pad = np.zeros((YY,XX))    
    pad[pad_y:pad_y+conv_fpa.shape[0],pad_x:pad_x+conv_fpa.shape[1]] = conv_fpa
    
    conv_fpa = pad
    
    
        
    pix_count = conv_fpa[int(osf/2)+pad_y: int(osf/2)+pad_y + y_pixels*osf: osf, 
                         int(osf/2)+pad_x: int(osf/2)+pad_x + x_pixels*osf: osf]
                         
                             
    plt.figure(999)
    plt.imshow(pix_count)
    
    pix_count *= 0
    N = 100
    pca = np.zeros((pix_count.shape[0],pix_count.shape[1],N))
    j = 0
    ct = 0
        
    
    for i in range(0,osf*N):
        ct += 1
        offset_y = ra_jitter[i]
        offset_x = dec_jitter[i]
        
#        print offset_y, offset_x
        
#        print int(osf/2)+pad_y, conv_fpa.shape[0]+pad_y, (conv_fpa.shape[0]+pad_y) - (int(osf/2)+pad_y)
#        
#        print int(osf/2)+pad_y+offset_y, conv_fpa.shape[0]+pad_y+offset_y, (conv_fpa.shape[0]+pad_y+offset_y) - (int(osf/2)+pad_y+offset_y)
        
#        print i, pix_count.shape, conv_fpa[int(osf/2)+pad_y+offset_y: int(osf/2)+pad_y + offset_y+ y_pixels*osf: osf, 
#                                           int(osf/2)+pad_x+offset_x: int(osf/2)+pad_x + offset_x+ x_pixels*osf: osf].shape
        pix_count +=              conv_fpa[int(osf/2)+pad_y+offset_y: int(osf/2)+pad_y + offset_y+ y_pixels*osf: osf, 
                                           int(osf/2)+pad_x+offset_x: int(osf/2)+pad_x + offset_x+ x_pixels*osf: osf]

        if ct == osf:
            ct = 0
            pca[...,j] = pix_count
            j += 1
            pix_count *= 0
#            print np.sum(pix_count)
    
    

#        
#    print "mean count", np.mean(jitter_list)
#    print "std count", np.std(jitter_list)

#    conv_fpa = conv_fpa0*1    

    return pca
     
    plt.figure(1)
    plt.clf()
    plt.subplot(2,1,1)
    plt.plot(time,ra_jitter,'ro',markersize=0.6)
    plt.subplot(2,1,2)
    plt.plot(time,dec_jitter,'ro',markersize=0.6)





    
    
    
