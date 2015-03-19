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
    
    key = 'SWIR'
      
    int_time = 12.0 #enter intergration time in sec for one "exposure" = NDR
    obs_time_hours = 6.0 #enter total observing time in hours
    obs_time = obs_time_hours*60*60.0
    time_osf = 100
 # no of osr samples per integration time
    rms = 4.0e-4
     
    #Call Jitter function
     
    ra_jitter, dec_jitter, time = exosim.lib.exolib.jitter(obs_time,int_time,time_osf,rms,mode=2)
    
    focal_length = opt.channel[key]['wfno'].val
    pixel_size = opt.channel[key]['pixel_size'].val
    total_osf = channel[key].osf*channel[key].ad_osf
    conv_fpa = channel[key].conv_fpa
    fpn = np.fromstring(opt.channel[key]['array_geometry'].val, 
	  sep=' ', dtype=np.float64)
   
    
    pix_count = conv_fpa[int(total_osf/2)::total_osf, int(total_osf/2)::total_osf]*1
            
    plate_scale =  206265/focal_length #arcsec per metre
    plate_scale /= 1e6 #arcsec per micron
    plate_scale /= 60*60 #degrees per micron
                
    ra_jitter /= plate_scale  # convert jitter from degrees to microns
    dec_jitter /= plate_scale

    print "std ra_jitter in microns", np.std(ra_jitter)
    print "std dec_jitter in microns", np.std(dec_jitter)
        
    ra_jitter *= total_osf/pixel_size  # convert jitter to no of osf 'units'
    dec_jitter *= total_osf/pixel_size 
    
    print "maximum jitter in osf units", abs(ra_jitter).max(), abs(dec_jitter).max()
    print "minimum jitter in osf units", abs(ra_jitter).min(), abs(dec_jitter).min()

    
    pad_y = total_osf*(np.round(ra_jitter.max()/total_osf))*2
    pad_x= total_osf*(np.round(dec_jitter.max()/total_osf))*2
        
#    conv_fpa0 = conv_fpa
    XX = pad_x*2 + conv_fpa.shape[1]  
    YY = pad_y*2 + conv_fpa.shape[0] 
    
    y_pixels = fpn[0]
    x_pixels = fpn[1]
    
    pad = np.zeros((YY,XX))    
    pad[pad_y:pad_y+conv_fpa.shape[0],pad_x:pad_x+conv_fpa.shape[1]] = conv_fpa
    
    conv_fpa = pad
    
    pix_count = conv_fpa[int(total_osf/2)::total_osf, int(total_osf/2)::total_osf]*1
           
    pix_count = np.zeros((y_pixels,x_pixels))
                                          
    N = 100
    pca = np.zeros((pix_count.shape[0],pix_count.shape[1],N))
    j = 0
    ct = 0
        
    
    for i in range(0,N*time_osf):
        ct += 1
        offset_y = ra_jitter[i]
        offset_x = dec_jitter[i]
        

        count = conv_fpa[int(total_osf/2)+pad_y+offset_y: int(total_osf/2)+pad_y + offset_y+ y_pixels*total_osf: total_osf, 
                int(total_osf/2)+pad_x+offset_x: int(total_osf/2)+pad_x + offset_x+ x_pixels*total_osf: total_osf]
         
#        print count.sum()                                  
        pix_count += count                                   
        if  count.sum() > 2000:    
       
             print count.sum(),i,"!!!!!"
    
                                   
        if ct == time_osf:
            ct = 0
            pca[...,j] = pix_count
            j += 1
            pix_count *= 0
            
        
    

    return pca
     
    plt.figure(1)
    plt.clf()
    plt.subplot(2,1,1)
    plt.plot(time,ra_jitter,'ro',markersize=0.6)
    plt.subplot(2,1,2)
    plt.plot(time,dec_jitter,'ro',markersize=0.6)





    
    
    
