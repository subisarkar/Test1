#Time line generator module
"""
Created on Wed Mar 11 12:51:06 2015

@author: Subi
"""

import numpy as np
from scipy import interpolate
from ..lib import exolib
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.nan)

#set observation parameters

def run(opt,channel):
    
    #  for WFC3 : ndr_number = 12, ndr_time = 7.3465 sec, exp_number = 67
                
    ndr_number = 12 # the number of ndrs per exposure
    exp_number = 67 # the total number of exposures
    ndr_time = 7.3465 #enter intergration time in sec for one "exposure" = NDR
    exp_time = ndr_time *ndr_number
    ndr_osf = 10  # no of osr samples per ndr

    N = exp_number * ndr_number  # number of ndrs in total
    N = 100
    
    obs_time_hours = 6.0 #enter total observing time in hours
    obs_time = obs_time_hours*60*60.0
    rms = 4.0e-5
    
 
     
    #Call Jitter function
    
    jitter_file = opt.common_exosym_path.val+"/data/instrument/herschel_long_pointing.fits"
    ra_jitter, dec_jitter, time = exolib.jitter(jitter_file,obs_time,ndr_time,ndr_osf,rms,mode=2)
        
        
    for key in opt.channel.keys():
        
        osf = opt.channel[key]['osf'].val
        ad_osf = opt.channel[key]['ad_osf'].val
        ad_osf = 7
        new_osf =  ad_osf * osf
               
        print""
        print"Jittered timeline in %s channel being created"%(key)
        
        pixel_size = opt.channel[key]['pixel_size'].val
        osf = channel[key].osf 
        fpn = np.fromstring(opt.channel[key]['array_geometry'].val, 
    	  sep=' ', dtype=np.float64)
        fp = channel[key].fp
                 
        plate_scale =opt.channel[key]['plate_scale'].val
        
        jitter_x = (ra_jitter/((plate_scale)/3600))/pixel_size
        jitter_y = (dec_jitter/((plate_scale)/3600))/pixel_size
            

     
        # Apply additional oversampling
     

        fp_count = fp[int(osf/2) :int(osf/2) + fp.shape[0]:osf, \
                     int(osf/2) :int(osf/2) + fp.shape[1]:osf]
                     
        
        xin = np.linspace(0,fp.shape[1]-1,fp.shape[1])
        yin = np.linspace(0,fp.shape[0]-1,fp.shape[0])
                        
        x_step =  abs(xin[1]) - abs(xin[0])
        y_step =  abs(yin[1]) - abs(yin[0])
                
        # calculates the equivalent step sizes for new grid and ad_osf        
        x_step_new = np.float(x_step/ad_osf)
        y_step_new = np.float(y_step/ad_osf)
        
        # new grid must start with an offset to produce correct number of new points
        x_start = -x_step_new * np.float((ad_osf-1)/2)
        y_start = -y_step_new * np.float((ad_osf-1)/2)
        
        # new grid points, correct start, end and spacing
        xout = np.arange(x_start, x_start + x_step_new*fp.shape[1]*ad_osf, x_step_new)
        yout = np.arange(y_start, y_start + y_step_new*fp.shape[0]*ad_osf, y_step_new)
        
        # interpolate fp onto new grid
        fn = interpolate.RectBivariateSpline(yin,xin, fp)
        new_fp = fn(yout,xout)
        
        # 'central subpixel' = the subpixel corresponding to the centre of the real pixel
        # count on central subpixel requires following offset in new fp
        start =(ad_osf-1)/2 + ad_osf
        
        new_fp_count = new_fp[start: start + fpn[0]*new_osf :new_osf, \
                              start: start + fpn[1]*new_osf :new_osf]


        # summed counts on central subpixels should be the same on fp and new_fp
        print "fp count and shape", fp_count.sum(), fp.shape
        print "new _fp count and shape", new_fp_count.sum(), new_fp.shape
        
                 
           
        # convert jitter from real pixel units to new subpixel units   
           
        jitter_x = jitter_x*new_osf
        jitter_y = jitter_y*new_osf
                
        # Apply a zeropad to allow for jitter beyond the size of the fp array    
    
        zeropad_x = round(jitter_x.max())*2
        zeropad_y = round(jitter_y.max())*2
        
        pad = np.zeros((zeropad_y*2 + new_fp.shape[0], zeropad_x*2 + new_fp.shape[1]))
    
        pad[zeropad_y:zeropad_y + new_fp.shape[0], zeropad_x:zeropad_x + new_fp.shape[1]] = new_fp
    
        new_fp = pad
        
        #check zeropadded new_fp produces the same count
        
        start_y = zeropad_y + (ad_osf-1)/2 + ad_osf
        start_x = zeropad_x + (ad_osf-1)/2 + ad_osf
        
        new_fp_count = new_fp[start_y : start_y + fpn[0]*new_osf :new_osf, \
                              start_x : start_x + fpn[1]*new_osf :new_osf]
                
        print "zero-padded new _fp count and shape", new_fp_count.sum(), new_fp.shape

    # Generate focal plane jitter
    #  - RPE simulated by ndr_osf loops where fp jitters and pixel count accumulates            
    #  - At end of each ndr, the accumulated count is stored as one 2d array in pca data cube

                
        pca = np.zeros((fpn[0],fpn[1],N))
                        
        for i in range(N):
            
            accum = np.zeros((fpn[0],fpn[1]))
    
            for j in range(ndr_osf):
                                
                ofx =  -jitter_x[i*ndr_osf +j]
                ofy =  -jitter_y[i*ndr_osf +j]
            
                start_x = zeropad_x + (ad_osf-1)/2 + ad_osf + round(ofx)
                start_y = zeropad_y + (ad_osf-1)/2 + ad_osf + round(ofy)
                            
                count = new_fp[start_y : start_y + fpn[0]*new_osf :new_osf, \
                                start_x : start_x + fpn[1]*new_osf :new_osf]
                accum += count * ndr_time/ndr_osf
                                                                           
            pca[...,i] = accum            

    channel[key].timeline = pca 
        
        #pca is an array of NDRs
    

    T = pca[...,0:12].sum(axis=2)
    C = T.sum(axis=0)
    
    print "NDR 1-12 total summed count in electrons", C.sum()

    plt.figure(111)
    plt.plot(C)
    
    return channel    

        



    
    
