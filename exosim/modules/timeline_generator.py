#Time line generator module
"""
Created on Wed Mar 11 12:51:06 2015

@author: Subi
"""

import numpy as np
from scipy import interpolate
from ..lib import exolib
import matplotlib.pyplot as plt

#set observation parameters

def run(opt,channel):
                
    ndr_number = 13 # the number of ndrs per exposure
    exp_number = 10 # the total number of exposures
    ndr_time = 7.0 #enter intergration time in sec for one "exposure" = NDR
    exp_time = ndr_time *ndr_number
    ndr_osf = 100  # no of osr samples per ndr

    N = exp_number * ndr_number  # number of ndrs in total
    
    obs_time_hours = 6.0 #enter total observing time in hours
    obs_time = obs_time_hours*60*60.0
    rms = 4.0e-5
    
    QE = 0.6
    QEsd = 0.05
     
    #Call Jitter function
     
    ra_jitter, dec_jitter, time = exolib.jitter(opt,obs_time,ndr_time,ndr_osf,rms,mode=2)
        
        
    for key in opt.channel.keys():
        
        osf = opt.channel[key]['osf'].val
        ad_osf = opt.channel[key]['ad_osf'].val
        ad_osf = 3
        new_osf =  ad_osf * osf
               
        print""
        print"Jittered timeline in %s channel being created"%(key)
        
        focal_length = opt.channel[key]['wfno'].val
        pixel_size = opt.channel[key]['pixel_size'].val
        osf = channel[key].osf 
        fpn = np.fromstring(opt.channel[key]['array_geometry'].val, 
    	  sep=' ', dtype=np.float64)
        fp = channel[key].fp
    
        plate_scale =  206265/ focal_length #arcsec per metre
        
        jitter_x = (ra_jitter/((plate_scale/1e6)/3600))/pixel_size
        jitter_y = (dec_jitter/((plate_scale/1e6)/3600))/pixel_size
    
        # Apply quantum efficiency variations 
       
        QE_array = np.random.normal(QE,QEsd*QE,fpn)
        QE_array = np.repeat(QE_array,osf,axis=0)
        QE_array = np.repeat(QE_array,osf,axis=1)  
       
        fp = fp*QE_array
     
        # Apply additional oversampling
     


        fp_count = fp[int(osf/2) :int(osf/2) + fp.shape[0]:osf, \
                     int(osf/2) :int(osf/2) + fp.shape[1]:osf]
                     
        
        xin = np.linspace(0,fp.shape[1]-1,fp.shape[1])
        yin = np.linspace(0,fp.shape[0]-1,fp.shape[0])
                
        x_step =  abs(xin[1]) - abs(xin[0])
        y_step =  abs(yin[1]) - abs(yin[0])
                
        x_step = np.float(x_step/ad_osf)
        y_step = np.float(y_step/ad_osf)
        
        xout = np.arange(0,fp.shape[1],x_step)
        yout = np.arange(0,fp.shape[0],y_step)


        fn = interpolate.RectBivariateSpline(yin,xin, fp)

        new_fp = fn(yout,xout)

        start = ad_osf

        new_fp_count = new_fp[start: start + fpn[0]*new_osf :new_osf, \
                              start: start + fpn[1]*new_osf :new_osf]


        print "fp", fp_count.sum(), fp.shape
        print "new _fp", new_fp_count.sum(), new_fp.shape, new_fp_count.shape
        
           
        # convert jitter from pixel units to new_osf units   
        
        
        jitter_x = jitter_x*new_osf
        jitter_y = jitter_y*new_osf
                
        # Apply a zeropad to allow for jitter beyond the size of the fp array    
    
        zeropad_x = round(jitter_x.max())*2
        zeropad_y = round(jitter_y.max())*2
        
        pad = np.zeros((zeropad_y*2 + new_fp.shape[0], zeropad_x*2 + new_fp.shape[1]))
    
        pad[zeropad_y:zeropad_y + new_fp.shape[0], zeropad_x:zeropad_x + new_fp.shape[1]] = new_fp
    
        new_fp = pad  
        
        start_y = zeropad_y + ad_osf
        start_x = zeropad_x + ad_osf
        
        new_fp_count = new_fp[start_y : start_y + fpn[0]*new_osf :new_osf, \
                              start_x : start_x + fpn[1]*new_osf :new_osf]
                
        print "new _fp", new_fp_count.sum(), new_fp.shape, new_fp_count.shape

    # Generate focal plane jitter
    #  - RPE simulated by ndr_osf loops where fp jitters and pixel count accumulates            
    #  - At end of each ndr, the accumulated count is stored as one 2d array in pca data cube
    #  - At end of each ndr, the RPE jitter standard deviation for each pixel 
    #         is stored in the pna data cube 
                
        pca = np.zeros((fpn[0],fpn[1],N))
                        
        for i in range(N):
            
            accum = np.zeros((fpn[0],fpn[1]))
    
            for j in range(ndr_osf):
                                
                ofx =  -jitter_x[i*ndr_osf +j]
                ofy =  -jitter_y[i*ndr_osf +j]
            
                start_x = ad_osf +zeropad_x + round(ofx)
                start_y = ad_osf +zeropad_y + round(ofy)
                            
                accum += new_fp[start_y : start_y + fpn[0]*new_osf :new_osf, \
                               start_x : start_x + fpn[1]*new_osf :new_osf]
                                                   
        
            pca[...,i] = accum            
    
        channel[key].timeline = pca 
        
        #pca is an array of NDRs
        
    return channel    

        



    
    
