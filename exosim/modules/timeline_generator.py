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
    
    for key in ['SWIR']:  #opt.channel.keys():
        
        focal_length = opt.channel[key]['wfno'].val
        pixel_size = opt.channel[key]['pixel_size'].val
        osf = channel[key].osf 
        fpn = np.fromstring(opt.channel[key]['array_geometry'].val, 
    	  sep=' ', dtype=np.float64)
        fp = channel[key].fp
    
             
        ndr_number = 13 # the number of ndrs per exposure
        exp_number = 10 # the total number of exposures
        ndr_time = 7.0 #enter intergration time in sec for one "exposure" = NDR
        exp_time = ndr_time *ndr_number
        ndr_osf = 100  # no of osr samples per ndr
    
        N = exp_number * ndr_number  # number of ndrs in total
        
        obs_time_hours = 6.0 #enter total observing time in hours
        obs_time = obs_time_hours*60*60.0
        rms = 4.0e-5
        
        ad_osf = 7
        QE = 0.6
        QEsd = 0.05
         
        #Call Jitter function
         
        ra_jitter, dec_jitter, time = exolib.jitter(obs_time,ndr_time,ndr_osf,rms,mode=2)
    
        plate_scale =  206265/ focal_length #arcsec per metre
        
        ra_jitter = (ra_jitter/((plate_scale/1e6)/3600))/pixel_size
        dec_jitter = (dec_jitter/((plate_scale/1e6)/3600))/pixel_size
    
        # Apply quantum efficiency variations 
       
        QE_array = np.random.normal(QE,QEsd*QE,fpn)
        QE_array = np.repeat(QE_array,osf,axis=0)
        QE_array = np.repeat(QE_array,osf,axis=1)  
       
        fp = fp*QE_array
     
        # Apply additional oversampling
       
        fp_count = fp[int(osf/2) :int(osf/2) + fpn[0]*osf :osf, \
                      int(osf/2) :int(osf/2) + fpn[1]*osf :osf]
            
        xin = np.linspace(-1,1,fp.shape[1])
        xout = np.linspace(-1,1,fp.shape[1]*ad_osf)
        yin = np.linspace(-1,1,fp.shape[0])
        yout = np.linspace(-1,1,fp.shape[0]*ad_osf)
    
        fn = interpolate.RectBivariateSpline(yin,xin,fp)
        
        new_fp = fn(yout, xout)
        new_osf = osf*ad_osf   
            
        count = new_fp[int(new_osf/2):int(new_osf/2)+fpn[0]*new_osf:new_osf, \
                       int(new_osf/2):int(new_osf/2) +fpn[1]*new_osf:new_osf]
        
        print""               
        print "fp_count, unmodified count", fp_count.sum(),count.sum() 
                       
        new_fp  = new_fp * fp_count.sum()/count.sum() 
        #this compensates for increased count caused by the oversampling               
        
        count = new_fp[int(new_osf/2):int(new_osf/2)+fpn[0]*new_osf:new_osf, \
                       int(new_osf/2):int(new_osf/2) +fpn[1]*new_osf:new_osf]  
        
        new_fp0 = new_fp*1
                      
        print "modified  count", count.sum()   
           
           
        # convert jitter from pixel units to new_osf units   
        
        
        jitter_x = ra_jitter*new_osf
        jitter_y = dec_jitter*new_osf
                
        # Apply a zeropad to allow for jitter beyond the size of the fp array    
    
        zeropad_x = round(jitter_x.max())*2
        zeropad_y = round(jitter_y.max())*2
        
        pad = np.zeros((zeropad_y*2 + new_fp.shape[0], zeropad_x*2 + new_fp.shape[1]))
    
        pad[zeropad_y:zeropad_y + new_fp.shape[0], zeropad_x:zeropad_x + new_fp.shape[1]] = new_fp
    
        new_fp = pad  
                
    # Generate focal plane jitter
    #  - RPE simulated by ndr_osf loops where fp jitters and pixel count accumulates            
    #  - At end of each ndr, the accumulated count is stored as one 2d array in pca data cube
    #  - At end of each ndr, the RPE jitter standard deviation for each pixel 
    #         is stored in the pna data cube 
                
        pca = np.zeros((fpn[0],fpn[1],N))
        pna = np.zeros((fpn[0],fpn[1],N))
        
        nda = np.zeros((fpn[0],fpn[1],ndr_osf))
        
        count = np.zeros((fpn[0],fpn[1]))
        
        for i in range(N):
            
            nda = np.zeros((fpn[0],fpn[1],ndr_osf))
    
            for j in range(ndr_osf):
                
                ofx = jitter_x[i*ndr_osf +j]
                ofy = jitter_y[i*ndr_osf +j]
            
                start_x = int(new_osf/2)+zeropad_x
                start_y = int(new_osf/2)+zeropad_y
                
                start_x = start_x + round(ofx)
                start_y = start_y + round(ofy)
            
                nda[...,j] = new_fp[start_y:start_y +fpn[0]*new_osf:new_osf, \
                                     start_x:start_x +fpn[1]*new_osf:new_osf] 
        
            pca[...,i] = np.sum(nda, axis =2, dtype=np.float64)
            pna[...,i] = np.std(nda, axis =2, dtype=np.float64)         
             
        return pca
        
            
            ###
            
    #        start_x = int(new_osf/2)
    #        start_y = int(new_osf/2)
    #        
    #        start_x = start_x + round(ofx)
    #        start_y = start_y + round(ofy)
    #    
    #        count0 = new_fp0[start_y:start_y +fpn[0]*new_osf:new_osf, \
    #                   start_x:start_x +fpn[1]*new_osf:new_osf]     
    #        
    #        print count.sum(), count0.sum(), ofx, ofy
            
            ####
            
    #        C.append(count.sum())
    #    print 'mean', np.mean(C), 'std', np.std(C), len(C)



    
    
