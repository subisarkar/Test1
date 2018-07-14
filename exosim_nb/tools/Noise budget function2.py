# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 16:54:18 2015

@author: c1341133
"""

import scipy.ndimage as nd
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import pyfits 
import exosim_nb as exosim
from exosim_nb.PointingDecorSubi import ExosimPointingDecorr
import os



use_hard_drive = 1
apply_cds = 1

decorr_make = True
#decorr_make = False

apply_decorr = True

decorr_method = 'xcorr'
decorr_method = 'pointing'
decorr_method = 'fft'
 

#ch_list= ['AIRS CH0', 'AIRS CH1', 'NIR Spec', 'NIR Phot', 'FGS Red', 'FGS Prime']
ch_list= ['FGS Red', 'FGS Prime']
 
 
pl = 'GJ 1214 b'
#pl = 'HD 209458 b'
#pl = 'HD 219134 b'
pl = 'Fake'  


No_real = 2
len_jitterless = 7
len_jitterful = 4
 
qe = np.load('/Users/c1341133/Desktop/ExoSim-stellar/exosim_s/data/qe_rms.npy') 
                   
for noise_group in ['jitterless','jitterful'] :
    
    if noise_group == 'jitterless':
        len_group = len_jitterless
        decorr = False
    elif noise_group == 'jitterful':
        len_group = len_jitterful
        if apply_decorr == True:
            decorr = True
        else:
            decorr = False
  
    for ch in ch_list:
        
        ct = -1
        noise_results = 0
        pixSize = 18.0
        if ch == 'AIRS CH0' or ch == 'AIRS CH1':
            pixSize = 15.0
        print "Channel...  ", ch    
        print "Pixel size...  ", pixSize
    
        for j in range(No_real):
                        
            for i in range(0,len_group):
                
                ct+=1
                
                sim_no = "%04d" % (j*len_group + i)
                print "Sim file number...  ", sim_no
                                
                if use_hard_drive == 1:
                    sim_file ="/Volumes/My Passport for Mac/ExoSimOutput/%s_%s--%s/sim_%s/%s_signal.fits"%(pl, ch, noise_group, sim_no, ch) 
                    txt_file ="/Volumes/My Passport for Mac/ExoSimOutput/%s_%s--%s/sim_%s/%s_info.txt"%(pl, ch, noise_group, sim_no, ch)  
                else:
                    sim_file ="~/ExoSimOutput/%s_%s--%s/sim_%s/%s_signal.fits"%(pl, ch, noise_group, sim_no, ch) 
                    txt_file ="~/ExoSimOutput/%s_%s--%s/sim_%s/%s_info.txt"%(pl, ch, noise_group, sim_no, ch)  


                if 'Flat field and bkg subtraction applied:  True' in open(txt_file).read():
                    print("File has been flat-fielded in ExoSim already...")
                    apply_flat = 0
                else:
                    apply_flat = 1
                    
                if 'QE grid applied:  False' in open(txt_file).read():
                    print("QE variations were not applied in ExoSim... no flat fielding needed")
                    apply_flat = 0
                    
                diff = 0
                if 'Diffuse:  True' in open(txt_file).read():
                    diff =1
                if 'Noise option:  3' in open(txt_file).read():
                    diff =1
                if 'Noise option:  4' in open(txt_file).read():
                    diff =1
                if 'Noise option:  5' in open(txt_file).read():
                    diff =1
                if 'Noise option:  6' in open(txt_file).read():
                    diff =1                                        
                if diff ==1:
                    print("File has diffuse source ...")
                    
                if decorr == True:
                    if ch == 'AIRS CH0' or ch == 'AIRS CH1' or ch == 'NIR Spec':
                        decorr_file = "/Volumes/My Passport for Mac/ExoSimOutput/%s_%s--%s/sim_%s/%s_signal_decorr.fits"%(pl, ch, noise_group, sim_no, ch) 
                        if os.path.isfile(decorr_file):
                            print "Decorr file already exists..."
                            sim_file = decorr_file
                            apply_flat = 0
                            apply_cds = 0
                        else:
                            print "Making decorr file..."
                            if 'JN(Spectral):  False'  in open(txt_file).read():
                                p_spec = False
                            else:
                                p_spec = True
                            if 'JN(Spatial):  False' in open(txt_file).read():
                                p_spat = False
                            else:
                                p_spat = True                                                     
                            if ch == 'AIRS CH0':
                                pscale = 6.11e-5
                            elif ch == 'AIRS CH1':
                                pscale = 1.23e-4
                            elif ch == 'NIR Spec':
                                pscale = 3.01e-5  # use for baseline
    #                            pscale = 6.11e-5
                            elif ch == 'FGS Red':
                                pscale = 2.23e-5
                            elif ch == 'FGS Prime':
                                pscale = 3.10e-5
                            elif ch == 'NIR Phot':
                                pscale = 2.43e-5                            
                                
                            ps_info = [p_spec, p_spat, pscale]
                            epd = ExosimPointingDecorr(sim_file, qe, ch, decorr_method, apply_flat, diff, ps_info)
                            hdu = epd.getHdu()
                            sim_file = decorr_file
                            apply_flat = 0
                            apply_cds = 0
                    elif ch == 'FGS Red'or ch == 'FGS Prime'or ch == 'NIR Phot':
                        print "Photometric channel... will decorrelate by aperture method"
                        
      
                if ch == 'AIRS CH0':
                    R = 100.
                    F = 13.2
                    wavrange = [1.8,4.2]
                    
                if ch == 'AIRS CH1':
                    R = 30.
                    F = 6.36
                    wavrange = [3.80,8.0]
  
                if ch == 'NIR Spec':
                    R = 20.0                   
                    F = 38.68
                    wavrange = [0.9,2.2]

                if ch == 'FGS Red':
                    Ap = 5.0
                    if pl == 'HD 219134 b':
                        Ap = 5.0
                    if pl == 'GJ 1214 b':
                        Ap = 5.0
                    if pl == 'HD 209458 b':
                        Ap = 5.0
                    F = 31.30*Ap
                    
                if ch == 'FGS Prime':
                    Ap = 7.0
                    if pl == 'HD 219134 b':
                        Ap = 6.0
                    if pl == 'GJ 1214 b':
                        Ap = 7.0
                    if pl == 'HD 209458 b':
                        Ap = 7.0
                    F =  24.62*Ap
  
                if ch == 'NIR Phot':
                    Ap = 9.0
                    if pl == 'HD 219134 b':
                        Ap = 9.0
                    if pl == 'GJ 1214 b':
                        Ap = 8.0
                    if pl == 'HD 209458 b':
                        Ap = 9.0
                    F = 39.46*Ap
    
                hdu = pyfits.open(sim_file)           
      
                SNR = exosim.tools.snrprism2.SnrTool()
                if ch == 'FGS Red' or ch == 'FGS Prime'  or ch == 'NIR Phot':
                    SNR = exosim.tools.snrphotometer.SnrTool()
                                     
                              
                if ch == 'FGS Red' or ch == 'FGS Prime' or ch == 'NIR Phot':
                    wl, sig, no = SNR.calc_exosim_SNR(sim_file, pixSize, apply_cds, 
                                                                         F, qe, decorr, 
                                                                         diff, ch, apply_flat) 
                    print wl, sig, no
                    if ct==0:
                        noise_results = np.hstack((wl,sig,no))
                    elif 'Noise option:  10' in open(txt_file).read():
                        noise_results[1] = sig
                    else:
                        noise_results = np.hstack((noise_results,no))
                    
                    print noise_results
                        
        
                                                                         
                else:
                    wl, sig, no = SNR.calc_exosim_SNR(sim_file, R, pixSize, apply_cds, 
                                                                         F, qe, decorr, 
                                                                         diff, wavrange, ch, apply_flat)  


                    for i in range (len(wl)):
                        print wl[i], sig[i], no[i]
 
                    if ct==0:
                        noise_results = np.hstack((wl.reshape(len(wl),1),sig.reshape(len(sig),1),no.reshape(len(no),1)))
                    elif 'Noise option:  10' in open(txt_file).read():
                        noise_results[:,1] = sig
                    else:
                        noise_results = np.hstack((noise_results,no.reshape(len(no),1)))
                        

        np.save('/Volumes/My Passport for Mac/ExoSimOutput/%s_%s--%s_%sx%s_noise_results'%(pl, ch, noise_group, len_group, No_real), noise_results) 
                   
                
                
                
                            
                     
 