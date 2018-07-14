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
#import exosim_t as exosim
import exosim_p as exosim
from exosim_p.PointingDecorSubi import ExosimPointingDecorr
#import jitter_ft



#a = '/Users/c1341133/ExoSimOutput/HD 219134 b_AIRS CH0/sim_0008/AIRS CH0_signal.fits'
#
#hdu = pyfits.open(a)
#
#print hdu[0].header['NDR_time']
#
#xx


sw_st = 0
mw_st = 0

qe_on = True # i.e. flat field
#qe_on = False # i.e. flat field


decorr_make = True
#decorr_make = False


# for simple version
decorr = True
decorr = False 
#

R_set = 20. # for NIR Spec

method = 'xcorr'
method = 'pointing'
method = 'fft'
p_spec = True
p_spat = True

#
#for ch in ['FGS Red', 'FGS Prime', 'NIR Phot']:

#for ch in ['FGS Prime']: 
#for ch in ['FGS Red']:  
#for ch in ['NIR Phot']:
#for ch in ['NIR Spec']:
    
#ch =  'NIR Spec'
#ch =  'AIRS CH1'
ch =  'AIRS CH0'

#for decorr in [False, True]:    
    
for decorr in [False]:    
#for decorr in [True]:
#for R_set in [10, 20]:
  
#    for pl in ['HD 219134 b', 'GJ 1214 b','HD 209458 b']:
#    for pl in ['GJ1214','HD209458']:
    for pl in ['GJ 1214 b']:
#    for pl in ['HD 209458 b']:
#    for pl in ['HD 219134 b']:
    
    
    

        sw_st0 =[]
        mw_st0 =[] 
        
        #for ResMag in [10,20,50,100,200]:
        for ResMag in [10]:
        
            
            if qe_on == True:
                qe = np.load('/Users/c1341133/Desktop/ExoSim-prism/exosim_p/data/qe_rms.npy') 
                qe_uncert = np.load('/Users/c1341133/Desktop/ExoSim-prism/exosim_p/data/qe_uncert.npy') 
            else:
                qe =  np.load('/Users/c1341133/Desktop/ExoSim-prism/exosim_p/data/qe_rms.npy')*0 + 1.0
            
            X =  ['0000','0001','0002', '0003','0004','0005','0006', '0007', '0008', '0009']
            diff =  [0, 0, 0, 1, 1, 1, 0, 0, 0, 0,]
            
#            X = ['0010', '0011','0012']
#            diff =  [0, 0, 0, 0]
        
#            X =  ['0002', '0003','0006', '0009']
#            diff =  [0, 1, 1, 0]   
#            
#            
#            X =  ['0000','0001']
#            diff =  [0, 0,]
            
            X =  ['0006']
            diff =  [0]
##
#            X =  ['0003','0004', '0005']
#            diff =  [1, 1, 1]  
#
#            
#            if ResMag ==10:
#                X =  ['0000','0001','0002', '0003']
#            elif ResMag ==20:
#                X =  ['0004','0005','0006', '0007']    
#            elif ResMag ==50:
#                X =  ['0008','0009','0010', '0011']
#            elif ResMag ==100:
#                X =  ['0012','0013','0014', '0015']
#            elif ResMag ==200:
#                X =  ['0016','0017','0018', '0019']
 
     
            if decorr == True:
                X = ['0000', '0007','0008', '0009']
                diff =  [0, 0, 0, 0]
                
#                X = ['0010', '0011','0012']
#                diff =  [0, 0, 0, 0]
##                
##                
##                
#                X = ['0010']
#                diff =  [0]
#        #        
        #        if ResMag ==10:
        #            X =  ['0000','0002', '0003']
        #        elif ResMag ==20:
        #            X =  ['0004','0006', '0007']    
        #        elif ResMag ==50:
        #            X =  ['0008','0010', '0011']
        #        elif ResMag ==100:
        #            X =  ['0012','0014', '0015']
        #        elif ResMag ==200:
        #            X =  ['0016','0018', '0019']
                
            XS= X[0]
            
            
##  -------------------------------------------------------------------------               
                 
#            #
#            #
            if decorr == True and decorr_make == True: 
                for j in range(len(X)):
                    i = X[j]
                    if method == 'pointing' or method == 'fft':
                        p_spec = True
                        p_spat = True
                      
        #    #            if i == '0007':
        #                if i == '0002':
        #                    p_spec = False
        #    #            if i == '0008':
                        if i == '0010':
                            p_spat = False
                            
            #            if i == '0007':
#                        if i == '0002' or i == '0006' or i =='0010' or i=='0014' or i=='0018':
#                            p_spec = False
#            #            if i == '0008':
#                        if i == '0003' or i == '0007' or i =='0011' or i=='0015' or i=='0019':
#                            p_spat = False
                            
                            
                        print "p_spec", p_spec
                        print "p_spat", p_spat
                        
                        
        #            for ch in ['SWIR', 'MWIR']:    
        #            for ch in ['AIRS CH0', 'AIRS CH1']:
        #            for ch in ['NIR Spec']:
#                    for ch in ['FGS Red']:
        #            for ch in ['FGS Prime']:
        #            for ch in ['NIR Phot']:
                    for k in [0]:
        
        
                        #fitsFileName = '/Users/c1341133/ExoSimOutput/Twinkle_test_2/HD 209458 b/airy_pri/%s_signal.fits'%(ch)
                #        fitsFileName = '/Users/c1341133/ExoSimOutput/GJ1214 LRS/sim_%s/%s_signal.fits'%(i,ch)
            #            fitsFileName = '/Users/c1341133/ExoSimOutput/generic_1m/sim_%s/%s_signal.fits'%(i,ch)
            #            if ch == 'AIRS CH0':
            #                pscale = 7.61e-5
            #            elif ch == 'AIRS CH1':
            #                pscale = 1.52e-4
            #            elif ch == 'NIR Spec':
            #                pscale = 3.20e-5
            #                
            #                
                        
                        if ch == 'AIRS CH0':
                            pscale = 6.11e-5
                        elif ch == 'AIRS CH1':
                            pscale = 1.23e-4
                        elif ch == 'NIR Spec':
                            pscale = 3.01e-5  # use for baseline
#                            pscale = 6.11e-5
                        elif ch == 'SWIR':
                            pscale =5.583333e-5
                        elif ch == 'MWIR':
                            pscale = 11.166666e-5
                        elif ch == 'FGS Red':
                            pscale = 2.23e-5
                        elif ch == 'FGS Prime':
                            pscale = 3.10e-5
                        elif ch == 'NIR Phot':
                            pscale = 2.43e-5
                            
                            
                            
                            
#                        fitsFileName = '/Users/c1341133/ExoSimOutput/sim_%s/%s_signal.fits'%(i,ch)
                        fitsFileName = '/Users/c1341133/ExoSimOutput/%s_%s/sim_%s/%s_signal.fits'%(pl,ch,i,ch)

                        if ch == 'SWIR' or ch == 'MWIR':
                            fitsFileName = '/Users/c1341133/ExoSimOutput/sim_%s/%s_signal.fits'%(i,ch)

                
                        hdu = pyfits.open(fitsFileName)
            
                        fileName = fitsFileName
                        epd = ExosimPointingDecorr(fileName, qe, ch, method, p_spec, p_spat, pscale, diff)
                        hdu = epd.getHdu()
#             
 #-------------------------------------------------------------------------               
             
            for j in range(len(X)):
                i= X[j]
                
        #        for ch in ['SWIR', 'MWIR']:        
        #        for ch in ['AIRS CH0', 'AIRS CH1']:
        #        for ch in ['NIR Spec']:
        #        for ch in ['FGS Red']:
        #        for ch in ['FGS Prime']:
#                for ch in ['NIR Phot']:
                for k in [0]:
        #
        #    
                    #fitsFileName = '/Users/c1341133/ExoSimOutput/Twinkle_test_2/HD 209458 b/airy_pri/%s_signal.fits'%(ch)
            #        fitsFileName = '/Users/c1341133/ExoSimOutput/GJ1214 LRS/sim_%s/%s_signal.fits'%(i,ch)
            #        fitsFileName = '/Users/c1341133/ExoSimOutput/generic_1m/sim_%s/%s_signal.fits'%(i,ch)
            
#                    fitsFileName = '/Users/c1341133/ExoSimOutput/sim_%s/%s_signal.fits'%(i,ch)
                    fitsFileName = '/Users/c1341133/ExoSimOutput/%s_%s/sim_%s/%s_signal.fits'%(pl,ch,i,ch)
                    
#                    fitsFileName = '/Volumes/My Passport for Mac/Ariel_mcr/%s/AIRS/sim_%s/%s_signal.fits'%(pl,i,ch)

                    if decorr == True:
#                        fitsFileName = '/Users/c1341133/ExoSimOutput/sim_%s/%s_signal_decorr.fits'%(i,ch)
                        fitsFileName = '/Users/c1341133/ExoSimOutput/%s_%s/sim_%s/%s_signal_decorr.fits'%(pl,ch,i,ch)
                        
#                        fitsFileName = '/Volumes/My Passport for Mac/Ariel_mcr/%s/AIRS/sim_%s/%s_signal_decorr.fits'%(pl,i,ch)
               
                    if ch == 'SWIR' or ch == 'MWIR':
                        fitsFileName = '/Users/c1341133/ExoSimOutput/sim_%s/%s_signal.fits'%(i,ch)
                        if decorr == True:
                            fitsFileName = '/Users/c1341133/ExoSimOutput/sim_%s/%s_signal_decorr.fits'%(i,ch)
               
                    fitsFileName = '/Users/c1341133/Downloads/AIRS CH0_signal (1).fits' 
#                    fitsFileName = '/Users/c1341133/Downloads/NIR Spec_signal (1).fits' 

            
                    hdu = pyfits.open(fitsFileName)
             
                    print fitsFileName
                    
            
                    #Use for Ariel prism
                    pixSize = 18.0
                    
                    pixSize = 15.0

                    
                    doCds = True
                    if decorr == True: 
                        doCds = False
            
                    # for normal F
              
        #            if ch == 'AIRS CH0':
        #                R = 100.
        #                F = 13.2
        #                wavrange = [1.8,4.0]
        #            elif ch == 'AIRS CH1':
        #                R = 30.
        #                F = 6.36
        #                wavrange = [3.8,8.0]
        #    
        #            if ch == 'NIR Spec':
        #                R = 10.
        #                F = 32.14
        #                wavrange = [0.9,2.2]
                    
                    # for modified F
        #            
        #            if ch == 'AIRS CH0':
        #                R = 100.
        #                F = 15.3
        #                wavrange = [1.8,4.0]
        #            elif ch == 'AIRS CH1':
        #                R = 30.
        #                F = 7.7
        #                wavrange = [3.8,8.0]
                        
                    #for modified model            
            
        #      
                    if ch == 'AIRS CH0':
                        R = 100.
        #                F = 13.2
        #                F = 20.5
        #                F = 30.0
                        F = 18.36
                        F = 13.2
                        wavrange = [1.8,4.0]
                        
                    elif ch == 'AIRS CH1':
                        R = 30.
        #                F = 6.36
        #                F = 10.3
        #                F = 15.0
                        F = 9.24
                        F = 6.36
                        wavrange = [3.8,8.0]
        #    
                    if ch == 'NIR Spec':
                        R = R_set
        #                F = 32.14
        #                F = 15.3
        #                F = 30.6
        #                F = 40.0
        #                F = 64.28
#                        F = 48.0
        #                F = 77.136
                        F = 63.408
                        F = 2*26.42
                        
        #                F=40.0
                        wavrange = [0.9,2.2]
        #
                    if ch == 'FGS Red':
                        F = 41.364
                        F = 41.364 * 2.35 # correction for WFE 1st min
                        F = 41.364 * 4 # correction for WFE 2nd min
                        F = 34.47 *4
                    if ch == 'FGS Prime':
                        F = 30.0
                        # correction for WFE
                        F = 2.75*30.0
                        F = 4.5 * 30.0 # correction for WFE 2nd min
                        F = 25*4.5   
                        
                    if ch == 'NIR Phot':
                        F = 38.568
                        # correction for WFE
                        F = 4.55* 38.568
                        F = 7.7 *38.568 # correction for WFE 2nd min
                        F = 32.14*7.7     
                        F = 32.14*9.0  
                        
        #            if ch == 'SWIR':
        #                R = 100
        #                F = 20.5
        ##                F = 13.2
        #                wavrange = [1.8,4.0]
        #                pixSize = 18.0
        #                ld = [1.95, 0.0004896, 360]
        #    
        #            if ch == 'MWIR':
        #                R = 30
        #                F = 10.3
        ##                F = 6.36
        #                wavrange = [3.8,8.0]
        #                pixSize = 18.0
        #                ld = [3.9, 0.0032883, 360]
            
                        
                    SNR = exosim.tools.snrprism.SnrTool()
        #            if ch == 'SWIR' or ch == 'MWIR':
        #                SNR = exosim.tools.snr4.SnrTool()
                    if ch == 'FGS Red' or ch == 'FGS Prime'  or ch == 'NIR Phot':
                        SNR = exosim.tools.snrphotometer.SnrTool()
                        
                    
                    
                    decorr_status = False
                    if decorr == True: 
                        decorr_status = True
                      
                    print "decorr status", decorr_status
                    print "do CDS?", doCds
                    print "qe", qe_on
            
                    if ch == 'FGS Red' or ch == 'FGS Prime' or ch == 'NIR Phot':
                        wl, sig, no = SNR.calc_exosim_SNR(fitsFileName, pixSize, doCds, 
                                                                             F, qe, decorr_status, 
                                                                             diff[j]) 
        #                                                                     
            
            
                    else:
                        wl, sig, no, binsize = SNR.calc_exosim_SNR(fitsFileName, R, pixSize, doCds, 
                                                                             F, qe, decorr_status, 
                                                                             diff[j], wavrange, ch, False)    
        #                                                                     
        #            if ch == 'SWIR' or ch == 'MWIR':
        #                wl, sig, no = SNR.calc_exosim_SNR(fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff, optimal =False, superbin=False,maskon=True)
        #                
        #            else:
        #                wl, sig, no = SNR.calc_exosim_SNR(fitsFileName, R, pixSize, doCds, 
        #                                                                     F, qe, decorr_status, 
        #                                                                     diff[j], wavrange, ch)
                                                                             
                                                                             
            #        
            #        wl, sig, no  = SNR.calc_exosim_SNR(fitsFileName, R, pixSize, doCds, F, wavrange, 0)
            
            
                    if ch == 'AIRS CH0' or ch == 'NIR Spec' or ch == 'SWIR':
                        if i == XS:
            #                sw_st = np.array((wl,binsize,binpix,sig,no1))
                            sw_st = np.array((wl,binsize, sig,no))
                        else:
                            sw_st = np.vstack((sw_st, no))
                            
                    elif ch == 'AIRS CH1' or ch == 'MWIR':
                        if i == XS :
            #                mw_st = np.array((wl,binsize,binpix,sig,no1))
                            sw_st = np.array((wl,binsize,sig,no))
                            
                        else:
                            sw_st = np.vstack((sw_st, no))
                    elif ch =='FGS Red' or ch == 'FGS Prime' or ch == 'NIR Phot':
                         if i == XS:
                            sw_st = np.array((wl, sig,no))
                         else:
                            sw_st = np.hstack((sw_st, no))
                         
                        
                            
            #        wl0, sig0, no0 = SNR.fit_lc_mandel_agol(fitsFileName, R, ld, pixSize, doCds, F)
                
                    hdu.close()
            if ch == 'FGS Red' or ch == 'FGS Prime'or ch == 'NIR Phot':
                pass
            else:
                sw_st = np.transpose(sw_st)
                mw_st = np.transpose(mw_st)
                print "hello"
               
                AA = np.zeros(sw_st.shape[1])
            
                if ResMag == 10:
                    sw_st0 = np.zeros(sw_st.shape[1])
            #        mw_st0 = np.zeros(mw_st.shape[1])
                    
                sw_st0 = np.vstack((sw_st0,AA,AA,sw_st))  # used when using multiple Rs
            #    mw_st0 = np.vstack((mw_st0,AA,mw_st))
        
        
        import csv
        import time
        cl = time.clock()
        
        if ch == 'FGS Red' or ch =='FGS Prime' or ch == 'NIR Phot':
            sw_st = np.vstack((sw_st, sw_st))
        #    with open("/Users/c1341133/Desktop/sw_st_%s_qe-%s_d-%s.csv"%(1, qe_on, decorr),"wb") as f:
            with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s.csv"%(ch, pl, decorr),"wb") as f:
        
                writer = csv.writer(f)
                writer.writerows(sw_st)

                
        if ch == 'NIR Spec':
            xx = np.vstack((np.zeros(sw_st.shape[1]),np.zeros(sw_st.shape[1])))
            
            with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_narrow.csv"%(ch, pl,decorr,R_set),"wb") as f:
                writer = csv.writer(f)
                writer.writerows(sw_st)
                writer.writerows(xx)

        if ch == 'AIRS CH0':
            xx = np.vstack((np.zeros(sw_st.shape[1]),np.zeros(sw_st.shape[1])))
            
            with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_narrow.csv"%(ch, pl,decorr,R),"wb") as f:
                writer = csv.writer(f)
                writer.writerows(sw_st)
                writer.writerows(xx)
        
        if ch == 'AIRS CH1':
            xx = np.vstack((np.zeros(sw_st.shape[1]),np.zeros(sw_st.shape[1])))
            
            with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_narrow.csv"%(ch, pl,decorr,R),"wb") as f:
                writer = csv.writer(f)
                writer.writerows(sw_st)
                writer.writerows(xx)

#        else:        
#            xx = np.vstack((np.zeros(sw_st.shape[1]),np.zeros(sw_st.shape[1])))
#            
#            with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s.csv"%(ch, pl,decorr),"wb") as f:
#                writer = csv.writer(f)
#                writer.writerows(sw_st)
#                writer.writerows(xx)
#                writer.writerows(mw_st)
                
#            with open("/Users/c1341133/Desktop/mw_st_%s_qe-%s_d-%s-%s.csv"%(1, qe_on, decorr, R), "wb") as f:
#                writer = csv.writer(f)
#                writer.writerows(mw_st)
            #    
            #with open("/Users/c1341133/Desktop/sw_st_%s_qe-%s_d-%s.csv"%(1, qe_on, decorr),"wb") as f:
            #    writer = csv.writer(f)
            #    writer.writerows(sw_st0)
            #    
            #with open("/Users/c1341133/Desktop/mw_st_%s_qe-%s_d-%s.csv"%(1, qe_on, decorr), "wb") as f:
            #    writer = csv.writer(f)
            #    writer.writerows(mw_st0)
        
            
            
        
        
        
  