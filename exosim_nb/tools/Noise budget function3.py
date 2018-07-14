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
#decorr = False 
#


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
#    
#ch =  'NIR Spec'
ch =  'AIRS CH0'
ch =  'AIRS CH1'
#ch = 'FGS Red'
#ch = 'FGS Prime'
#ch = 'NIR Phot'

#jitter_batch = True
#jitter_batch = False

pixSize = 18.0

if ch == 'AIRS CH0' or ch == 'AIRS CH1':
    pixSize = 15.0
    
print "pix size", pixSize

#for decorr in [False, True]:    
#    
#for decorr in [True]:   
#for decorr in [False]: 

#for jitter_batch in [True, False]:
for jitter_batch in [True]:
#for jitter_batch in [False]:
#
   
    
    if jitter_batch == False:
        decorr = False
    elif jitter_batch == True:
        decorr = True
#        decorr = False # to get undecorrelated jitter noise

        if ch == 'FGS Red' or ch=='FGS Prime' or ch =='NIR Phot':
            decorr = False
            decorr = True

    no_exp = 30
    
    step = 6
    if jitter_batch == True:
#        no_exp = 1
#        step = 4
        no_exp = 1
        step = 4        
        
        
        if ch == 'FGS Red' or ch=='FGS Prime' or ch =='NIR Phot':
            no_exp = 40
            step = 2   
            
        
    print 'jitter_batch', jitter_batch, 'decorr', decorr, 'pixSize', pixSize
        
#for decorr in [True]:
#for R_set in [10, 20]:
  
#    for pl in ['HD 219134 b', 'GJ 1214 b','HD 209458 b']:
#    for pl in ['GJ 1214 b','HD 209458 b']:
#    for pl in ['GJ 1214 b']:
#    for pl in ['HD 209458 b']:
    for pl in ['HD 219134 b']:
        
#        for kk in range (0,no_exp,step):
        for kk in range (0,no_exp,step):


    
    
    

            sw_st0 =[]
            mw_st0 =[] 
            
            #for ResMag in [10,20,50,100,200]:
            for ResMag in [10]:
            
                
                if qe_on == True:
                    qe = np.load('/Users/c1341133/Desktop/ExoSim-prism/exosim_p/data/qe_rms.npy') 
                    qe_uncert = np.load('/Users/c1341133/Desktop/ExoSim-prism/exosim_p/data/qe_uncert.npy') 
                else:
                    qe =  np.load('/Users/c1341133/Desktop/ExoSim-prism/exosim_p/data/qe_rms.npy')*0 + 1.0
                
    #            X =  ['0000','0001','0002', '0003','0004','0005','0006', '0007', '0008', '0009']
    #            diff =  [0, 0, 0, 1, 1, 1, 0, 0, 0, 0,]
                
                a0 ='000%s'%(kk)
                a1 ='000%s'%(kk+1)
                a2 ='000%s'%(kk+2)
                a3= '000%s'%(kk+3)
                a4= '000%s'%(kk+4)
                a5= '000%s'%(kk+5)
                
                if len(a0) > 4:
                    a0 = a0[1:]
                if len(a1) > 4:
                    a1 = a1[1:]
                if len(a2) > 4:
                    a2 = a2[1:]
                if len(a3) > 4:
                    a3 = a3[1:]
                if len(a4) > 4:
                    a4 = a4[1:]                
                if len(a5) > 4:
                    a5 = a5[1:]
                
                X = [a0,a1,a2,a3,a4,a5]    
                diff =  [0, 0, 1, 1, 1,0] 
                
                if jitter_batch == True:
                    X = [a0,a1,a2,a3]    
                    diff =  [0, 0, 0,0] 
                    if ch == 'FGS Red' or ch=='FGS Prime' or ch =='NIR Phot':
                        X = [a0,a1]    
                        diff =  [0, 0] 
 
    
                    
                
                
#                X = ['000%s'%(kk)]
#                diff =  [0]
                
                print X
                
                
    #            X = ['0010', '0011','0012']
    #            diff =  [0, 0, 0, 0]
            
    #            X =  ['0002', '0003','0006', '0009']
    #            diff =  [0, 1, 1, 0]   
    #            
    #            
    #            X =  ['0000','0001']
    #            diff =  [0, 0,]
                
    #            X =  ['0007']
    #            diff =  [0]
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
    #                X = ['0000', '0007','0008', '0009']
    #                diff =  [0, 0, 0, 0]
    
#                    X =  ['0000', '0001', '0002']
#                    diff =  [0, 0, 0,]
                                
                    a0 ='000%s'%(kk)
                    a1 ='000%s'%(kk+1)
                    a2 ='000%s'%(kk+2)
                    a3= '000%s'%(kk+3)
                    
                    if len(a0) > 4:
                        a0 = a0[1:]
                    if len(a1) > 4:
                        a1 = a1[1:]
                    if len(a2) > 4:
                        a2 = a2[1:]
                    if len(a3) > 4:
                        a3 = a3[1:]
                    
                    X = [a0,a1,a2,a3]                     
                    diff =  [0, 0, 0, 0]
                    
                    if ch == 'FGS Red' or ch=='FGS Prime' or ch =='NIR Phot':
                        X = [a0,a1]    
                        diff =  [0, 0] 
 
                    
                    
                    
#                    X = ['000%s'%(kk)]
#                    diff =  [0]
                    print X                
                    
    #                X = ['0010', '0011','0012']
    #                diff =  [0, 0, 0, 0]
    ##                
    #                X =  ['0010','0011','0012', '0013','0014','0015','0016', '0017', '0018', '0019']
    #                diff =  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,]    
    ###                
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
                          
                            if i == '0001':
            #                if i == '0002':
                                p_spec = False
            #    #            if i == '0008':
                            if i == '002':
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
                                
                                
                                
                                
#                            fitsFileName = '/Users/c1341133/ExoSimOutput/%s_%s/sim_%s/%s_signal.fits'%(pl,ch,i,ch)                             
                            fitsFileName = '/Volumes/My Passport for Mac/ExoSimOutput_x_oldsims/%s_%s2/sim_%s/%s_signal.fits'%(pl,ch,i,ch)
#
#                            fitsFileName = '/Volumes/My Passport for Mac/ExoSimOutput_x_oldsims//%s_%s-esa/sim_%s/%s_signal.fits'%(pl,ch,i,ch)


                    
                            hdu = pyfits.open(fitsFileName)
                
                            fileName = fitsFileName
                            epd = ExosimPointingDecorr(fileName, qe, ch, method, p_spec, p_spat, pscale, diff[j])
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

#                        fitsFileName = '/Users/c1341133/ExoSimOutput/%s_%s/sim_%s/%s_signal.fits'%(pl,ch,i,ch)
#                      fitsFileName = '/Volumes/My Passport for Mac/Ariel_mcr/%s/AIRS/sim_%s/%s_signal.fits'%(pl,i,ch)
                        
                        if jitter_batch == False:
                            fitsFileName = '/Volumes/My Passport for Mac/ExoSimOutput_x_oldsims/%s_%s/sim_%s/%s_signal.fits'%(pl,ch,i,ch)
                        else:
                            fitsFileName = '/Volumes/My Passport for Mac/ExoSimOutput_x_oldsims/%s_%s2/sim_%s/%s_signal.fits'%(pl,ch,i,ch)

                        if decorr == True:
#                            fitsFileName = '/Users/c1341133/ExoSimOutput/%s_%s/sim_%s/%s_signal_decorr.fits'%(pl,ch,i,ch)
                            
#                            fitsFileName = '/Volumes/My Passport for Mac/Ariel_mcr/%s/AIRS/sim_%s/%s_signal_decorr.fits'%(pl,i,ch)
        
                            fitsFileName = '/Volumes/My Passport for Mac/ExoSimOutput_x_oldsims/%s_%s2/sim_%s/%s_signal_decorr.fits'%(pl,ch,i,ch)
                            
#                            fitsFileName = '/Volumes/My Passport for Mac/ExoSimOutput_x_oldsims/%s_%s-esa/sim_%s/%s_signal_decorr.fits'%(pl,ch,i,ch)

                
                        hdu = pyfits.open(fitsFileName)
                 
                        print fitsFileName
                        
                
                        #Use for Ariel prism

                        
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
#                            F = 18.36
                            F = 13.2
                            wavrange = [1.8,4.2]
                            
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
                            R = 20
            # 
#                            F = 50.0
                            F = 38.68
                            
            #                F=40.0
                            wavrange = [0.9,2.2]
                            


                        if ch == 'FGS Red':
                            if pl == 'HD 219134 b':
#                                Ap = 5.0
#                                Ap = 5.0
                                Ap = 5.0
                            if pl == 'GJ 1214 b':
#                                Ap = 5.0
#                                Ap = 6.0
                                Ap = 5.0
                            if pl == 'HD 209458 b':
#                                Ap = 5.0
#                                Ap = 5.0
                                Ap = 5.0
                            F = 31.30*Ap
                            
                        if ch == 'FGS Prime':
                            if pl == 'HD 219134 b':
#                                Ap = 6.0
#                                Ap = 5.0
                                Ap = 6.0
                            if pl == 'GJ 1214 b':
#                                Ap = 6.0
#                                Ap=7.0
                                Ap = 7.0
                            if pl == 'HD 209458 b':
#                                Ap = 5.0
#                                Ap = 5.0
                                Ap = 7.0
                            F =  24.62*Ap
      
                        if ch == 'NIR Phot':
                            if pl == 'HD 219134 b':
#                                Ap = 11.0
#                                Ap = 10.0
                                Ap = 9.0
                            if pl == 'GJ 1214 b':
#                                Ap = 8.0
#                                Ap = 8.0
                                Ap = 8.0

                            if pl == 'HD 209458 b':
#                                Ap = 9.0   
#                                Ap = 8.0
                                Ap = 9.0
                            F = 39.46*Ap
        
            
            
                            
                        print 'F number used', F
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
                
                            
                        SNR = exosim.tools.snrprism2.SnrTool()
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
                  
            
            import csv
            import time
            cl = time.clock()
            
            if kk ==0:
                Arr = sw_st
            else:
                Arr = np.hstack(((Arr,sw_st)))
                
        if ch == 'FGS Red' or ch == 'FGS Prime'or ch == 'NIR Phot':
            
            Arr = Arr.transpose()
            Arr1 = np.vstack((Arr,Arr))
            with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr.csv"%(ch, pl,decorr,jitter_batch),"wb") as f:
                writer = csv.writer(f)
                writer.writerows(Arr1)
            if jitter_batch == True: 
                to_no = []
                comb=[]
                for ii in range (0,len(Arr),4):
                    to_no.append(Arr[ii+2])
                    comb.append(Arr[ii+3])
                to_no0 = np.array(to_no).mean()
                comb0 = np.array(comb).mean()
            
                to_no_sd = np.array(to_no).std()
                comb_sd = np.array(comb).std()  
                Arr0 = np.hstack((to_no0,comb0, to_no_sd, comb_sd))
                Arr0 = np.vstack((Arr0,Arr0))
                with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr0.csv"%(ch, pl,decorr,jitter_batch),"wb") as f:
                    writer = csv.writer(f)
                    writer.writerows(Arr0)
            
            elif jitter_batch == False: 
                sig = []
                sn =[]; pn=[]; dn =[]; zn=[]; en=[]; rn=[]
                for ii in range (0,len(Arr),8):
                    print ii
                    sig.append(Arr[ii+1])
                    sn.append(Arr[ii+2])
                    pn.append(Arr[ii+3])
                    dn.append(Arr[ii+4])
                    zn.append(Arr[ii+5])
                    en.append(Arr[ii+6])
                    rn.append(Arr[ii+7])
                sig0 = np.array(sig).mean()
                sn0 = np.array(sn).mean()
                pn0 = np.array(pn).mean()
                dn0 = np.array(dn).mean()
                zn0 = np.array(zn).mean()
                en0 = np.array(en).mean()
                rn0 = np.array(rn).mean()
                
                sig_sd = np.array(sig).std()
                sn_sd = np.array(sn).std()
                pn_sd = np.array(pn).std()
                dn_sd = np.array(dn).std()
                zn_sd = np.array(zn).std()
                en_sd = np.array(en).std()
                rn_sd = np.array(rn).std()
            
                Arr0 = np.hstack((sig0, sn0,pn0,dn0,zn0, en0, rn0, sig_sd, sn_sd, pn_sd, dn_sd, zn_sd, en_sd, rn_sd))
                Arr0 = np.vstack((Arr0,Arr0))
                with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr0.csv"%(ch, pl,decorr,jitter_batch),"wb") as f:
                    writer = csv.writer(f)
                    writer.writerows(Arr0)
                    
        else:   
        
            if jitter_batch == True:        
                print Arr.shape
                xx = np.vstack((np.zeros(Arr.shape[1]),np.zeros(Arr.shape[1])))
                with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr-%s-2.csv"%(ch, pl,decorr,jitter_batch,method),"wb") as f:
                            writer = csv.writer(f)
                            writer.writerows(Arr)
                            writer.writerows(xx)       
                to_no = []
                spat=[]
                spec=[]
                comb=[]
                for ii in range (0,Arr.shape[1],7):
                    to_no.append(Arr[:,ii+3])
                    spat.append(Arr[:,ii+4])
                    spec.append(Arr[:,ii+5])
                    comb.append(Arr[:,ii+6])
                to_no0 = np.array(to_no).transpose().mean(axis=1)
                spat0 = np.array(spat).transpose().mean(axis=1)
                spec0 = np.array(spec).transpose().mean(axis=1)
                comb0 = np.array(comb).transpose().mean(axis=1)
                
                to_no_sd = np.array(to_no).transpose().std(axis=1)
                spat_sd = np.array(spat).transpose().std(axis=1)
                spec_sd = np.array(spec).transpose().std(axis=1)
                comb_sd = np.array(comb).transpose().std(axis=1)
        #        
                Arr0 = np.vstack((wl,to_no0,spat0,spec0,comb0, to_no_sd, spat_sd, spec_sd, comb_sd)).transpose()
                xx = np.vstack((np.zeros(Arr0.shape[1]),np.zeros(Arr0.shape[1])))
                with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr0-%s-2.csv"%(ch, pl,decorr,jitter_batch,method),"wb") as f:
                            writer = csv.writer(f)
                            writer.writerows(Arr0)
                            writer.writerows(xx)
    ###########################################################
    #        for ii in range (0,Arr.shape[1],4):
    #            to_no.append(Arr[:,ii+2])
    #            spat.append(Arr[:,ii+3])
    #        to_no0 = np.array(to_no).transpose().mean(axis=1)
    #        spat0 = np.array(spat).transpose().mean(axis=1)
    #    
    #        to_no_sd = np.array(to_no).transpose().std(axis=1)
    #        spat_sd = np.array(spat).transpose().std(axis=1)
    #        
    #        Arr0 = np.vstack((wl,to_no0,spat0,to_no_sd, spat_sd)).transpose()
    #        xx = np.vstack((np.zeros(Arr0.shape[1]),np.zeros(Arr0.shape[1])))
    #        with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr0.csv"%(ch, pl,decorr,R),"wb") as f:
    #                    writer = csv.writer(f)
    #                    writer.writerows(Arr0)
    #                    writer.writerows(xx)
    ###########################################################
            else:
                sig = []
                sn =[]; pn=[]; dn =[]; zn=[]; en=[]; rn=[]
                for ii in range (0,Arr.shape[1],9):
                    sig.append(Arr[:,ii+2])
                    sn.append(Arr[:,ii+3])
                    pn.append(Arr[:,ii+4])
                    dn.append(Arr[:,ii+5])
                    zn.append(Arr[:,ii+6])
                    en.append(Arr[:,ii+7])
                    rn.append(Arr[:,ii+8])
                sig0 = np.array(sig).transpose().mean(axis=1)
                sn0 = np.array(sn).transpose().mean(axis=1)
                pn0 = np.array(pn).transpose().mean(axis=1)
                dn0 = np.array(dn).transpose().mean(axis=1)
                zn0 = np.array(zn).transpose().mean(axis=1)
                en0 = np.array(en).transpose().mean(axis=1)
                rn0 = np.array(rn).transpose().mean(axis=1)
                
                sig_sd = np.array(sig).transpose().std(axis=1)
                sn_sd = np.array(sn).transpose().std(axis=1)
                pn_sd = np.array(pn).transpose().std(axis=1)
                dn_sd = np.array(dn).transpose().std(axis=1)
                zn_sd = np.array(zn).transpose().std(axis=1)
                en_sd = np.array(en).transpose().std(axis=1)
                rn_sd = np.array(rn).transpose().std(axis=1)
                
                Arr0 = np.vstack((wl,sig0, sn0,pn0,dn0,zn0, en0, rn0, sig_sd, sn_sd, pn_sd, dn_sd, zn_sd, en_sd, rn_sd)).transpose()
                xx = np.vstack((np.zeros(Arr0.shape[1]),np.zeros(Arr0.shape[1])))
                with open("/Users/c1341133/Desktop/sw_st_%s_%s_%s_%s_Arr02.csv"%(ch, pl,decorr,jitter_batch),"wb") as f:
                            writer = csv.writer(f)
                            writer.writerows(Arr0)
                            writer.writerows(xx)
    #         
    #      