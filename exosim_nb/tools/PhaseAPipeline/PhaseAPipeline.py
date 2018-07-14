# -*- coding: utf-8 -*-
"""
Created on Tue May 15 10:23:25 2018

@author: user1

Phase A ExoSim pipeline
"""

import numpy as np
from calibration import*
from jitterdecorr import*
from binning import*
from signal_noise import*



#ch = 'AIRS CH0' 
##
##ch = 'AIRS CH1' 
##
##ch = 'NIR Spec'  
##ch = 'FGS Red'
##ch = 'FGS Prime'
###
##ch = 'NIR Phot'
##
##fileName = '/Users/user1/ExoSimOutput_JN_5/GJ 1214b_%s--jitterful/sim_0000/%s_signal.fits'%(ch, ch)
##
#fileName = '/Users/user1/ExoSimOutput_JN_6_highPSD/GJ 1214b_%s--jitterful/sim_0000/%s_signal.fits'%(ch, ch)
##
##
#aa, info = loadFits(fileName)
##aa = aa[...,0:10000]
##
##fileName = '/Users/user1/ExoSimOutput_JN_6/GJ 1214b_%s--jitterful/sim_0000/%s_signal.fits'%(ch, ch)
##
##
#fileName = '/Users/user1/ExoSimOutput_JN_6_highPSD_noRPE/GJ 1214b_%s--jitterful/sim_0000/%s_signal.fits'%(ch, ch)
##
#bb, info = loadFits(fileName)
##
##
##if ch == 'AIRS CH1' :
##    T= 5.83648039021     #ch1 
##
##elif ch == 'AIRS CH0'  :   
##    T= 5.52851157725   #ch0                      
##
##
##elif ch == 'NIR Spec'   :  
##
##    T = 1.49184895196 #ns
##
##elif ch == 'FGS Red'   :  
##
##    T = 1.52305165888  #fgs2
##    
##elif ch == 'FGS Prime'     :
##
##    T = 0.689974611433 #fgs1
##    
##elif ch == 'NIR Phot'    : 
##          
##    T = 4.19670980005
#
#
#
#
#if ch == 'AIRS CH1' :
#    T= 37.7833314049     #ch1 
#
#elif ch == 'AIRS CH0'  :   
#    T= 71.2560978389   #ch0                      
#
#
#elif ch == 'NIR Spec'   :  
#
#    T = 23.7516087501 #ns
#
#elif ch == 'FGS Red'   :  
#
#    T = 30.173994731  #fgs2
#    
#elif ch == 'FGS Prime'     :
#
#    T = 23.49142178533 #fgs1
#    
#elif ch == 'NIR Phot'    : 
#          
#    T = 300.0
#    
#    
#    
#
#N = int(np.ceil(90/T))
#
#print N
#
#tot = int(aa.shape[2]/(N))
#if tot%2 !=0:
#    tot = tot-1
#
#print tot*N*2
#
#
#cc = np.zeros((aa.shape[0], aa.shape[1], 2*N*tot ))
#
#ct=0
#for i in range(0, cc.shape[2], 2*2*N):
#    
#    cc[..., i:i+(N*2)] = aa[..., ct:ct+(N*2)]
#    cc[..., i+(N*2) : i+ (2*2*N)] = bb[..., ct:ct+(N*2)]
#    
#    ct = ct + N*2
#  
#
#s=[]
#for i in range(1,cc.shape[2],2):
#    s.append(cc[...,i].sum())
#  
#plt.plot(s, 'rx-') 
#
#np.save('/Users/user1/ExoSimOutput_JN_6_highPSD/SvRuns/GJ 1214b_%s-interleaved_sequence.npy'%(ch), cc)
#
#xxx 
#
#xxx
##
##bb = np.load('/Users/user1/Desktop/ExoSimOutput-255A-PN/SvRuns/GJ1214_NIR Spec_sim1.npy')
#
#
#aa =    np.load('/Users/user1/ExoSimOutput_pointing_JN_3/SvRuns/GJ 1214 b_NIR Spec_sim0000.npy')
#
#aa = aa[...,0::2] 
#
#bb = aa[:,6]
#
#plt.plot(bb)
#
#no = []
#for i in range(2, 10, 1):
##    idx = np.arange(0, aa.shape[0], i)
##    
##    dd = np.add.reduceat(aa, idx, axis =0)[:-1] / i
##    no =   dd.std(axis = 0)/ dd.mean(axis=0)  
##    
##    if i ==2:
##        no_stack = no
##    else: no_stack = np.vstack((no_stack, no))
#    
#    idx = np.arange(0, len(bb), i)
#    
#    dd = np.add.reduceat(bb, idx)[:-1] 
#    no.append(dd.std()/ dd.mean()  )
#  
#
#    
#    
#    
#    
#
#
#
#xxxx


#inputs
pl = 'GJ 1214 b'
#pl = 'HD 209458 b'

#pl = 'GJ1214'

#ch = 'AIRS CH0' 
#
#ch = 'AIRS CH1' 

ch = 'NIR Spec'  
#ch = 'FGS Red'
#ch = 'FGS Prime'
#
#ch = 'NIR Phot'

 


#root = '/Users/user1/ExoSimOutput_pointing_PN'

#root = '/Users/user1/Desktop/ExoSimOutput-255A-PN' 

data_stack =[]

#for ii in ['jitterful' , 'jitterless']:
for ii in ['jitterful']:
#for ii in ['jitterless']:

#    for kk in [0,1]:
#    
#           
#        diff = 0
#        fullTransit = False #make False if doing OOT simulation with no light curve fit
#            
#        if ii == 'jitterless':
#            
#            if kk==0:
#                root = '/Users/user1/ExoSimOutput_pointing_PN_2'      
#            else:
#                root = '/Users/user1/ExoSimOutput_pointing_PN_3'      
#    
#        
#    #        sim_list=[ '0000', '0001', '0002', '0003', '0004',
#    #        '0005', '0006', '0007', '0008', '0009',
#    #        '0010', '0011', '0012', '0013', '0014',
#    #        '0015', '0016', '0017', '0018', '0019']
#    
#            sim_list=[ '0000']    
#            
#            
#        elif ii == 'jitterful':
#            
#            if kk==0:
#                root = '/Users/user1/ExoSimOutput_pointing_JN_2'      
#            else:
#                root = '/Users/user1/ExoSimOutput_pointing_JN_3' 
#            
#    #        sim_list=[ '0000', '0001', '0002', '0003', '0004',
#    #        '0005', '0006', '0007', '0008', '0009',
#    #        '0010', '0011', '0012', '0013', '0014',
#    #        '0015', '0016', '0017', '0018', '0019']
#    
#            sim_list=[ '0000']
#    
#    
#    
#        elif ii == 'full_transit':
#            sim_list=['0001']  
#            
#            
#        for jj in sim_list:
#            print jj
#            
#            if ii == 'jitterless':
#    #            if jj == '0001' or jj == '0003' or jj== '0004':
#    #                diff =1
#                doJitterDecorr = False
#            elif ii == 'jitterful':
#                doJitterDecorr = True
#            elif ii == 'full_transit':
#                doJitterDecorr = False 
#                fullTransit = True #make False if doing OOT simulation with no light curve fit
#                
#                
#                
#    
#            fileName = '%s/%s_%s--%s/sim_%s/%s_signal.fits'%(root,pl, ch,ii,jj, ch)
#            
#    #        fileName = '%s/%s_%s/sim_%s/%s_signal.fits'%(root,pl, ch,jj, ch)
#            
#            
#    #        fileName = '/Users/user1/Desktop/ExoSimOutput-255-PN/GJ1214_AIRS CH1/sim_0001/AIRS CH1_signal.fits'
#            
#            
#    #        fileName = '/Users/user1/ExoSimOutput/GJ 1214 b_AIRS CH0--jitterless/sim_0000/AIRS CH0_signal.fits'
#            
#            exosim_path = '/Users/user1/Desktop/ExoSim-nb/exosim_nb'
#            QE_rms_file = '/Users/user1/Desktop/ExoSim-nb/exosim_nb/data/qe_rms.npy'  #IMPORTANT! This must be the same QE rms file used in the simulation
#                 
#            
#    
#          
#    #        ch = 'AIRS CH1'
#            ICF = '/Users/user1/Desktop/ExoSim-nb/exosim_nb/xml_files/exosim_ariel_mcr_Euro.xml'
#            binSize = 5
#            binning = 'R-bin'
#    #        binning = 'fixed-bin'
#            
#                   
#            print fileName
#            data, info = loadFits(fileName)
#            
#    
#    #        data = subtractDark(data, info, ICF, ch)
#    
#    #        data = flatField(data, ICF, ch, QE_rms_file)
#    
#    #        if diff == 0:
#    #            data = backSub(data)  
#          
#      
#            data = doCDS(data,info) 
#            
#            data_stack.append(data)
#  
#            
#            
#    print data_stack[0].shape
#    d1 =  np.zeros((data_stack[0].shape[0],data_stack[0].shape[1],data_stack[0].shape[2]*2 ))
#    
#    
#
#    
#    ct = 0
#    for cc in range (0, data_stack[0].shape[2]*2, 2):
#        d1[...,cc] = data_stack[0][...,ct]
#        d1[...,cc+1] = data_stack[1][...,ct]
#        ct = ct+ 1
#        
#        
#        
#    
#
#    data = d1 
    




    diff = 0
    fullTransit = False #make False if doing OOT simulation with no light curve fit
        
    if ii == 'jitterless':

        
        

        root = '/Users/user1/ExoSimOutput_pointing_PN_4'      
        root = '/Users/user1/ExoSimOutput_PN_5'      
#        root = '/Users/user1/ExoSimOutput_NB'      
        
 
        sim_list=[ '0000']    
#        sim_list=[ '0000', '0001', '0002', '0003', '0004']
#        label_list = ['source', 'dark', 'read', 'zodiac', 'emission']
        
        
    elif ii == 'jitterful':
        
        root = '/Users/user1/ExoSimOutput_pointing_JN_4'   
        root = '/Users/user1/ExoSimOutput_JN_5'      
#        root = '/Users/user1/ExoSimOutput_NB'      
#        root = '/Users/user1/ExoSimOutput_JN_7'
#        root = '/Users/user1/ExoSimOutput_JN_6_highPSD'
       
        

        sim_list=[ '0000']
#        sim_list=['0000', '0001','0002']  
#        label_list = ['total', 'spatial jitter', 'spectral jitter']        
                
        
        

    elif ii == 'full_transit':
        sim_list=['0001']  
        
        
    for jj in sim_list:
        print jj
        
        if ii == 'jitterless':
            if jj == '0001' or jj == '0003' or jj== '0004':
                diff =1
            doJitterDecorr = False
        elif ii == 'jitterful':
            doJitterDecorr = True
        elif ii == 'full_transit':
            doJitterDecorr = False 
            fullTransit = True #make False if doing OOT simulation with no light curve fit
            
            
            

#        fileName = '%s/%s_%s--%s/sim_%s/%s_signal.fits'%(root,pl, ch,ii,jj, ch)      
        fileName = '/Users/user1/ExoSimOutput_JN_6_highPSD/NIR Spec_signal.fits'

   
#        fileName = '%s/%s_%s/sim_%s/%s_signal.fits'%(root,pl, ch,jj, ch)
        
        
#        fileName = '/Users/user1/Desktop/ExoSimOutput-255-PN/GJ1214_AIRS CH1/sim_0001/AIRS CH1_signal.fits'
        
        
#        fileName = '/Users/user1/ExoSimOutput/GJ 1214 b_AIRS CH0--jitterless/sim_0000/AIRS CH0_signal.fits'
        
        exosim_path = '/Users/user1/Desktop/ExoSim-nb/exosim_nb'
        QE_rms_file = '/Users/user1/Desktop/ExoSim-nb/exosim_nb/data/qe_rms.npy'  #IMPORTANT! This must be the same QE rms file used in the simulation
             
        

      
#        ch = 'AIRS CH1'
        ICF = '/Users/user1/Desktop/ExoSim-nb/exosim_nb/xml_files/exosim_ariel_mcr_Euro.xml'
        binSize = 5
        binning = 'R-bin'
#        binning = 'fixed-bin'
        
               
        print fileName
        data, info = loadFits(fileName)
        
        data = np.load('/Users/user1/ExoSimOutput_JN_5/SvRuns/GJ 1214b_%s-interleaved_sequence.npy'%(ch))
        
#        print info['NEXP'], data.shape[2]/2
        
#        data = data[...,0:5000]
##        
        info['NEXP'] = data.shape[2]/2
        
#        print info['NEXP'], data.shape[2]/2     

#        data = subtractDark(data, info, ICF, ch)

#        data = flatField(data, ICF, ch, QE_rms_file)

#        if diff == 0:
#            data = backSub(data)  
      
  
        data = doCDS(data,info) 

    
        if ch == 'AIRS CH0' or ch== 'AIRS CH1' or ch == 'NIR Spec':
            
           print "spectroscopic channel...", ch
           
           if doJitterDecorr == True:
                print "doing jitter decorrelation"
                method ='fft'  
                epd = ExosimPointingDecorr_(data, info, method)
                data = epd.getHdu(data)       
           
           print "doing masking" 
        
           data = spectroscopic_mask(data, info, ICF, ch, diff)
    
           print "extracting 1D spectra"
    
           data = extract1D(data)
    
           print "binning spectra"
              
           data, wl, binsizes  = binSpectrum(data, info, ch, binning, binSize)
            
        else:
    
           print "photometric channel...", ch
           
           print "doing masking / decorrelation" 
         
           data, wl = photometric_mask(data, info, ICF, ch, diff)
        
        print data.shape
        
        import os
        if not os.path.exists('%s/SvRuns'%(root)):
           os.makedirs('%s/SvRuns'%(root))
    
        np.save('%s/SvRuns/%s_%s_sim%s-new-interleaved.npy'%(root,pl, ch, jj), data)  
#        np.save('%s/SvRuns/%s_%s_sim%s_wl.npy'%(root,pl, ch, jj), wl)  
        
        xxx
    

      
 
         
        if fullTransit == False:
            
            signal, noise = getSNR(data)   
    #        idx1 = np.argwhere(wl<1.95)[0].item()
    #        idx2 = np.argwhere(wl<3.9)[0].item()
    
    #        idx1 = np.argwhere(wl<3.9)[0].item()
    #        idx2 = np.argwhere(wl<7.8)[0].item()        
    #        
    #        
    #        wl = wl[idx2:idx1]
    #        noise = noise[idx2:idx1]
            plt.figure('noise variance/ time')
            plt.semilogy(wl,noise**2/info['CDS_time'],label = '%s-%s'%(ii,jj))
            plt.legend(loc='best'); plt.grid()
            
            plt.figure('noise/sinal')
            plt.semilogy(wl,noise/signal,label = '%s-%s'%(ii,jj))
            plt.legend(loc='best'); plt.grid()            
            
    
    
            plt.figure('signal')
            plt.semilogy(wl,signal/info['CDS_time'],label = '%s-%s'%(ii,jj))
            plt.legend(loc='best'); plt.grid()            
                
              
                  
        else:
            print 'hello'
            spectrum, noise = getSpectrum(data, info)

            plt.figure('spectrum')
            plt.semilogy(wl,spectrum,label = '%s-%s'%(ii,jj))
            plt.legend(loc='best'); plt.grid()                 
            

