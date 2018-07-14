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
import os

import pyfits


#a = pyfits.open('/Users/user1/ExoSimOutput_FGS2_JN_test/HD 209458 b_FGS Red--jitterful/sim_0000/FGS Red_signal.fits' )
#
#
#
#
##
#x_offsets = a[-3].data[:,1]
#y_offsets = a[-3].data[:,2]
#NDRs = a[-3].data[:,0]
#
#x_av_list = []
#x_list = []
#y_av_list = []
#y_list = []
#
#for i in range (len(NDRs)):
#    if i>1:
#        if NDRs[i] != NDRs[i-1]:
#            x_av = np.mean(x_list)
#            x_av_list.append(x_av)  
#            x_list =[]
#
#            y_av = np.mean(y_list)
#            y_av_list.append(y_av)  
#            y_list =[]            
#            
#            
#    x_list.append(x_offsets[i])
#    y_list.append(y_offsets[i])
#    
# 
#    if i == len(NDRs)-1:
#        x_av = np.mean(x_list)
#        x_av_list.append(x_av)
#            
#        y_av = np.mean(y_list)
#        y_av_list.append(y_av)
#    
#pscale = 3.60163e-05
#pointing_x =    np.array(x_av_list)/pscale
#pointing_y =    np.array(y_av_list)/pscale
#pointing_x = pointing_x - pointing_x[0]
#pointing_y = pointing_y - pointing_y[0]
#
#
#np.save('/Users/user1/ExoSimOutput_FGS2_JN_test/pointing_x.npy', pointing_x)      
#np.save('/Users/user1/ExoSimOutput_FGS2_JN_test/pointing_y.npy', pointing_y)      
#       
#xxx




#inputs
pl = 'GJ 1214 b'
#pl = 'HD 209458 b'
 
ch = 'AIRS CH0' 
#ch = 'NIR Spec'  
#ch = 'FGS Red'
#ch = 'FGS Prime'
####
#ch = 'NIR Phot'
#ch = 'AIRS CH1' 
# 
diff = 0
fullTransit = False


RPE_var = 0
Inc_PSD = 0
no_RPE = 0
  
incAp = 0
  
PN = 0

a = 'no_RPE_var'
b = 'normal_PSD'
c ='full_RPE'
  
if RPE_var == 1:
  a = 'with_RPE_var'
if Inc_PSD == 1:
  b = 'Inc_PSD'
if no_RPE == 1:
  c = 'no_RPE'
  
label = '%s_%s_%s'%(a,b,c)

if PN == 1:
    label = ''



if PN ==0:
    root = '/Users/user1/ExoSimOutput_JN_166/%s/%s_%s--jitterful'%(label, pl, ch)
else: 
    root = '/Users/user1/ExoSimOutput_PN_166/%s/%s_%s--jitterless'%(label, pl, ch)
    


if ch == 'NIR Spec' or ch == 'FGS Red' or ch == 'FGS Prime' or ch == 'NIR Phot':
    if pl == 'HD 209458 b':
        sim_list = ['0000' , '0001' , '0002', '0003', '0004']
    else:
        sim_list = ['0000']
else:
    sim_list = ['0000']   
    
    
print "...", pl, ch, sim_list

#for ii in ['jitterful' , 'jitterless']:
for ii in ['jitterful']:
#for ii in ['jitterless']:

#for ii in ['full_transit']:

    
    if ii == 'jitterless' or PN ==1:
#            if jj == '0001' or jj == '0003' or jj== '0004':
#                diff =1
        doJitterDecorr = False
    elif ii == 'jitterful':
        doJitterDecorr = True
    elif ii == 'full_transit':
        doJitterDecorr = False 
        fullTransit = True #make False if doing OOT simulation with no light curve fit
        

    exosim_path = '/Users/user1/Desktop/ExoSim-nb/exosim_nb'
    QE_rms_file = '/Users/user1/Desktop/ExoSim-nb/exosim_nb/data/qe_rms.npy'  #IMPORTANT! This must be the same QE rms file used in the simulation
             
    ICF = '/Users/user1/Desktop/ExoSim-nb/exosim_nb/xml_files/exosim_ariel_mcr_Euro.xml'
    binSize = 1
    binning = 'R-bin'
#    binning = 'fixed-bin'
    
    
    
    for sim in sim_list:
    
        aa = '%s/sim_%s/%s_signal.fits'%(root, sim, ch)
                 
        data, info = loadFits(aa)
        
        print "doing CDS, sim..", sim, data.shape
        data = doCDS(data,info)

        if sim == sim_list[0]:
            data_stack = data
        else:
            data_stack = np.dstack((data_stack, data))
            
    data = data_stack
        
    info['NEXP'] = data.shape[2]
    print data.shape[2]
        
        

#        data = subtractDark(data, info, ICF, ch)

#        data = flatField(data, ICF, ch, QE_rms_file)

#        if diff == 0:
#            data = backSub(data)  
      
    


 
    if ch == 'AIRS CH0' or ch== 'AIRS CH1' or ch == 'NIR Spec':
        
       print "spectroscopic channel...", ch
       
       if doJitterDecorr == True:
            print "doing jitter decorrelation"
            method ='fft'  
            epd = ExosimPointingDecorr_(data, info, method)
            data = epd.getHdu(data)       
       
       print "doing masking" 
    
#       data = spectroscopic_mask(data, info, ICF, ch, diff)
#
#       print "extracting 1D spectra"
#
#       data = extract1D(data)
#       
       
       data = spectroscopic_mask_plus_extraction(data, info, ICF, ch, diff)

       print "binning spectra"
          
       data, wl, binsizes  = binSpectrum(data, info, ch, binning, binSize)
        
    else:

       print "photometric channel...", ch
       
       print "doing masking / decorrelation" 
     
       data, wl = photometric_mask(data, info, ICF, ch, diff, incAp, pl)
       

       
#       d_av = obtain_average_Gaussian(data, info, ICF, ch, diff)
#       data, wl = use_average_Gaussian(data, info, ICF, ch, diff, d_av)      
       
#       pointing_x = np.load('/Users/user1/ExoSimOutput_FGS2_JN_test/pointing_x.npy')
#       pointing_y = np.load('/Users/user1/ExoSimOutput_FGS2_JN_test/pointing_y.npy')
#       
#       data, wl =use_pointing(data, info, ICF, ch, diff, pointing_x, pointing_y)

    
    print data.shape
 

#    np.save('%s/SvRuns/%s_%s_sim%s_RPEvar.npy'%(root_2,pl, ch, jj), data)  
#    np.save('%s/SvRuns/%s_%s_sim%s_incAp.npy'%(root_2,pl, ch, jj), data)  
    
    sim = sim_list[0]
    if ch == 'FGS Red' or  ch == 'FGS Prime' or ch == 'NIR Phot':
        if incAp ==1:
            if PN ==0:
                root2 = '/Users/user1/ExoSimOutput_JN_166/%s/%s_%s--jitterful'%(label+'_incAp', pl, ch)
                if not os.path.exists(root2):
                    os.makedirs(root2)
                np.save('%s/%s.npy'%(root2,sim), data)  
            elif PN ==1:
                root2 = '/Users/user1/ExoSimOutput_PN_166/%s/%s_%s--jitterless'%(label+'_incAp', pl, ch)
                if not os.path.exists(root2):
                    os.makedirs(root2)
                np.save('%s/%s.npy'%(root2,sim), data)                  
        else:
            np.save('%s/%s.npy'%(root,sim), data)  
    else:
            np.save('%s/%s.npy'%(root,sim), data)        
                
    
    
