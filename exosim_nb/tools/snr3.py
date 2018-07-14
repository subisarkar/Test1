"""
Created on Thu May 14 09:21:01 2015

@author: Subhajit Sarkar
Cardiff University

A Tool that reads in Exosim-created fits files 
and retrieves spectra and SNR from the 'data' output of ExoSim
The spectra are binned on a wavelength grid defined by R-based binning

Basic usage:

from exosim.tools.snr import SnrTool
snrTool = SnrTool()
wl, signal, noise = snrTool.calc_exolib_SNR(fitsFilePath, R, ld, pix_size)
wl, fittedDepth = snrTool.fit_lc_mandel_agol(fitsFilePath, R, ld, pix_size, doCds)

where:
fitsFilePath : pathe to fits file created by ExoSim
R            : R-based binning parameter 
ld           : Wavelength solution 
               (q, m, x_ref): wl = q + m(x - x_ref) where x is the pixel coordinate 
               in units of micron"
pix_size     : pixel size in microns              
"""

import pyfits
import numpy as np
import quantities  as aq
import pytransit
from scipy.optimize import fmin

class SnrTool():
    def __init__(self):
        self.lcf = LightcurveFitting()
        self.fdp = FitsDataProcessing()
    def calc_exosim_SNR(self, fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff):
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin, noiseBin, bin_size_r, bin_pix = self.fdp.bin_spec_data_R(fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff)
        signal = np.zeros(signalBin.shape[1])
        noise = np.zeros(signalBin.shape[1])
        noise_corr = np.zeros(signalBin.shape[1])
        nExp = signalBin.shape[0]
        for i in range (signalBin.shape[1]):
            #signal[i] = (signalBin[:,i].max() - signalBin[:,i].min()) /signalBin[:,i].max()
            #noise[i]  =  (np.std(noiseBin[:,i])    / signalBin[:,i].max()   )/ np.sqrt(nExp/2.)
#            signal[i] = signalBin[:,i].max()
            signal[i] = np.mean(noiseBin[:,i])

            noise[i]  =  np.std(noiseBin[:,i])
            noise_corr[i] = np.std(noiseBin[:,i])/ np.sqrt(nExp/2.)
        
        return wl, signal, noise, noise_corr, bin_size_r, bin_pix
    
    def fit_lc_mandel_agol(self, fitsFileName, R, ld, pixSize, doCds, F, qe, decorr):
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin, noiseBin, bin_size_r, bin_pix  = self.fdp.bin_spec_data_R(fitsFileName, R, ld, pixSize, doCds,F, qe, decorr, diff)
        hdul = pyfits.open(fitsFileName)
        ##z = hdul[-1].data['z'][1::2]
        multiaccum = hdul[0].header['MACCUM']
        z = hdul[-2].data['z'][multiaccum-1::multiaccum]
        hdul.close()
        u = [0,0] # modify by using u from fit files later
        sig = signalBin*1
        no2 = noiseBin*1
        nExp = signalBin.shape[0]

        fittedDepth =  np.zeros(signalBin.shape[1])
        for i in range (signalBin.shape[1]):
            no2[:,i] = no2[:,i] /  sig[:,i].max()
            sig[:,i] = sig[:,i] / sig[:,i].max()
            sig[:,i] = sig[:,i] + no2[:,i]
            no2[:,i] = no2[:,i] / np.sqrt(nExp/2.)

            fittedDepth[i] = self.lcf.fit_curve_fmin(sig[:,i],z,u)                      
        return wl, fittedDepth, no2




class LightcurveFitting():
    def __init__(self):
        self.modelMandelAgol =  pytransit.MandelAgol(eclipse=False)

    
    
    def fit_curve_fmin(self, lc,z, u):  
        # fits a Mandel & Agol lightcurve using simplex downhill minimisation
        #' returns the fitted Depth
        p_init = 0.01  # planet/star radius ratio initial estimate
        data_std = np.std(lc[0:len(lc)/4])+ 1e-8  # data variance estimate
#        data_std = 0.0000000000000001
        p_fit, err, out3, out4, out5 = fmin(self.chi_sq, p_init, args=(data_std, lc, z, u), disp=0, full_output=1)
        depth = p_fit[0]**2 
        return depth

    def chi_sq (self, p, sigma, lc, z, u): 
        # calculates the sum of squares or RESIDUAL = DATA - MODEL
        model = self.modelMandelAgol(z, p, u)
        chi_sq  = np.sum((lc - model) / sigma)**2
        return chi_sq



    def fit_curve_fmin(self, lc,z, u):  
        # fits a Mandel & Agol lightcurve using simplex downhill minimisation
        #' returns the fitted Depth
        p_init = 0.01  # planet/star radius ratio initial estimate
        data_std = np.std(lc[0:len(lc)/4])+ 1e-8  # data variance estimate
        p_fit, err, out3, out4, out5 = fmin(self.chi_sq, p_init, args=(data_std, lc, z, u), disp=0, full_output=1)
        depth = p_fit[0]**2 
        return depth

    def chi_sq (self, p, sigma, lc, z, u): 
        # calculates the sum of squares or RESIDUAL = DATA - MODEL
        model = self.modelMandelAgol(z, p, u)
        chi_sq  = np.sum((lc - model) / sigma)**2
        return chi_sq


class FitsDataProcessing():
    def bin_spec_data_R(self, fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff):
        R = float(R) ## If R is not float, this can lead to wrong results.
        hdul = pyfits.open(fitsFileName)
        wl0 = hdul[-3].data['Wavelength 1.0 um'][0]
        wl00 = hdul[-3].data['Wavelength 1.0 um']
        signal_stack, noise_stack, mask_x = self.extract_exosim_spectrum(fitsFileName, doCds, wl00, F, pixSize, qe, decorr, diff)
        signalBin, wl, bin_size_r = self.div_in2_wav_2(signal_stack, wl0, R, ld, pixSize)
        noiseBin , wl, bin_size_r = self.div_in2_wav_2(noise_stack,  wl0, R, ld, pixSize)
        
#        bin_size_r = [5]* (len(wl00)/5)
#        signalBin, wl = self.div_in2_wav(signal_stack, bin_size_r[0], wl00)
#        noiseBin , wln = self.div_in2_wav(noise_stack, bin_size_r[0], wl00)

        ct = 1
        q= 1
        pix = 0
        bin_size_pix = []
        print len(mask_x)

        for i in range (np.int(bin_size_r[0]/2),len(mask_x)):
            L = bin_size_r[q]
            
            if ct > L:
                q = q + 1
                ct = 1
                bin_size_pix.append(pix)
                pix = 0
                
            pix = pix + mask_x[i]
            ct = ct+1
        print "KKKK", len(bin_size_pix)
        hdul.close()
        
        print "HELLLLLLO", signalBin[0].sum()
        return wl, signalBin, noiseBin, bin_size_r, bin_size_pix

    def extract_exosim_spectrum(self, fitsFilePath, doCds, wl, F, pixSize, qe, decorr, diff):
        ## Assumes the fits file contains:
        ## 0th HDU: PrimaryHDU
        ## 1th-Nth: ImageHDUs of the same shape
        ## Perhaps a CREATOR:Exolib in fits header meta key
        ## would resolve potential issues
        hdul = pyfits.open(fitsFilePath)
        multiaccum = hdul[0].header['MACCUM']
        nExp = hdul[0].header['NEXP']
        
        if decorr == False:
        
            for i in range (2,2+ nExp*multiaccum):
                
                 a = hdul[i].data
              
                 background_1 = a[5:10,10:502]
                 background_2 = a[55:60,10:502]   
                 background = np.mean(np.vstack ( (background_1, background_2) ) )
                 a = a - background
                 
                 a = a/qe
                 
                 hdul[i].data = a       
                       
#        b = hdul[-1].data
#        lc = np.zeros((len(b[0]),len(b)))
#        for xx in range(len(b[0])):
#            lc[xx] = b['light curve %s'%(xx)]
#        
        sig_frame = hdul[1].data
#        
#        ndr_duration = []
#        for i in range(multiaccum): 
#            ndr_duration.append(hdul[1].header['T_NDR%s'%(i)])
#    
#        sig_tl = np.ones((sig_frame.shape[0],sig_frame.shape[1],nExp*multiaccum))
#        
#        
#        
        sig_noise_tl = np.ones((sig_frame.shape[0],sig_frame.shape[1],nExp*multiaccum))


        for i in range(0, nExp*multiaccum, multiaccum):
            for j in range(multiaccum):
#                sig_tl[...,i+j] = (sig_frame-DC)*ndr_duration[j]
                sig_noise_tl[...,i+j] = hdul[i+j+2].data

        noise_tl = sig_noise_tl


        noise = np.zeros((hdul[1].shape[0], hdul[1].shape[1], nExp))
        
                
        
        for i in range (0, nExp*multiaccum, multiaccum):
          
            if multiaccum>1:
                if doCds:
                    #if CDS being used.... 
                    cds_noise = noise_tl[...,i+multiaccum-1] - noise_tl[...,i]
                else:
                    #if CDS not used - use last NDR only
                    cds_noise = noise_tl[...,i+multiaccum-1]
        
            else:
                #use last NDR only
                 cds_noise = noise_tl[...,i+multiaccum-1]
        
            # perform background subtraction
            #background_1 = cds_sig[10:16,10:502]
            #background_2 = cds_sig[48:54,10:502]
            
#            background_1 = cds_sig[5:10,10:502]
#            background_2 = cds_sig[55:60,10:502]
#            
#            background = np.mean(np.vstack ( (background_1, background_2) ) )
#            cds_sig = cds_sig - background
#
#            background_1 = cds_noise[5:10,10:502]
#            background_2 = cds_noise[55:60,10:502]
#            
#            background = np.mean(np.vstack ( (background_1, background_2) ) )
#                                                
#            cds_noise = cds_noise - background
#            
            mask, mask_x = self.calc_signal_mask(cds_noise, wl, F, pixSize, diff)
            
            
            cds_noise = cds_noise * mask
#            
    
            # collect signal and noise into array
            
            
            noise[...,i*1.0/multiaccum] = cds_noise 
            
            signal = noise
            
            
        
        # extract spectrum (produces 1 D spectrum from 2 D array)
        signal_stack = np.zeros((nExp, hdul[1].shape[1]))
        noise_stack = np.zeros((nExp, hdul[1].shape[1]))
        for i in range(nExp):
            signal_stack[i] = self.extract_spectrum(signal[...,i]) 
            noise_stack[i] = self.extract_spectrum(noise[...,i])
        
        print "oot wlc signal", np.sum(signal_stack, axis = 1)[0]
        print "noise on wlc", np.std(np.sum(noise_stack, axis = 1))
        
        import matplotlib.pyplot as plt
        plt.figure('kkk')
        xx1 = []
        xx2 = []        
        for i in range (signal.shape[2]):
            xx1.append(np.sum(signal[...,i]))
            xx2.append(np.sum(noise[...,i]))
#        plt.plot(xx1)
        plt.plot(xx2)
        plt.plot(xx1)

        print np.std(xx2)

#            
##        plt.plot(np.sum(signal_stack, axis = 1))
##        plt.plot(np.sum(noise_stack, axis = 1), 'ro')
#        print (np.sum(signal_stack, axis = 1).max() - np.sum(signal_stack, axis = 1).min())  / np.sum(signal_stack, axis = 1).max() 
#        print (np.sum(noise_stack, axis = 1).max() - np.sum(noise_stack, axis = 1).min())  / np.sum(noise_stack, axis = 1).max() 
    
        hdul.close()


        return signal_stack, noise_stack, mask_x

    def extract_spectrum(self, f): # extracts 1 D spectrum from 2 D image
        # this can become a more sophisticated extraction in the future using Hornes' algorithm
        count = f.sum(axis = 0)
        return count
    
    
    def calc_signal_mask(self, inImage, wl,F, pixSize, diff):
        ## Requires Background Subtracted images (?) AP
        # apply mask    
        # 1) find max position of image
        if diff == 1:
            y_max = inImage.shape[0] /2.
        else:

            indices = np.where(inImage == inImage.max())            
            y_max = indices[0].max() # max excludes rare situation of having 2 or more indices

        # create and apply mask
        mask =  np.zeros(inImage.shape)
    
##        # variable mask
        mask_x = 2* 1.22 * wl * F / pixSize      
        mask_x = np.ceil(mask_x)
        mask_x = np.where(mask_x%2 != 0, mask_x, mask_x + 1)
##        print "MASKLIST", len(mask_x), mask_x[0], mask_x[-1], wl[0], wl[-1]

        for ii in range(mask.shape[1]):
            mask[y_max- mask_x[ii]//2 : 1 + y_max + mask_x[ii]//2, ii] = 1.
###        
        import matplotlib.pyplot as plt
        plt.figure('mask..')
#        plt.plot(mask[:,50],'ro')
        plt.imshow(mask, interpolation = 'none')
        
        return mask, mask_x
    
    def div_in2_wav_2(self, spec_stack, wl0, R, ld, pix_size): 
            # divides input 1 D array into bins based on R parameter
            x_pix = spec_stack.shape[1]
            m = ld[1]
            bin_size_0 = wl0 / (R* pix_size *m)   # initial bin size in pixel units
            bin_size =[] # bin size in pixel units list
            bin_pos = [] # bin edges position in pixel units list
            tot_pix = 0        
            # bin sizes calculated across array per R starting with bin_size_0
            i = 0
            while tot_pix < x_pix :
                if len(bin_size) == 0:
                    b = (bin_size_0)
                    i += 1
                else:
                    b = ((1+(1/R))* bin_size[i-1])
                    i += 1
                bin_size.append(b)
                tot_pix += b
            bin_size = np.array(bin_size) 
            bin_size_r = np.round(bin_size)        
            bin_pos0 = np.cumsum(bin_size_r) - np.round(bin_size_r[0]/2)
            bin_pos =  np.cumsum(bin_size_r) - np.round(bin_size_r[0]/2)  
            bin_size0 =[]
            for i in range (len(bin_pos)-1):
                bin_size0.append( bin_pos[i+1] - bin_pos[i] )
                
            x_pos = bin_pos*pix_size 
            # find dist of centre of each bin to assign wl
            ##x_cen = []
            x_cen =  np.zeros(len(x_pos)-1)
            for i in range (len(x_cen)):
                x_cen[i] =  (x_pos[i].item() + x_pos[i+1].item() )/ 2
            # find wl at centre of each bin
            wl = ld[0] + ld[1]*(x_cen-ld[2])      
            bin_pos = bin_pos[0:-1]
            bin_pos = bin_pos.tolist()
            # apply binning
            spec = np.add.reduceat(spec_stack, bin_pos, axis = 1)
            return spec, wl, bin_size_r
        
        
        
    def div_in2_wav(self, spec_stack, bin_size, wl):
        
        spec = np.add.reduceat(spec_stack, np.arange(spec_stack.shape[1])[::bin_size], axis = 1)
        wav = np.add.reduceat(wl, np.arange(len(wl))[::bin_size]) / bin_size
        if wav[-1] < wav [-2]:
            wav = wav[0:-2]
            spec = spec[:,0:-2]
            
        
        return spec, wav