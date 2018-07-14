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
import scipy

class SnrTool():
    def __init__(self):
        self.lcf = LightcurveFitting()
        self.fdp = FitsDataProcessing()
    def calc_exosim_SNR(self, fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff,optimal,superbin, maskon):
        
      
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin = self.fdp.bin_spec_data_R(fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff,optimal, superbin, maskon)
        signal = np.zeros(signalBin.shape[1])
        noise = np.zeros(signalBin.shape[1])
        for i in range (signalBin.shape[1]):
#            signal[i] = (signalBin[:,i].max() - signalBin[:,i].min()) /signalBin[:,i].max()
#            noise[i]  =  (np.std(noiseBin[:,i])    / signalBin[:,i].max()   )/ np.sqrt(nExp/2.)
            signal[i] = signalBin[:,i].mean()
            noise[i]  = np.std(signalBin[:,i])
        return wl, signal, noise
    
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
    def bin_spec_data_R(self,fitsFileName, R, ld, pixSize, doCds, F, qe, decorr, diff, optimal, superbin, maskon):
        R = float(R) ## If R is not float, this can lead to wrong results.
        hdul = pyfits.open(fitsFileName)
        wl0 = hdul[-4].data['Wavelength 1.0 um'][0]
        wl00 = hdul[-4].data['Wavelength 1.0 um']
        signal_stack = self.extract_exosim_spectrum(fitsFileName, doCds, wl00, F, pixSize, qe, decorr, diff, optimal, maskon)
        if superbin == True:        
            signalBin, wl, bin_size_r = self.div_in2_wav_3(signal_stack, wl0, R, ld, pixSize)
            noiseBin, wl, bin_size_r = self.div_in2_wav_3(signal_stack, wl0, R, ld, pixSize)

        else:
            signalBin, wl = self.div_in2_wav_2(signal_stack, wl0, R, ld, pixSize)
            noiseBin, wl = self.div_in2_wav_2(signal_stack, wl0, R, ld, pixSize)

        hdul.close()
 
        return wl, signalBin

    def extract_exosim_spectrum(self, fitsFilePath, doCds, wl, F, pixSize, qe, decorr, diff, optimal, maskon):
        ## Assumes the fits file contains:
        ## 0th HDU: PrimaryHDU
        ## 1th-Nth: ImageHDUs of the same shape
        ## Perhaps a CREATOR:Exolib in fits header meta key
        ## would resolve potential issues
        hdul = pyfits.open(fitsFilePath)
        multiaccum = hdul[0].header['MACCUM']
        nExp = hdul[0].header['NEXP']
        qe = qe[0:hdul[1].data.shape[0],0:hdul[1].data.shape[1]]
        print "NEXP", nExp

        
        if decorr == False:
        
            for i in range (1,1+ nExp*multiaccum):
                
                 a = hdul[i].data
                 background_1 = a[5:10,10:-10]
                 background_2 = a[-10:-5,10:-10]   
                 background = np.mean(np.vstack ( (background_1, background_2) ) )
                 a = a - background
                 a = a/qe
                 hdul[i].data = a    
                 
        signal = np.zeros((hdul[1].shape[0], hdul[1].shape[1], nExp))
        for i in range (nExp):
            ##
            idx_1 =  i*multiaccum   + 1
            idx_2 =  idx_1 + (multiaccum - 1)
            ## FIXME: Its hould be "if doCDS" and if "multiaccum"
            
            if multiaccum>1:
                if doCds:
                    #if CDS being used.... 
                    cds_sig = hdul[idx_2].data - hdul[idx_1].data 
                else:
                    #if CDS not used - use last NDR only
                    cds_sig = hdul[idx_2].data 
            else:
                #use last NDR only
                 cds_sig = hdul[idx_2].data 

            
            mask = self.calc_signal_mask(cds_sig, wl, F, pixSize, diff)
            cds_sig = cds_sig * mask
            
            # collect signal and noise into array
            signal[...,i] = cds_sig
        
        # extract spectrum (produces 1 D spectrum from 2 D array)
        signal_stack = np.zeros((nExp, hdul[1].shape[1]))
        for i in range(signal_stack.shape[0]):
            signal_stack[i] = self.extract_spectrum(signal[...,i]) 
        hdul.close()
        return signal_stack

    def extract_spectrum(self, f): # extracts 1 D spectrum from 2 D image
        # this can become a more sophisticated extraction in the future using Hornes' algorithm
        count = f.sum(axis = 0)
        return count
    
    def extract_spectrum_optimal(self, f): # extracts 1 D spectrum from 2 D image
        # this can become a more sophisticated extraction in the future using Hornes' algorithm
        ybox = 30
        xbox = 256
        f  =  optimal_extract3.extraction_box(f,xbox,ybox)
 
        gain = 1
        readnoise = 22
        var = np.abs(f/gain + (readnoise/gain)**2) # instead read the variance
 
        count = optimal_extract3.optimalExtract(f,var,gain,readnoise,verbose=1,
                                    bsigma=5,psigma=15,csigma=15, 
                                    dispaxis=1, extract_radius=11 , 
                                    bord = 0, pord = 13, bkg_radii = [0,0],nreject = 25)[0]
        count = count.transpose()                           
    
                               
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
        
        return mask






    def div_in2_wav_3(self, spec_stack, wl0, R, ld, pix_size): 
        
            x_pix = spec_stack.shape[1]
            
            osf = 100.
            a = pix_size/osf
            x0 =  a/2.
            new_grid = np.arange(x0, x0 + a*osf*x_pix, a)
            old_grid = np.arange(pix_size/2.,  pix_size/2. + pix_size*x_pix, pix_size)
            
            print len(old_grid), len(new_grid)
            
            spec_stack_new = np.zeros((spec_stack.shape[0],spec_stack.shape[1]*osf))
            print spec_stack_new.shape
            
            
#            for i in range (spec_stack.shape[0]):
            y = spec_stack
            x = old_grid
            f = scipy.interpolate.interp1d(x, y, axis = 1, kind='cubic', bounds_error=False)
            ynew = f(new_grid)
            spec_stack_new = ynew
#            import matplotlib.pyplot as plt 
#            plt.figure('testttt')
#            plt.plot(old_grid, y[110], 'ro')
#            plt.plot(new_grid, ynew[110], 'bo')
#
#
#            xx
            # divides input 1 D array into bins based on R parameter
                        

#            bin_size_r = np.round(bin_size)        
            x_pix = spec_stack_new.shape[1]
            pix_size = pix_size/osf
            spec_stack = spec_stack_new / osf
        
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
            spec = np.add.reduceat(spec_stack_new, bin_pos, axis = 1)
            
            return spec, wl, bin_size_r

    
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

            return spec, wl
        
        
        
    def div_in2_wav(self, spec_stack, bin_size, wl):
        
        spec = np.add.reduceat(spec_stack, np.arange(spec_stack.shape[1])[::bin_size], axis = 1)
        wav = np.add.reduceat(wl, np.arange(len(wl))[::bin_size]) / bin_size
        if wav[-1] < wav [-2]:
            wav = wav[0:-2]
            spec = spec[:,0:-2]
            
        
        return spec, wav