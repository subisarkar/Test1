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
wl, signal, noise = snrTool.calc_exolib_SNR(fitsFileName, R, ld, pix_size)
wl, fittedDepth = snrTool.fit_lc_mandel_agol(fitsFileName, R, ld, pix_size, doCds)

where:
fitsFileName : pathe to fits file created by ExoSim
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
import matplotlib.pyplot as plt
import matplotlib
import gc

class SnrTool():
    def __init__(self):
        self.lcf = LightcurveFitting()
        self.fdp = FitsDataProcessing()
    def calc_exosim_SNR(self, fitsFileName, R, pixSize, doCds, F, qe, decorr, diff, wavrange, ch, apply_flat):
                
        print "processing..", fitsFileName      
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin, binsize= self.fdp.bin_spec_data_R(fitsFileName, R, pixSize, doCds, F, qe, decorr, diff, wavrange, ch, apply_flat)
            
        signal = np.zeros(signalBin.shape[1])
        noise = np.zeros(signalBin.shape[1])
        for i in range (signalBin.shape[1]):
            signal[i] = signalBin[:,i].mean()
            noise[i]  = np.std(signalBin[:,i])
  
        return wl, signal, noise
        
    def fit_lc_mandel_agol(self, fitsFileName, R, ld, pixSize, doCds):
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin, noiseBin, bin_pos = self.fdp.bin_spec_data_R(fitsFileName, R, ld, pixSize, doCds)
        hdul = pyfits.open(fitsFileName)
        ##z = hdul[-1].data['z'][1::2]
        multiaccum = hdul[0].header['MACCUM']
        z = hdul[-1].data['z'][multiaccum-1::multiaccum]
        hdul.close()
        u = [0,0] # modify by using u from fit files later
        sig = signalBin*1
        no2 = noiseBin*1
        fittedDepth =  np.zeros(signalBin.shape[1])
        for i in range (signalBin.shape[1]):
            no2[:,i] = no2[:,i] /  sig[:,i].max()
            sig[:,i] = sig[:,i] / sig[:,i].max()
            sig[:,i] = sig[:,i] + no2[:,i]
            fittedDepth[i] = self.lcf.fit_curve_fmin(sig[:,i],z,u)   
        
                   
        return wl, fittedDepth




class LightcurveFitting():
    def __init__(self):
        self.modelMandelAgol =  pytransit.MandelAgol(eclipse=False)

    
    
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
    def bin_spec_data_R(self, fitsFileName, R, pixSize, doCds, F, qe, decorr, diff, wavrange, ch, apply_flat):
        R = float(R) ## If R is not float, this can lead to wrong results.
        hdul = pyfits.open(fitsFileName, memmap=False)
        wl0 = hdul[-4].data['Wavelength 1.0 um']
        hdul.close()
        gc.collect()
 
        signal_stack = self.extract_exosim_spectrum(fitsFileName, doCds, wl0, F, pixSize, qe, decorr, diff,ch, apply_flat)
  
        signalBin, wl,binsize = self.div_in2_wav_prism(signal_stack, wl0, R, wavrange, ch)

        return wl, signalBin,binsize

    def extract_exosim_spectrum(self, fitsFileName, doCds, wl0, F, pixSize, qe, decorr, diff,ch, apply_flat):
        ## Assumes the fits file contains:
        ## 0th HDU: PrimaryHDU
        ## 1th-Nth: ImageHDUs of the same shape
        ## Perhaps a CREATOR:Exolib in fits header meta key
        ## would resolve potential issues
        hdul = pyfits.open(fitsFileName, memmap=False)
        multiaccum = hdul[0].header['MACCUM']
        nExp = hdul[0].header['NEXP']
        print "NEXP", nExp
        qe = qe[0:hdul[1].data.shape[0],0:hdul[1].data.shape[1]]

        if apply_flat==1:
            print "Applying flat field... in pipeline"
            for i in range (1,1+ nExp*multiaccum):
                
                 a = hdul[i].data
              
                 background_1 = a[5:10,10:-10]
                 background_2 = a[-10:-5,10:-10]   
                 background = np.mean(np.vstack ( (background_1, background_2) ) )
                 if diff == 0:
                     a = a - background
                 
                 a = a/qe
                 
                 hdul[i].data = a  
        else:
            print "NOT applying flat field... in pipeline"
     
        signal_stack = np.zeros((nExp, hdul[1].shape[1]))        
        
        if doCds:
            print "Applying CDS... in pipeline"
        else:
            print "Not applying CDS... in pipeline"  
        
      
        for i in range (nExp):
            print i
        
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
                 

            signal_stack[i] = self.calc_signal_mask2(cds_sig, wl0, F, pixSize, diff, ch)
     
        hdul.close()
        return signal_stack

    def extract_spectrum(self, f): # extracts 1 D spectrum from 2 D image
        # this can become a more sophisticated extraction in the future using Hornes' algorithm
        count = f.sum(axis = 0)
        return count
    
    def calc_signal_mask2(self, inImage, wl,F, pixSize, diff, ch):
        
        import scipy
        import matplotlib.pyplot as plt
        ## Requires Background Subtracted images (?) AP
        # apply mask    
        # 1) find max position of image
        fact = 100
        # oversampling factor
        inImage2 = np.repeat(inImage, fact, axis = 0)/fact
        #oversampled image
        
        if diff == 1:
            y_max = inImage.shape[0] /2.
            y_max2 = inImage2.shape[0] /2.
            Cen = y_max
            Cen2 = y_max*fact
        else:
            indices = np.where(inImage == inImage.max())            
            y_max = indices[0].max() # max excludes rare situation of having 2 or more indices
            x_max =indices[1].max()
#            print y_max, x_max
#            print inImage[y_max, x_max]
            ydata = inImage[:,x_max]
            xdata = np.arange(0,inImage.shape[0])
            xdata2 =np.arange(0,inImage.shape[0],0.1)
        
            fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
            errfunc  = lambda p, x, y: (y - fitfunc(p, x))
            init  = [inImage[y_max,x_max], y_max , 2.0]
            out   = scipy.optimize.leastsq(errfunc, init, args=(xdata, ydata))[0]
            ymodel = out[0]*np.exp(-0.5*((xdata2-out[1])/out[2])**2)
          
            Cen = out[1]
#            print Cen
            Cen2 = np.round(Cen*fact+fact/2.) # centre of mask in inImage2
 
        # create and apply mask
        mask =  np.zeros(inImage.shape)
        mask2 = np.zeros(inImage2.shape)
        
        # variable mask
        mask_x = 2* 1.22 * wl * F / pixSize      
        mask_x = np.ceil(mask_x)
        mask_x = np.where(mask_x%2 != 0, mask_x, mask_x + 1)

#        m = -2; c= 6
#        m = -1.111; c= 4.222 # Ch0
#        m = -0.5; c=3.9 # Ch1
    
        if ch =='AIRS CH0':
#            print 'Ch0 mask'
            start_wav = 1.95
            end_wav = 3.9
            extra_pix = 0.5
            
        if ch =='AIRS CH1':
#            print 'Ch1 mask'
            start_wav = 3.9 #Ch1
            end_wav = 7.8
            extra_pix = 0.5
            
        if ch =='NIR Spec':
#            print 'Nir spec mask'
            start_wav = 1.25
            end_wav = 1.9 
            extra_pix = 0.0
            
#        extra_pix = 0.5 # double for total width        
        m = (0-extra_pix)/(end_wav- start_wav)
        c = -end_wav*m
        
#        mask_x2 = (2* 1.22 * wl * F / pixSize)*fact   # causes higher spatial noise  
#        mask_x2 = (2* 2.23 * wl * F / pixSize)*fact     # second minimum        
        mask_x2 = 2*((m*wl + c ) + (1.22 * wl * F / pixSize))*fact   # to widen mask at blue end
    
        mask_x2 = np.ceil(mask_x2)
        mask_x2 = np.where(mask_x2%2 != 0, mask_x2, mask_x2 + 1) 
                
        Cen2 = int(Cen2)     
        for ii in range(mask2.shape[1]):
            mask2[Cen2- int(mask_x2[ii])/2 : 1 + Cen2 +int( mask_x2[ii])/2, ii] = 1.

 
        # apply mask

        inImage2 = inImage2 * mask2
 
        sig_ct2 = inImage2.sum(axis = 0)
 
        return sig_ct2

     
    def div_in2_wav_prism(self, arr, wl, R, wavrange, ch):
        
        x  = np.arange(len(wl))
        
        # oversampling factor to apply for smoother R
        fact = 100
        arr2 = np.repeat(arr,fact, axis = 1)/fact
        wl2 = np.repeat(wl,fact)
        
#        plt.figure(99)
#        plt.plot(x,wl)

        
#        ch = 'NIR Spec'
        
        NPT = 1 # nirspec type 0 = enzo, 1= piotr
#        NPT = 1
        if ch == 'NIR Spec' and NPT == 0:
#              use for Enxo design
            lim = np.argwhere(wl<0.6)[0].item()
            x0 = x[0:lim]
            wl0 = wl[0:lim]
            
            r = 7 # degree 7 works well
            z = np.polyfit(wl0, x0, r)
            zz = np.polyfit(x0, wl0, r)
              
        elif ch == 'NIR Spec' and NPT == 1:          
#            use for Piotr design 
            lim = 107
            lim0 = 60
            x0 = x[lim0:lim]
            wl0 = wl[lim0:lim]
                        
            r = 7 # degree 7 works well
            z = np.polyfit(wl0, x0, r)
            zz = np.polyfit(x0, wl0, r)
            
        else:
            r = 7 # degree 7 works well
            z = np.polyfit(wl, x, r)
            zz = np.polyfit(x, wl, r)       
        
 
        # bin sizes from wavmin to wavmax
        
        idx = np.argwhere(wl < wavrange[0])[0]
        if ch == 'MWIR':
            idx = [169]
        
        w0 = wl[idx][0]
        if ch == 'NIR Spec' and NPT ==1:
            w0 = wl0[-1] #adjust starting wavelength if Piotr design
        dw = w0/R
        
        
        binlist=[dw]
        for i in range(200):
            dw2 = (1+1/R)*dw
            binlist.append(dw2)
            dw = dw2
            if np.sum(binlist) > wavrange[1]-w0:
                break
        
        
        binlist[0] /= 2.
        binpos = np.cumsum(binlist)    
        wavpos = w0+binpos 
        #positions in wavelength of bin edges
        
        # convert to pixels    
        pixpos = 0
        for i in range (0,r+1):
            pixpos = pixpos + z[i]*wavpos**(r-i)    
        #pixpos = binedges in exact pixel units (but not whole pixels) 
        
        pixpos = pixpos[::-1]
        # reverse pix positions for prism : gives exact positions of bin edges in pixel units
        pixpos2 = pixpos*fact
        # exact pixel positions of expanded array

        BinSizes = pixpos[1:]-pixpos[:-1]
        # exact sizes of bins in pixel units
        BinSizes2 = pixpos2[1:]-pixpos2[:-1]

#        plt.figure(88)
#        plt.plot(BinSizes, 'bo-') 

#        PixStart = np.round(pixpos[0])
        PixStart = np.int(pixpos[0])+1
        PixEnd = np.round(pixpos[-1])
        PixTot = PixEnd-PixStart

        PixStart2 = np.int(pixpos2[0])+1
        PixEnd2 = np.round(pixpos2[-1])
        PixTot2 = PixEnd2-PixStart2        
        

        BinSizes0 = np.round(BinSizes)
        # bin sizes are rounded to nearest whole pixel
        BinSizes02 = np.round(BinSizes2)

        # adjust binsizes so that total pixel range is unchanged
#        if np.sum(BinSizes0) > PixTot:   
#            for i in range(200):
#                if np.sum(BinSizes0) > PixTot:
#                    BinSizes = BinSizes*0.9999
#                    BinSizes0 = np.round(BinSizes)
#                    print np.sum(BinSizes0)
#                else:
#                    break            
#        elif np.sum(BinSizes0) < PixTot:   
#            for i in range(200):
#                if np.sum(BinSizes0) < PixTot:
#                    BinSizes = BinSizes/0.9999
#                    BinSizes0 = np.round(BinSizes)
#                    print np.sum(BinSizes0)
#                else:
#                    break
#            
        if np.sum(BinSizes02) > PixTot2:   
            for i in range(200):
                if np.sum(BinSizes02) > PixTot2:
                    BinSizes2 = BinSizes2*0.99999
                    BinSizes02 = np.round(BinSizes2)
#                    print np.sum(BinSizes02)
                else:
                    break            
        elif np.sum(BinSizes02) < PixTot2:   
            for i in range(200):
                if np.sum(BinSizes02) < PixTot2:
                    BinSizes2 = BinSizes2*1.000001
                    BinSizes02 = np.round(BinSizes2)
#                    print np.sum(BinSizes02)
                else:
                    break               
        print PixTot2, np.sum(BinSizes02)
        
        # starting pixel for rounded bins = first bin edge in pixel units
        CumBinSizes = np.cumsum(BinSizes0)    
        pixpos = PixStart + CumBinSizes
        pixpos = np.hstack((PixStart,pixpos))
        # new bin edges based on rounded bins and rounded starting pixel (in whole pixel units)
        
        # starting pixel for rounded bins = first bin edge in pixel units
        CumBinSizes2 = np.cumsum(BinSizes02)    
        pixpos2 = PixStart2 + CumBinSizes2
        pixpos2 = np.hstack((PixStart2,pixpos2))
        # new bin edges based on rounded bins and rounded starting pixel (in whole pixel units)
  
#        plt.figure(89)
#        plt.plot(BinSizes0, 'ro-') 
#        plt.plot(BinSizes02/fact, 'go-') 
#        
  
        #positions in pixels of bin edges rounded to integers +1 (gives best rounding) to prep for binning by add.reduceat
#                
#        
#        # find wavelength at bin centre 
#        
#        cenpos = (pixpos[1:] + pixpos[:-1] )/ 2
#        cenpos = cenpos - 0.5
#        #positions in pixels of bin centre to find wavelength at bin centre
#
#        cenpos2 = (pixpos2[1:] + pixpos2[:-1] )/ 2
#        cenpos2 = cenpos2 - 0.5
#        cenpos2 /=fact
#        #positions in pixels of bin centre to find wavelength at bin centre
#
#        # find central wavelength of each bin
#        cenwav = 0
#        for i in range (0,r+1):
#            cenwav = cenwav + zz[i]*cenpos**(r-i) 
#
#        cenwav2 = 0
#        for i in range (0,r+1):
##            cenwav2 = cenwav2 + zz[i]*cenpos2**(r-i)
#
#        
#        plt.figure(99)
#        plt.plot(x,wl,'rx-')
#        plt.plot(cenpos2,cenwav2,'go-')
#        plt.plot(cenpos,cenwav,'bo')
        
        #safety
        if pixpos[-1] >= arr.shape[1]:
            pixpos[-1] = arr.shape[1]-1
        pixpos = np.where(pixpos >= 0, pixpos, 0.0)
        
        if pixpos2[-1] >= arr2.shape[1]:
            pixpos2[-1] = arr2.shape[1]-1
            
            
        pixpos2 = np.where(pixpos2 >= 0, pixpos2, 0.0)

        # binsizes for return
        binsize = pixpos[1:]- pixpos[:-1]          
        binsize2 = (pixpos2[1:]- pixpos2[:-1])/fact     
        binsize0 = (pixpos2[1:]- pixpos2[:-1])
           
        # r-bin
        bin_arr = np.add.reduceat(arr,pixpos.tolist(), axis = 1)[:,:-1]    
        bin_arr2 = np.add.reduceat(arr2,pixpos2.tolist(), axis = 1)[:,:-1]
        cenwav2 = np.add.reduceat(wl2,pixpos2.tolist())[:-1]/binsize0
   
         # one way to get cenpos
#        x2  = np.repeat(x,fact)
#        cenpos = np.add.reduceat(x2,pixpos2.tolist())[:-1]/binsize0 

         #this way tells us we have correctly assigned the wavelength to the bin centre: 
         #bin centre = (start pixel + binsize/2 ) -0.5
        cenpos =[]
        bb = PixStart2-1
 
        for i in range(len(cenwav2)):
            aa = bb + binsize0[i]/2
            cenpos.append(aa)
            bb = bb+binsize0[i]
        cenpos = np.array(cenpos)/fact
        cenpos -= 0.5  # subtract the 0.5 pixel to get the centre in exact pixel coords
 
   
#        return bin_arr, cenwav, binsize
        
        return bin_arr2, cenwav2, binsize2 
        
        
     
    def div_in2_wav(self, spec_stack, wl, bin_size):
        
        spec = np.add.reduceat(spec_stack, np.arange(spec_stack.shape[1])[::bin_size], axis = 1)
        wav = np.add.reduceat(wl, np.arange(len(wl))[::bin_size]) / bin_size
        if wav[-1] < wav [-2]:
            wav = wav[0:-2]
            spec = spec[:,0:-2]
            
        binsizes = np.array([bin_size]*spec.shape[1])
        
        return spec, wav, binsizes
