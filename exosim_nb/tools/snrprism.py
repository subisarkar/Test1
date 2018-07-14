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

class SnrTool():
    def __init__(self):
        self.lcf = LightcurveFitting()
        self.fdp = FitsDataProcessing()
    def calc_exosim_SNR(self, fitsFileName, R, pixSize, doCds, F, qe, decorr, diff, wavrange, ch, ap_phot_mask):
                
        print "processing..", fitsFileName      
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin, binsize= self.fdp.bin_spec_data_R(fitsFileName, R, pixSize, doCds, F, qe, decorr, diff, wavrange, ch, ap_phot_mask)
        signal = np.zeros(signalBin.shape[1])
        noise = np.zeros(signalBin.shape[1])
        for i in range (signalBin.shape[1]):
#            signal[i] = (signalBin[:,i].max() - signalBin[:,i].min()) /signalBin[:,i].max()
#            noise[i]  =  (np.std(noiseBin[:,i])    / signalBin[:,i].max()   )/ np.sqrt(nExp/2.)
            signal[i] = signalBin[:,i].mean()
            noise[i]  = np.std(signalBin[:,i])
        return wl, signal, noise, binsize
    
    def fit_lc_mandel_agol(self, fitsFileName, R, ld, pixSize, doCds):
        R = float(R) ## If R is not float, this can lead to wrong results.
        wl, signalBin, noiseBin, bin_pos = self.fdp.bin_spec_data_R(fitsFileName, R, ld, pixSize, doCds)
        print signalBin, noiseBin
    
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
    def bin_spec_data_R(self, fitsFileName, R, pixSize, doCds, F, qe, decorr, diff, wavrange, ch, ap_phot_mask):
        R = float(R) ## If R is not float, this can lead to wrong results.
        hdul = pyfits.open(fitsFileName)
        wl0 = hdul[-5].data['Wavelength 1.0 um']
        signal_stack = self.extract_exosim_spectrum(fitsFileName, doCds, wl0, F, pixSize, qe, decorr, diff, ap_phot_mask)

        
#        import matplotlib.pyplot as plt
#        plt.figure(99)
#        qq = []
#        ww = []
#        for i in range (signal_stack.shape[1]):
#            qq.append(np.mean(signal_stack[:,i]))
#            ww.append(np.var(signal_stack[:,i]))
#        plt.plot(qq)
#        plt.plot(ww)

# (signal_stack[:,i].std())**2


        
        
        signalBin, wl,binsize = self.div_in2_wav_prism(signal_stack, wl0, R, wavrange, ch)
        hdul.close()
 
        return wl, signalBin,binsize

    def extract_exosim_spectrum(self, fitsFileName, doCds, wl0, F, pixSize, qe, decorr, diff, ap_phot_mask):
        ## Assumes the fits file contains:
        ## 0th HDU: PrimaryHDU
        ## 1th-Nth: ImageHDUs of the same shape
        ## Perhaps a CREATOR:Exolib in fits header meta key
        ## would resolve potential issues
        hdul = pyfits.open(fitsFileName)
        multiaccum = hdul[0].header['MACCUM']
        nExp = hdul[0].header['NEXP']
        print "NEXP", nExp
        qe = qe[0:hdul[1].data.shape[0],0:hdul[1].data.shape[1]]

        if decorr == False:
            for i in range (1,1+ nExp*multiaccum):
                
                 a = hdul[i].data
              
                 background_1 = a[5:10,10:-10]
                 background_2 = a[-10:-5,10:-10]   
                 background = np.mean(np.vstack ( (background_1, background_2) ) )
                 if diff == 0:
                     a = a - background
                 
                 a = a/qe
                 
                 hdul[i].data = a  
#     
#        
        signal = np.zeros((hdul[1].shape[0], hdul[1].shape[1], nExp))
        signal_stack = np.zeros((nExp, hdul[1].shape[1]))        
        
        signal_stack = np.zeros(((nExp/4)-10, hdul[1].shape[1]))        

        for i in range ((nExp/4)-10):
            
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
                 
            if ap_phot_mask == True:
                 signal_stack[i] = self.calc_signal_mask_ap(cds_sig, wl0, F, pixSize, diff)
                
            else:
                signal_stack[i] = self.calc_signal_mask2(cds_sig, wl0, F, pixSize, diff)
     
        hdul.close()
        return signal_stack

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
        
        return mask


    def calc_signal_mask2(self, inImage, wl,F, pixSize, diff):
        
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
            Cen2 = np.round(Cen*fact) # centre of mask in inImage2
            
        # create and apply mask
#        mask =  np.zeros(inImage.shape)
        mask2 = np.zeros(inImage2.shape)
        
    
        # variable mask
#        mask_x = 2* 1.22 * wl * F / pixSize      
#        mask_x = np.ceil(mask_x)
#        mask_x = np.where(mask_x%2 != 0, mask_x, mask_x + 1)
#        
#        
#        mask_x2 = (2* 1.22 * wl * F / pixSize)*fact   # causes higher spatial noise  
#        mask_x2 = (2* 2.23 * wl * F / pixSize)*fact     # second minimum 
#        m = -2; c= 6
        m = -1.111; c= 4.222 # Ch0
#        m = -0.5; c=3.9 # Ch1
        
        mask_x2 = ((m*wl + c ) + (2* 1.22 * wl * F / pixSize))*fact   # causes higher spatial noise  


        mask_x2 = np.ceil(mask_x2)
        mask_x2 = np.where(mask_x2%2 != 0, mask_x2, mask_x2 + 1) 
        
#        print mask_x, mask_x2/fact
#        plt.figure('mask width')
#        plt.plot(mask_x)
#        plt.plot(mask_x2/fact)
#        
#        for ii in range(mask.shape[1]):
#            mask[y_max- mask_x[ii]//2 : 1 + y_max + mask_x[ii]//2, ii] = 1.
#            
        for ii in range(mask2.shape[1]):
            mask2[Cen2- mask_x2[ii]//2 : 1 + Cen2 + mask_x2[ii]//2, ii] = 1.
###        
#        import matplotlib.pyplot as plt
#        plt.figure('mask..')
##        plt.plot(mask[:,50],'ro')
#        plt.imshow(mask, interpolation = 'none')
#
        
#        plt.figure('mask2')
#        plt.imshow(mask2, interpolation = 'none')

        # apply mask
#        inImage = inImage * mask
        inImage2 = inImage2 * mask2
        
        # extract spectrum
#        sig_ct = inImage.sum(axis = 0)
        sig_ct2 = inImage2.sum(axis = 0)
        
#        plt.figure('sig_ct')
#        plt.plot(sig_ct)
#        plt.plot(sig_ct2, 'r-')        
                
        return sig_ct2



    def calc_signal_mask_ap(self, sig, wl,F, pixSize, diff):
        
        from photutils import RectangularAperture as RectAp
        from photutils import aperture_photometry as ApPhot
        import scipy
        import matplotlib.pyplot as plt
        print sig.shape

        if diff == 1:
            y_max = inImage.shape[0] /2. 
            Cen = y_max
        # 1) find max position of image
            
        
        else:
            indices = np.where(sig == sig.max())            
            y_max = indices[0].max() # max excludes rare situation of having 2 or more indices
            x_max =indices[1].max()
            print y_max, x_max
            print sig[y_max, x_max]
            ydata = sig[:,x_max]
            xdata = np.arange(0,sig.shape[0])
            xdata2 =np.arange(0,sig.shape[0],0.1)
        
            fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
            errfunc  = lambda p, x, y: (y - fitfunc(p, x))
            init  = [sig[y_max,x_max], y_max , 2.0]
            out   = scipy.optimize.leastsq(errfunc, init, args=(xdata, ydata))[0]
            ymodel = out[0]*np.exp(-0.5*((xdata2-out[1])/out[2])**2)
          
            Cen = out[1]
            
            Cen2 = y_max
#            plt.figure(7777)
#            plt.plot(xdata,ydata)
#            plt.plot(xdata2,ymodel)
#            
         

        # create and apply mask  
##        # variable mask
        mask_x = 2* 1.22 * wl * F / pixSize      

        sig_ct = []
        no_ct = []
#        for hh in range(0,sig.shape[1]):
        for hh in range(0,100):

            p = RectAp((hh,Cen),1,mask_x[hh],0)
            SS = ApPhot(sig,p, method='subpixel', subpixels=32)[0][0]
            sig_ct.append(SS)
#        print sig_ct, len(sig_ct)
        plt.figure('sig_ct')
        plt.plot(sig_ct)
        plt.figure('cen2')
        plt.plot(cenn2)
  
        
        sig_ct = np.sum(sig, axis = 0)
#        no_ct = np.sum(no, axis = 0)

        return sig_ct



       
        
    def div_in2_wav_prism(self, arr, wl, R, wavrange, ch):
        
        x  = np.arange(len(wl))
        
        # oversampling factor to apply for smoother R
        fact = 100
        arr2 = np.repeat(arr,fact, axis = 1)/fact

                
#        ch = 'NIR Spec'
         
        if ch == 'NIR Spec':
            
            lim = np.argwhere(wl<0.6)[0].item()
            print "limit on interp", lim
            x0 = x[0:lim]
            wl0 = wl[0:lim] 
            r = 7 # degree 7 works well
            z = np.polyfit(wl0, x0, r)
            zz = np.polyfit(x0, wl0, r)
        
        elif ch == 'MWIR':
            print wavrange
            lim0 = 180
            lim1 = 40
            print "limit on interp", lim0,lim1
            x0 = x[lim1:lim0]
            wl0 = wl[lim1:lim0] 
        #            plt.figure(99)
        #            plt.plot(x0,wl0,'bo')
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
        dw = w0/R
        
        print w0, dw
        
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


        plt.figure(88)
        plt.plot(BinSizes, 'bo-') 
  


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
  
        plt.figure(88)
        plt.plot(BinSizes0, 'ro-') 
        plt.plot(BinSizes02/fact, 'go-') 
  
        #positions in pixels of bin edges rounded to integers +1 (gives best rounding) to prep for binning by add.reduceat
                
        
        # find wavelength at bin centre 
        
        cenpos = (pixpos[1:] + pixpos[:-1] )/ 2
        cenpos = cenpos - 0.5
        #positions in pixels of bin centre to find wavelength at bin centre

        cenpos2 = (pixpos2[1:] + pixpos2[:-1] )/ 2
        cenpos2 = cenpos2 - 0.5
        cenpos2 /=fact
        #positions in pixels of bin centre to find wavelength at bin centre

        # find central wavelength of each bin
        cenwav = 0
        for i in range (0,r+1):
            cenwav = cenwav + zz[i]*cenpos**(r-i) 

        cenwav2 = 0
        for i in range (0,r+1):
            cenwav2 = cenwav2 + zz[i]*cenpos2**(r-i)

        
        plt.figure(99)
        plt.plot(x,wl,'rx-')
        plt.plot(cenpos2,cenwav2,'go-')
        plt.plot(cenpos,cenwav,'bo')
        
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

        # r-bin
        bin_arr = np.add.reduceat(arr,pixpos.tolist(), axis = 1)[:,:-1]    
        bin_arr2 = np.add.reduceat(arr2,pixpos2.tolist(), axis = 1)[:,:-1]
        
        print bin_arr[0]
        print bin_arr2[0]
        
        plt.figure(999)
        plt.plot(bin_arr[0])
        plt.plot(bin_arr2[0])

#        return bin_arr, cenwav, binsize
        
        return bin_arr2, cenwav2, binsize2
