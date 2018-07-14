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
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.optimize as opt

from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as ApPhot


class SnrTool():
    def __init__(self):
        self.fdp = FitsDataProcessing()
    def calc_exosim_SNR(self, fitsFileName, pixSize, doCds, F, qe, decorr, diff, ch, apply_flat):
                
        print "processing..", fitsFileName      
        signal, wl= self.fdp.extract_signal(fitsFileName, doCds, F, pixSize, qe, decorr, diff, apply_flat)
        noise = np.std(signal)
        signal_ = np.mean(signal)
        
        return wl, signal_, noise
 



class FitsDataProcessing():
 
    def extract_signal(self, fitsFileName, doCds, F, pixSize, qe, decorr, diff, apply_flat):
   
        hdul = pyfits.open(fitsFileName)
        multiaccum = hdul[0].header['MACCUM']
        nExp = hdul[0].header['NEXP']
        wl0 = hdul[-4].data['Wavelength 1.0 um'][0]       
        
        print "NEXP", nExp
        print "WL", wl0
        print "F", F
         
        qe = qe[0:hdul[1].data.shape[0],0:hdul[1].data.shape[1]]
        
    
        if apply_flat==1:
            for i in range (1,1+ nExp*multiaccum):                
                 a = hdul[i].data                 
                 background_1 = a[5:10,10:-10]
                 background_2 = a[-10:-5,10:-10]    
                 background = np.mean(np.vstack ( (background_1, background_2) ) )
                 if diff == 0:
                     a = a - background
                 a = a/qe      
                 hdul[i].data = a  
  

        
        signal = []
        subymaxlist = []
        subxmaxlist = []
        cds_sig_list =[]
  
        if doCds:
            print "Applying CDS... in pipeline"
        else:
            print "Not applying CDS... in pipeline"  
  
        for i in range (nExp):        
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
                 

            # for varying centering
            if decorr == True:
                Ap_Count = self.calc_aperture2(cds_sig, wl0, F, pixSize, diff)
                signal.append(Ap_Count)
            else:
                subymax, subxmax = self.calc_aperture3a(cds_sig, diff)
                subymaxlist.append(subymax)
                subxmaxlist.append(subxmax)
                cds_sig_list.append(cds_sig)   
        
        if decorr == False:
            print "Using fixed aperture method..."
            subymax = np.mean(subymaxlist)
            subxmax = np.mean(subxmaxlist)
            for i in range (0,nExp):
                cds_sig = cds_sig_list[i]
                subymax = subymaxlist[i]
                subxmax = subxmaxlist[i]            
                Ap_Count = self.calc_aperture3b(cds_sig, wl0, F, pixSize, subymax, subxmax)
                signal.append(Ap_Count)
        else:
            print "Using moving aperture method... will decorrelate"
   
                    
        hdul.close()
        
        return np.array(signal), wl0
        
        
 
  
    def calc_aperture2(self, inImage, wl, F, pixSize, diff): 
        
        if diff == 0:            
            b = inImage
            def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
                xo = float(xo)
                yo = float(yo)    
                a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
                g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                        + c*((y-yo)**2)))
                return g.ravel()
            
            # Create x and y indices
            x = np.linspace(0, b.shape[1]-1, b.shape[1])
            y = np.linspace(0, b.shape[0]-1, b.shape[0])
            x, y = np.meshgrid(x, y)
            
            data = b.ravel()
            Y,X =  np.unravel_index(b.argmax(),b.shape)
            
            initial_guess = (b.max(),X,Y,10,10,0,0)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data, p0=initial_guess)
            
            data_fitted = twoD_Gaussian((x, y), *popt)
            data_fitted = data_fitted.reshape(b.shape[0], b.shape[1])
            
#            fig, ax = plt.subplots(1, 1)
#            ax.hold(True)
#            ax.imshow(b, cmap=plt.cm.jet, origin='bottom',
#                extent=(x.min(), x.max(), y.min(), y.max()))
#            ax.contour(x, y, data_fitted, 8, colors='w')
#            plt.show()
            
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10. # 30 is too slow           
            y = np.arange(0,b.shape[0],1)
            y2 = np.arange(0,b.shape[0], 1/fact)
            x = np.arange(0,b.shape[1],1)
            x2 = np.arange(0,b.shape[0], 1/fact)
            
            f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
            data_fitted_new = f(x2,y2)
            
            Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
            X /= fact
            Y /= fact
            
            subymax = Y
            subxmax = X
                 
        
        elif diff == 1:
            subymax = inImage.shape[0]/2
            subxmax = inImage.shape[1]/2
                    
        # apply circular aperture centred at ymax, xmax
        
#        radius = 1.22 * wl * F / pixSize  
        radius = wl * F / pixSize  # F = F * position of 1st or 2nd minimum
        p = CircAp((subxmax,subymax), radius)
        Ap_Count = ApPhot(inImage,p)[0][3]

#        print "radius of aperture in pixels", radius                
##                
#        apertures = p
#        plt.figure('aperture')
#        plt.imshow(inImage, cmap='gray_r', origin='lower', interpolation = 'None')
#        apertures.plot(color='blue', lw=1.5, alpha=0.5)
#        xx
#        
        return Ap_Count
    
  
    def calc_aperture3a(self, inImage, diff): 
        
        if diff == 0:
            
            b = inImage
            def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
                xo = float(xo)
                yo = float(yo)    
                a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
                g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                        + c*((y-yo)**2)))
                return g.ravel()
            
            # Create x and y indices
            x = np.linspace(0, b.shape[1]-1, b.shape[1])
            y = np.linspace(0, b.shape[0]-1, b.shape[0])
            x, y = np.meshgrid(x, y)
            
            data = b.ravel()
            Y,X =  np.unravel_index(b.argmax(),b.shape)
            
            initial_guess = (b.max(),X,Y,10,10,0,0)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data, p0=initial_guess)
            
            data_fitted = twoD_Gaussian((x, y), *popt)
            data_fitted = data_fitted.reshape(b.shape[0], b.shape[1])
            
#            fig, ax = plt.subplots(1, 1)
#            ax.hold(True)
#            ax.imshow(b, cmap=plt.cm.jet, origin='bottom',
#                extent=(x.min(), x.max(), y.min(), y.max()))
#            ax.contour(x, y, data_fitted, 8, colors='w')
#            plt.show()
    
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10.           
            y = np.arange(0,b.shape[0],1)
            y2 = np.arange(0,b.shape[0], 1/fact)
            x = np.arange(0,b.shape[1],1)
            x2 = np.arange(0,b.shape[0], 1/fact)
            
            f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
            data_fitted_new = f(x2,y2)
            
            Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
            X /= fact
            Y /= fact
            
            subymax = Y
            subxmax = X
                 
        
        elif diff == 1:
            subymax = inImage.shape[0]/2
            subxmax = inImage.shape[1]/2
            
        return subymax, subxmax


    def calc_aperture3b(self, inImage, wl, F, pixSize, subymax, subxmax):         
        # apply circular aperture centred at ymax, xmax
        
#        radius = 1.22 * wl * F / pixSize  
        radius = wl * F / pixSize  # F = F * position of 1st or 2nd minimum
        p = CircAp((subxmax,subymax), radius)
        Ap_Count = ApPhot(inImage,p)[0][3]

#        print "radius of aperture in pixels", radius                
#                
#        apertures = p
#        plt.figure('aperture')
#        plt.imshow(inImage, cmap='gray_r', origin='lower', interpolation = 'None')
#        apertures.plot(color='blue', lw=1.5, alpha=0.5)
#        xx
#        
        return Ap_Count 
