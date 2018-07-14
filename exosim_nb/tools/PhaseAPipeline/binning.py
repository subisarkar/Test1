# -*- coding: utf-8 -*-
"""
Created on Tue May 15 16:16:45 2018

@author: c1341133

Masking and binning

"""

import numpy as np
from scipy import optimize
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from scipy import interpolate
import scipy.optimize as opt
from photutils import CircularAperture as CircAp
from photutils import aperture_photometry as ApPhot
import sys  # for progress bar (sys.stdout)


def spectroscopic_mask(data, info, ICF, ch, diff):
     
    wl = info['WL']
    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            pixSize = np.float(child.find('detector_pixel').find('pixel_size').get('val'))
            F = np.float(child.find('wfno').get('val'))

    fact = 100    # oversampling factor
    data2 = np.zeros((data.shape[0]*100,data.shape[1],data.shape[2]))
    
    if ch =='AIRS CH0':
        start_wav = 1.95
        end_wav = 3.9
        extra_pix = 0.5
        
    if ch =='AIRS CH1':
        start_wav = 3.9  
        end_wav = 7.8
        extra_pix = 0.5
        
    if ch =='NIR Spec':
        start_wav = 1.25
        end_wav = 1.9 
        extra_pix = 0.0
        F = F*2 # since PSF simulated using 2*F to account for WFE    


    nWLs = data.shape[2]  # how many steps in the loop

    # Progress Bar setup:
    ProgMax = 100    # number of dots in progress bar
    if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
    print "|" +    ProgMax*"-"    + "|     MyFunction() progress"
    sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
    nProg = 0   # fraction of progress    
    

    for i in range(data.shape[2]):
        
        if ( i >= nProg*nWLs/ProgMax ):
                '''Print dot at some fraction of the loop.'''
                sys.stdout.write('*'); sys.stdout.flush();  
                nProg = nProg+1
        if ( i >= nWLs-1 ):
                '''If done, write the end of the progress bar'''
                sys.stdout.write('|     done  \n'); sys.stdout.flush(); 


        inImage = data[...,i]
         # 1) find max position of image
        inImage2 = np.repeat(inImage, fact, axis = 0)/fact #oversampled image
        indices = np.where(inImage == inImage.max())  
        y_max = indices[0].max() # max excludes rare situation of having 2 or more indices
        x_max =indices[1].max()
        ydata = inImage[:,x_max]
        xdata = np.arange(0,inImage.shape[0])
        xdata2 =np.arange(0,inImage.shape[0],0.1)
        
        if diff ==0:
        
            fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
            errfunc  = lambda p, x, y: (y - fitfunc(p, x))
            init  = [inImage[y_max,x_max], y_max , 2.0]
            out   = optimize.leastsq(errfunc, init, args=(xdata, ydata))[0]
            ymodel = out[0]*np.exp(-0.5*((xdata2-out[1])/out[2])**2) 
            Cen = out[1]
            Cen2 = np.int(np.round(Cen*fact+fact/2.)) # centre of mask in inImage2
        elif diff ==1:
            y_max = inImage.shape[0] /2.
            y_max2 = inImage2.shape[0] /2.
            Cen = y_max
            Cen2 = y_max*fact
            
        # create and apply mask
        mask2 = np.zeros(inImage2.shape)

        # adjustments for blue end widening to reduce spatial jitter noise

      
        m = (0-extra_pix)/(end_wav- start_wav)
        c = -end_wav*m
        
#        m=0;c=0  #comment out to apply blue-end widening
        
#        mask_x2 = (2* 1.22 * wl * F / pixSize)*fact   # 1st min  higher spatial noise               
#        mask_x2 = (2* 2.23 * wl * F / pixSize)*fact   # second minimum        
                   
        mask_x2 = 2*(((m*wl + c ) + (1.22 * wl * F) / pixSize))*fact   # to widen mask at blue end
        
#        mask_x2 = (2* 1.22 * wl * F / pixSize)*fact   # causes higher spatial noise  
  
        
        mask_x2 = np.ceil(mask_x2)
        mask_x2 = np.where(mask_x2%2 != 0, mask_x2, mask_x2 + 1) 
                
        Cen2 = int(Cen2)   
        for ii in range(mask2.shape[1]):
            mask2[Cen2- int(mask_x2[ii])/2 : 1 + Cen2 +int( mask_x2[ii])/2, ii] = 1.

        # apply mask
        inImage2 = inImage2 * mask2
        
        data2[...,i] = inImage2
        
        
    return data2 #expanded array
    
    
def extract1D(data2):

    spectra = np.zeros((data2.shape[2],data2.shape[1]))  
    
    # extract spectrum
    for i in range(data2.shape[2]):
        
        print 'extracting spectrum', i
    
        spectra[i] = data2[...,i].sum(axis = 0)

    return spectra
    
    
def spectroscopic_mask_plus_extraction(data, info, ICF, ch, diff):
     
    wl = info['WL']
    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            pixSize = np.float(child.find('detector_pixel').find('pixel_size').get('val'))
            F = np.float(child.find('wfno').get('val'))

    fact = 100    # oversampling factor
    
    if ch =='AIRS CH0':
        start_wav = 1.95
        end_wav = 3.9
        extra_pix = 0.5
        
    if ch =='AIRS CH1':
        start_wav = 3.9  
        end_wav = 7.8
        extra_pix = 0.5
        
    if ch =='NIR Spec':
        start_wav = 1.25
        end_wav = 1.9 
        extra_pix = 0.0
        F = F*2 # since PSF simulated using 2*F to account for WFE    

    spectra = np.zeros((data.shape[2],data.shape[1]))  
    
    nWLs = data.shape[2]  # how many steps in the loop

    # Progress Bar setup:
    ProgMax = 100    # number of dots in progress bar
    if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
    print "|" +    ProgMax*"-"    + "|     MyFunction() progress"
    sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
    nProg = 0   # fraction of progress    
    

    for i in range(data.shape[2]):
        
        if ( i >= nProg*nWLs/ProgMax ):
                '''Print dot at some fraction of the loop.'''
                sys.stdout.write('*'); sys.stdout.flush();  
                nProg = nProg+1
        if ( i >= nWLs-1 ):
                '''If done, write the end of the progress bar'''
                sys.stdout.write('|     done  \n'); sys.stdout.flush(); 


        inImage = data[...,i]
         # 1) find max position of image
        inImage2 = np.repeat(inImage, fact, axis = 0)/fact #oversampled image
        indices = np.where(inImage == inImage.max())  
        y_max = indices[0].max() # max excludes rare situation of having 2 or more indices
        x_max =indices[1].max()
        ydata = inImage[:,x_max]
        xdata = np.arange(0,inImage.shape[0])
        xdata2 =np.arange(0,inImage.shape[0],0.1)
        
        if diff ==0:
        
            fitfunc  = lambda p, x: p[0]*np.exp(-0.5*((x-p[1])/p[2])**2)
            errfunc  = lambda p, x, y: (y - fitfunc(p, x))
            init  = [inImage[y_max,x_max], y_max , 2.0]
            out   = optimize.leastsq(errfunc, init, args=(xdata, ydata))[0]
            ymodel = out[0]*np.exp(-0.5*((xdata2-out[1])/out[2])**2) 
            Cen = out[1]
            Cen2 = np.int(np.round(Cen*fact+fact/2.)) # centre of mask in inImage2
        elif diff ==1:
            y_max = inImage.shape[0] /2.
            y_max2 = inImage2.shape[0] /2.
            Cen = y_max
            Cen2 = y_max*fact
            
        # create and apply mask
        mask2 = np.zeros(inImage2.shape)

        # adjustments for blue end widening to reduce spatial jitter noise

      
        m = (0-extra_pix)/(end_wav- start_wav)
        c = -end_wav*m
        
#        m=0;c=0  #comment out to apply blue-end widening
        
#        mask_x2 = (2* 1.22 * wl * F / pixSize)*fact   # 1st min  higher spatial noise               
#        mask_x2 = (2* 2.23 * wl * F / pixSize)*fact   # second minimum        
                   
        mask_x2 = 2*(((m*wl + c ) + (1.22 * wl * F) / pixSize))*fact   # to widen mask at blue end
        
#        mask_x2 = (2* 1.22 * wl * F / pixSize)*fact   # causes higher spatial noise  
  
        
        mask_x2 = np.ceil(mask_x2)
        mask_x2 = np.where(mask_x2%2 != 0, mask_x2, mask_x2 + 1) 
                
        Cen2 = int(Cen2)   
        for ii in range(mask2.shape[1]):
            mask2[Cen2- int(mask_x2[ii])/2 : 1 + Cen2 +int( mask_x2[ii])/2, ii] = 1.

        # apply mask
        inImage2 = inImage2 * mask2
               
        spectra[i] = inImage2.sum(axis = 0)
 
    return spectra        
 
    
    
    
    


def binSpectrum(spectra, info, ch, binning, binSize):
    
    wl = info['WL']
    
    if ch == 'AIRS CH0':
        R = 100.
        wavrange = [1.8,3.89]
        
    elif ch == 'AIRS CH1':
        R = 30.
        wavrange = [3.80,8.0]
   
    if ch == 'NIR Spec':
        R = 20.0                                     
        wavrange = [0.9,2.2]
     
    if binning == 'R-bin':
            
        # oversampling factor to apply for smoother R
        fact = 100
        x  = np.arange(len(wl))
        spectra2 = np.repeat(spectra,fact, axis = 1)/fact
        wl2 = np.repeat(wl,fact)
   
        if ch == 'NIR Spec' :          
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
        w0 = wl[idx][0]
        if ch == 'NIR Spec':
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
 
        PixStart = np.int(pixpos[0])+1
        PixEnd = np.round(pixpos[-1])
        PixTot = PixEnd-PixStart

        PixStart2 = np.int(pixpos2[0])+1
        PixEnd2 = np.round(pixpos2[-1])
        PixTot2 = PixEnd2-PixStart2        
        

        BinSizes0 = np.round(BinSizes)
        # bin sizes are rounded to nearest whole pixel
        BinSizes02 = np.round(BinSizes2)
 
        if np.sum(BinSizes02) > PixTot2:   
            for i in range(200):
                if np.sum(BinSizes02) > PixTot2:
                    BinSizes2 = BinSizes2*0.99999
                    BinSizes02 = np.round(BinSizes2)
                else:
                    break            
        elif np.sum(BinSizes02) < PixTot2:   
            for i in range(200):
                if np.sum(BinSizes02) < PixTot2:
                    BinSizes2 = BinSizes2*1.000001
                    BinSizes02 = np.round(BinSizes2)
                else:
                    break               
        
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
 
        #safety
        if pixpos[-1] >= spectra.shape[1]:
            pixpos[-1] = spectra.shape[1]-1
        pixpos = np.where(pixpos >= 0, pixpos, 0.0)
        
        if pixpos2[-1] >= spectra2.shape[1]:
            pixpos2[-1] = spectra2.shape[1]-1
            
            
        pixpos2 = np.where(pixpos2 >= 0, pixpos2, 0.0)

        # binsizes for return
        binsize = pixpos[1:]- pixpos[:-1]          
        binsize2 = (pixpos2[1:]- pixpos2[:-1])/fact     
        binsize0 = (pixpos2[1:]- pixpos2[:-1])
  
        # r-bin
        bin_spectra = np.add.reduceat(spectra,pixpos.tolist(), axis = 1)[:,:-1]    
        bin_spectra2 = np.add.reduceat(spectra2,pixpos2.tolist(), axis = 1)[:,:-1]
        cenwav2 = np.add.reduceat(wl2,pixpos2.tolist())[:-1]/binsize0
 
        cenpos =[]
        bb = PixStart2-1
 
        for i in range(len(cenwav2)):
            aa = bb + binsize0[i]/2
            cenpos.append(aa)
            bb = bb+binsize0[i]
        cenpos = np.array(cenpos)/fact
        cenpos -= 0.5  # subtract the 0.5 pixel to get the centre in exact pixel coords
         
        return bin_spectra2, cenwav2, binsize2
        
    elif binning =='fixed-bin': # to bin to a fixed bin size given by number of pixels

        spec = np.add.reduceat(spectra, np.arange(spectra.shape[1])[::binSize], axis = 1)
        wav = np.add.reduceat(wl, np.arange(len(wl))[::binSize]) / binSize
        if wav[-1] < wav [-2]:
            wav = wav[0:-2]
            spec = spec[:,0:-2]
            
        binsizes = np.array([binSize]*spec.shape[1])

        return spec, wav, binsizes
        


def photometric_mask(data, info, ICF, ch, diff, incAp, pl):
    
    wl = info['WL'][0]

    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            pixSize = np.float(child.find('detector_pixel').find('pixel_size').get('val'))
            F_ = np.float(child.find('wfno').get('val'))  
   
    if ch == 'FGS Red':
        Ap = 5.0
#        Ap = 7.910455
        if incAp ==1 :
            Ap = Ap*1.5                 
        F = F_*Ap                
    elif ch == 'FGS Prime':
        Ap = 7.0
        if incAp ==1 :
            Ap = Ap*1.5         
        
        F =  F_*Ap
    elif ch == 'NIR Phot':
        Ap = 9.0
        if pl == 'GJ 1214 b':
            Ap=8.0    
        if incAp ==1 :
            Ap = Ap*1.5 
        
        F = F_*Ap
    else:
        F = 1.22*F_
       
    print "aperture of radius", Ap, "x F lambda, where F = ", F_, "and lambda =", wl
    
    signal = []
    subymaxlist = []
    subxmaxlist = []  
    
    import sys  # for progress bar (sys.stdout)
    
    nWLs = data.shape[2]  # how many steps in the loop

    # Progress Bar setup:
    ProgMax = 100    # number of dots in progress bar
    if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
    print "|" +    ProgMax*"-"    + "|     MyFunction() progress"
    sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
    nProg = 0   # fraction of progress    
    

    for i in range(data.shape[2]):
        
        if ( i >= nProg*nWLs/ProgMax ):
                '''Print dot at some fraction of the loop.'''
                sys.stdout.write('*'); sys.stdout.flush();  
                nProg = nProg+1
        if ( i >= nWLs-1 ):
                '''If done, write the end of the progress bar'''
                sys.stdout.write('|     done  \n'); sys.stdout.flush(); 

        inImage = data[...,i]  
              
        if diff == 0:
            
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

            x = np.linspace(0, inImage.shape[1]-1, inImage.shape[1])
            y = np.linspace(0, inImage.shape[0]-1, inImage.shape[0])
            x, y = np.meshgrid(x, y)
            
            data_ravel = inImage.ravel()
            Y,X =  np.unravel_index(inImage.argmax(),inImage.shape)
            
            initial_guess = (inImage.max(),X,Y,10,10,0,0)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_ravel, p0=initial_guess)
            
#            print popt
            
            data_fitted = twoD_Gaussian((x, y), *popt)
            data_fitted = data_fitted.reshape(inImage.shape[0], inImage.shape[1])
            
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10. # 30 is too slow           
            y = np.arange(0,inImage.shape[0],1)
            y2 = np.arange(0,inImage.shape[0], 1/fact)
            x = np.arange(0,inImage.shape[1],1)
            x2 = np.arange(0,inImage.shape[0], 1/fact)
            
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
        signal.append(Ap_Count)
        
    return np.array(signal),wl
    
 




def use_pointing(data, info, ICF, ch, diff, pointing_x, pointing_y):
    
    wl = info['WL'][0]

    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            pixSize = np.float(child.find('detector_pixel').find('pixel_size').get('val'))
            F_ = np.float(child.find('wfno').get('val'))  
   
    if ch == 'FGS Red':
        Ap = 5.0
#        Ap = 7.910455
        
        
        F = F_*Ap                
    elif ch == 'FGS Prime':
        Ap = 7.0
        F =  F_*Ap
    elif ch == 'NIR Phot':
        Ap = 9.0        
        
        F = F_*Ap
    else:
        F = 1.22*F_
       
    print "aperture of radius", Ap, "x F lambda, where F = ", F_, "and lambda =", wl
    
    signal = []
    subymaxlist = []
    subxmaxlist = []  
    
    import sys  # for progress bar (sys.stdout)
    
    nWLs = data.shape[2]  # how many steps in the loop

    # Progress Bar setup:
    ProgMax = 100    # number of dots in progress bar
    if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
    print "|" +    ProgMax*"-"    + "|     MyFunction() progress"
    sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
    nProg = 0   # fraction of progress    
    

    for i in range(1):
        
        if ( i >= nProg*nWLs/ProgMax ):
                '''Print dot at some fraction of the loop.'''
                sys.stdout.write('*'); sys.stdout.flush();  
                nProg = nProg+1
        if ( i >= nWLs-1 ):
                '''If done, write the end of the progress bar'''
                sys.stdout.write('|     done  \n'); sys.stdout.flush(); 

        inImage = data[...,i]  
              
        if diff == 0:
            
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

            x = np.linspace(0, inImage.shape[1]-1, inImage.shape[1])
            y = np.linspace(0, inImage.shape[0]-1, inImage.shape[0])
            x, y = np.meshgrid(x, y)
            
            data_ravel = inImage.ravel()
            Y,X =  np.unravel_index(inImage.argmax(),inImage.shape)
            
            initial_guess = (inImage.max(),X,Y,10,10,0,0)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_ravel, p0=initial_guess)
            
#            print popt
            
            data_fitted = twoD_Gaussian((x, y), *popt)
            data_fitted = data_fitted.reshape(inImage.shape[0], inImage.shape[1])
            
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10. # 30 is too slow           
            y = np.arange(0,inImage.shape[0],1)
            y2 = np.arange(0,inImage.shape[0], 1/fact)
            x = np.arange(0,inImage.shape[1],1)
            x2 = np.arange(0,inImage.shape[0], 1/fact)
            
            f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
            data_fitted_new = f(x2,y2)
            
            Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
            X /= fact
            Y /= fact
            
            subymax = Y
            subxmax = X
            
#            print subymax, subxmax
                 
        
        elif diff == 1:
            subymax = inImage.shape[0]/2
            subxmax = inImage.shape[1]/2

    for i in range(data.shape[2]):
        # apply circular aperture centred at ymax, xmax      
        radius = wl * F / pixSize  # F = F * position of 1st or 2nd minimum
        p = CircAp((subxmax + pointing_x[i],subymax +pointing_y[i]), radius)
        Ap_Count = ApPhot(inImage,p)[0][3]
        
        
        print subxmax + pointing_x[i],subymax +pointing_y[i]

        print "radius of aperture in pixels", radius                
#                
        apertures = p
        plt.figure('aperture%s'%(i))
        plt.imshow(inImage, cmap='gray_r', origin='lower', interpolation = 'None')
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
       
    
        signal.append(Ap_Count)
        
    return np.array(signal),wl




   
        
        
def obtain_average_Gaussian(data, info, ICF, ch, diff):
    
    sys.stdout.write('obtaining average Guassian...')
    
    wl = info['WL'][0]

    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            pixSize = np.float(child.find('detector_pixel').find('pixel_size').get('val'))
            F_ = np.float(child.find('wfno').get('val'))  
   
    if ch == 'FGS Red':
        Ap = 5.0
        F = F_*Ap                
    elif ch == 'FGS Prime':
        Ap = 7.0
        F =  F_*Ap
    elif ch == 'NIR Phot':
        Ap = 9.0        
        
        F = F_*Ap
    else:
        F = 1.22*F_
       
    print "aperture of radius", Ap, "x F lambda, where F = ", F_, "and lambda =", wl
    
    signal = []
    subymaxlist = []
    subxmaxlist = []  
    
    
    nWLs = data.shape[2]  # how many steps in the loop

    # Progress Bar setup:
    ProgMax = 100    # number of dots in progress bar
    if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
    print "|" +    ProgMax*"-"    + "|     MyFunction() progress"
    sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
    nProg = 0   # fraction of progress    
    

    for i in range(data.shape[2]):
        
        if ( i >= nProg*nWLs/ProgMax ):
                '''Print dot at some fraction of the loop.'''
                sys.stdout.write('*'); sys.stdout.flush();  
                nProg = nProg+1
        if ( i >= nWLs-1 ):
                '''If done, write the end of the progress bar'''
                sys.stdout.write('|     done  \n'); sys.stdout.flush(); 

        inImage = data[...,i]  
              
        if diff == 0:
            
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

            x = np.linspace(0, inImage.shape[1]-1, inImage.shape[1])
            y = np.linspace(0, inImage.shape[0]-1, inImage.shape[0])
            x, y = np.meshgrid(x, y)
            
            data_ravel = inImage.ravel()
            Y,X =  np.unravel_index(inImage.argmax(),inImage.shape)
            
            initial_guess = (inImage.max(),X,Y,10,10,0,0)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_ravel, p0=initial_guess)
            
#            print popt
            
            data_fitted = twoD_Gaussian((x, y), *popt)
            data_fitted = data_fitted.reshape(inImage.shape[0], inImage.shape[1])
            
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10. # 30 is too slow           
            y = np.arange(0,inImage.shape[0],1)
            y2 = np.arange(0,inImage.shape[0], 1/fact)
            x = np.arange(0,inImage.shape[1],1)
            x2 = np.arange(0,inImage.shape[0], 1/fact)
            
            f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
            data_fitted_new = f(x2,y2)
            
            Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
            X /= fact
            Y /= fact
            
            subymax = Y
            subxmax = X
            
#            print subymax, subxmax
            
            if i == 0:
                d_stack =  popt
            else:
                d_stack = np.vstack((d_stack, popt))
                
            
    d_av = d_stack.mean(axis =0 )
                 
       
        
    return d_av
        

def use_average_Gaussian(data, info, ICF, ch, diff, d_av):
    
    sys.stdout.write('applying average Guassian...')
    
    sigma_x = d_av[3]
    sigma_y = d_av[4]
    theta = d_av[5]
    
    offset_init = d_av[6]
    amplitude_init = d_av[0]
    
    
    wl = info['WL'][0]

    root = ET.parse(ICF).getroot()
    for child in root.findall('channel'):
        if child.get('name') == ch:
            pixSize = np.float(child.find('detector_pixel').find('pixel_size').get('val'))
            F_ = np.float(child.find('wfno').get('val'))  
   
    if ch == 'FGS Red':
        Ap = 5.0
        F = F_*Ap                
    elif ch == 'FGS Prime':
        Ap = 7.0
        F =  F_*Ap
    elif ch == 'NIR Phot':
        Ap = 9.0        
        
        F = F_*Ap
    else:
        F = 1.22*F_
       
    print "aperture of radius", Ap, "x F lambda, where F = ", F_, "and lambda =", wl
    
    signal = []
    subymaxlist = []
    subxmaxlist = []  
    
    
    nWLs = data.shape[2]  # how many steps in the loop

    # Progress Bar setup:
    ProgMax = 100    # number of dots in progress bar
    if nWLs<ProgMax:   ProgMax = nWLs   # if less than 20 points in scan, shorten bar
    print "|" +    ProgMax*"-"    + "|     MyFunction() progress"
    sys.stdout.write('|'); sys.stdout.flush();  # print start of progress bar
    nProg = 0   # fraction of progress    
    

    for i in range(data.shape[2]):
        
        if ( i >= nProg*nWLs/ProgMax ):
                '''Print dot at some fraction of the loop.'''
                sys.stdout.write('*'); sys.stdout.flush();  
                nProg = nProg+1
        if ( i >= nWLs-1 ):
                '''If done, write the end of the progress bar'''
                sys.stdout.write('|     done  \n'); sys.stdout.flush(); 

        inImage = data[...,i]  
              
        if diff == 0:

            def twoD_Gaussian((x, y, sigma_x, sigma_y, theta), amplitude, xo, yo, offset):
                xo = float(xo)
                yo = float(yo)    
                a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
                g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                                        + c*((y-yo)**2)))
                return g.ravel()
            
            # Create x and y indices

            x = np.linspace(0, inImage.shape[1]-1, inImage.shape[1])
            y = np.linspace(0, inImage.shape[0]-1, inImage.shape[0])
            x, y = np.meshgrid(x, y)
            
            data_ravel = inImage.ravel()
            Y,X =  np.unravel_index(inImage.argmax(),inImage.shape)
            
            initial_guess = (amplitude_init, X,Y, offset_init)
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y, sigma_x, sigma_y, theta), data_ravel, p0=initial_guess)
            
#            print popt
            
            data_fitted = twoD_Gaussian((x, y, sigma_x, sigma_y, theta), *popt)
            data_fitted = data_fitted.reshape(inImage.shape[0], inImage.shape[1])
            
            Y,X =  np.unravel_index(data_fitted.argmax(),data_fitted.shape)
             
            fact = 10. # 30 is too slow           
            y = np.arange(0,inImage.shape[0],1)
            y2 = np.arange(0,inImage.shape[0], 1/fact)
            x = np.arange(0,inImage.shape[1],1)
            x2 = np.arange(0,inImage.shape[0], 1/fact)
            
            f = interpolate.interp2d(x, y, data_fitted, kind='cubic')
            data_fitted_new = f(x2,y2)
            
            Y,X =  np.unravel_index(data_fitted_new.argmax(),data_fitted_new.shape)
            X /= fact
            Y /= fact
            
            subymax = Y
            subxmax = X
            
#            print subymax, subxmax
            
        elif diff == 1:
            subymax = inImage.shape[0]/2
            subxmax = inImage.shape[1]/2
             
        # apply circular aperture centred at ymax, xmax      
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
        signal.append(Ap_Count)
        
    return np.array(signal),wl
                 
