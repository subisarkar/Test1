from ..classes.sed import Sed
from ..classes.channel import Channel
from ..lib         import exolib

from exosim_nb.lib.exolib import exosim_msg
import time
import numpy           as np
import quantities      as pq
import scipy.constants as spc
import scipy.interpolate
import matplotlib.pyplot as plt

def run(opt, star, planet, zodi):
  
  exosim_msg('Run instrument model ... ')
  st = time.time()
  instrument_emission  = Sed(star.sed.wl, 
                             np.zeros(star.sed.wl.size, dtype=np.float64)* \
                             pq.W/pq.m**2/pq.um/pq.sr)
  instrument_transmission = Sed(star.sed.wl, np.ones(star.sed.wl.size, dtype=np.float64))

  # use for prism xml file
  

  TR =np.array([1.]*len(opt.common.common_wl))
  
#  
  for op in opt.common_optics.optical_surface:
    dtmp=np.loadtxt(op.transmission.replace('__path__', opt.__path__), delimiter=',')
    tr = Sed(dtmp[:,0]*pq.um,dtmp[:,1]*pq.dimensionless)
    tr.rebin(opt.common.common_wl)
    TR *=tr.sed
#    plt.figure('transmission')
#    plt.plot(opt.common.common_wl,TR)

    em = Sed(dtmp[:,0]*pq.um,dtmp[:,2]*pq.dimensionless)
    em.rebin(opt.common.common_wl)
    
    exolib.sed_propagation(star.sed, tr)
    exolib.sed_propagation(zodi.sed, tr)
    exolib.sed_propagation(instrument_emission, tr, emissivity=em, temperature=op())
    instrument_transmission.sed = instrument_transmission.sed*tr.sed

  # use for LRS xml file

#  for op in opt.common_optics.optical_surface:
#    dtmp=np.loadtxt(op.transmission.replace('__path__', opt.__path__), delimiter=',')
#    tr = Sed(dtmp[:,0]*pq.um,dtmp[:,1]*pq.dimensionless)
#    tr.rebin(opt.common.common_wl)
#
#    dtmp=np.loadtxt(op.emissivity.replace('__path__', opt.__path__), delimiter=',')
#    em = Sed(dtmp[:,0]*pq.um,dtmp[:,1]*pq.dimensionless)
#    em.rebin(opt.common.common_wl)
#    
#    exolib.sed_propagation(star.sed, tr)
#    exolib.sed_propagation(zodi.sed, tr)
#    exolib.sed_propagation(instrument_emission, tr, emissivity=em, temperature=op())
#    instrument_transmission.sed = instrument_transmission.sed*tr.sed
###

    
  channel = {}
  FPCOUNTLIST =[]
  for ch in opt.channel:
    #if ch.name != 'NIR Spec': continue
    channel[ch.name] = Channel(star.sed, planet.cr, 
			       zodi.sed, 
			       instrument_emission,   
			       instrument_transmission,
			       options=ch)
    print "CHANNEL", ch.name
    
    ch_optical_surface = ch.optical_surface if isinstance(ch.optical_surface, list) else \
      [ch.optical_surface]
    for op in ch.optical_surface:
      dtmp=np.loadtxt(op.transmission.replace(
          '__path__', opt.__path__), delimiter=',')
      tr = Sed(dtmp[:,0]*pq.um, \
              dtmp[:,1]*pq.dimensionless)
      tr.rebin(opt.common.common_wl)
      TR *=tr.sed
#      plt.figure('transmission')
#      plt.plot(opt.common.common_wl,TR)    
#      plt.figure('transmission2')
#      plt.plot(tr.wl, tr.sed) 
      
      em = Sed(dtmp[:,0]*pq.um, \
               dtmp[:,2]*pq.dimensionless)
      em.rebin(opt.common.common_wl)
      exolib.sed_propagation(channel[ch.name].star, tr)
      exolib.sed_propagation(channel[ch.name].zodi, tr)
      exolib.sed_propagation(channel[ch.name].emission, \
              tr, emissivity=em,temperature=op())
      channel[ch.name].transmission.sed *= tr.sed
   
    # BUG workaround. There is a bug in the binning function. If transmission is zero,
    # it is rebiined to a finite, very small value. This needs to be fixed!
    # For now, I set to zero all transmission smaller than an arbitrary value
    #idx = np.where(channel[ch.name].transmission.sed < 1.0e-5)
    #channel[ch.name].star.sed[idx] = 0.0*channel[ch.name].star.sed.units
    #channel[ch.name].zodi.sed[idx] = 0.0*channel[ch.name].zodi.sed.units
    #channel[ch.name].emission.sed[idx] = 0.0*channel[ch.name].emission.sed.units
    #channel[ch.name].transmission.sed[idx] = 0.0*channel[ch.name].transmission.sed.units
#    channel[ch.name].star.sed *=0.80
    
    # Convert spectral signals
    dtmp=np.loadtxt(ch.qe().replace(
	    '__path__', opt.__path__), delimiter=',')
    qe = Sed(dtmp[:,0]*pq.um, \
		 dtmp[:,1]*pq.dimensionless)
   

    Responsivity = qe.sed * qe.wl.rescale(pq.m)/(spc.c * spc.h * pq.m * pq.J)*pq.UnitQuantity('electron', symbol='e-')
    
    Re = scipy.interpolate.interp1d(qe.wl, Responsivity)
    
    Aeff = 0.25*np.pi*opt.common_optics.TelescopeEffectiveDiameter()**2
    Omega_pix = 2.0*np.pi*(1.0-np.cos(np.arctan(0.5/ch.wfno())))*pq.sr
    Apix = ch.detector_pixel.pixel_size()**2
    channel[ch.name].star.sed     *= Aeff             * \
      Re(channel[ch.name].star.wl)*pq.UnitQuantity('electron', 1*pq.counts, symbol='e-')/pq.J
    channel[ch.name].zodi.sed     *= Apix * Omega_pix * \
      Re(channel[ch.name].zodi.wl)*pq.UnitQuantity('electron', 1*pq.counts, symbol='e-')/pq.J
    channel[ch.name].emission.sed *= Apix * Omega_pix * \
      Re(channel[ch.name].emission.wl)*pq.UnitQuantity('electron', 1*pq.counts, symbol='e-')/pq.J
    
    ### create focal plane 

    #1# allocate focal plane with pixel oversampling such that Nyquist sampling is done correctly 
    fpn = ch.array_geometry()

#    print ch.dispersion() 
    print "Focal plane array geometry:  ",fpn
    
    fp  = np.zeros( (fpn*ch.osf()).astype(np.int) )

    #2# This is the current sampling interval in the focal plane.  
    fp_delta = ch.detector_pixel.pixel_size() / ch.osf()
    
    #3# Load dispersion law 
    if ch.type == 'spectrometer':
      if hasattr(ch, "dispersion"):
          ch.dispersion.val = (ch.detector_pixel.pixel_size()*fpn[1]  )/2.
          dtmp=np.loadtxt(ch.dispersion.path.replace(
	  '__path__', opt.__path__), delimiter=',')
          ld = scipy.interpolate.interp1d(dtmp[...,2]*pq.um + ch.dispersion().rescale(pq.um), 
					dtmp[...,0],
					bounds_error=False, 
					fill_value=0.0)
      elif hasattr(ch, "ld"):
          # wl = ld[0] + ld[1](x - ld[2]) = ld[1]*x + ld[0]-ldp[1]*ld[2]
          ld = np.poly1d( (ch.ld()[1], ch.ld()[0]-ch.ld()[1]*ch.ld()[2]) )
      else:
          exolib.exosim_error("Dispersion law not defined.")
      
      #4a# Estimate pixel and wavelength coordinates
      x_pix_osr = np.arange(fp.shape[1]) * fp_delta  
#      print len(x_pix_osr)
      x_wav_osr = ld(x_pix_osr.rescale(pq.um))*pq.um # walength on each x pixel
      channel[ch.name].wl_solution = x_wav_osr
      
#      print x_wav_osr[1::3]
   
    elif ch.type == 'photometer':
      #4b# Estimate pixel and wavelength coordinates
#      print channel[ch.name].transmission.sed.max()
      idx = np.where(channel[ch.name].transmission.sed > channel[ch.name].transmission.sed.max()/np.e)
      
#      print channel[ch.name].transmission.wl[idx].min().item()
      x_wav_osr = np.linspace(channel[ch.name].transmission.wl[idx].min().item(),
			      channel[ch.name].transmission.wl[idx].max().item(),
			      8 * ch.osf()) * channel[ch.name].transmission.wl.units
      x_wav_center = (channel[ch.name].transmission.wl[idx]*channel[ch.name].transmission.sed[idx]).sum() / \
	channel[ch.name].transmission.sed[idx].sum()
 
      print "Central wavelength:  ", x_wav_center
     
      channel[ch.name].wl_solution = np.repeat(x_wav_center, fp.shape[1])
     
    else:
      exolib.exosim_error("Channel should be either photometer or spectrometer.")
      
    d_x_wav_osr = np.zeros_like (x_wav_osr)
    idx = np.where(x_wav_osr > 0.0)
    d_x_wav_osr[idx] = np.gradient(x_wav_osr[idx])
    if np.any(d_x_wav_osr < 0): d_x_wav_osr *= -1.0
        
  
  
    
#     for defined psf in photometric channels: 4 is latest F numbers. 3 are previous F numbers
    if ch.name == 'FGS Red':
        zfile = '/Users/user1/Desktop/zfile_fgs_red.fits'
        zfile = '/Users/user1/Desktop/zfile_fgs_red_Euro4.fits'     
        print "PSF file chosen:  ", zfile
    if ch.name == 'FGS Prime':  
        zfile = '/Users/user1/Desktop/zfile_fgs_prime.fits'
        zfile = '/Users/user1/Desktop/zfile_fgs_prime_Euro4.fits'
        print "PSF file chosen:  ", zfile
    if ch.name == 'NIR Phot':  
        zfile = '/Users/user1/Desktop/zfile_fgs_nirphot.fits'
        zfile = '/Users/user1/Desktop/zfile_fgs_nirphot_Euro4.fits'
        print "PSF file chosen:  ", zfile

    if ch.name == 'FGS Red' or ch.name == 'FGS Prime' or ch.name == 'NIR Phot':
        psf = exolib.Psf_photometer(zfile, fp_delta, x_wav_osr)
        print "WFE error PSF chosen..."
    else:
        psf = exolib.Psf(x_wav_osr, ch.wfno(), fp_delta, shape='airy') 

    #6# Save results in Channel class
    channel[ch.name].fp_delta    = fp_delta
    channel[ch.name].psf         = psf
    channel[ch.name].fp          = fp
    channel[ch.name].osf         = np.int(ch.osf())
    channel[ch.name].offs        = np.int(ch.pix_offs())
    
    channel[ch.name].planet.sed  *= channel[ch.name].star.sed

    channel[ch.name].star.rebin(x_wav_osr)
    channel[ch.name].planet.rebin(x_wav_osr)
    channel[ch.name].zodi.rebin(x_wav_osr)
    channel[ch.name].emission.rebin(x_wav_osr)
    channel[ch.name].transmission.rebin(x_wav_osr)
    
#    wav_mid = 1.25+(1.9-1.25)/2
##        
#    wav_mid = 1.95+(3.9-1.95)/2
#    wav_mid = 3.9+(7.8-3.9)/2
#
##
#    wav_mid =x_wav_center
#    print "wav mid", wav_mid
#    idx = np.argwhere(x_wav_osr.magnitude>wav_mid)[-1]
#    print idx
#    print x_wav_osr[idx]
#    print "mid channel transmission", channel[ch.name].transmission.sed[idx]
#    xx
    channel[ch.name].star.sed     *= d_x_wav_osr
    channel[ch.name].planet.sed   *= d_x_wav_osr
    channel[ch.name].zodi.sed     *= d_x_wav_osr 
    channel[ch.name].emission.sed *= d_x_wav_osr
    
    #7# Populate focal plane with monochromatic PSFs
    if ch.type == 'spectrometer':
      j0 = np.round(np.arange(fp.shape[1]) - psf.shape[1]/2).astype(np.int) 
    
    elif ch.type == 'photometer':
      j0 = np.repeat(fp.shape[1]//2 - psf.shape[1]/2, x_wav_osr.size)
    else:
      exolib.exosim_error("Channel should be either photometer or spectrometer.")
      
 
    j1 = j0 + psf.shape[1]
    idx = np.where((j0>=0) & (j1 < fp.shape[1]))[0]
    i0 = fp.shape[0]/2 - psf.shape[0]/2 + channel[ch.name].offs 
    i1 = i0 + psf.shape[0]
    
    import copy
    FPCOPY = copy.deepcopy(channel[ch.name].fp)
    for k in idx:FPCOPY[i0:i1, j0[k]:j1[k]] += psf[...,k] * \
    channel[ch.name].star.sed[k] 
    
    if opt.source_switch ==1 :
        for k in idx: channel[ch.name].fp[i0:i1, j0[k]:j1[k]] += psf[...,k] * \
        channel[ch.name].star.sed[k] 
        
    
    #9# Now deal with the planet
    planet_response = np.zeros(fp.shape[1])
    i0p = np.unravel_index(np.argmax(channel[ch.name].psf.sum(axis=2)), channel[ch.name].psf[...,0].shape)[0]
    for k in idx: planet_response[j0[k]:j1[k]] += psf[i0p,:,k] * channel[ch.name].planet.sed[k] 
    
    #9# Allocate pixel response function
    kernel, kernel_delta = exolib.PixelResponseFunction(
        channel[ch.name].psf.shape[0:2],
        7*ch.osf(),   # NEED TO CHANGE FACTOR OF 7 
        ch.detector_pixel.pixel_size(),
        lx = ch.detector_pixel.pixel_diffusion_length())


    channel[ch.name].fp = exolib.fast_convolution(
        channel[ch.name].fp, 
        channel[ch.name].fp_delta,
        kernel, kernel_delta)
        
    FPCOPY = exolib.fast_convolution(
        FPCOPY, 
        channel[ch.name].fp_delta,
        kernel, kernel_delta) 
        
    ## TODO CHANGE THIS: need to convolve planet with pixel response function
    channel[ch.name].planet = Sed(channel[ch.name].wl_solution, planet_response/(1e-30+fp[(i0+i1)//2, ...]))
    ## Fix units
    channel[ch.name].fp = channel[ch.name].fp*channel[ch.name].star.sed.units   
    channel[ch.name].planet.sed = channel[ch.name].planet.sed*pq.dimensionless
  
    
    ## Deal with diffuse radiation
    if ch.type == 'spectrometer':
      channel[ch.name].zodi.sed     = scipy.convolve(channel[ch.name].zodi.sed, 
		      np.ones(np.int(ch.slit_width()*channel[ch.name].opt.osf())), 
		      'same') * channel[ch.name].zodi.sed.units
      channel[ch.name].emission.sed = scipy.convolve(channel[ch.name].emission.sed, 
		      np.ones(np.int(ch.slit_width()*channel[ch.name].opt.osf())), 
		      'same') * channel[ch.name].emission.sed.units
    elif ch.type == 'photometer':
      channel[ch.name].zodi.sed = np.repeat(channel[ch.name].zodi.sed.sum(),
					    channel[ch.name].wl_solution.size)
      channel[ch.name].zodi.wl = channel[ch.name].wl_solution
      channel[ch.name].emission.sed = np.repeat(channel[ch.name].emission.sed.sum(),
						channel[ch.name].wl_solution.size)
      channel[ch.name].emission.wl = channel[ch.name].wl_solution
      
      

      
    FPCOPY += channel[ch.name].zodi.sed.magnitude
    FPCOPY += channel[ch.name].emission.sed.magnitude
    FPCOUNT = FPCOPY[1::3,1::3] + ch.detector_pixel.Idc.val.magnitude
    FW1 = 50e3
    FW2 = 100e3
    FW3 = 70e3
    A,B = np.unravel_index(FPCOUNT.argmax(), FPCOUNT.shape)
    print "maximum index", A,B
    print "calculcated int time (50 K)", FW1 / FPCOUNT.max()
    print "calculcated int time (100 K)", FW2 / FPCOUNT.max()
    print "calculcated int time (70 K)", FW3 / FPCOUNT.max()
    
    if ch.name == 'FGS Red' or ch.name == 'FGS Prime' or ch.name == 'NIR Phot' or ch.name == 'NIR Spec': 
        FPCOUNTLIST.append(FW2 / FPCOUNT.max())
        print "FULL WELL chosen = 100 k"
    else:
        FPCOUNTLIST.append(FW3 / FPCOUNT.max())
        print "FULL WELL chosen = 70 k"
        
    if len(FPCOUNTLIST) == len(opt.channel):
        t_min = np.min(FPCOUNTLIST)
        t1 = (opt.timeline.nGND.val / opt.timeline.frame_rate.val).magnitude
        t2 = (opt.timeline.nRST.val / opt.timeline.frame_rate.val).magnitude
        
        print "CALCULATED EXPOSURE TIME", t_min
        # min time to read frame = 0.1 sec - estimate of NDR time
        if t_min < 0.1:
            print "CALCULATED EXPOSURE TIME < 0.1 sec, therefore set to 0.1 sec"
            t_min = 0.1
        if t_min > 300.0:
            print "CALCULATED EXPOSURE TIME > 300 sec, therefore set to 300 sec"
            t_min = 300.0  

        if opt.use_sat == False:
            t_min = opt.int_time             
            
        print "SELECTED EXP TIME", t_min     
        opt.timeline.exposure_time.val = (t_min + t1 + t2) *pq.s
        print "TOTAL CYCLE TIME", t_min + t1 + t2
        print "CDS time estimate", t_min - (opt.timeline.nNDR0.val / opt.timeline.frame_rate.val).magnitude
   

    if opt.dc_switch !=1: 
        ch.detector_pixel.Idc.val = 0*1/pq.s
        
    print "EMISSION SWITCH:  ",bool(opt.emm_switch)
    print "ZODI SWITCH:  ",bool(opt.zodi_switch)
    print "DARK CURRENT SWITCH:  ",bool(opt.dc_switch)    
    print "DARK CURRENT", ch.detector_pixel.Idc.val
   
        
    print "FP max", channel[ch.name].fp[1::3,1::3].max()
    print "FP max x int time", channel[ch.name].fp[1::3,1::3].max()*t_min
    print FPCOUNT.max()*t_min
    

 
#    plt.figure('focal plane%s'%(ch.name))
#    x = np.arange(0.5*ch.detector_pixel.pixel_size(),0.5*ch.detector_pixel.pixel_size()+fpn[1]*ch.detector_pixel.pixel_size(),ch.detector_pixel.pixel_size()) 
#    plt.plot(channel[ch.name].fp[1::3,1::3].sum(axis=1))
#    plt.xlabel('focal plane x axis distance (m) %s'%(ch.name))
#    plt.ylabel('count per second per pixel column')
#    np.save('/Users/user1/Desktop/p_mod_fcounts_%s.npy'%(ch.name),np.vstack((x,channel[ch.name].fp[1::3,1::3].sum(axis=0))))
#    np.save('/Users/user1/Desktop/wl_%s'%(ch.name),channel[ch.name].wl_solution[1::3])
#    plt.figure('fp2')
#    plt.imshow(channel[ch.name].fp[1::3,1::3].magnitude)
#    plt.figure('wavelength solution')
#    plt.plot(x-0.0012,x_wav_osr[1::3],'bx-')
#    plt.plot(x00/1e6,w00, 'g-')
#    plt.figure('focal plane%s2'%(ch.name))
#    x = np.arange(0.5*ch.detector_pixel.pixel_size(),0.5*ch.detector_pixel.pixel_size()+fpn[1]*ch.detector_pixel.pixel_size(),ch.detector_pixel.pixel_size()) 
#    plt.plot(t_min*channel[ch.name].fp[1::3,1::3].sum(axis=0))
#    
#    xxx
    
  exosim_msg(' - execution time: {:.0f} msec.\n'.format((time.time()-st)*1000.0))
  return channel
 
  pass
