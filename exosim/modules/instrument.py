from ..classes.sed import Sed
from ..lib         import exolib
import numpy     as     np
import copy
from scipy import interpolate
import matplotlib.pyplot as plt


class Channel(object):
  planet       = None
  star         = None
  zodi         = None
  emission     = None
  transmission = None
  psf          = None

  
  def __init__(self, star, planet, zodi, emission, transmission):
    self.star         = copy.deepcopy(star)
    self.planet       = copy.deepcopy(planet)
    self.zodi         = copy.deepcopy(zodi)
    self.emission     = copy.deepcopy(emission)
    self.transmission = copy.deepcopy(transmission)
   
  def rebin(self):
    self.star.rebin(self.wl_solution)
    self.star.planet(self.wl_solution)
    self.zodi.rebin(self.wl_solution)
    self.emission.rebin(self.wl_solution)
    self.transmission.rebin(self.wl_solution)

def run(opt, star, planet, zodi):
  
  instrument_emission       = Sed(star.sed.wl, np.zeros(star.sed.wl.size, dtype=np.float64))
  instrument_transmission   = Sed(star.sed.wl, np.ones(star.sed.wl.size, dtype=np.float64))
  
  for op in opt.common_optic['optical_surface']:
     
    dtmp=np.loadtxt(op.transmission.replace('$root$', opt.common_exosym_path.val), delimiter=',')
    tr = Sed(dtmp[:,0],dtmp[:,1])
    tr.rebin(opt.common_wl.val)
    dtmp=np.loadtxt(op.emissivity.replace('$root$', opt.common_exosym_path.val), delimiter=',')
    em = Sed(dtmp[:,0],dtmp[:,1])
    em.rebin(opt.common_wl.val)
        
    exolib.sed_propagation(star.sed, tr)
    exolib.sed_propagation(zodi.sed, tr)
    exolib.sed_propagation(instrument_emission, tr, emissivity=em, temperature=op.temperature)
    
    instrument_transmission.sed = instrument_transmission.sed*tr.sed
    
    channel = {}
    for key in opt.channel.keys():
        
      channel[key] = Channel(star.sed, planet.cr, zodi.sed, instrument_emission, instrument_transmission)
      
      for op in opt.channel[key]['optical_surface']:

          dtmp=np.loadtxt(op.transmission.replace('$root$', opt.common_exosym_path.val), delimiter=',')
          tr = Sed(dtmp[:,0],dtmp[:,1])
          tr.rebin(opt.common_wl.val)
          em = Sed(dtmp[:,0],dtmp[:,2])	  
          em.rebin(opt.common_wl.val)
          exolib.sed_propagation(channel[key].star, tr)
          exolib.sed_propagation(channel[key].zodi, tr)
          exolib.sed_propagation(channel[key].emission, tr, emissivity=em, temperature=op.temperature)
          channel[key].transmission.sed = channel[key].transmission.sed*tr.sed
     
      ### create focal plane
           

      channel[key].osf          = opt.channel[key]['osf'].val    # Oversample psf by this factor to #   ensure Nyquist
      channel[key].kernel_osf   = opt.channel[key]['kernel_osf'].val   # Oversample kernel to ensure Nyquist sampling
      channel[key].psf_osf      = opt.channel[key]['psf_osf'].val
           
      #1# Obtain wavelength dispertion relation 
      ld  = np.fromstring(opt.channel[key]['ld'].val, 
	  sep=' ', dtype=np.float64)
   
      #2# allocate focal plane with pixel oversampling such that Nyquist sampling is done correctly 
      fpn = np.fromstring(opt.channel[key]['array_geometry'].val, 
	  sep=' ', dtype=np.float64)
      fp = np.zeros(fpn*channel[key].osf)
      
      #3# This is the current sampling interval in the focal plane.  
      fp_delta = opt.channel[key]['pixel_size'].val/channel[key].osf
      
      x_dist = np.arange(opt.channel[key]['pixel_size'].val/2,  
                         opt.channel[key]['pixel_size'].val/2 + fpn[1] *opt.channel[key]['pixel_size'].val, 
                         opt.channel[key]['pixel_size'].val) 
      x_wav =  ld[0] + ld[1]*(x_dist-ld[2])
      
      x_dist_osf =  np.arange(fp_delta /2,  
                         fp_delta/2 + fp.shape[1] * fp_delta, 
                         fp_delta) 
                         
      x_wav_osf = ld[0] + ld[1]*(x_dist_osf-ld[2]) # wvalength on each x pixel
      
                                  
      #5# Generate PSFs, one for each detector pixel along spectral axis
      psf = exolib.Psf(x_wav, opt.channel[key]['wfno'].val, fp_delta, shape='airy') 
      
      #6# Repeate same PSF on each pixel
      psf = np.repeat(psf, channel[key].osf, axis=2)
      
      #7# Save results in Channel class
      channel[key].fp_delta    = fp_delta
      channel[key].psf         = psf
      channel[key].fp          = fp
      channel[key].wl_solution = x_wav_osf
      
      channel[key].planet.sed  *= channel[key].star.sed
      channel[key].planet.units = channel[key].star.units
      
      channel[key].star.rebin(channel[key].wl_solution)
      channel[key].planet.rebin(channel[key].wl_solution)
#      print channel[key].wl_solution.size, channel[key].planet.sed.size
      
      
      #parameters needed to obtain photon count per wavelength from sed
       
      delta_wav = x_wav_osf[1]- x_wav_osf[0]
      ap_area = np.pi*(opt.common_optic['diameter'].val/2)**2
      photon_E = 6.62606957e-34* 299792458/(x_wav_osf*1e-6)
                        
      #8# Populate the point source focal plane
 
      psf_delta = (1/channel[key].psf_osf) *channel[key].osf
 
      j0 = np.round(np.arange(-psf.shape[1]/2,  fp.shape[1] - psf.shape[1]/2, psf_delta)).astype(np.int) 
      j1 = j0 + psf.shape[1]
      idx = np.where((j0>=0) & (j1 < fp.shape[1]))[0]
      i0 = fp.shape[0]/2 - psf.shape[0]/2 
      i1 = i0 + psf.shape[0]
      for k in idx: 
         fp[i0:i1, j0[k]:j1[k]] += psf[...,k]*channel[key].star.sed[k]*ap_area*delta_wav/photon_E[k]


      # apply quantum effiicency and variations
        
      QE = 0.6  # to be replaced with file
      QEsd = 0.05
            
      QE_array = np.random.normal(QE,QEsd*QE,fpn)
      QE_array = np.repeat(QE_array,channel[key].osf,axis=0)
      QE_array = np.repeat(QE_array,channel[key].osf,axis=1)  
       
      fp *= QE_array         

        
                
      #10# Allocate pixel response function
       
      kernel, kernel_delta = exolib.PixelResponseFunction(
                    channel[key].psf.shape[0:2],
                    channel[key].fp_delta,
                    channel[key].kernel_osf, 
                    opt.channel[key]['pixel_size'].val)
                    
      
      channel[key].fp = exolib.fast_convolution(
                                   fp, 
                                   channel[key].fp_delta,
                                   kernel, kernel_delta)
                                         
      
    return channel

  pass