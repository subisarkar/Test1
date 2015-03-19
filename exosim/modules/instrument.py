from ..classes.sed import Sed
from ..lib         import exolib
import numpy     as     np
import copy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


import numpy     as     np
import copy

class Channel(object):
  planet       = None
  star         = None
  zodi         = None
  emission     = None
  transmission = None
  psf          = None
  psf_osf      = 4    # Oversample psf by this factor to ensure Nyquist- use even
  osf          = 21   # Oversample each pixel by this factor per axis
  ad_ovs     = 5    # Oversample each pixel by this factor after osf and convolution with prf kernal
  fp           = 0    # focal plane
  fp_delta     = 0    # Focal plane sampling interval (equal for spatial -y- and spectral -x- directions)
  wl_solution  = 0    # this is the wavelength solution for the focal plane at its current sampling. 
  def __init__(self, star, planet, zodi, emission, transmission):
    self.star         = copy.deepcopy(star)
    self.planet       = copy.deepcopy(planet)
    self.zodi         = copy.deepcopy(zodi)
    self.emission     = copy.deepcopy(emission)
    self.transmission = copy.deepcopy(transmission)
   

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
    for key in ['SWIR']:
#    for key in opt.channel.keys():
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
      
      
      channel[key].psf_osf      = 4    # Oversample psf by this factor to ensure Nyquist- use even
      channel[key].osf          = 21   # Oversample each pixel by this factor per axis prior to pixel convolution
      channel[key].ad_osf     = 5    # Oversample each pixel again
      channel[key].quantum_efficiency = 0.6
      channel[key].quantum_efficiency_sd = 0.005

      osf = channel[key].osf 
      ad_osf = channel[key].ad_osf   

      ### create focal plane
      #1# Obtain wavelength dispersion relation
      ld  = np.fromstring(opt.channel[key]['ld'].val, 
	  sep=' ', dtype=np.float64)
      #2# allocate focal plane with pixel oversampling such that Nyquist sampling is done correctly 
      fpn = np.fromstring(opt.channel[key]['array_geometry'].val, 
	  sep=' ', dtype=np.float64)
      fp = np.zeros(fpn*channel[key].osf)  
      #3# This is the sampling interval in the oversampled focal plane.
      fp_delta = opt.channel[key]['pixel_size'].val/channel[key].osf
      
      pix_size = opt.channel[key]['pixel_size'].val
      
      x_pix     = np.arange(fpn[1]) * pix_size
      x_pix_osr = np.arange(fp.shape[1])  * fp_delta
      x_wav     = ld[0] + ld[1]*(x_pix-ld[2]) # wvalength on each x pixel
      x_wav_osr = ld[0] + ld[1]*(x_pix_osr-ld[2]) # wvalength on each x pixel
                             
      psf = exolib.psf(x_wav, opt.channel[key]['wfno'].val,  fp_delta, shape='airy')      
      psf = np.repeat(psf, channel[key].psf_osf, axis=2)
                  
      psf_delta = (1./channel[key].psf_osf)*channel[key].osf  
          
      channel[key].fp_delta    = fp_delta
      channel[key].psf         = psf
      channel[key].fp          = fp
      channel[key].wl_solution = x_wav_osr

        
      QE = np.random.normal(channel[key].quantum_efficiency,
                            channel[key].quantum_efficiency_sd*channel[key].quantum_efficiency
                            ,fpn)
      QE = np.repeat(QE,osf,axis=0)
      QE = np.repeat(QE,osf,axis=1)
     
      kernal = exolib.kernal(pix_size, channel[key].osf)
#      sed = np.ones((psf.shape[2]))
      sed = channel[key].star.sed[0:psf.shape[2]]

      
      fpa, conv_fpa= exolib.fpa(fp,psf,psf_delta,sed,kernal,ad_osf,pix_size,QE)
            
      channel[key].fpa = fpa
      channel[key].conv_fpa = conv_fpa
         
    return channel

  pass