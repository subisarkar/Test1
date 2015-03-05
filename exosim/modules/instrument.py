from ..classes.sed import Sed
from ..lib         import exolib
import numpy     as     np
import copy

class Channel(object):
  planet       = None
  star         = None
  zodi         = None
  emission     = None
  transmission = None
  psf          = None
  osf          = 3    # Oversample psf by this factor to ensure Nyquist
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
      #1# Obtain wavelength dispertion relation 
      ld  = np.fromstring(opt.channel[key]['ld'].val, 
	  sep=' ', dtype=np.float64)
      #2# allocate focal plane with pixel oversampling such that Nyquist sampling is done correctly 
      fpn = np.fromstring(opt.channel['SWIR']['array_geometry'].val, 
	  sep=' ', dtype=np.float64)
      fp = np.zeros(fpn*channel[key].osf)
      #2# This is the current sampling interval in the focal plane.  
      fp_delta = opt.channel[key]['pixel_size'].val/channel[key].osf
      
      x_pix     = np.arange(fpn[1])  * opt.channel[key]['pixel_size'].val
      x_pix_osr = np.arange(fp.shape[1])  * fp_delta
      x_wav     = ld[0] + ld[1]*(x_pix-ld[2]) # wvalength on each x pixel
      x_wav_osr = ld[0] + ld[1]*(x_pix_osr-ld[2]) # wvalength on each x pixel
      
      psf = exolib.Psf(x_wav, opt.channel[key]['wfno'].val,  fp_delta, shape='airy') 
      psf = np.repeat(psf, channel[key].osf, axis=2)
      
      channel[key].fp_delta    = fp_delta
      channel[key].psf         = psf
      channel[key].fp          = fp
      channel[key].wl_solution = x_wav_osr
      
      j0 = np.round(np.arange(fp.shape[1]) - psf.shape[1]/2).astype(np.int) 
      j1 = j0 + psf.shape[1]
      idx = np.where((j0>=0) & (j1 < fp.shape[1]))[0]
      i0 = fp.shape[0]/2 - psf.shape[0]/2 
      i1 = i0 + psf.shape[0]  
      for k in idx: fp[i0:i1, j0[k]:j1[k]] += psf[...,k]
      
    
    
    
    #return star, planet, zodi, instrument_emission, instrument_transmission
    return channel

  pass