import numpy   as np
import os
from   ..classes import sed
from   ..lib     import exolib

class zodiacal_light(object):
  def __init__(self, wl, level=1.0):
    
    spectrum = level*(3.5e-14*exolib.planck(wl, 5500) + 
			    exolib.planck(wl, 270) * 3.58e-8)
    
    self.sed 		= sed.Sed(wl, spectrum )
    self.transmission 	= sed.Sed(wl, np.ones(wl.size))
    self.units         = 'W m**-2 sr**-1 micron**-1'
    