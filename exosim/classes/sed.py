import numpy as np
import sys 
from lib import exolib


class Sed(object):
  def __init__(self, wl=None, sed=None, units=None):
    self.sed        = sed
    self.wl         = wl
    self.units      = units
  
  
  def rebin_(self, R, wl=None, sed=None, wl_min=None, wl_max=None):
    """
      Rebins the SED to the given spectral resolution
      
      Parameters
      __________
	R: scalar
	    spectral resolution
    """
    if wl != None: self.wl = np.copy(wl)
    if sed != None: self.sed = np.copy(sed)
    if wl_min == None: wl_min = self.wl.min()
    if wl_max == None: wl_max = self.wl.max()
    
    self.wl, self.sed = exolib.logbin(self.wl, self.sed, R, xmin=wl_min, xmax=wl_max)
    
  def rebin(self, wl):
    self.wl, self.sed = exolib.rebin(wl, self.wl, self.sed)
    
    