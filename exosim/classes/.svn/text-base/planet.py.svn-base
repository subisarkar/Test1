import os
import numpy as np
from classes import sed
from lib import occultquad 
from lib.exolib import exosim_error

class Planet(object):
  """
    Instantiate a Planet class 
    
    Attributes
    ----------
    wl				array
				wavelength [micron]
    sed 			array
				contains planet-star contrast ratio
  """
  t14 = 0.0
  phi = None
  z   = None
  lc  = None

  def __init__(self, contrast_path):
    """
    Initialize a new Planet class 
    
    Parameters
    __________
    contrast_path		string
				path name of file containg the planet-star contrast spectrum
    """

    pl_wl, pl_sed = self.read_planet_spectrum(contrast_path)
    idx = np.argsort(pl_wl)      
    self.cr = sed.Sed(pl_wl[idx], pl_sed[idx], units='')  
    
  def get_t14(self, inc, a, period, planet_radius, star_radius):
    """ t14
    Calculates the transit time 
    
    Parameters
    __________
    inc:			scalar
				Planet oprbital inclination [rad]
    a:				scalar
				Semimajor axis [meters]
    period:			scalar
				Orbital period [seconds]
    planet_radius	:	scalar
				Planet radius [meters]
    star_radius	:		scalar
				Star radius [meters]
    
    Returns
    __________
    transit duration : float
	Returns the transit duration [seconds]
    
    Notes
    _____
    Seager, S., & Mallen-Ornelas, G. 2003, ApJ, 585, 1038
      
    
    """
    impact_parameter = np.cos(inc)*a/star_radius
    dtmp = 1+planet_radius/star_radius
    if impact_parameter < dtmp:
      self.t14 = period/np.pi * \
	    star_radius/a * \
	    np.sqrt(dtmp**2 - impact_parameter**2)
    else:
      print "WARNING: planet not transiting"
      self.t14 = 0.0
    return self.t14
    
  def read_planet_spectrum(self, contrast_path):
    try:
      _d = np.loadtxt(contrast_path)
    except IOError:
      exosim_error( 'problem reading '+ contrast_path) 

    return _d[:,0], _d[:,1]
  
  def get_orbital_phase(self, t14, period, N=1000):
    f = t14/period
    self.phi = np.linspace(-f, f, N)
    return self.phi
    
  def eccentric(self, phi, inc, ecc, omega, a, period, star_radius):
    """ eccentric
    Implements an eccentric orbit and calculates the projectedseparation (z)
    
    Parameters
    __________
    phi : 	array
      orbital phase 
    inc :	float
      orbital inclination [radiants]
    ecc : 	float
      orbital eccentricity [radiants]
    omega :	float
      argument of periastron [radiants]
    a	:	float
      semimajor axis [meters]
    star_radius : 	float
      star radius [meters]
    
    Returns
    _______
    z	: 	array
      orbital separations
    
    Notes
    _____
    This implementation written by Ingo Waldman 2012
    """
    theta = 2.0 * np.pi * phi 
    aR    = a / star_radius
    
    if ecc < 1e-5:
      # calculating z for a circular orbit
      self.z = aR * np.sqrt(1-((np.cos(theta))**2.*(np.sin(inc))**2.))
      return self.z
    
    # calculating z for eccentric orbits
    n = len(theta)
    E = np.zeros(n)
    ecc2 = np.sqrt((1.+ecc)/(1.-ecc))
    fref = (np.pi / 2.) - omega #setting reference point for the true anomaly
    
    Eref = 2. * np.arctan(1./ecc2 * np.tan(fref/2.)) #reference point for eccentric anomaly
    if Eref < (-np.pi/2.):
      Eref = Eref + 2. * np.pi
      
    Mref = Eref - (ecc * np.sin(Eref)) #reference point for mean anomaly
   
    for i in range(n): 
      # calculating eccentric anomaly using Newton-Rhapson method
      Mtmp = theta[i] + Mref #calculating mean anomaly for phase point
      Etmp = Mtmp
      for j in range(10):
        Etmp = Etmp + ((Mtmp + ecc*np.sin(Etmp) - Etmp) / (1.-ecc*np.cos(Etmp)))
      E[i] = Etmp
    
      # calculating true anomaly
      f = 2.*np.arctan(ecc2*np.tan(E/2.))
      # calculating distance from true anomaly as fraction
      r_frac = (1.-ecc**2.)/(1. + ecc*np.cos(f))
      # computing z
      self.z = 1. - (np.sin(inc)**2.*np.sin(f+omega)**2.)
      self.z = aR*r_frac*np.sqrt(self.z)

    return self.z
    
  def get_light_curve(self, z, u1, u2, p0, primary_transit = True):
    self.lc = occultquad.occultquad(z, u1, u2, p0)[0 if primary_transit else 1] 
    return self.lc
