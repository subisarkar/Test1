import numpy   as np
import os, time, tables, pyfits
from   ..classes import sed
from   ..lib import exolib

class Star(object):
  """
    Instantiate a Stellar class using Phenix Stellar Models
    
    Attributes
    ----------
    lumiosity : 		float
				Stellar bolometric luminosity computed from Phenix stellar modle. Units [W]
    wl				array
				wavelength [micron]
    sed 			array
				Spectral energy density [[W m**-2 micron**-1]]
    ph_wl			array
				phoenix wavelength [micron]. Phenix native resolution
    ph_sed 			array
				phenix spectral energy density [W m**-2 micron**-1]. Phenix native resolution
  
  """
  def __init__(self, star_sed_path, star_distance, star_temperature, star_logg, star_f_h, star_radius):
    """
    Parameters
    __________
      star_sed_path : 		string
				path to the file containing the stellar SED [erg/s/cm^2/cm]
      star_temperature : 	scalar
				Stellar temperature [K]
      star_logg : 		scalar
				Stellar log_10(g), where g is the surface gravity
      star_F_H : 		scalar
				Stellar metallicity [F/H]	
      star_radius:      Stellar radius [m]
      
    """
    self.distance, self.temperature, self.logg, self.f_h, self.radius = star_distance, star_temperature, star_logg, star_f_h, star_radius

    ph_wl, ph_sed, ph_L = self.read_phenix_spectrum(star_sed_path, star_distance, star_temperature, star_logg, star_f_h, star_radius)
    self.ph_luminosity = ph_L
    self.ph_sed   = sed.Sed(ph_wl, ph_sed, units='W m**-2 micron**-1')
    self.sed = sed.Sed(ph_wl, ph_sed, units='W m**-2 micron**-1')
    self.luminosity = ph_L
          
  def read_phenix_spectrum(self, path, star_distance, star_temperature, star_logg, star_f_h, star_radius):
    """Read a PHENIX Stellar Spectrum. 
    
    Parameters
    __________
      path : 			string
				path to the file containing the stellar SED [erg/s/cm^2/cm]
      star_temperature : 	scalar
				Stellar temperature [K]
      star_logg : 		scalar
				Stellar log_10(g), where g is the surface gravity
      star_F_H : 		scalar
				Stellar metallicity [F/H]
      spot_temperature : 	scalar, optional
				Spot temperature [K]
    Returns
    -------
      wl:			array
				The Wavelength at which the SED is sampled. Units are [micron]
      sed :			array
				The SED of the star. Units are [W m**-2 micron**-1]
      L :			scalar
				The bolometric luminosity of the star. Units are [W]
      
				
    """
    
    ### This block reads BT-Settl text files 
    
    #wl_filename  = os.path.join(path,"lte{:03.0f}-{:01.1f}-{:01.1f}a+0.0.BT-Settl.spec.7.bz2".format(
      #np.round(star_temperature/100.0), star_logg, star_f_h))
    #import string
    #rule = string.maketrans('D', 'E')
    #raw_data = np.loadtxt(wl_filename, usecols=(0,1),
			  #converters={1:lambda val: float(val.translate(rule))})
    #wl  = raw_data[...,0]
    #sed = raw_data[...,1]
    ## Unit conversion
    #wl *= 1.0e-4 # [um]
    #sed = 10**(sed + 8.0 - 15.0) # [W m^-2 mu^-1]
    
    


####################### USING PHOENIX FITS FILES 
    
#    wl_filename  = os.path.join(path,"WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")
#    wl           = pyfits.getdata(wl_filename) # [angstrom]
#    sed_filename = os.path.join(path, 
#      "lte{:05.0f}-{:03.2f}-{:02.1f}.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits".format( np.round(star_temperature/100.0)*100.0, star_logg, star_f_h))
# 
#    sed = pyfits.getdata(sed_filename) # [erg/s/cm^2/cm]
#
#    hdulist      = pyfits.open(sed_filename)
#    radius       = hdulist[0].header['PHXREFF'] # [cm]
#    hdulist.close()
#
#    sed = pyfits.getdata(sed_filename) # [erg/s/cm^2/cm
#      
#    # Unit conversion
#    wl         *= 1.0e-4 # [um]
#    sed        *= 1.0e-7 # [W m^-2 mu^-1]
#    radius     *= 1.0e-2 # [m]
      
####################### USING PHOENIX BIN_SPECTRA BINARY FILES (h5)

    wl_filename = os.path.join(path,"lte{:03.0f}-{:01.1f}-{:01.1f}a+0.0.BT-Settl.h5".format(np.round(star_temperature/100.0), star_logg, star_f_h))
    
    f = tables.openFile(wl_filename, mode = "r")
    
    w1 = np.array(f.root.Spectrum.cmtdis)
    w2 = np.array(f.root.Spectrum.cmtber)

    wl = []
    for i in range(0,len(w2)):    
        if w1[i]>0:
            aa  = np.arange(w2[i],w2[i+1],w1[i]).tolist()
            wl = wl+aa
            
    wl = np.array(wl)
    sed = f.root.Spectrum.flux[0:len(wl)]
    f.close()
    
    # Unit conversion
    wl *= 1.0e-4 # [um]
    sed *= 1e-7 # [W m^-2 mu^-1]        

########################    
    
    # Calculate Luinosity for consistency check
    bolometric_flux        =  np.trapz(sed, x = wl) 			      # [W m**-2]
    bolometric_luminosity  =  4*np.pi * star_radius**2 * bolometric_flux 	# [W]
    sed                   *=  (star_radius/star_distance)**2 		      # [W/m^2/mu]
    
    return wl, sed, bolometric_luminosity
    
  def get_limbdarkening(self, filename):    
    lddata = np.loadtxt(filename)
    return lddata
    
def find_nearest(arr, value):
    arr = np.array(arr)
    # find nearest value in array
    idx = (abs(arr-value)).argmin()
    return arr[idx], idx

    