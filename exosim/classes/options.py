import xml.etree.ElementTree as ET
import numpy as np
import sys, os
from   scipy import constants

# Constants definition
jupiter_radius2m	= 6.955e7
sun_radius2m		= 6.955e8
pc2m			= constants.parsec
au2m			= constants.au
days2s			= constants.day

EXOSIM_DEFAULTS        = 'exosim_defaults_2.xml'
WRONGUNITS = -1
FIELDNOTFOUND = -2

def error_msg(err_code):
    sys.stderr.write("Non valid option or physical units given - error code {:d}\n".format(err_code))
    sys.exit(0)
      

class Token(object):
  val          = None
  units        = None
  comment      = None
  transmission = None
  emission     = None
  emissivity   = None
  temperature  = None
  tag_type     = None
  tag          = None
  
  def __init__(self, parent=None, child_name=None, units=None, val=None):
    
    if parent != None : 
      self.get(parent=parent, child_name=child_name, units=units)
    elif val != None:
      self.val = val
      self.units = units
      
    
  def get(self, parent, child_name = None, units=None):
    child = parent.find(child_name) if child_name else parent 
    if child == None: error_msg(FIELDNOTFOUND)
    self.tag = child.tag
 
    try:
      self.val   = np.float64(child.get("val"))
    except ValueError:
      self.val = child.get("val")
    
    self.units = child.get("units")
    self.comment = child.get("comment")
    self.transmission = child.get("transmission")
    self.emissivity = child.get("emissivity")
    try:
      self.temperature = np.float64(child.get("T"))
    except ValueError:
      self.temperature = np.nan
  
    self.tag_type = child.get("type")
    if units != None : self.check_units(units)
    
  def units_conversion(self, from_units, to_units, from2to):
    if self.units == from_units:
      self.units = to_units
      self.val *= from2to
    else:
      error_msg(WRONGUNITS)
  
  def check_units(self, units):
    if self.units != units: error_msg(WRONGUNITS)
    
  
class Options(object):
  def __init__(self, filename=None):
    """
    
    """
    if filename == None: filename = EXOSIM_DEFAULTS
    tree = ET.parse(filename)
    root = tree.getroot()
    self.root = root
    
    self.root_directory = os.path.dirname(os.path.split(__file__)[0])
    self.initialize_common_block(root)
    self.initialize_planet(root)
    self.initialize_star(root)  
    self.initialize_common_optics(root)
    self.initialize_channels(root)
  
  def initialize_common_block(self, root):
    common_blk = root.find("common")
    self.common_exosym_path 	= Token(common_blk, 'exosym_path')
    logbinres  			= Token(common_blk, "logbinres")
    wl_min          		= Token(common_blk, 'wl_min', units='micron')
    wl_max          		= Token(common_blk, 'wl_max', units='micron')
    wl_delta			= wl_min.val/logbinres.val
    n_wl                   	= np.round(1.0+(wl_max.val-wl_min.val)/wl_delta).astype(np.int)
    self.common_wl         	= Token(val=np.linspace(wl_min.val,
							 wl_max.val,
						         n_wl),
						         units = wl_min.units)
    
      
  def initialize_planet(self, root):
    # Get planet parameters
    planet = root.find("planet")
    self.planet_radius = Token(planet, 'radius')
    self.planet_radius.units_conversion('Rjup', 'm', jupiter_radius2m)
    self.planet_period = Token(planet, 'period')
    self.planet_period.units_conversion('days', 's', days2s)
   
    self.planet_a = Token(planet, 'semimajor_axis')
    self.planet_a.units_conversion('AU', 'm', au2m)
    self.planet_inclination = Token(planet, 'inclination')
    self.planet_inclination.units_conversion('deg', 'rad', np.deg2rad(1))
    self.planet_eccentricity = Token(planet, 'eccentricity')
    self.planet_eccentricity.units_conversion('deg', 'rad', np.deg2rad(1))
    self.planet_omega = Token(planet, 'omega')
    self.planet_omega.units_conversion('deg', 'rad', np.deg2rad(1))
    self.planet_name = Token(planet, 'name') 
    self.planet_contrast_path = Token(planet, "contrast")
    self.planet_transit = Token(planet, "transit")
    
  def initialize_star(self, root):
    # Get planet parameters
    star = root.find("star")
    #
    self.star_radius = Token(star, 'radius')
    self.star_radius.units_conversion('Rsun', 'm', sun_radius2m)
    self.star_period = Token(star, 'period')
    self.star_period.units_conversion('days', 's', days2s)
    self.star_distance = Token(star, 'distance')
    self.star_distance.units_conversion('pc', 'm', pc2m)
    self.star_temperature = Token(star, 'temperature')
    self.star_temperature.units_conversion('K', 'K', 1.0)
    self.star_logg = Token(star, 'logg')
    self.star_logg.units_conversion('', '', 1.0)
    self.star_f_h = Token(star, 'F_H')
    self.star_f_h.units_conversion('', '', 1.0)
    self.star_sed_path = Token(star, 'data_path')
    self.star_limb_darkening_path = Token(star, 'limb_darkening')
    
  def initialize_common_optics(self, root):
    parent = root.find("common_optics")
    self.common_optic = {'optical_surface': [Token(ch) for ch in parent.findall('optical_surface')]} 
    self.common_optic['diameter'] = Token(parent, 'effective_diameter')
    self.common_optic['diameter'].units_conversion('m', 'm', 1.0)
    
  def initialize_channels(self, root):
  
    self.channel = {}
    for ch in root.findall('channel'):
      tk = Token(ch)
      self.channel[tk.val] = {'optical_surface': [Token(i) for i in ch.findall('optical_surface')]}
      self.channel[tk.val]['ld'] = Token(ch, 'ld')
      self.channel[tk.val]['array_geometry'] = Token(ch, 'array_geometry')
      self.channel[tk.val]['pixel_size'] = Token(ch, 'pixel_size')
      self.channel[tk.val]['pixel_size'].units_conversion('micron', 'micron', 1.0)
      self.channel[tk.val]['wfno'] = Token(ch, 'wfno')
      self.channel[tk.val]['plate_scale'] = Token(ch, 'plate_scale')
      self.channel[tk.val]['osf'] = Token(ch, 'osf')
      self.channel[tk.val]['kernel_osf'] = Token(ch, 'kernel_osf')
      self.channel[tk.val]['psf_osf'] = Token(ch, 'psf_osf')
      self.channel[tk.val]['ad_osf'] = Token(ch, 'ad_osf')

      
 	
  def error(self, err_code):
    sys.stderr.write("Non valid option or physical units given - error code {:d}\n".format(err_code))
    sys.exit(0)
    


