from   ..lib             import exolib
from   ..classes.star    import Star
from   ..classes.planet  import Planet

def run(opt):
  star = Star(opt.star_sed_path.val.replace('$root$', opt.common_exosym_path.val),
              opt.star_distance.val,
              opt.star_temperature.val,
              opt.star_logg.val, 
              opt.star_f_h.val,
	      opt.star_radius.val)     
  
  star.sed.rebin(opt.common_wl.val) 
  #star.get_limbdarkening(opt.star_limb_darkening_path.val.replace('$root$', opt.common_exosym_path.val))
  
  
  
  
  planet = Planet(opt.planet_contrast_path.val.replace('$root$', opt.common_exosym_path.val))
  planet.cr.rebin(opt.common_wl.val) 
  
  planet.get_t14(opt.planet_inclination.val, 
  opt.planet_a.val, 
  opt.planet_period.val, 
  opt.planet_radius.val, 
  opt.star_radius.val)
  
  planet.get_orbital_phase(planet.t14, opt.planet_period.val)
  planet.eccentric(planet.phi,
  opt.planet_inclination.val, 
  opt.planet_eccentricity.val, 
  opt.planet_omega.val, 
  opt.planet_a.val, 
  opt.planet_period.val, 
  opt.star_radius.val)
    
  return star, planet
  
if __name__ == "__main__":
  
  exolib.exosim_error("This module not made to run stand alone")
    
    