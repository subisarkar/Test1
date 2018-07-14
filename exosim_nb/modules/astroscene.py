from   ..lib             import exolib
from   ..classes.star    import Star
from   ..classes.planet  import Planet

from   ..lib.exolib import exosim_msg

import os, time

import exodata
import exodata.astroquantities as aq


class fake_star_class():
    def __init__(self, fake_star):
        self.R = fake_star['R'] 
        self.M = fake_star['M'] 
        self.d = fake_star['d']
        self.T = fake_star['T'] 
        self.Z = fake_star['Z'] 
        self.logg = fake_star['logg'] 
        self.name = fake_star['name'] 
        
    def calcLogg(self):
        return self.logg
        
class fake_exosystem():
    def __init__(self, fake_planet, fake_star):
        self.i = fake_planet['i']
        self.P = fake_planet['P']      
        self.a = fake_planet['a']
        self.e = fake_planet['e'] 
        self.R = fake_planet['R'] 
        self.M = fake_planet['M'] 
        self.T = fake_planet['T'] 
        self.name = fake_planet['name'] 
        self.albedo = fake_planet['albedo']
        self.star = fake_star_class(fake_star)

def run(fake_system, opt):
  
  # Load the openexoplanet catalogue
  exosim_msg('Run astroscene ...')
  st = time.time()
  
  if opt.astroscene.OpenExoplanetCatalogue():
    
    exocat_path = opt.astroscene.OpenExoplanetCatalogue().replace('__path__', opt.__path__)
    exocat = exodata.OECDatabase(os.path.expanduser(exocat_path))
  else:
    exocat = exodata.load_db_from_url()
  
  fake_system= [0,0,False]

  if fake_system[2] == False:  
      exosystem = exocat.searchPlanet(opt.astroscene.planet())
  else:
      exosystem = fake_exosystem(fake_system[0],fake_system[1])
  

  star = Star(exosystem.star, 
	      opt.astroscene.StarSEDPath().replace('__path__', opt.__path__),
	      use_planck_spectrum=opt.astroscene.use_planck_spectrum().lower()=='true')

#  star = Star(exosystem.star, 
#	      opt.astroscene.StarSEDPath().replace('__path__', opt.__path__),
#	      use_planck_spectrum=True)


  star.sed.rebin(opt.common.common_wl)
  #star.get_limbdarkening(opt.star_limb_darkening_path.val.replace('$root$', opt.common_exosym_path.val))  
  cr_path = opt.astroscene.planetCR().replace('__path__', opt.__path__)
  limb_darkenning_path = opt.astroscene.StarLimbDarkening().replace('__path__', opt.__path__)
  planet = Planet(exosystem, cr_path, limb_darkenning_path)
  planet.cr.rebin(opt.common.common_wl) 
  
  planet.get_t14(planet.planet.i.rescale(aq.rad),
		 planet.planet.a.rescale(aq.m), 
		 planet.planet.P.rescale(aq.s), 
		 planet.planet.R.rescale(aq.m), 
		 planet.planet.star.R.rescale(aq.m))
  
  planet.get_orbital_phase(planet.t14, planet.planet.P.rescale(aq.s))
  planet.eccentric(planet.phi,
		   planet.planet.i.rescale(aq.rad), 
		   0.0, #eccetricity 
		   0.0, # argument of periastron
                   planet.planet.a.rescale(aq.m), 
                   planet.planet.P.rescale(aq.s), 
                   planet.planet.star.R.rescale(aq.m))
  
  exosim_msg(' - execution time: {:.0f} msec.\n'.format((time.time()-st)*1000.0))
  return star, planet
  
if __name__ == "__main__":
  
  exolib.exosim_error("This module not made to run stand alone")
    
    