#from   lib             import exolib
#from   classes.options import Options
#from   classes.zodi    import zodiacal_light
#from   modules         import astroscene, instrument
import numpy           as     np
import pylab           as     pl
import sys, time
import exosim
# hello
data = {}
opt = 0 
def run_exosim(parameters=None):
  global data 
  global opt
  

  exosim.exolib.exosim_msg('Reading options from file ... \n')
  opt = exosim.Options(parameters)

  exosim.exolib.exosim_msg('Run astroscene ... ')
  st = time.clock()
  star, planet = exosim.modules.astroscene.run(opt)
  exosim.exolib.exosim_msg(' - execution time: {:.0f} msec.\n'.format((time.clock()-st)*1000.0))
  
  exosim.exolib.exosim_msg('Instantiate Zodi ... ')
  st = time.clock()
  zodi = exosim.classes.zodiacal_light(opt.common_wl.val, level=1.0)
  exosim.exolib.exosim_msg(' - execution time: {:.0f} msec.\n'.format((time.clock()-st)*1000.0))

  exosim.exolib.sed_propagation(star.sed, zodi.transmission)
  
  exosim.exolib.exosim_msg('Run instrument model ... ')
  st = time.clock()
  channel = exosim.modules.instrument.run(opt, star, planet, zodi)
  exosim.exolib.exosim_msg(' - execution time: {:.0f} msec.\n'.format((time.clock()-st)*1000.0))
  
  data['qstar'], data['qplanet'], data['qzodi'], data['channel'] = star, planet, zodi, channel
  
if __name__ == "__main__":
  
  run_exosim()
  
  #########################################
  print 'Star luminosity {:.2e} [w]'.format(data['qstar'].luminosity)
  
  data['qplanet'].get_light_curve(data['qplanet'].z, 0.0, 0.0, 0.1, primary_transit = False)
  
  pl.figure(1)
  pl.subplot(3,3,1)
  #pl.plot(data['qstar'].ph_sed.wl,data['qstar'].ph_sed.sed)
  pl.plot(data['qstar'].sed.wl, data['qstar'].sed.sed)
  pl.plot(data['channel']['SWIR'].star.wl, data['channel']['SWIR'].star.sed)
  pl.plot(data['channel']['MWIR'].star.wl, data['channel']['MWIR'].star.sed)
  pl.subplot(3,3,2)
  pl.plot(data['qplanet'].cr.wl, data['qplanet'].cr.sed)
  pl.subplot(3,3,3)
  pl.plot(data['qplanet'].phi, data['qplanet'].lc, 'r')
  pl.subplot(3,3,4)
  pl.plot(data['qzodi'].sed.wl, data['qzodi'].sed.sed, 'r')
  pl.plot(data['channel']['SWIR'].zodi.wl, data['channel']['SWIR'].zodi.sed)
  pl.plot(data['channel']['MWIR'].zodi.wl, data['channel']['MWIR'].zodi.sed)
  pl.subplot(3,3,5)
  pl.plot(data['channel']['SWIR'].emission.wl, data['channel']['SWIR'].emission.sed)
  pl.plot(data['channel']['MWIR'].emission.wl, data['channel']['MWIR'].emission.sed)
  pl.subplot(3,3,6)
  pl.plot(data['channel']['SWIR'].transmission.wl, data['channel']['SWIR'].transmission.sed)
  pl.plot(data['channel']['MWIR'].transmission.wl, data['channel']['MWIR'].transmission.sed)
  pl.subplot(3,3,7)
  pl.imshow(data['channel']['SWIR'].fp, 
    extent=[data['channel']['SWIR'].wl_solution.min(),data['channel']['SWIR'].wl_solution.max(),0,1])
  pl.subplot(3,3,8)
  pl.imshow(data['channel']['MWIR'].fp, 
    extent=[data['channel']['MWIR'].wl_solution.min(),data['channel']['MWIR'].wl_solution.max(),0,1])

  
  pl.show()