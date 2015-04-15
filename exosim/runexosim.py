#from   lib             import exolib
#from   classes.options import Options
#from   classes.zodi    import zodiacal_light
#from   modules         import astroscene, instrument
import numpy           as     np
import pylab           as     pl
import sys, time
import exosim
import matplotlib.pyplot as plt


# hello
data = {}
opt = 0 
channel = 0
def run_exosim(parameters=None):
  global data 
  global opt
  global channel
#  timeline = exosim.modules.timeline_generator.run()

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
  
#####
      
  exosim.lib.exolib.exosim_msg('Create jittered timeline ... ')
  st = time.clock()
  
  exosim.modules.timeline_generator.run(opt, channel)

  exosim.lib.exolib.exosim_msg(' - execution time: {:.0f} msec.\n'.format((time.clock()-st)*1000.0))

  
if __name__ == "__main__":
      
  run_exosim()
  
  #########################################
  print 'Star luminosity {:.2e} [w]'.format(data['qstar'].luminosity)
  
  data['qplanet'].get_light_curve(data['qplanet'].z, 0.0, 0.0, 0.1, primary_transit = False)
    

  
  pl.figure(2)
  pl.subplot(3,3,1)
#  pl.plot(data['qstar'].ph_sed.wl,data['qstar'].ph_sed.sed)
  pl.plot(data['qstar'].sed.wl, data['qstar'].sed.sed)
  pl.subplot(3,3,2)
  pl.plot(data['qplanet'].cr.wl, data['qplanet'].cr.sed)
  pl.subplot(3,3,3)
  pl.plot(data['qplanet'].phi, data['qplanet'].lc, 'r') 
  pl.subplot(3,3,4)
  pl.plot(data['qzodi'].sed.wl, data['qzodi'].sed.sed, 'r')

  
  i=0
  for key in opt.channel.keys():
      i += 1
  
      pl.figure(2)
      pl.subplot(3,3,1)
      pl.plot(data['channel'][key].star.wl, data['channel'][key].star.sed)
      pl.subplot(3,3,4)
      pl.plot(data['channel'][key].zodi.wl, data['channel'][key].zodi.sed)
      pl.subplot(3,3,5)
      pl.plot(data['channel'][key].emission.wl, data['channel'][key].emission.sed)
      pl.subplot(3,3,6)
      pl.plot(data['channel'][key].transmission.wl, data['channel'][key].transmission.sed)
      pl.subplot(3,3,6+i)
      pl.imshow(data['channel'][key].fp, 
        extent=[data['channel'][key].wl_solution.min(),data['channel'][key].wl_solution.max(),0,1])
    
 
#  plt.figure(45)
#  plt.plot(channel['WFC3IR'].fp[32*channel['WFC3IR'].osf])
      
  pl.show()

#tl = []
#plt.figure(878)
#for i in range (0,data['timeline'].shape[2]):
#    plt.plot(data['timeline'][...,100,i])
#    tl.append(data['timeline'][...,i].sum())
#print "mean count", np.mean(tl)
#print "sd of count", np.std(tl)


#exosim.lib.exolib.animate(data['channel'][opt.channel.keys()[0]].timeline)


#  need to animate, check why only short section illumimated, generalize for MWIR