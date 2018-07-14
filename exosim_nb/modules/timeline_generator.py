
#Time line generator module
"""
Created on Wed Mar 11 12:51:06 2015

@author: Subi
"""
import time
import numpy as np
import quantities as pq
from ..lib import exolib
from exosim_nb.lib.exolib import exosim_msg

def run(opt, channel, planet):
  exosim_msg('Create signal-only timelines ... ')
  print "hello"
  st = time.time()
  # Estimate simulation length. Needs to be in units of hours.
  
  
  T14 =   planet.get_t14(planet.planet.i.rescale(pq.rad),
		 planet.planet.a.rescale(pq.m), 
		 planet.planet.P.rescale(pq.s), 
		 planet.planet.R.rescale(pq.m), 
		 planet.planet.star.R.rescale(pq.m)).rescale(pq.hour)
  
  total_observing_time = T14*(1.0+opt.timeline.before_transit()+opt.timeline.after_transit())
  time_at_transit      = T14*(0.5+opt.timeline.before_transit())
  
  NDR_time_estimate  = opt.timeline.exposure_time()/(opt.timeline.multiaccum()-1)
#  if NDR_time_estimate >= 0.01:
#      opt.timeline.frame_rate.val = 10000* 1/pq.s   
#  if NDR_time_estimate >= 0.1:
#      opt.timeline.frame_rate.val = 1000* 1/pq.s  
#  if NDR_time_estimate < 1.0:
#      opt.timeline.frame_rate.val = 100* 1/pq.s
#  if NDR_time_estimate >= 10.0:
#      opt.timeline.frame_rate.val = 10* 1/pq.s
#  opt.timeline.frame_rate.val = np.int(400./NDR_time_estimate) *  1/pq.s
  
  if opt.var_frame_rate == 1:
      opt.timeline.frame_rate.val = (400./NDR_time_estimate) 
      print "using frame rate of 400 frames/NDR"
  
#  if opt.timeline.frame_rate.val < 10 *1/pq.s:
#      opt.timeline.frame_rate.val = 10 * 1/pq.s
  
#   USE ONLY FOR GJ1214 AND VERY LONG INT TIME IN NIR PHOT  
#  if opt.timeline.frame_rate.val.magnitude < 100:
#        if opt.timeline.frame_rate.val.magnitude < 10:
#            print "frame rate < 10, therefore adopting min frame rate of 10"
#            opt.timeline.frame_rate.val = 10 *1/pq.s
#        else:
#            opt.timeline.frame_rate.val = 100 *1/pq.s
#            "frame rate < 100, therefore adopting min frame rate of 100"      

#  if opt.timeline.frame_rate.val.magnitude < 100:
#      opt.timeline.frame_rate.val = 100 *1/pq.s
#      print "frame rate < 100, therefore adopting min frame rate of 100" 

#  opt.timeline.frame_rate.val = 1 *1/pq.s
 

#  if channel.keys()[0]=='NIR Phot' or channel.keys()[0]=='AIRS CH0':
#      
#  if channel.keys()[0]=='NIR Phot':
#      if opt.astroscene.planet.val == 'GJ 1214 b':
#          opt.timeline.frame_rate.val = 10 *1/pq.s
      
#  opt.timeline.frame_rate.val = 10*1/pq.s
  print "NDR time estimate", NDR_time_estimate,
  print "selected frame rate", opt.timeline.frame_rate.val
  
  frame_time           = 1.0/opt.timeline.frame_rate.val   # Frame exposure, CLK
  print "T14", T14
      
  for key in channel.keys():
    
    # Having exposure_time here will allow to have different integration times 
    # for different focal planes.
    exposure_time   = opt.timeline.exposure_time() # Exposure time
    # Estimate NDR rates
    multiaccum     = opt.timeline.multiaccum()    # Number of NDRs per exposure
    allocated_time = (opt.timeline.nGND()+
		      opt.timeline.nNDR0()+
		      opt.timeline.nRST()) * frame_time
        
    print exposure_time
    print allocated_time
    print multiaccum
    NDR_time       = (exposure_time-allocated_time)/(multiaccum-1)
    
 
    # to make allocated time a 100th of NDR time
#    NDR_time += allocated_time     
    print "initial NDR time", NDR_time 
#    if NDR_time < 0.1 *pq.s:
#        print "Since true NDR time cannot be < 0.1 correcting to 0.1 sec"
#        NDR_time = 0.1*pq.s    
    print "overheads", opt.timeline.nGND()*frame_time, opt.timeline.nNDR0()*frame_time, opt.timeline.nRST()*frame_time
    
#    nNDR           = np.ceil(NDR_time/frame_time).astype(np.int).take(0)
    
    nNDR = np.round(NDR_time/frame_time).astype(np.int).take(0)
    
    # Estimate the base block of CLK cycles
#    if opt.timeline.nGND().magnitude > 0 or opt.timeline.nRST().magnitude > 0:
#        base = [opt.timeline.nGND().take(0), opt.timeline.nNDR0().take(0)]
#        for x in xrange(multiaccum-1): base.append(nNDR)
#        base.append(opt.timeline.nRST().take(0))
#    
    base = [opt.timeline.nGND().take(0), opt.timeline.nNDR0().take(0)]
    for x in xrange(multiaccum-1): base.append(nNDR)
    base.append(opt.timeline.nRST().take(0))  
    
    # Recalculate exposure time and estimates how many exposures are needed
    exposure_time = sum(base)*frame_time

    number_of_exposures = np.ceil(
      (total_observing_time/exposure_time).simplified.take(0)).astype(np.int)

    if opt.apply_LC ==0:
#        number_of_exposures = np.int(2000./exposure_time.magnitude)# use for esa only
        number_of_exposures = int(opt.n_exp)
    
#        if number_of_exposures  < 250:
#            number_of_exposures = 250
#        if number_of_exposures  > 500:
#            number_of_exposures = 500
        
    
    print "number of exposures projected", number_of_exposures

#    if number_of_exposures > 3000:
#        number_of_exposures  = 2000
#    if number_of_exposures < 250:
#        number_of_exposures  = 250
        
#    number_of_exposures  = 100
    print "number of exposures used", number_of_exposures

   
    total_observing_time = exposure_time*number_of_exposures
    frame_sequence=np.tile(base, number_of_exposures) # This is Nij
    time_sequence = frame_time * frame_sequence.cumsum() # This is Tij
    
    
    # Physical time of each NDR
    ndr_time = np.dstack([time_sequence[1+i::len(base)] \
      for i in xrange(multiaccum)]).flatten()*time_sequence.units
          
          
    print "actual NDR time", ndr_time[1]-ndr_time[0]
    
 
    opt.NDR_time = ndr_time[1]-ndr_time[0]
  
    # Number of frames contributing to each NDR
    ndr_sequence = np.dstack([frame_sequence[1+i::len(base)] \
      for i in xrange(multiaccum)]).flatten()
    # CLK counter of each NDR
    ndr_cumulative_sequence = (ndr_time/frame_time).astype(np.int).magnitude

    # Create the noise-less timeline
    channel[key].set_timeline(exposure_time,
			      frame_time,
			      ndr_time, 
			      ndr_sequence,
			      ndr_cumulative_sequence)
         
    # Apply lightcurve model
    cr    =  channel[key].planet.sed[channel[key].offs::channel[key].osf]
    cr_wl =  channel[key].planet.wl[channel[key].offs::channel[key].osf]
    
    
    isPrimaryTransit = True if opt.astroscene.transit_is_primary()=='True' else False
    
    channel[key].lc, channel[key].z, aa = planet.get_light_curve(channel, cr, cr_wl, 
							channel[key].ndr_time, 
							time_at_transit, 
							isPrimaryTransit)
    channel[key].ld1 = aa[1]
    channel[key].ld2 = aa[2] 

        
  exosim_msg(' - execution time: {:.0f} msec.\n'.format(
  (time.time()-st)*1000.0))
  return frame_time, total_observing_time, exposure_time

 
 
      
    