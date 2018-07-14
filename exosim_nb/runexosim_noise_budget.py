import numpy           as     np
import sys, time, os
import exosim_nb as exosim
from exosim_nb.lib.exolib import exosim_msg
import matplotlib.pyplot as plt
import quantities as pq
from exodata import astroquantities as aq
import gc

def run_exosim(seed, fake_system, opt=None):
  print "RANDOM SEED", seed
  np.random.seed(seed)


  star, planet = exosim.modules.astroscene.run(fake_system, opt)
  
  exosim_msg(' Stellar SED: {:s}\n'.format(os.path.basename(star.ph_filename)))
  exosim_msg(' Star luminosity {:s}\n'.format(star.luminosity))
  
  #Instanciate Zodi
  zodi = exosim.classes.zodiacal_light(opt.common.common_wl, level=1.0)
  
  exosim.exolib.sed_propagation(star.sed, zodi.transmission)
  #Run Instrument Model
  channel = exosim.modules.instrument.run(opt, star, planet, zodi)
  #Create Signal timelines
  frame_time, total_observing_time, exposure_time=  exosim.modules.timeline_generator.run(opt, channel, planet)
  #Generate noise timelines
  exosim.modules.noise.run(opt, channel, frame_time, total_observing_time, exposure_time)
   
  exosim.modules.output.run(opt, channel, planet) 
  
def run_exosim2(seed, fake_system, pl, opt=None):
    
  print "RANDOM SEED", seed
  np.random.seed(seed)
  
  opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_bright.csv"     
  
  if pl == 'GJ 1214 b': 
      opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_dim.csv"        
  
  star, planet = exosim.modules.astroscene.run(fake_system, opt)
  
  exosim_msg(' Stellar SED: {:s}\n'.format(os.path.basename(star.ph_filename)))
  exosim_msg(' Star luminosity {:s}\n'.format(star.luminosity))
  
  #Instanciate Zodi
  zodi = exosim.classes.zodiacal_light(opt.common.common_wl, level=1.0)
  
  exosim.exolib.sed_propagation(star.sed, zodi.transmission)
  #Run Instrument Model
  channel = exosim.modules.instrument.run(opt, star, planet, zodi)
  #Create Signal timelines
  frame_time, total_observing_time, exposure_time=  exosim.modules.timeline_generator.run(opt, channel, planet)
  #Generate noise timelines
  noise1, TL1 = exosim.modules.noise.run(opt, channel, frame_time, total_observing_time, exposure_time)

  ######### repeat

  
  opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_bright_no_RPE.csv"  
  if pl == 'GJ 1214 b': 
      opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_dim_no_RPE.csv"    
                 
  star, planet = exosim.modules.astroscene.run(fake_system, opt)
  
  exosim_msg(' Stellar SED: {:s}\n'.format(os.path.basename(star.ph_filename)))
  exosim_msg(' Star luminosity {:s}\n'.format(star.luminosity))
  
  #Instanciate Zodi
  zodi = exosim.classes.zodiacal_light(opt.common.common_wl, level=1.0)
  
  exosim.exolib.sed_propagation(star.sed, zodi.transmission)
  #Run Instrument Model
  channel = exosim.modules.instrument.run(opt, star, planet, zodi)
  #Create Signal timelines
  frame_time, total_observing_time, exposure_time=  exosim.modules.timeline_generator.run(opt, channel, planet)
  #Generate noise timelines
  noise2, TL2 = exosim.modules.noise.run(opt, channel, frame_time, total_observing_time, exposure_time)

  ######### feed output

  exosim.modules.noise.run2(opt, channel, noise1, noise2, TL1, TL2)
  
  exosim.modules.output.run(opt, channel, planet) 
  
  

if __name__ == "__main__":
  
  exosim_msg('Reading options from file ... \n')
 
  seed= [35644, 88008, 75816,  6250, 78083, 93567,   990, 25692, 78722,
        6374, 94116, 21669,  3042, 35373, 58928, 69149, 64405, 26163,
       28229, 15492]
       
  fake_planet = {'R': 0.5*aq.R_e, 'M': 0.5*aq.M_e, 'T': 1000.0*aq.K, 'i': 90.0*aq.deg, 
  'e': 0.0, 'a': 1.0*aq.au, 'P': 365.0*aq.day, 'albedo':0.3, 'name': 'Fake'}
  fake_star = {'R': 1.0*aq.R_s, 'M': 1.0*aq.M_s, 'T': 5800.0*aq.K, 'd': 12.0*aq.pc, 
  'logg':4.4, 'Z': 0.0, 'name':'Fake'}
  
  fake_planet = {'R': 0.5*aq.R_e, 'M': 0.5*aq.M_e, 'T': 1000.0*aq.K, 'i': 90.0*aq.deg, 
  'e': 0.0, 'a': 1.0*aq.au, 'P': 365.0*aq.day, 'albedo':0.3, 'name': 'Fake'}
  fake_star = {'R': 0.978486*aq.R_s, 'M': 0.974764*aq.M_s, 'T': 5760.0*aq.K, 'd': 578.72*aq.pc, 
  'logg':4.4, 'Z': 0.0, 'name':'Fake'}  
  
  
  fake_switch = 0 # 1 to overide exodata choices and planet
  
  
  fake_system = [fake_planet, fake_star, fake_switch]  
  
  apply_LC = 0

  # apply QE grid in ExoSim
  apply_QE = 1
  use_random_QE = 0
  # apply flat and bkg subtraction in ExoSim  
  apply_flat = 1
  
  No_real = 1
  
  
  
  
  
  
  
  
  use_hard_drive = 1
  
  int_time = 120.
  
  var_frame_rate = 0
  
#  n_exp = 2250.*2

  n_exp = 20000 # hd 209 x 5 sims for fgs 2, fgs 1, nir phot, nirspec = hd209
#  
#  n_exp = 10000  # for GJ 1214
  
#  n_exp = 5000  # nirphot - GJ 1214
 



  RPE_var = 0
  Inc_PSD = 0
  no_RPE = 0
  
  PN = 0
    
  a = 'no_RPE_var'
  b = 'normal_PSD'
  c ='full_RPE'
  
  if RPE_var == 1:
      a = 'with_RPE_var'
  if Inc_PSD == 1:
      b = 'Inc_PSD'
  if no_RPE == 1:
      c = 'no_RPE'
  
  label = '%s_%s_%s'%(a,b,c)
  

  use_sat = True
  use_T14 = False
  
  
  if PN == 1:
      no_group = ['jitterless']
      label = ''
      RPE_var = 0
      
  else:
      no_group = ['jitterful']
      
                        
#  
#  ch_list= ['AIRS CH0']
#  ch_list= ['FGS Red']
#  ch_list= ['NIR Spec']
#  ch_list= ['NIR Phot']
  ch_list= ['AIRS CH0', 'AIRS CH1']
#  ch_list= ['AIRS CH1']
#  ch_list= [ 'FGS Prime']
                        
#  ch_list= [ 'NIR Spec' , 'FGS Red',  'FGS Prime', 'NIR Phot']          

#  ch_list= ['NIR Phot']
#  ch_list= ['AIRS CH0', 'AIRS CH1','NIR Spec' , 'FGS Red',  'FGS Prime'] 
    
#  for pl in ['GJ 1214 b']:
      
  for pl in ['HD 209458 b']:
#  for pl in ['HD 219134 b']:
#  for pl in ['GJ 1214 b']:
#  for pl in ['Fake']:

      
      seed = np.random.uniform(0,100000000)
    
                  
#              
#      for noise_group in ['jitterful', 'jitterless']:
#      for noise_group in ['jitterless']:
      for noise_group in no_group:
              
#      for noise_group in ['full_transit']:

          
          if noise_group == 'jitterless':
#                    i_list = [2, 3, 4, 5, 6]
                i_list = [2]
#                root = '/Users/user1/ExoSimOutput_NB'                
#                root = '/Users/user1/ExoSimOutput_PN_5'                
#                root = '/Users/user1/ExoSimOutput_PN_5'                
#
#                root = '/Volumes/My Passport for Mac/ExoSimOutput_76'    
#                root = '/Users/user1/ExoSimOutput_FGS2_PN'

                root = '/Users/user1/ExoSimOutput_PN_166/%s'%(label)
                
              
          elif noise_group == 'jitterful':
#                  i_list = [0, 7, 8]
              i_list = [9]
#              root = '/Users/user1/ExoSimOutput_NB' 
              
#              root = '/Users/user1/ExoSimOutput_Ch1_JN_test'   
#              root = '/Users/user1/ExoSimOutput_NP_JN_noRPE'   
              
##                  root = '/Users/user1/ExoSimOutput_JN_5_noRPE'   
#              root = '/Users/user1/ExoSimOutput_JN_6_highPSD_noRPE'                
              root = '/Users/user1/ExoSimOutput_JN_166/%s'%(label)
              
 
          elif noise_group == 'full_transit':
                i_list = [0]
 
#                root = '/Users/user1/ExoSimOutput_nb' 
#                root = '/Volumes/My Passport for Mac/ExoSimOutput_GJ1214-new'                 
                    


          for j in range(No_real):    


              seed = int(np.random.uniform(0,100000000))
              print 'random seed', seed

                 
              for ch in ch_list :
            
                  if ch == 'AIRS CH0':
                      a = 'exosim_ariel_mcr_modified_AIRS_Ch0_Euro.xml'
                  if ch == 'AIRS CH1':
                      a = 'exosim_ariel_mcr_modified_AIRS_Ch1_Euro.xml'
                  if ch == 'NIR Spec':
                      a = 'exosim_ariel_mcr_modified_NIR_Euro_18.xml'
                  if ch == 'NIR Phot':
                      a = 'exosim_ariel_mcr_modified_NIR_phot_Euro_18.xml'
                  if ch == 'FGS Red':
                      a = 'exosim_ariel_mcr_modified_FGS_red_Euro_18.xml'
                  if ch == 'FGS Prime':              
                      a = 'exosim_ariel_mcr_modified_FGS_prime_Euro_18.xml'    

                                
                  for i in i_list:
                      gc.collect()
                      print "Noise Group:  ", noise_group
                      print "Realization:  ", j, " of ", No_real
                      print "Noise option:  ", i
                      
                      opt = exosim.Options(filename=a).opt #, default_path = exosim.__path__[0]).opt
                                  
                      opt.astroscene.planet.val =  pl
                      
                      if no_RPE == 1:
                          opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_bright_no_RPE.csv"     
                      else:
                          opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_bright.csv"
                      
                      opt.timeline.frame_rate.val = 10*1/pq.s    
                      
        #              if opt.astroscene.planet.val == 'GJ 1214 b':
        #                            opt.timeline.frame_rate.val = 1*1/pq.s     #okay for long NDR since supersupaled      
        #              if opt.astroscene.planet.val == '55 Cnc e':
        #        #          opt.timeline.exposure_time.val = 0.79*pq.s
        #                  opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_bright.csv"       
                      if opt.astroscene.planet.val == 'GJ 1214 b':
                #          opt.timeline.exposure_time.val = 45.03*pq.s
                          if no_RPE == 1:
                              opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_dim_no_RPE.csv"    
                          else:
                              opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd_dim.csv"   
                 
        #                  opt.timeline.frame_rate.val = 10 *1/pq.s
                #      opt.aocs.PointingModel.val="__path__/data/ariel/pointing_model_psd.csv"      
                
                      opt.timeline.nGND.val = 0.0 *pq.dimensionless
                      opt.timeline.nNDR0.val = 0.0*pq.dimensionless
                      opt.timeline.nRST.val = 0.0*pq.dimensionless
                      
        #              opt.astroscene.planetCR.val="__path__/data/planetary/test/Gliese 1214 b_pri.dat"
        #              print "dispersion file", opt.channel[0].dispersion.path
                      
                      
                      opt.common.ExoSimOutputPath.val="%s/%s_%s--%s/"%(root,pl, ch, noise_group)
                      if not os.path.exists(opt.common.ExoSimOutputPath.val):
                          os.makedirs(opt.common.ExoSimOutputPath.val)
                          
                  
                
                      if i==0 : # all noise 
                          opt.noise.EnableReadoutNoise.val = 'True'
                          opt.noise.EnableShotNoise.val = 'True'
                          opt.noise.EnableSpatialJitter.val= 'True'
                          opt.noise.EnableSpectralJitter.val= 'True'   
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1
                          process_switch = 1
                          diff = 0
                          jitter_switch = 1 
               
                      elif i == 1: # All shot
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'True'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'                  
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1
                          process_switch = 1
                          diff = 0
                          jitter_switch = 0                 
                          
                      elif i == 2: # source noise only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'True'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'    
                          emm_switch = 0
                          zodi_switch = 0
                          dc_switch = 0
                          source_switch = 1
                          process_switch = 1
                          diff = 0
                          jitter_switch = 0
                                  
                      elif i == 3: # dark current noise only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'True'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'    
                          emm_switch = 0
                          zodi_switch = 0
                          dc_switch = 1
                          source_switch = 0
                          process_switch = 1
                          diff = 1
                          jitter_switch = 0        
                #                  
                      elif i == 4: # zodi noise only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'True'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'    
                          emm_switch = 0
                          zodi_switch = 1
                          dc_switch = 0
                          source_switch = 0
                          process_switch = 1
                          diff = 1
                          jitter_switch = 0 
                         
                      elif i == 5:# emission noise only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'True'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'    
                          emm_switch = 1
                          zodi_switch = 0
                          dc_switch = 0
                          source_switch = 0
                          process_switch = 1
                          diff = 1
                          jitter_switch = 0
                          
                      elif i==6 : # RN only
                          opt.noise.EnableReadoutNoise.val = 'True'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'    
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1
                          process_switch = 1
                          diff = 0
                          jitter_switch = 0
                          
                      elif i==7 : # spatial jitter only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'True'
                          opt.noise.EnableSpectralJitter.val= 'False'
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1  
                          process_switch = 1
                          diff = 0                  
                          jitter_switch = 1                  
                                    
                      elif i==8 : # spectral jitter only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'True'         
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1  
                          process_switch = 1
                          diff = 0
                          jitter_switch = 1
                          
                      elif i==9 : # combined jitter only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'True'
                          opt.noise.EnableSpectralJitter.val= 'True' 
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1  
                          process_switch = 1
                          diff = 0
                          jitter_switch = 1
                          

                      elif i==9 : # combined jitter only
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'True'
                          opt.noise.EnableSpectralJitter.val= 'True' 
                          emm_switch = 0
                          zodi_switch = 0
                          dc_switch = 1
                          source_switch = 1  
                          process_switch = 1
                          diff = 0
                          jitter_switch = 1                          
                          
        
                          
                      elif i==10 : # no noise- no background
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'   
                          emm_switch = 0
                          zodi_switch = 0
                          dc_switch = 0
                          source_switch = 1
                          process_switch = 1
                          diff = 0
                          jitter_switch = 0  
                          apply_QE = 0
                          apply_flat = 0                       
                          
                          
                      elif i==11 : # no noise - all background
                          opt.noise.EnableReadoutNoise.val = 'False'
                          opt.noise.EnableShotNoise.val = 'False'
                          opt.noise.EnableSpatialJitter.val= 'False'
                          opt.noise.EnableSpectralJitter.val= 'False'   
                          emm_switch = 1
                          zodi_switch = 1
                          dc_switch = 1
                          source_switch = 1
                          process_switch = 1
                          diff = 0
                          jitter_switch = 0                  
                          
                       
                      opt.emm_switch = emm_switch
                      opt.zodi_switch = zodi_switch
                      opt.dc_switch = dc_switch
                      opt.source_switch = source_switch
                      opt.diff = diff
                      opt.jitter_switch = jitter_switch
                      
                      opt.apply_QE = apply_QE
                      opt.use_random_QE = use_random_QE
                      opt.apply_flat = apply_flat
                      
                      opt.noise_option = i
                      opt.noise_group = noise_group
                      opt.apply_LC = apply_LC
                      
                      opt.n_exp = n_exp
                      opt.int_time = int_time
                      opt.use_sat = use_sat
                      opt.use_T14 = use_T14
                      opt.var_frame_rate = var_frame_rate
                      
                      opt.Inc_PSD = Inc_PSD
                      
                      
                      if RPE_var ==0:
                     
                          run_exosim(seed, fake_system, opt)
                          
                      elif RPE_var ==1:
                          
                          run_exosim2(seed, fake_system, pl, opt)
                          
                      

                   
         
