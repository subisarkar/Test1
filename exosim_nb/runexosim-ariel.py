import numpy           as     np
import sys, time, os, glob
import exosim



if __name__ == "__main__":

  xml_filename = 'exosim_ariel_mcr.xml'
  # Run simulation
  exosim.exolib.exosim_msg('Reading options from file ... \n')
  opt = exosim.Options(filename=xml_filename).opt #, default_path = exosim.__path__[0]).opt

  # modify_opt(opt)

  star, planet, zodi, channel = exosim.run_exosim(opt)
  
  
  
  
  