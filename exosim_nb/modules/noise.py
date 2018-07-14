import time
import numpy as np
import quantities as pq
import scipy.weave
import matplotlib.pyplot as plt

from exosim_nb.lib.exolib import exosim_msg
from ..lib import exolib

def c_create_jitter_noise(fp, osf, 
			  ndr_sequence, ndr_cumulative_sequence, frame_osf,
			  x_jit, y_jit, 
			  x_offs = 0, y_offs = 0):
  """
  c_create_timeline: create a timeline cube
  
  Parameters
  ----------
  fp     : 2D Array
	     The oversampled focal plane
  osf    : scalar
	     The oversample factor applied. Relative to physical focal plane
  ndr_sequence: Array
	   number of frames in each NDR
  ndr_cumulative_sequence: Array
	   number of frames since simulation start corresponding to each NDR 
  frame_osf: scalar
	  number of samples in jitter timeline for each frame.
  x_jit  : Array like 
	    Jitter in the spectral direction
	    Has to be a ndarray of type np.int64
	    with shape (nndr, nsndr) where nndr is the number of desired ndr_offs
	    and sndr is the number of desired sub_ndr for RPE sampling.
  y_jit  : Array like
            Similar to x_jit, but contains jitter in the spatial direction
  x_offs : Scalar
	    Offset in fp in the spectral direction
  y_offs : Scalar
	    Offset in fp in the spatial direction
  
  Returns
  -------
  time_line: 3D Array
	      The timeline cube of shape (nndr, Nx, Ny) where Nx and Ny are
	      fp.shape/osf. The average is returned in each NDR.
  """
 
  # Here time is first axis for memory access efficiency. Will be transposed before return
  time_line = np.zeros( (ndr_sequence.size, fp.shape[0]/osf, fp.shape[1]/osf ) ).astype(np.float32)
  index =   frame_osf*ndr_cumulative_sequence 
  nclk  = frame_osf*ndr_sequence
  print "Frames per NDR", nclk[0:4], len(nclk)


  osf = np.int(osf)
  code2 = r"""
      int32_t ndr, idx, x, y, i, j,i_jit, j_jit;
      float *time_line_[Ntime_line[0]][Ntime_line[1]];
      float *fp_[Nfp[0]];
      float *dtmp;

      // Allocate memory space
      for (ndr = 0; ndr < Ntime_line[0]; ndr++) {
	dtmp = time_line + ndr*Ntime_line[2]*Ntime_line[1];
	for (y = 0; y < Ntime_line[1]; y++, dtmp += Ntime_line[2]) {
	  time_line_[ndr][y] = dtmp; 
	}
      }

      for (j = 0, dtmp = fp; j < Nfp[0]; j++, dtmp+=Nfp[1]) fp_[j] = dtmp; 


      for (ndr = 0; ndr < Nindex[0]; ndr++) {
	for (y = 0, j = y_offs; y < Ntime_line[1]; y++, j += osf) {
	  for (x = 0, i = x_offs; x < Ntime_line[2]; x++, i+=osf) {    
	    for (idx = index[ndr]-(nclk[ndr]-1); idx <= index[ndr]; idx++) {
	      
	      j_jit = y_jit[idx]+j;
	      i_jit = x_jit[idx]+i;
	      
	      if(i_jit < 0) i_jit += Nfp[1];
	      else if(i_jit >= Nfp[1]) i_jit -= Nfp[1];
	      if(j_jit < 0) j_jit += Nfp[0];
	      else if(j_jit >= Nfp[0]) j_jit -= Nfp[0];
	      
	      time_line_[ndr][y][x] += fp_[j_jit][i_jit];
	    }
	    time_line_[ndr][y][x] /= nclk[ndr];
	  }
	}
      }
  """
  
   
  scipy.weave.inline(code2, ['time_line', 'fp', \
    'nclk', 'index', \
    'x_jit', 'y_jit', 'osf', 'x_offs', 'y_offs'], verbose=1)
  
  ## Disabled since ExoSim  Issue#42
  #######time_line -= fp[x_offs::osf, y_offs::osf]
  
  return time_line.transpose( (1,2,0) )

def create_jitter_noise(channel, x_jit, y_jit, frame_osf, frame_time, key, opt):
  
  outputPointingTL = create_output_pointing_timeline(x_jit, y_jit, frame_osf, 
                                ndrCumSeq = channel.ndr_cumulative_sequence )
  
  jitter_x = channel.osf*(x_jit/channel.opt.plate_scale()).simplified
  jitter_y = channel.osf*(y_jit/channel.opt.plate_scale()).simplified

  fp_units = channel.fp.units
  fp  	   = channel.fp.magnitude
  osf      = np.int32(channel.osf)
  offs     = np.int32(channel.offs)

  magnification_factor = np.ceil( max(3.0/jitter_x.std(), 3.0/jitter_y.std()) )
  
  if (magnification_factor > 1):
    try:
      mag = np.int(magnification_factor.item()) | 1
    except:
      mag = np.int(magnification_factor) | 1
      
      
    print mag  
    fp = exolib.oversample(fp, mag)
    print mag  
 
    
    #### See ExoSim Issue 42, for following. 
#    fp = np.where(fp >= 0.0, fp, 1e-10)
    
    osf *= mag
    offs = mag*offs + mag//2
    jitter_x *= mag
    jitter_y *= mag
  

  if opt.noise.EnableSpatialJitter() != 'True': jitter_y *= 0.0
  if opt.noise.EnableSpectralJitter() != 'True': jitter_x *= 0.0

  jitter_x = np.round(jitter_x)
  jitter_y = np.round(jitter_y)
  noise = np.zeros((fp.shape[0]/osf, fp.shape[1]/osf, 0)).astype(np.float32)
  
  print channel.tl_shape[2], channel.tl_shape[2]/5.
  indxRanges = np.arange(0,6)*channel.tl_shape[2]/5.
#  indxRanges = np.arange(0,10)*channel.tl_shape[2]/9

  for i in range(len(indxRanges)-1):
    startIdx = int(indxRanges[i])
    print startIdx
    endIdx = int(indxRanges[i+1])
    print endIdx
    print i
    noise = np.append(noise , c_create_jitter_noise(fp.astype(np.float32), 
  				osf.astype(np.int32),
  				channel.ndr_sequence[startIdx:endIdx].astype(np.int32),
  				channel.ndr_cumulative_sequence[startIdx:endIdx].astype(np.int32), 
  				frame_osf.astype(np.int32),
  				jitter_x.magnitude.astype(np.int32), 
  				jitter_y.magnitude.astype(np.int32), 
  				x_offs = offs.astype(np.int32), 
  				y_offs = offs.astype(np.int32)).astype(np.float32), 
  				axis=2)
#  noise = c_create_jitter_noise(fp.astype(np.float32), 
#  				osf.astype(np.int32),
#  				channel.ndr_sequence.astype(np.int32),
#  				channel.ndr_cumulative_sequence.astype(np.int32), 
#  				frame_osf.astype(np.int32),
#  				jitter_x.magnitude.astype(np.int32), 
#  				jitter_y.magnitude.astype(np.int32), 
#  				x_offs = offs.astype(np.int32), 
#  				y_offs = offs.astype(np.int32))

  ## Multiply units to noise in 2 steps, to avoid 
  ##        Quantities memmory inefficiency 
  qq = channel.ndr_sequence* fp_units*frame_time
  noise = noise*qq
  print noise.shape, noise[...,0].max(), noise[...,1].max()
  print noise[...,1].max()/channel.fp[1::3,1::3].max()
  
  
  
  return  noise, outputPointingTL
  

def create_output_pointing_timeline(jitter_x, jitter_y, frame_osf, ndrCumSeq ):
  """
  Returns  an array  containing the high time resolution pointing offset and the corresponding NDR number.  
  Inputs:
    jitter_x  : array containing spectral pointing jitter time-line
    jitter_x  : array containing spatial pointing jitter time-line
    frameOsf  : Oversampling factor of the 2 pointing arrays
    ndrCumSeq : array associating index number of the above arrays  (before 
                oversampling)where a new NDR takes place.
  Returns:
     A 2d array of 3 columns, where 1st column is NDR number, 
     2nd Column is Oversampled Spectral Pointing jitter
     3rd column is Oversampled Spatial Pointing Jitter
     
  """
  
  pointingArray = np.zeros((ndrCumSeq[-1]*frame_osf, 3))
  

  indxPrev = 0
  for i, indx in enumerate(ndrCumSeq):
    pointingArray[indxPrev*frame_osf:indx*frame_osf, 0] = i
    pointingArray[indxPrev*frame_osf:indx*frame_osf, 1] = jitter_x[indxPrev*frame_osf:indx*frame_osf]
    pointingArray[indxPrev*frame_osf:indx*frame_osf, 2] = jitter_y[indxPrev*frame_osf:indx*frame_osf]
    indxPrev = indx
  
  return  pointingArray
  
  

def channel_noise_estimator(channel, key, yaw_jitter, pitch_jitter, frame_osf, frame_time, opt):

  # Jitter Noise
  
  fp = channel.fp[1::3,1::3]
  jitterless = np.ones((fp.shape[0], fp.shape[1], len(channel.ndr_sequence)))
  jitterless =  np.rollaxis(jitterless,2,0)
  jitterless = jitterless*fp
  jitterless =  np.rollaxis(jitterless,0,3)
#  plt.figure('fp')
#  plt.imshow(jitterless[...,10].magnitude)
  jitterless = jitterless*channel.ndr_sequence*frame_time

  
  if opt.jitter_switch == 0: 
     noise = jitterless
     outputPointingTL= np.zeros((channel.ndr_cumulative_sequence[-1], 3))
     print "using jitterless array"
  else:
     print "using jittered arry"
     jitNoise, outputPointingTL = create_jitter_noise(channel, yaw_jitter, pitch_jitter, frame_osf, frame_time, key, opt)
     noise = jitNoise

 
  #### activate for stellar bart work-----
  if opt.apply_LC ==0:
      print "OMITTING LIGHT CURVE..."
  else:
      print "APPLYING LIGHT CURVE..."
      noise *= channel.lc  

  if opt.zodi_switch ==1:
        noise += channel.zodi.sed[channel.offs::channel.osf].reshape(-1,1) *\
    frame_time * channel.ndr_sequence * channel.tl_units
  
  if opt.emm_switch ==1:   
      noise += channel.emission.sed[channel.offs::channel.osf].reshape(-1,1) *\
        frame_time * channel.ndr_sequence * channel.tl_units
  
  qe = np.load('%s/data/qe_rms.npy'%(opt.__path__))
  qe_uncert = np.load('%s/data/qe_uncert.npy'%(opt.__path__))
  qe = qe[0:noise.shape[0],0:noise.shape[1]]
  qe_uncert = qe_uncert[0:noise.shape[0],0:noise.shape[1]]  
  channel.qe = qe
  
  if opt.use_random_QE == 1:
      qe = np.random.normal(1,0.05, qe.shape) # for random uncertainty
      qe_uncert = np.random.normal(1,0.005, qe_uncert.shape) # for random uncertainty
      
      print "RANDOM QE GRID SELECTED..." 
  else:
      print "STANDARD QE GRID SELECTED..."
  channel.qe = qe 

 
  print "QE AND RESIDUAL QE STD", qe.std(), qe_uncert.std()

  if opt.apply_QE == 1:
      print "APPLYING QE GRID..."
      noise =  np.rollaxis(noise,2,0)
      noise = noise*qe*qe_uncert 
 
      noise =  np.rollaxis(noise,0,3)
  else:
      print "QE GRID NOT APPLIED..."
      
      
  noise  += channel.opt.detector_pixel.Idc() * frame_time * channel.ndr_sequence * channel.tl_units

  noise = np.where(noise >= 0.0, noise, 1e-30) 
  
  
  # Shot Noise
  if opt.noise.EnableShotNoise() == 'True':
    indxRanges = np.arange(0,10)*channel.tl_shape[2]/9
    for i in range(len(indxRanges)-1):
      startIdx = indxRanges[i]
      endIdx = indxRanges[i+1]
      noise[:,:,startIdx:endIdx] = np.random.poisson(noise[:,:,startIdx:endIdx])
#      noise[:,:,startIdx:endIdx] = np.random.normal(noise[:,:,startIdx:endIdx], np.sqrt(noise[:,:,startIdx:endIdx]))

 
  # Create ramps
  tl_shape = channel.tl_shape
  ####noise = noise.reshape(tl_shape[0], tl_shape[1], -1, opt.timeline.multiaccum()).cumsum(axis=3).reshape(tl_shape[0], tl_shape[1], -1)
  #### The above refactored to the bellow, for memmory efficiency. 
  for i in range(1, 1+ tl_shape[2]/opt.timeline.multiaccum()):
    noise[:,:, int(opt.timeline.multiaccum())*i-1] = np.sum(noise[:,:, int(opt.timeline.multiaccum())*(i-1):int(opt.timeline.multiaccum())*i], 2)


#  
  if opt.noise.EnableReadoutNoise() == "True":
    indxRanges = np.arange(0,20)*channel.tl_shape[2]/19
    for i in range(len(indxRanges)-1):
      startIdx = indxRanges[i]
      endIdx = indxRanges[i+1]
      useShape = (channel.tl_shape[0], channel.tl_shape[1],  endIdx-startIdx)
      noise[:,:,startIdx:endIdx] += np.random.normal(scale=channel.opt.detector_pixel.sigma_ro(), size=useShape)*channel.tl_units


  return key, noise,  outputPointingTL

def run(opt, channel, frame_time, total_observing_time, exposure_time):
  exosim_msg('Create noise timelines ... ')
  st = time.time()
  yaw_jitter = pitch_jitter = frame_osf = None
  
  jitter_file = opt.aocs.PointingModel().replace("__path__", opt.__path__) 
  if hasattr(opt.aocs, 'pointing_rms'): 
    jit_rms = opt.aocs.pointing_rms()
  else:
    jit_rms = None
    
  if opt.jitter_switch ==1:
    
      yaw_jitter, pitch_jitter, frame_osf = exolib.pointing_jitter(jitter_file,
                         total_observing_time, frame_time, opt, rms = jit_rms)
      print "RMS jitter", np.std(yaw_jitter), np.std(pitch_jitter)
      
      
          

      if opt.aocs.pointing_scan_throw()>0:       
           pitch_jitter = exolib.pointing_add_scan(pitch_jitter,
                                scan_throw_arcsec=opt.aocs.pointing_scan_throw(), 
                                frame_time = frame_time, frame_osf = frame_osf, 
                                exposure_time = exposure_time)
  
  
  for key in channel.keys():
    
    key, noise,  outputPointingTl = channel_noise_estimator(channel[key], key, yaw_jitter, 
                                                            pitch_jitter, frame_osf, frame_time, opt)
    channel[key].noise = noise
    channel[key].outputPointingTl = outputPointingTl
    
    if opt.apply_flat == 1:
        print "APPLYING FLAT FIELD AND BACKGROUND SUBTRACTION IN EXOSIM..."
        qe = channel[key].qe
        for kk in range (channel[key].noise.shape[2]):         
                 a = channel[key].noise[...,kk]    
                 a = a - channel[key].opt.detector_pixel.Idc() * frame_time * channel[key].ndr_sequence[kk] 
                              
                 a = a/qe
                 
                 background_1 = a[5:10,10:-10]
                 background_2 = a[-10:-5,10:-10]   
                 background = np.mean(np.vstack ( (background_1, background_2) ) )     
                 
                 if opt.diff == 0:
                     a = a - background                                
                 
                 channel[key].noise[...,kk] = a
    else:
        print "FLAT FIELD NOT APPLIED IN EXOSIM..."
  
  
  exosim_msg(' - execution time: {:.0f} msec.\n'.format(
  (time.time()-st)*1000.0))
  
  return noise, outputPointingTl


def run2(opt, channel, noise1, noise2, TL1, TL2):
  int_time = opt.NDR_time
  
  print "int time", int_time
  for key in channel.keys():
      
        bn = int(90./int_time)
        if bn < 1: 
            bn =1
            
        print bn
        ss =  int(noise1.shape[2]*2 / (bn*4)) *bn*4
        

        ct =0   
        noise = np.zeros((noise1.shape[0], noise1.shape[1], ss))
        for i in range (0, ss, bn*4):
            for j in range (2*bn):
                noise[...,i+j] = noise1[...,ct]
                noise[...,i+2*bn+j] = noise2[...,ct]
                ct = ct+1
       
        s3 = []
        for i in range (1,noise.shape[2],2) :
            s3.append(noise[...,i].sum())
        s3 = np.array(s3)
        
        plt.figure('s3')
        plt.plot(s3, 'rx-') 
        print noise.shape
         
        
  channel[key].noise = noise
  channel[key].outputPointingTl = TL1
 
  return noise
  
  
if __name__ == "__main__":
  
  exolib.exosim_error("This module not made to run stand alone")
