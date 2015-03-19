import numpy as np
from scipy import signal,interpolate,special,linalg
import sys, os, pyfits

fd = sys.stderr.fileno()

def exosim_error(error_msg):
    sys.stderr.write("Non valid option or physical units given - error code {:s}\n".format(error_msg))
    sys.exit(0)
    
def exosim_msg(msg):
    os.write(fd, msg)
    #os.fsync(fd)
      
    
def logbin(x, a,  R, xmin=None, xmax=None):
  n = a.size
  imin = 0
  imax = n-1
  
  if xmin == None or xmin < x.min(): xmin = x.min()
  if xmax == None or xmax > x.max(): xmax = x.max()
  
  idx = np.argsort(x)
  xp = x[idx]
  yp = a[idx]
  
  delta_x = xmax/R
  N = 20.0 * (xmax-xmin)/delta_x
  _x = np.linspace(xmin,xmax, N)
  _y = np.interp(_x, xp, yp)
  
  nbins = 1+np.round( np.log(xmax/xmin)/np.log(1.0 + 1.0/R) ).astype(np.int)
  bins = xmin*np.power( (1.0+1.0/R), np.arange(nbins))
  
  slices  = np.searchsorted(_x, bins)
  counts = np.ediff1d(slices)
  
  mean = np.add.reduceat(_y, slices[:-1])/(counts)
  bins = 0.5*(bins[:-1] + bins[1:])
  return bins[:-1], mean[:-1]


def rebin(x, xp, fp):
  ''' Resample a function fp(xp) over the new grid x, rebinning if necessary, 
    otherwise interpolates
    Parameters
    ----------
    x	: 	array like
	New coordinates
    fp 	:	array like
	y-coordinates to be resampled
    xp 	:	array like
	x-coordinates at which fp are sampled
	
    Returns
    -------
    out	: 	array like
	new samples
  
  '''
  if (x.min() < xp.min() or x.max() > xp.max()): exosim_error('array can not be rebinned')
  bin_width = xp[1]-xp[0]
  idx = np.argsort(xp)
  xp = xp[idx]
  fp = fp[idx]
  
  if x.size >= xp.size :
    return x, np.interp(x, xp, fp)
 
  elif x.size < xp.size :
    xbin = np.append(x-bin_width/2.0, x[-1]+bin_width/2.0)
    
    
    slices = np.searchsorted(xp, xbin)
    slices[slices == fp.size] = fp.size - 1
    s = np.add.reduceat(fp, slices)/np.add.reduceat(np.ones(fp.size), slices)
    
    #test_integral_of_input= np.trapz(fp, x = xp) 
    #test_integral_of_output= np.trapz(s[:-1], x = x) 
    #print test_integral_of_input, test_integral_of_output
    
    return x, s[:-1]
    
   
def planck(wl, T):
  """ Planck function. 
    
    Parameters
    __________
      wl : 			array
				wavelength [micron]
      T : 			scalar
				Temperature [K]
				Spot temperature [K]
    Returns
    -------
      spectrum:			array
				The Planck spectrum  [W m^-2 sr^-2 micron^-1]
  """
    
  a = np.float64(1.191042768e8)
  b = np.float64(14387.7516)
  try:
    x = b/(wl*T)
    bb = a/wl**5 / (np.exp(x) - 1.0)
  except ArithmeticError:
    bb = np.zeros(np.size(wl))
  return bb
 
 
def sed_propagation(sed, transmission, emissivity=None, temperature = None):
  sed.sed = sed.sed*transmission.sed
  if emissivity and temperature:
    sed.sed = sed.sed + emissivity.sed*planck(sed.wl, temperature)

  return sed


  
def psf_stack(zfile, delta_pix, wav_range):
    '''
        PSF Interpolation
        Parametes
        ---------
        zfile : string
        input PSF fits file
        Delta : scalar
        Sampling interval in micron
        WavRange : ndarray
        array of wavelengths in micron
        
        Returns
        -------
        PSF interpolated data cube. Area normalised to unity.
        
        '''
    zfile = "/Users/c1341133/Desktop/zgrid7.fits" # file containing core psfs
    
    hdulist = pyfits.open(zfile)
    NAXIS1, NAXIS2 = hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2']
    
    pix_size0 = hdulist[0].header['CDELT1'] # size of the original pixel in micron
    ad_ovs = pix_size0/delta_pix
    
    inwl   = np.zeros(len(hdulist))
    redata = np.zeros((NAXIS1*ad_ovs, NAXIS2*ad_ovs, len(hdulist)))
    
    xin = np.linspace(-1.0, 1.0, NAXIS1)
    yin = np.linspace(-1.0, 1.0, NAXIS2)
    
    xout = np.linspace(-1.0, 1.0, NAXIS1*ad_ovs)
    yout = np.linspace(-1.0, 1.0, NAXIS2*ad_ovs)
    
    for i, hdu in enumerate(hdulist):
        inwl[i]   = np.float64(hdu.header['WAV'])
        f = interpolate.RectBivariateSpline(xin, yin, hdu.data)
        redata[..., i] = f(xout,yout)
        redata[..., i] /= redata[..., i].sum()
    
    outdata = interpolate.interp1d(inwl, redata, axis=2, bounds_error=False, fill_value=0.0, kind='quadratic')(wav_range)
    psf_stack = outdata

    return psf_stack


def psf(wl, fnum, delta, nzero = 4, shape='airy'):
    '''
        Calculates an Airy Point Spread Function arranged as a data-cube. The spatial axies are
        0 and 1. The wavelength axis is 2. Each PSF area is normalised to unity.
        
        Parameters
        ----------
        wl	: ndarray [physical dimension of length]
        array of wavelengths at which to calculate the PSF
        fnum : scalar
        Instrument f/number
        delta : scalar
        the increment to use [physical units of length]
        nzero : scalar
        number of Airy zeros. The PSF kernel will be this big. Calculated at wl.max()
        shape : string
        Set to 'airy' for a Airy function,to 'gauss' for a Gaussian
        
        Returns
        ------
        Psf : ndarray
        three dimensional array. Each PSF normalised to unity
        '''
    print "generating psf stack"    
    
    Nx = np.round(special.jn_zeros(1, nzero)[-1]/(2.0*np.pi) * fnum*wl.max()/delta).astype(np.int)
  
    Ny = Nx
  
    if shape=='airy':
      d = 1.0/(fnum*wl)
    elif shape=='gauss':
      sigma = 1.029*fnum*wl/np.sqrt(8.0*np.log(2.0))
      d     = 0.5/sigma**2
    
    x = np.linspace(-Nx*delta, Nx*delta, 2*Nx+1)
    y = np.linspace(-Ny*delta, Ny*delta, 2*Ny+1)
  
    yy, xx = np.meshgrid(y, x)
  
    if shape=='airy':
      arg = 1.0e-20+np.pi*np.multiply.outer(np.sqrt(yy**2 + xx**2), d)
      img   = (special.j1(arg)/arg)**2
    elif shape=='gauss':
      arg = np.multiply.outer(yy**2 + xx**2, d)
      img = np.exp(-arg)
    
    norm = img.sum(axis=0).sum(axis=0)
    
    print "psf stack generated"

    return img/norm
                                                                                


def fpa(fp,psf_stack,delta_psf,sed,kernal,ad_osf,pix_size,QE): # to be placed in exolib

    print "generating oversampled fpa"
    
    fpa = fp*1 #nullify array copying effect
    
    psf_position = np.arange(0,fpa.shape[1]-psf_stack.shape[1],delta_psf) # psf spacings along x axis
    
    j0 =  psf_position #array of psf 'tiles' left most positions along x axis
    j1 = j0 + psf_stack.shape[1] #right position of psf 'tiles'
    idx = range(0,len(psf_position)) #index = no of psfs being coadded
    i0 = (fpa.shape[0]/2 - psf_stack.shape[0]/2) #upper edge of psf 'tiles'
    i1 = i0 + psf_stack.shape[0] #lower edge of psf 'tiles'
        
    for k in idx: # coadds psfs onto fp array and multiplies by sed
        fpa[i0:i1, j0[k]:j1[k]] += sed[k]*psf_stack[...,k]  #sed needs to be oversampled at rate ovs
    
    print "applying QE grid with variations"
    
    fpa = fpa*QE
                   
    print "convolving fpa with pixel kernal"

    fpa_crop = fpa[i0-100:i1+100]
        
    [U,S,V] = linalg.svd(kernal)
    s = S[0]
    v = U[:,0]
    v = v.reshape(len(v),1)*np.sqrt(s)
    h = V[0]*np.sqrt(s)
    h= np.array([h])
    aa = signal.convolve2d(fpa_crop,v,'same')
    cc1 = signal.convolve2d(aa,h,'same')
    
    print "convolution done"
    
    corr = convolution_normalization(kernal) 
    
    cc1=cc1*corr
    
    print "kernal normalization correction applied:", corr
    
    print "applying additional oversampling factor to convolved fpa"
    
    xin = np.linspace(-1.0, 1.0, cc1.shape[0])
    yin = np.linspace(-1.0, 1.0, cc1.shape[1])
    xout = np.linspace(-1.0, 1.0, cc1.shape[0]*ad_osf)
    yout = np.linspace(-1.0, 1.0, cc1.shape[1]*ad_osf)
    
    f = interpolate.RectBivariateSpline(xin, yin, cc1)
    cc2 = f(xout,yout)
    
            
    conv_fpa = np.zeros((fpa.shape[0]*ad_osf,fpa.shape[1]*ad_osf))
        
    conv_fpa = np.zeros((fpa.shape[0]*ad_osf,fpa.shape[1]*ad_osf)) 
    conv_fpa[(i0-100)*ad_osf:(i0-100)*ad_osf+cc2.shape[0],(j0[0])*ad_osf:(j0[0])*ad_osf+cc2.shape[1]] = cc2
    
    osf = kernal.shape[0]*ad_osf
    pix_count = conv_fpa[int(osf/2)::osf, int(osf/2)::osf]*1
                
    conv_fpa = conv_fpa*(fpa.sum()/pix_count.sum())
    pix_count = conv_fpa[int(osf/2)::osf, int(osf/2)::osf]*1
    
    
    print "addition oversampling done"
    


    
    return fpa,conv_fpa

def convolution_normalization(kernal):
    osf = kernal.shape[0]
    corr = osf**2/np.sum(kernal)
    return corr   

def kernal(pix_size,ovs):  #need to add edge
    ld = 1.7
    
    ax = np.linspace(-pix_size/2,pix_size/2,ovs)
    ay = np.linspace(-pix_size/2,pix_size/2,ovs)
    ax,ay = np.meshgrid(ax,ay)
    resp = (np.arctan(np.tanh(0.5*( 0.5*pix_size-ax)/ld))-np.arctan(np.tanh(0.5*(-0.5*pix_size-ax)/ld)))*(np.arctan(np.tanh(0.5*( 0.5*pix_size-ay)/ld))-np.arctan(np.tanh(0.5*(-0.5*pix_size-ay)/ld)))
    kernal = resp
    
    return kernal


def jitter(obs_time,int_time,osr,rms,mode=2):
    
    """
        Jitter
        
        Simulates 2 d jitter (pointing variation) as a timeline of positional offsets
        Uses Herschel jitter timeline as basis.
        
        Inputs:
        
        1)  obs_time : total observation time in seconds
        
        2)  int_time : time for one integration (NDR)
        
        3)  osr : oversampling rate of the integration time (e.g 100 = time step of int_time/100)
        
        4)  rms : rms of the desired jitter in degrees
        
        5)  mode = 1 : one PSD used to obtain jitter in 2 dimensions
        mode = 2 : two PSDs used - one for each orthogonal dimeension
        
        Output:
        
        1) RA jitter time series in radians (xt) in degrees
        
        2) Dec jitter time series in radians (yt) in degrees
        
        Requirements:
        
        1) jitter_file : file with known jitter timelines (e.g. Herschel data)
        
        """
    
    int_N = np.int(obs_time/int_time) # total number of full NDRs in the total obs time
    new_N = int_N*int_time*osr
    x = int(1+np.log(new_N)/np.log(2))
    new_N = 2**x
    new_dt = int_time/osr
    new_time = np.arange(0,new_N*new_dt,new_dt)
    new_fs=1/new_dt
    new_df = new_fs/new_N
    
    jitter_file = "/Users/c1341133/Desktop/herschel_long_pointing.fits"
    f = pyfits.open(jitter_file)  # open a FITS file
    tbdata = f[1].data  # assume the first extension is a table
    time = tbdata['Time']
    ra_t = tbdata['RA']
    dec_t = tbdata['Dec']
    
    if len(time)%2 != 0:  # timeline needs to be even number for real fft
        time = time[0:-1]
        ra_t = ra_t[0:-1]
        dec_t = dec_t[0:-1]
    
    ra_t = (ra_t-np.mean(ra_t))*(np.cos(dec_t*np.pi/180))
    dec_t = dec_t-np.mean(dec_t)

    N = np.float(len(time))
    dt = time[1]-time[0]
    fs = 1.0/dt
    df = fs/N
    
    freq = np.fft.rfftfreq(np.int(N),d=dt)
    new_freq = np.fft.rfftfreq(np.int(new_N),d=new_dt)
    new_df = new_freq[1]-new_freq[0]
    
    ra_f = np.fft.rfft(ra_t)/N
    dec_f = np.fft.rfft(dec_t)/N
    
    ra_psd= 2*abs(ra_f)**2/df
    dec_psd = 2*abs(dec_f)**2/df
    
    ra_psd[0]=1e-30
    ra_psd[-1]=ra_psd[-1]/2
    dec_psd[0]=1e-30
    dec_psd[-1]=dec_psd[-1]/2
    
    # smooth the psd
    
    window_size = 10
    window = np.ones(int(window_size))/float(window_size)
    ra_psd = np.convolve(ra_psd, window, 'same')
    dec_psd = np.convolve(dec_psd, window, 'same')
    
    # resample and 'zero pad' to new frequency grid and N
    
    f1 = interpolate.interp1d(freq, ra_psd,bounds_error=False,fill_value=1e-30, kind='linear')
    f2 = interpolate.interp1d(freq, dec_psd,bounds_error=False,fill_value=1e-30, kind='linear')
    
    new_ra_psd = f1(new_freq)
    new_dec_psd = f2(new_freq)
    
    psd = [new_ra_psd,new_dec_psd]
    N = new_N
    fs = new_fs
    df = fs/N
    
    if mode == 1:
        #new to work on how to scale this
        
        comb_t = np.sqrt(ra_t**2 + dec_t**2)
        comb_f = np.fft.rfft(comb_t)/N
        comb_psd = 2*abs(comb_f)**2/df
        comb_psd[0]=1e-30
        comb_psd[-1]=comb_psd[-1]/2
        
        phi_t = np.arctan2(dec_t,ra_t)
        phi_f = np.fft.rfft(phi_t)/N
        phi_psd = 2*abs(phi_f)**2/df
        phi_psd[0]=1e-30
        phi_psd[-1]=phi_psd[-1]/2
        
        psd1 = psd[0]
        ps1 = psd1*df/2
        amp = np.random.normal(0,np.sqrt(ps1),len(ps1))
        phi = np.random.uniform(0,np.pi,len(ps1))
        zf = amp*np.exp(phi*1j)
        zt = np.fft.irfft(zf)*N
        
        psd2 = psd[1]
        ps2 = psd2*df/2
        amp = np.random.normal(0,np.sqrt(ps2),len(ps2))
        phi = np.random.uniform(0,np.pi,len(ps2))
        anglef = amp*np.exp(phi*1j)
        anglet = np.fft.irfft(anglef)*N
        
        xt = zt*np.cos(anglet)
        yt = zt*np.sin(anglet)


    elif mode == 2:
        
        psd1 = psd[0]
        ps = psd1*df/2
        amp = np.random.normal(0,np.sqrt(ps),len(ps))
        phi = np.random.uniform(0,np.pi,len(ps))
        xf = amp*np.exp(phi*1j)
        xt = np.fft.irfft(xf)*N
        
        psd2 = psd[1]
        ps = psd2*df/2
        amp = np.random.normal(0,np.sqrt(ps),len(ps))
        phi = np.random.uniform(0,np.pi,len(ps))
        yf = amp*np.exp(phi*1j)
        yt = np.fft.irfft(yf)*N


    else:
        print "error: maximum of 2 psds can be used"

#    xt = xt   
#    yt = yt 
    xt = xt*(rms*(np.sqrt(2)/2)/np.std(xt))    
    yt = yt*(rms*(np.sqrt(2)/2)/np.std(yt))
    return xt,yt,new_time
    
def animate(data):
    import matplotlib.pyplot as plt  

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    a = data['timeline']
    wframe = None
    
    def anim(X,Y,Z):
        return ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5, color = 'blue')
  
    for i in range(0,a.shape[2]):
        
#            oldcol = wframe
            ax.clear()
            
            
            Z = a[...,i]
            
            
            X = np.arange(0, Z.shape[1])
            Y = np.arange(0, Z.shape[0])
            X, Y = np.meshgrid(X, Y)
            
            anim(X,Y,Z)
            
            ax.set_zlim(0,a[...,0].max())
#            ax.set_title(j)
#            
    
            
    
        
            # Remove old line collection before drawing
#            if oldcol is not None:
#                ax.collections.remove(oldcol)
        
            plt.pause(.01)
            
def animate2(data):
    from matplotlib import animation    
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d.axes3d import Axes3D  #-- voor 3d assen
    import time

    a = data['timeline']
    Z = a[...,0]

    X = np.arange(0, Z.shape[1])
    Y = np.arange(0, Z.shape[0])
    X, Y = np.meshgrid(X, Y)

    def init_wave_profile(): 
        Z = a[...,0]


  # ----------------------------------------------
    def update_wave_profile(i,Z,p1):
  # ----------------------------------------------
        Z = a[...,i]



  #-- define plot and initialize profile (flat)
    fig = plt.figure(figsize=(6,6))
    init_wave_profile()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    p1 = ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10, color = 'blue')
#  ax.set_xlim3d(-0.1,1.1);
#  ax.set_ylim3d(-0.1, 1.1);
#  ax.set_zlim3d(-0.5, 1.5);

  #-- animate the profile  
    i= np.arange(0,100) 
    anim = animation.FuncAnimation(
            fig         = fig,
            func       = update_wave_profile,
            init_func = init_wave_profile,
            frames   = 10,
            interval  = 200,
            blit        = False,
            fargs     = (i, Z, p1)
            )
    anim()
    plt.show()

    
