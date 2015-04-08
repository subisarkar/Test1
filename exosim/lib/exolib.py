import numpy as np
from scipy import signal
from scipy import interpolate
import sys, os, pyfits

fd = sys.stderr.fileno()

def exosim_error(error_msg):
    sys.stderr.write("Error code: {:s}\n".format(error_msg))
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


def rebin_old(x, xp, fp):
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
  idx = np.where(np.logical_and(xp > 0.9*x.min(), xp < 1.1*x.max()))[0]
  xp = xp[idx]
  fp = fp[idx]
  new_x = np.arange(xp.min(), xp.max(), np.diff(xp).min())
  new_f = np.interp(new_x, xp, fp, left=0.0, right=0.0)
  
  c = np.cumsum(new_f, dtype=np.float64)[1:] *np.diff(new_x)
  x = np.hstack( (np.array(x[0] - (x[1]-x[0])), x) )
  fr = np.interp(x, new_x[1:], c)
  
  fr = np.diff(fr) /np.diff(x)
    
  return x[1:], fr
  
def fast_convolution(im, delta_im, ker, delta_ker):
  """ fast_convolution.
    Convolve an image with a kernel. Image and kernel can be sampled on different
      grids defined.
    
    Parameters
    __________
      im : 			array like
				the image to be convolved
      delta_im :		scalar
				image sampling interval
      ker : 			array like
				the convolution kernel
      delta_ker :		scalar
				Kernel sampling interval
    Returns
    -------
      spectrum:			array like
				the image convolved with the kernel.
  """
  fc_debug = False
  # Fourier transform the kernel
  kerf = (np.fft.rfft2(ker))
  ker_k = [ np.fft.fftfreq(ker.shape[0], d=delta_ker),
	   np.fft.rfftfreq(ker.shape[1], d=delta_ker) ]
  ker_k[0] = np.fft.fftshift(ker_k[0])
  kerf     = np.fft.fftshift(kerf, axes=0)
  
  # Fourier transform the image
  imf  = np.fft.rfft2(im)
  im_k = [ np.fft.fftfreq(im.shape[0], d=delta_im),
	   np.fft.rfftfreq(im.shape[1], d=delta_im) ]
  im_k[0] = np.fft.fftshift(im_k[0])
  imf     = np.fft.fftshift(imf, axes=0)
  
  # Interpolate kernel 
  kerf_r = interpolate.RectBivariateSpline(ker_k[0], ker_k[1],
					   kerf.real)
  kerf_i = interpolate.RectBivariateSpline(ker_k[0], ker_k[1],
					   kerf.imag)
  if (fc_debug):
    pl.plot(ker_k[0], kerf[:, 0].real,'.r')
    pl.plot(ker_k[0], kerf[:, 0].imag,'.g')
    pl.plot(im_k[0], kerf_r(im_k[0], im_k[1])[:, 0],'-r')
    pl.plot(im_k[0], np.abs(imf[:, 0]),'-b')

  # Convolve
  imf = imf * (kerf_r(im_k[0], im_k[1]) + 1j*kerf_i(im_k[0], im_k[1])) 
  
  if (fc_debug):
    pl.plot(im_k[0], np.abs(imf[:, 0]),'-y')

  imf = np.fft.ifftshift(imf, axes=0)
  
  return np.fft.irfft2(imf)*(delta_ker/delta_im)**2

   
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
  
def Psf_Interp(zfile, delta_pix, WavRange):
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
    hdulist = pyfits.open(zfile)    
    NAXIS1, NAXIS2 = hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2']
    in_ph_size_x, in_ph_size_y = hdulist[0].header['CDELT1']*NAXIS1, hdulist[0].header['CDELT2']*NAXIS2
    num_pix_x, num_pix_y = np.trunc(in_ph_size_x/delta_pix).astype(np.int), np.trunc(in_ph_size_y/delta_pix).astype(np.int)
               
    inwl   = np.zeros(len(hdulist))
    redata = np.zeros((num_pix_y, num_pix_x, len(hdulist)))

    xin = np.linspace(-1.0, 1.0, NAXIS1)
    yin = np.linspace(-1.0, 1.0, NAXIS2)

    xout = np.linspace(-1.0, 1.0, num_pix_x)
    yout = np.linspace(-1.0, 1.0, num_pix_y)

    for i, hdu in enumerate(hdulist):
        inwl[i]   = np.float64(hdu.header['WAV'])        
        f = interpolate.RectBivariateSpline(xin, yin, hdu.data)
	redata[..., i] = f(xout,yout)

        redata[..., i] /= redata[..., i].sum()
    return interpolate.interp1d(inwl, redata, axis=2, bounds_error=False, fill_value=0.0, kind='quadratic')(WavRange)


def Psf(wl, fnum, delta, nzero = 4, shape='airy'):
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
  Nx = np.round(scipy.special.jn_zeros(1, nzero)[-1]/(2.0*np.pi) * fnum*wl.max()/delta).astype(np.int)
  
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
    img   = (scipy.special.j1(arg)/arg)**2
  elif shape=='gauss':
    arg = np.multiply.outer(yy**2 + xx**2, d)
    img = np.exp(-arg)
    
  norm = img.sum(axis=0).sum(axis=0)
  return img/norm

  
  
  
def PixelResponseFunction(psf_shape, psf_delta, osf, delta, lx = 1.7, ipd = 0.0):
  '''
  Estimate the detector pixel response function with the prescription of 
  Barron et al., PASP, 119, 466-475 (2007).
  
  Parameters
  ----------
  psf_shape	: touple of scalars 
		  (ny, nx) defining the PSF size	
  psf_delta     : scalar
		  plate scale of PSF in microns/sample
  osf		: scalar
		  number of samples in each resolving element. The 
		  final shape of the response function would be shape*osf
  delta 	: scalar
		  Phisical size of the detector pixel in microns
  lx		: scalar
		  diffusion length in microns
  ipd           : scalar
		  distance between two adjacent detector pixels 
		  in microns
		 
  Returns
  -------
  kernel	: 2D array
		  the kernel image
  kernel_delta  : scalar
                  the kernel sampling interval in microns
  '''
  kernel = np.zeros( (psf_shape[0]*osf, psf_shape[1]*osf) )
  kernel_delta = delta/osf
  yc, xc = np.array(kernel.shape)//2
  yy = (np.arange(kernel.shape[0]) - yc) * kernel_delta
  xx = (np.arange(kernel.shape[1]) - xc) * kernel_delta
  mask_xx = np.where(np.abs(xx) > 0.5*(delta-ipd))
  mask_yy = np.where(np.abs(yy) > 0.5*(delta-ipd))
  xx, yy = np.meshgrid(xx, yy)
  
  kernel = np.arctan(np.tanh( 0.5*( 0.5*delta - xx)/lx )) - \
	   np.arctan(np.tanh( 0.5*(-0.5*delta - xx)/lx ))
	 
	 
  kernel*= np.arctan(np.tanh( 0.5*( 0.5*delta - yy)/lx )) - \
  	   np.arctan(np.tanh( 0.5*(-0.5*delta - yy)/lx )) 
  
  kernel[mask_yy, ...] = 0.0
  kernel[..., mask_xx] = 0.0

  # Normalise the kernel such that the pixel has QE=1
  kernel *= osf**2/kernel.sum()
  kernel = np.roll(kernel, -xc, axis=1)
  kernel = np.roll(kernel, -yc, axis=0)
  
  return kernel, kernel_delta
	  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  import numpy as n
import scipy.interpolate
import scipy.ndimage








def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None







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



def animate(Data):
    
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

#    Data = data['channel']['SWIR'].timeline

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    wframe = None
    
    for j in range(0,Data.shape[2]):
    
            oldcol = wframe

            X = np.arange(0, Data.shape[1])
            Y = np.arange(0, Data.shape[0])
            X, Y = np.meshgrid(X, Y)
            
           
        
            Z = Data[...,j]
            print Z.sum()
            
            wframe = ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5)
            
#            ax.set_zlim(0,20000)
#            ax.set_title(j)
#            
     
        

    
        # Remove old line collection before drawing
            if oldcol is not None:
                ax.collections.remove(oldcol)
    
            plt.pause(.01)