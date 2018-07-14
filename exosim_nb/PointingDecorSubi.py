import pyfits
import numpy as np
import scipy.interpolate
from sys import stdout
from scipy.interpolate import interp1d, interp2d

class ExoSimFits:
    exposureIndxDict = None
    hdulist = None
    modelBeam = None
    planetName = None
    pointingNdr = None
    pointingJittX = None
    pointingJittY = None
    modelJitterOffsetSpec = None
    modelJitterOffsetSpat = None
    modelJitterVarianceSpec  = None
    modelJitterVarianceSpat  = None
    wlArr = None
    crArr = None
    zArr = None
    pixArr = None
    nExp = None
    mAccum = None

    def __init__(self, fileName, qe, apply_flat, diff):
        self.hdulist = pyfits.open(fileName)
        self.planetName =  self.hdulist[0].header['PLANET']
        self.exposureIndxDict = self.getExposureNumDict()
        self.pixArr = np.arange(0, self.hdulist[1].data.shape[0])
        self.nExp =  self.hdulist[0].header['NEXP']
        self.mAccum = self.hdulist[0].header['MACCUM']
        
        qe = qe[0:self.hdulist[1].data.shape[0],0:self.hdulist[1].data.shape[1]]
#        print "QE STD", qe.std(), qe.shape
        
        print "Applying CDS... in decorr file"
        
        if apply_flat == 1:
            print "APPLYING FLAT FIELD AND BKG SUB IN DECORR FILE..."
            for i in range (1,1+ self.nExp*self.mAccum):          
                 a = self.hdulist[i].data 
    
                 background_1 = a[5:10,10:-10]
                 background_2 = a[-10:-5,10:-10]    
                 background = np.mean(np.vstack ( (background_1, background_2) ) )
                 if diff == 0:
                     a = a - background              
                 a = a/qe
                 self.hdulist[i].data = a

##             
        self.modelBeam = self.hdulist[2].data
        self.calcModelPointingOffsets()
        
    def getExposureNumDict(self):
		expNumDict = {}
		for i, hdu in enumerate(self.hdulist):
			if 'EXTNAME' in hdu.header.keys():
				if hdu.header['EXTNAME']=='NOISE':
					expNum = hdu.header['EXPNUM']
					ndrNum = hdu.header['ENDRNUM']
					if expNum not in expNumDict.keys():
						expNumDict[expNum] = {}
					expNumDict[expNum][ndrNum] = i
				if hdu.header['EXTNAME']=='SIM_POINTING':
					self.pointingNdr = hdu.data[:,0]
					self.pointingJittX = hdu.data[:,1]
					self.pointingJittY = hdu.data[:,2]
				if hdu.header['EXTNAME']=='INPUTS':
					tempInputs = hdu.data
					self.wlArr =  np.array(map(lambda x:x[0], hdu.data))
					self.crArr =  np.array(map(lambda x:x[1], hdu.data))
				if hdu.header['EXTNAME']=='TIME':
					self.timeArr =  np.array(map(lambda x:x[0], hdu.data))
					## HACK, Append to extra time
					self.timeArr = np.append(self.timeArr,(self.timeArr[-1] + self.timeArr[2] - self.timeArr[1]))
					self.timeArr = np.append(self.timeArr,(self.timeArr[-1] + self.timeArr[3] - self.timeArr[2]))
					self.zArr =  np.array(map(lambda x:x[1], hdu.data))
					
		return expNumDict
  
    def calcModelPointingOffsets(self):
		totalNdrs = self.mAccum* self.nExp#;print totalNdrs

		self.modelJitterOffsetSpec = np.zeros(totalNdrs)
		self.modelJitterOffsetSpat = np.zeros(totalNdrs)
		self.modelJitterVarianceSpec  = np.zeros(totalNdrs)
		self.modelJitterVarianceSpat  = np.zeros(totalNdrs)
		
#		for i in range(totalNdrs): 
#			self.modelJitterOffsetSpec[i] = np.mean(self.pointingJittX[self.pointingNdr==i]) - np.mean(self.pointingJittX[self.pointingNdr==1]) 
#			self.modelJitterOffsetSpat[i] = np.mean(self.pointingJittY[self.pointingNdr==i]) - np.mean(self.pointingJittY[self.pointingNdr==1]) 
#			self.modelJitterVarianceSpec[i] = np.std(self.pointingJittX[self.pointingNdr==i])
#			self.modelJitterVarianceSpat[i] = np.std(self.pointingJittY[self.pointingNdr==i])


		for i in range(totalNdrs): 
			self.modelJitterOffsetSpec[i] = np.median(self.pointingJittX[self.pointingNdr==i]) - np.median(self.pointingJittX[self.pointingNdr==1]) 
			self.modelJitterOffsetSpat[i] = np.median(self.pointingJittY[self.pointingNdr==i]) - np.median(self.pointingJittY[self.pointingNdr==1]) 
			self.modelJitterVarianceSpec[i] = np.std(self.pointingJittX[self.pointingNdr==i])
			self.modelJitterVarianceSpat[i] = np.std(self.pointingJittY[self.pointingNdr==i])
#	
    def getPoindingOffsets(self, pscale, p_spec, p_spat, mAccum):
        if mAccum >= self.mAccum:
            raise ValueError("Maximim mAccum value is %d"%(self.mAccum-1))
        XX = {'spec':self.modelJitterOffsetSpec[mAccum::self.mAccum]/pscale, 'spat':self.modelJitterOffsetSpat[mAccum::self.mAccum]/pscale}
        if p_spec == False:
            XX = {'spec':np.zeros((len(self.modelJitterOffsetSpat[mAccum::self.mAccum]))), 'spat':self.modelJitterOffsetSpat[mAccum::self.mAccum]/pscale}
        if p_spat == False:
            XX = {'spec':self.modelJitterOffsetSpec[mAccum::self.mAccum]/pscale,'spat':np.zeros((len(self.modelJitterOffsetSpat[mAccum::self.mAccum])))}

        return XX
	                
                
    def getPoindingVariance(self, mAccum):
			if mAccum >= self.mAccum:
				raise ValueError("Maximim mAccum value is %d"%(self.mAccum-1))
			return {'spec':self.modelJitterVarianceSpec[mAccum::self.mAccum], 'spat':self.modelJitterVarianceSpat[mAccum::self.mAccum]}
	
	
    def getExposure(self, expNum, ndrNum):
		indx = self.exposureIndxDict[expNum][ndrNum]
		intTime = self.timeArr[indx] - self.timeArr[indx - 1]
		return self.hdulist[indx].data, intTime
	
    def getBinEdges(self, R):
		binSize = [self.wlArr[1] - self.wlArr[0]]
		binEdges = [self.wlArr[0] , self.wlArr[1]]
		for i in range(1, len(self.wlArr)):
			binSize.append((1+1./R)*binSize[i-1])
			binEdges.append(binEdges[-1]+binSize[-1])
			if binEdges[-1]>np.max(self.wlArr):
				break
		binEdges = np.array(binEdges[:-1])
		return binEdges


class FitTools:
	gauss1d_mu = None
	gauss1d_sigma = None
	modelBeamWhite = None
	modelBeamWhiteXval = None
	f2 = None

	
	def __init__(self):
		pass
	
	def func_gauss1d(self, x, *p):
		A, mu, sigma, offset = p
		return offset + A*np.exp(-(x-mu)**2/(2.*sigma**2))
	
	def func_gauss1d_pixOff_sigma(self, x, *p):
		A, offset = p
		return offset + A*np.exp(-(x-self.gauss1d_mu)**2/(2.*self.gauss1d_sigma **2))
	
	def func_modelBeamWhite(self, x, *p):
		A, yOffset, xOffset = p
		x*=xOffset
		#return np.interp(x, , , 'cubic') + yOffset
		return A* self.f2(x) + yOffset
	
	def fitGauss1d(self, xVals, yVals):
		p0 = [np.max(yVals) - np.min(yVals), np.argmax(yVals), 3, np.min(yVals)]
		coeff, var_matrix = curve_fit(self.func_gauss1d, xVals, yVals, p0=p0)
		return coeff
	def fitGauss1d_pixOff_sigma(self, xVals, yVals, mu, sigma):
		self.gauss1d_mu = mu
		self.gauss1d_sigma = sigma
		p0 = [np.max(yVals) - np.min(yVals), np.min(yVals)]
		coeff, var_matrix = curve_fit(self.func_gauss1d_pixOff_sigma, xVals, yVals, p0=p0)
		return coeff
	
	def fitModelBeamWhite(self,modelBeamWhite,  xVals, yVals ):
		self.modelBeamWhite = modelBeamWhite
		self.modelBeamWhiteXval = np.arange(len(modelBeamWhite))
		self.f2 = interp1d(self.modelBeamWhiteXval, self.modelBeamWhite, kind='cubic')
		p0 = [np.max(yVals) - np.min(yVals), np.median(yVals) , 1]
		coeff, var_matrix = curve_fit(self.func_modelBeamWhite, xVals[5:-5], yVals[5:-5], p0=p0)
		return coeff



class JitterRemoval:
	@staticmethod
	def cubicIterp(esf, jiggOffsetMeasure):
		upsampleFactor =1#;print jiggOffsetMeasure
		shiftedMaps =[]
		for i in range(esf.nExp):
#			stdout.write("\r%d/%d" % (i,esf.nExp))
#			stdout.flush()
			map2, _ = esf.getExposure(i,1)#;print "hello",i,jiggOffsetMeasure['spec'][i],jiggOffsetMeasure['spec'][i]
   			map2_0, _ = esf.getExposure(i,0)
			map2 = map2-map2_0
			f = interp2d(np.arange(map2.shape[1]), np.arange(map2.shape[0]), map2, kind='cubic')
			tempMap2Cubic = f(np.linspace(0, map2.shape[1]-1, map2.shape[1]*upsampleFactor) -jiggOffsetMeasure['spec'][i], np.linspace(0, map2.shape[0]-1, map2.shape[0]*upsampleFactor) -jiggOffsetMeasure['spat'][i])
			tempMap2CubicLR = binData(tempMap2Cubic, upsampleFactor)
			shiftedMaps.append(tempMap2CubicLR)
#		print '\n'
		return shiftedMaps
	@staticmethod
	def subgridShift(esf, jiggOffsetMeasure, upsampleFactor):
		shiftedMaps =[]
		bufferPx = 2
		for i in range(esf.nExp):
			stdout.write("\r%d/%d" % (i,esf.nExp))
			stdout.flush()
			map2, _ = esf.getExposure(i,1)
			tempMap2 = map2.repeat(upsampleFactor, axis=0).repeat(upsampleFactor, axis=1)
			offsetSpec = int(np.round(jiggOffsetMeasure['spec'][i]*upsampleFactor))
			offsetSpatial = int(np.round(jiggOffsetMeasure['spat'][i]*upsampleFactor))
			tempMap2Shift = np.zeros_like(tempMap2)
			
			tempMap2Shift[bufferPx*upsampleFactor-offsetSpec:-bufferPx*upsampleFactor-offsetSpec, 
								bufferPx*upsampleFactor-offsetSpatial:-bufferPx*upsampleFactor-offsetSpatial] = tempMap2[bufferPx*upsampleFactor:-bufferPx*upsampleFactor, bufferPx*upsampleFactor:-bufferPx*upsampleFactor]
			shiftedMaps.append(tempMap2Shift)
#		print '\n'
		return shiftedMaps
        @staticmethod
        def fftShift(esf,jiggOffsetMeasure):
		shiftedMaps =[]
		for i in range(esf.nExp):
			map2, _ = esf.getExposure(i,1)#;print "hello",i
   			map2_0, _ = esf.getExposure(i,0)
			im = map2-map2_0
			dy = jiggOffsetMeasure['spec'][i]
			dx = jiggOffsetMeasure['spat'][i]
			imFft = np.fft.fftshift(np.fft.fft2(im))
			xF,yF = np.meshgrid(np.arange(im.shape[0]) - im.shape[0]/2,np.arange(im.shape[1]) - im.shape[1]/2)
			imFft=imFft*np.exp(-1j*2*np.pi*(xF.T*dx/im.shape[0]+yF.T*dy/im.shape[1]))
			shiftedIm = np.fft.ifft2(np.fft.ifftshift(imFft))
			shiftedMaps.append(np.abs(shiftedIm))  
		return shiftedMaps           

 
def crossCorr1d(refData, testData, xCorWin):
	xCor = np.zeros(2*xCorWin) 
	#plt.plot(refData[xCorWin:-xCorWin])
	for j in range(2*xCorWin):
		xCor[j] = np.std((refData[xCorWin:-xCorWin] - testData[j:-2*xCorWin+j])**2)
		####plt.plot(refData[xCorWin:-xCorWin] - testData[j:-2*xCorWin+j], color)
		#plt.plot(testData[j:-2*xCorWin+j], 'r')
	xInterp = np.linspace(0,len(xCor)-1, 3000)
	#plt.show()
	f2 =  scipy.interpolate.interp1d(np.arange(len(xCor)), xCor , kind='cubic')
#	if 0:
#		plt.plot(xCor)
#		plt.plot(xInterp, f2(xInterp))
#		plt.show()
	offsetPx = xInterp[np.argmin(f2(xInterp))]
	return offsetPx

def getRelativeOffsets(imageRef, imgeTest):
	# Colapse Spatial Axis
	tempRef = np.sum(imageRef,0)
	tempTest = np.sum(imgeTest, 0)
	tempRef/=np.sum(tempRef)
	tempTest/=np.sum(tempTest)
	xCorWin = 26
 	xCorWin = 10; 

	offsetSpec = -crossCorr1d(tempRef, tempTest, xCorWin)+xCorWin
 
	
	tempRef = np.sum(imageRef,1)
	tempTest = np.sum(imgeTest, 1)
	tempRef/=np.sum(tempRef)
	tempTest/=np.sum(tempTest)
	xCorWin = 26
 	xCorWin = 10; 

	offsetSpatial = -crossCorr1d(tempRef, tempTest, xCorWin)+xCorWin
 

	
	return {'spec':offsetSpec, 'spat':offsetSpatial}

def binData(hiresMap, binSize):
	sumMap = np.zeros((hiresMap.shape[0]/binSize, hiresMap.shape[1]/binSize))
	for i in range(sumMap.shape[0]):
		for j in range(sumMap.shape[1]):
			sumMap[i,j] = np.sum(hiresMap[i*binSize:(i+1)*binSize, j*binSize:(j+1)*binSize])
	return sumMap

def clear_folder():
  import os
  import glob

  files = glob.glob('/Users/c1341133/Desktop/TLdata/*')
  for f in files:
      os.remove(f)

    


class ExosimPointingDecorr():
    
    def __init__(self, fileName, qe, ch, method, apply_flat, diff):
                
#        p_spec = ps_info[0]
#        p_spat = ps_info[1]
#        pscale = ps_info[2]
        
        self.fileName = fileName
        ft = FitTools()
        self.esf = ExoSimFits(fileName, qe, apply_flat, diff)
        
        jiggOffsetMeasure = {'spec':np.zeros(self.esf.nExp), 'spat':np.zeros(self.esf.nExp)} 
        for i in range(self.esf.nExp):
            im, _ = self.esf.getExposure(i, 1)
            offset = getRelativeOffsets(self.esf.modelBeam, im)
            jiggOffsetMeasure['spec'][i] = offset['spec']		
            jiggOffsetMeasure['spat'][i] = offset['spat']

        if method =='xcorr':
            print "method xcorr", method
            self.shiftedMaps = JitterRemoval.cubicIterp(self.esf, jiggOffsetMeasure) 
        elif method == 'fft':
            print "method fft", method
            
            self.shiftedMaps = JitterRemoval.fftShift(self.esf, jiggOffsetMeasure) # does not need pointing info or pscale

#            self.shiftedMaps = JitterRemoval.fftShift(self.esf, self.esf.getPoindingOffsets(pscale, p_spec, p_spat,mAccum=1))
        else: #method == 'pointing'
            self.shiftedMaps = JitterRemoval.cubicIterp(self.esf, self.esf.getPoindingOffsets(pscale, p_spec, p_spat,mAccum=1))
#            self.shiftedMaps = JitterRemoval.fftShift(self.esf, self.esf.getPoindingOffsets(pscale, p_spec, p_spat,mAccum=1))
            
            
            print "method pointing", method
   

#        if ch == 'SWIR':     
#            self.shiftedMaps = JitterRemoval.cubicIterp(self.esf, jiggOffsetMeasure) 
##            self.shiftedMaps = JitterRemoval.cubicIterp(self.esf, self.esf.getPoindingOffsets(mAccum=1))
#        elif ch == 'MWIR':
#            aa= np.load('/Users/c1341133/Desktop/TLdata/specTL.npy') *  5.583333e-5/11.166666e-5
#            bb = np.load('/Users/c1341133/Desktop/TLdata/spatTL.npy')*  5.583333e-5/11.166666e-5
#            cc = {'spec':aa, 'spat':bb}
#            self.shiftedMaps = JitterRemoval.cubicIterp(self.esf, cc)
        
#        import matplotlib.pyplot as plt
#        plt.figure('%s'%(fileName))
#        plt.plot(jiggOffsetMeasure['spat'], 'b-')
#        plt.plot(self.esf.getPoindingOffsets(pscale, p_spec, p_spat, mAccum=1)['spat'], 'r-')
#        print "STD" ,np.std(jiggOffsetMeasure['spat']), np.std(self.esf.getPoindingOffsets(pscale, p_spec,p_spat, mAccum=1)['spat'])
#        print np.corrcoef(jiggOffsetMeasure['spat'],self.esf.getPoindingOffsets(pscale, p_spec, p_spat, mAccum=1)['spat'])
#        for i in range (len(jiggOffsetMeasure['spat'])):
#            print jiggOffsetMeasure['spat'][i], self.esf.getPoindingOffsets(mAccum=1)['spat'][i]
#        if ch == 'SWIR':       
#            np.save('/Users/c1341133/Desktop/TLdata/specTL.npy',np.array(jiggOffsetMeasure['spec']))
#            np.save('/Users/c1341133/Desktop/TLdata/spatTL.npy',np.array(jiggOffsetMeasure['spat']))
#        elif ch == 'MWIR':
#            clear_folder()

		
    def getHdu(self):
        for i in range(self.esf.nExp):
			accumNum = 1
			hduIndx =  self.esf.exposureIndxDict[i][accumNum]
			self.esf.hdulist[hduIndx].data = self.shiftedMaps[i]
        newFileName = "%s_decorr.fits"%self.fileName.split('/')[-1].split('.')[0]
        newFullFileName =  str.join('/', self.fileName.split('/')[:-1])+'/'+newFileName
        self.esf.hdulist.writeto(newFullFileName, clobber=True)
        hdulist0 = pyfits.open(newFullFileName)
        return hdulist0

if __name__ == "__main__":
	fileName = '/home/alp/ExoSimOutput/sim_0000/SWIR_signal.fits'
	epd = ExosimPointingDecorr(fileName)
	hdulist = epd.getHdu()











