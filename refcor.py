from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.stats import mad_std
import pdb

slpFiles = glob.glob('/data1/Local/AZLabStability/*I???.fits')

outDir = '/data1/Local/AZLabStability/refcor'

ampStarts = [0,512,1024,1536]
ampEnds = [512,1024,1536,2048]
namp = len(ampStarts)
nfile = len(slpFiles)

doplot=True

for fileInd, oneFile in enumerate(slpFiles):
    HDUList = fits.open(oneFile)
    
    newimg = np.zeros_like(HDUList[0].data,dtype=np.float)
    HDUList[0].header['REFPCORA']= ('T','Are reference pixels applied?')
    ngroup = HDUList[0].header['NGROUP']
    nint = HDUList[0].header['TOT_NINT']
    
    if fileInd == 0:
        correctionArray = np.zeros((ngroup,namp,nfile))
    
    for oneGroup in range(ngroup):
        img = np.array(HDUList[0].data[oneGroup,:,:],dtype=np.float)
        for ind in range(namp):
            ## Go through each amplifier and subtract average reference
            oneAmpRef = img[0:4,ampStarts[ind]:ampEnds[ind]]
            madstd = mad_std(oneAmpRef)
            goodp = np.abs(oneAmpRef - np.median(oneAmpRef)) < 3. * madstd
            avgRef = np.mean(oneAmpRef[goodp])
            newimg[oneGroup,:,ampStarts[ind]:ampEnds[ind]] = img[:,ampStarts[ind]:ampEnds[ind]] - avgRef
            #HDUList[0].header['REFPCOR'+str(ind)] = (avgRef,'Amplifier '+str(ind)+' correction flux')
            correctionArray[oneGroup,ind,fileInd] = avgRef

            if doplot == True and ind == 3:
                plt.plot(img[2,1084:1115],label='Raw')
                plt.plot(newimg[oneGroup,2,1084:1115],label='Ref Sub')
                plt.legend()
                plt.show()
                pdb.set_trace()

                    
    HDUList[0].data = newimg
    
    HDUList.writeto(os.path.join(outDir,'REFC_'+os.path.basename(oneFile)),overwrite=True)
        
    if np.mod(fileInd,10) == 0:
        print('Finished file '+str(fileInd)+' of '+str(len(slpFiles)))
    
corHDU = fits.PrimaryHDU(correctionArray)
corHDU.name = 'REFCOR'
corHDUList = fits.HDUList([corHDU])
corHDUList.writeto(os.path.join(outDir,'refcor_values.fits'),overwrite=True)

