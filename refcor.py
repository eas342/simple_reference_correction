from astropy.io import fits
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.stats import mad_std

slpFiles = glob.glob('/data1/Local/AZLabStability/*.slp.fits')

outDir = '/data1/Local/AZLabStability/refcor'

ampStarts = [0,512,1024,1536]
ampEnds = [512,1024,1536,2048]
namp = len(ampStarts)


for fileInd, oneFile in enumerate(slpFiles):
    HDUList = fits.open(oneFile)
    
    newimg = np.zeros_like(HDUList[0].data)
    HDUList[0].header['REFPCORA']= ('T','Are reference pixels applied?')
    img = HDUList[0].data[0,:,:]
    for ind in range(namp):
        ## Go through each amplifier and subtract average reference
        oneAmpRef = img[0:4,ampStarts[ind]:ampEnds[ind]]
        madstd = mad_std(oneAmpRef)
        goodp = np.abs(oneAmpRef - np.median(oneAmpRef)) < 3. * madstd
        avgRef = np.mean(oneAmpRef[goodp])
        newimg[0,:,ampStarts[ind]:ampEnds[ind]] = img[:,ampStarts[ind]:ampEnds[ind]] - avgRef
        HDUList[0].header['REFPCOR'+str(ind)] = (avgRef,'Amplifier '+str(ind)+' correction flux')
        
    HDUList[0].data = newimg
    
    HDUList.writeto(os.path.join(outDir,'REFC_'+os.path.basename(oneFile)),overwrite=True)
    
    if np.mod(fileInd,10) == 0:
        print('Finished file '+str(fileInd)+' of '+str(len(slpFiles)))
