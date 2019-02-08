#!/usr/bin/env python

"""
 Plotting the power spectra
"""

import os,sys
import pylab
import numpy as np
import optparse
import COSMO_constants as CNST1
import aipy as a

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1420405751.77 # FREQUENCY OF 21 CM HYDROGEN LINE

if '__name__==__main__':
  o = optparse.OptionParser()
  o.set_usage('plotPowerSpectra.py [options] *.npz')
  o.set_description(__doc__)
  a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
  o.add_option('--f0',dest='f0', type=float, default=150e6, help='Center Frequency of the bandwidth in Hz')
  o.add_option('--buff',dest='buff', type=float, default=50, help='Buffer in ns from the horizon limit')
  o.add_option('--plot', dest='plot', action='store_true', help='plot option enables the plotting of the spectra and a pop-up window of the spectrum will appear')
  o.add_option('--save', dest='save',action='store_true',help='Saves the 2D power spectra')
  opts, args = o.parse_args(sys.argv[1:])

  z = CNST1.fq2z(opts.f0*1e-9) # redshift
  LAMBDA = C/F21 # wavelength in meters
  buffr = opts.buff/3.3356
  aa = a.cal.get_aa(opts.cal, np.array([opts.f0*1e-9]))
  for npzfile in args: 
     print'Reading ',npzfile
     try:
        npz = np.load(npzfile)# loading the npz file
        outfile = 'wedgef%.3gMHz-'%(opts.f0*1e-6) + npzfile
     except IOError:
        print npzfile, ' does not exist'

     if os.path.exists(outfile):
        print outfile, 'exists.  Skipping...'
        continue

     kpl = npz['kpl']
     keys = npz.keys() # dict keys
     keys.remove('kpl') # removing kpl from dict keys

     Tspec , weights = {}, {}
     bls,mags,conj = [],[],[]
     for bl in keys:
        if bl[0] == 'w': continue
        wbl = 'w'+str(bl)
        i,j = a.miriad.bl2ij(bl)
        crd = aa.get_baseline(i,j)*(opts.f0*1e-9)
        mag = np.sqrt(crd[0]**2 + crd[1]**2) # calculating magnitude
        if crd[0] < 0.:
           conj.append(1)
        else:
           conj.append(0)
        mags.append(mag)
        bls.append(bl)
        Tspec[mag] = Tspec.get(mag,0) + npz[bl] # getting the power into an array
        weights[mag] = weights.get(mag,0) + npz[wbl] # getting the weights into an array
       
     mags,bls = np.unique(mags), np.unique(bls) # selecting unique baselines
     inds = np.argsort(mags)
     mags, bls, conj = np.take(mags,inds), np.take(bls,inds), np.take(conj,inds)
     kprs = np.zeros((len(mags)),)
     half = len(kpl)/2. 
     wfall = np.zeros((len(mags),half),dtype=complex)
     hors = np.zeros((len(mags)),)
     phors = np.zeros((len(mags)),)
    
     for ii, mg in enumerate(mags):
        print ii
        spec = Tspec[mg]/weights[mg] # averaging the power by divided the sum by the number of timestamps
        if conj[ii] ==1:
           spec = spec[::-1]
        inds = np.argsort(kpl)
        spec = np.take(spec,inds)
        #mag_m = mg/3.3356  # magnitude in metres
        kpr = CNST1.k_perp(z) * mg
        kprs[ii] = kpr
        pkpr = CNST1.k_perp(z)*(mg + opts.buff*opts.f0*1e-9)
        #pkpr = CNST1.k_perp(mag_m + buffr,opts.f0)
        hors[ii] = CNST1.horizon_limit(z) * kpr  # calculating the horizon limit
        phors[ii] = CNST1.horizon_limit(z) * pkpr # calculating the super-horizon limit
        #wfall[ii,0] = spec[0]
        wfall[ii,:] = (spec[0:half] + spec[half:][::-1])/2. # folding over k_parallel
     kpl.sort() 
     if opts.plot:
       print 'Plotting ...'
       fig = pylab.figure()
       axis = fig.add_subplot(111)
       if len(kprs)==1:
          kpr_min=0; kpr_max=2*kprs
       else:
          kpr_min=kprs[0]; kpr_max=kprs[-1]
       im=pylab.imshow(np.log10(np.abs(wfall.T)),aspect='auto',interpolation='nearest',cmap='jet',extent=(kpr_max,kpr_min,0,np.abs(kpl[0])))
       cb=pylab.colorbar(im)#,ticks=[0,2,4,6,8,10,12,14])
       #pylab.plot(kprs,hors,'.',lw=1,color='white')
       #pylab.plot(kprs,phors,'.',lw=1,color='orange')
       #pylab.title(r'${\rm log}_{10}[P(k)]$',fontsize = 15)
       #pylab.ylabel(r'$k_{\parallel}\ [h{\rm Mpc}^{-1}]$',fontsize = 12, labelpad=5)
       #pylab.xlabel(r'$k_{\perp}\ [h{\rm Mpc}^{-1}]$',fontsize = 12,labelpad=5)
       #pylab.tick_params(labelsize= 12)
       #pylab.xlim(0,0.12)
       #pylab.ylim(0,0.42)
       pylab.show()

     if opts.save:
       print 'Saving to ',outfile
       if os.path.exists(outfile):
          print outfile, 'exists.  Skipping...'
          continue
       np.savez(open(outfile, 'wb'),kprs=kprs, kpl=kpl, mags=mags*(opts.f0*1e-9) , hors=hors, phors=phors, wedge=wfall.T)    
