#! /usr/bin/env python

"""
  Generating power spectra using delay-filtering approach
"""

import numpy as np
import optparse
import scipy.constants as CNST1
import COSMO_constants as CNST2
import os,sys
import aipy as a

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1420405751.77 # FREQUENCY OF 21 CM HYDROGEN LINE
curtime, zen = None, None

#===============================================================
def genWindow(size):
    """
    Implements Blackmann-Harris filter

    size : Size/Lenth of frequency channel to which filter is applied; type:int
    """
    
    window = np.zeros((size),)
    alpha = 0.16
    _w = lambda n: (1-alpha)/2. - 0.5*np.cos((2*np.pi*n)/(size-1)) + alpha/2.*np.cos((4*np.pi*n)/(size-1))
    window[:] = _w(np.arange(size))
    return window

#===============================================================
def frequency_channels(freqs,f0,B):
    """
    Evaluates frequency boundaries given the center frequency and bandwidth

    freqs: Frequencies in GHz; type:numpy.ndarray
    f0   : Center frequency in GHz; type:float
    B    : Bandwidth in GHz; type:float
    """
    ch0 = np.argmin(np.abs(freqs-f0))
    df = freqs[1]-freqs[0]
    dCH = int(B/df)
    ch1, ch2 =  ch0-dCH/2, ch0+dCH/2 
    return ch1, ch2

#===============================================================
def f2etas(freqs):
    """
    Evaluates geometric delay (fourier conjugate of frequency)
   
    -freqs: Frequencies in GHz; type:numpy.ndarray
    """
    df = freqs[1] - freqs[0]
    etas = np.fft.fftfreq(freqs.size,df)
    return etas

#=============================================================== 
def delay_transform(data,fqs,convert=None):
    """
    Fourier transforms visibility along frequency axis

    - data: per baseline visibility; type:numpy.ndarray
    - fqs:  slected frequencies in GHz; dtypw:numpy.ndarray
    """
    N = fqs.size
    df = fqs[1] - fqs[0]
    window = genWindow(N)
    delaySpec = np.fft.ifft(data) * N * df
    return delaySpec 

#===============================================================
def compute_omegaB(beampath,jd,freqs,freq_wgts,nside=None,healpix=None):
    """
    Evaluates 3D vloume normlaization from Thyagarajan 2016
    
    - beampath : Path where the beam are stored; type:str
    - jd       : Julian date at which the beam is generated
    - freqs    : Frequency at which the beam is generated
    - freq_wgts: Weights applied to the visibility
    - healpix  : Enable healpix normalization
    """
    A0 = np.load(beampath+'/jones-f%.4g_j%.5f.npz'%(freqs[0]*1e3,np.float(t)))['mueller']
    df = freqs[1] - freqs[0] # frequency resolution in Hz
    A = np.zeros((4,A0.shape[-1],len(freqs)),dtype=complex) # initializing A marix
    for ii, f in enumerate(freqs):
      for p in range(4):
         jones = np.load(beampath+'/jones-f%.4g_j%.5f.npz'%(f*1e3,np.float(jd)))['mueller']
         A[p,:,ii] = jones[p,p,:]
    if healpix:
       domega = (4*np.pi)/(12*nside*nside)
    else:
       domega = 1.0/A0.shape[-1] # dividing by the number of sources
    Aw = A*freq_wgts
    omegaB = {}
    omegaB = np.nansum(np.nansum(Aw**2,axis=1),axis=1)
    OmegaB['xx'] = omegaB[0] 
    OmegaB['xy'] = omegaB[1]
    OmegaB['yx'] = omegaB[2]
    OmegaB['yy'] = omegaB[3]
    return omegaB     

#===============================================================
if '__name__==__main__':
   o = optparse.OptionParser()
   o.set_usage('genPowerSpectra.py [options] *.UV')
   o.set_description(__doc__)
   a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
   o.add_option('--f0', dest='f0', type=float, default=150e6, help='Center Frequency in Hz')
   o.add_option('--bw', dest='bw', type=float, default=8e6, help='Bandwidth in Hz')
   o.add_option('--convert', dest='convert', action='store_true', help='if convert option is set, temperature conversion is carried out')
   o.add_option('--jd', dest='jd', default='julian_dates.npz', help='npz file where the julian dates are stored')
   o.add_option('--beampath',dest='beampath',default=os.getcwd(), help='path where the beams are stored. Default is %s'%os.getcwd())
   o.add_option('--healpix',dest='healpix',action='store_true',help='Enable healpix normailization') 
   o.add_option('--nside',dest='nside',default=128,help='Nside of the healpiz map, needed if healpix option is enabled.')  
   o.add_option('-n',dest='normalization',default=None,help='Normalization value for power spectrum')
   opts,args = o.parse_args(sys.argv[1:])
   
   uv = a.miriad.UV(args[0])
   sfreq = uv['sfreq'] # Starting Frequency in GHz, aipy saves the freuqnency in GHz
   nchan = uv['nchan'] # Number of frequency channels
   sdf = uv['sdf'] # Frewquency resolution in GHz
   stop = sfreq + (nchan-1)*sdf
   freqs = np.linspace(sfreq,stop,nchan) # frequency range in GHz
   aa = a.cal.get_aa(opts.cal, sdf, sfreq, nchan)   
   jy2mK = 1e-26*(C/opts.f0)**2/(2*CNST1.k)*1e3 # conversion from Jy to mK
   z = CNST2.fq2z(opts.f0*1e-9) # corresponding redshift
   # choosing the frequency channels:
   ch1, ch2 = frequency_channels(freqs*1e9,opts.f0,opts.bw)
   fqs = freqs[ch1:ch2] # frequency in  GHz
   etas = np.fft.fftshift(f2etas(fqs))
   df = (fqs[1] - fqs[0])*1e9 # frequency resolution in Hz
   N = fqs.size #frequency channels in the bandwidth
   fqs_wgts = genWindow(N) # generating windowing function
   kpl =  CNST2.k_parallel(etas,z) # k_parallel
   cos_scalar = CNST2.transverse_comoving_distance(z)**2 * CNST2.comoving_depth(opts.bw,z) # cosmological conversions   
   
   jds = np.load(opts.jd)['jd']
   # computing normalization factor
   OmegaB = {}
   if opts.healpix:
      nside=opts.nside;healpix=True
   else:
      nside=None;healpix=None
   for ii, t in enumerate(jds):
      omegaB = compute_omegaB(opts.beampath,t,fqs,fqs_wgts,nside=nside,healpix=healpix) 
      if not t in OmegaB: OmegaB[t] = {}
      OmegaB[t]['xx'] = OmegaB[t].get('I',0) + omegaB[0]
      OmegaB[t]['xy'] = OmegaB[t].get('Q',0) + omegaB[1]
      OmegaB[t]['yx'] = OmegaB[t].get('U',0) + omegaB[2]
      OmegaB[t]['yy'] = OmegaB[t].get('V',0) + omegaB[3]

   tdat = {}
   Tspec = {}
   for uvfile in args:
     print 'Reading ',uvfile
     outfile = uvfile+'.npz'
     if os.path.exists(outfile):
        print outfile, 'exists.  Skipping...'
        continue
     uv = a.miriad.UV(uvfile) # reading uvfile
     a.scripting.uv_selector(uv,opts.ant,opts.pol) # selecting polarization
     for (uvw,t,(i,j)), data,f in uv.all(raw=True):
       bl = str(a.miriad.ij2bl(i,j))
       wbl = 'w'+ bl
       if t!= curtime:
         old_zen = zen
         aa.set_jultime(t) # setting the julian data
         lst = aa.sidereal_time() # local sidereal time
         zen = a.phs.RadioFixedBody(lst, aa.lat) #sets zenith
         zen.compute(aa)
         curtime = t
         pol = a.miriad.pol2str[uv['pol']] # extracting the polarization
       aa.set_active_pol(opts.pol) # activating the polarization
       if old_zen != None:
         data *= np.conj(aa.gen_phs(zen,i,j)) * aa.gen_phs(old_zen,i,j)

       if opts.convert:
           # convert to temperature
           d = jy2mK(fqs) * data[ch1:ch2]
       else:
           d = data[ch1:ch2]
           # delay transform one smaple at a time
       _d = np.fft.fftshift(np.fft.ifft(d*fqs_wgts)) # fourier transforming the visibilitier per baseline for the specified bandwidth
       if opts.normalization==None:
          omega_bw = OmegaB.get(t) # get normalization factor
          norm = omega_bw.get(opts.pol)
       else:
          norm = opts.normalization
       scalar = cos_scalar/(opts.bw * norm)
       try:        
         Tspec[bl] = Tspec.get(bl,0) + scalar * tdat[bl]*_d.conj() # evaluating the power
         Tspec[wbl] = Tspec.get(wbl,0) + 1 # calculating the weight/timestamps
       except(KeyError): pass
       tdat[bl] = _d

     Tspec['kpl'] = kpl
     print 'Writing', outfile
     np.savez(outfile,**Tspec)

         
