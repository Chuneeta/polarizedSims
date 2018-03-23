"""
Generated unpolarized or polarized foregrounds given a beam model and source catalgue/healpix map
The catalgoue shouls be in the following format:
                  SrcID  RA  DEC STOKES_I RM POLARIZATION_FRACTION SPECTRAL_INDEX
"""

from joblib import Parallel, delayed
import multiprocessing
import numpy as np
from math import *
import healpy as hp
import genBeamModel as bm
import os,sys
import random

c = 3e8 #speed of light

#=====================================================================
def ga2equ(ga):
    """
    Convert Galactic to Equatorial coordinates (J2000.0)
    The coordinates should be given in degrees
    
    -ga : Galactic coordinated (l,b); type:tuple of float
    """

    l = radians(ga[0])
    b = radians(ga[1])
    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = radians(192.859508)
    pole_dec = radians(27.128336)
    posangle = radians(122.932-90.0)

    ra = atan2( (cos(b)*cos(l-posangle)), (sin(b)*cos(pole_dec) - cos(b)*sin(pole_dec)*sin(l-posangle)) ) + pole_ra
    dec = asin( cos(b)*cos(pole_dec)*sin(l-posangle) + sin(b)*sin(pole_dec) )
    return np.array([degrees(ra), degrees(dec)])
#=====================================================================
def genForegrounds(stokesI,stokesQ,stokesU,stokesV,beam,nu,jd,filename):
    """
    Evaluates the sky as seen the telescope (I_xx,I_xy,I_yx,I_yy)

    -stokes I : Total intensity; type:numpy.ndarray
    -stokes Q : Linear polarization; type:numpy.ndarray
    -stokes U : Linear polarization; type:numpy.ndarray
    -stokes V : Circular polarization; type:numpy.ndarray
    -beam     : Beam Model; type:numpy.ndarray 
    -nu     : Frequency in Hz; type:float
    -jd       : Julian Date; type:float
    filename  : Output filename; type:string
    """

    jones = np.zeros(shape=(2,2,stokesI.shape[0]),dtype=complex)
    for i in range(0,stokesI.shape[0]):
       jones[0,0,i] = beam[i,0]
       jones[0,1,i] = beam[i,1]
       jones[1,0,i] = beam[i,2]
       jones[1,1,i] = beam[i,3]

    b_s = np.zeros(shape=(4,stokesI.shape[0]),dtype=float)
    b_s[0,:] = stokesI
    b_s[1,:] = stokesQ
    b_s[2,:] = stokesU
    b_s[3,:] = stokesV

    stokes=.5*np.matrix([[1.,1.,0.,0.],[0.,0.,1.,1.j],[0.,0.,1.,-1.j],[1.,-1.,0.,0.]])# matrix that converts from telescope X-Y coordinates to proper Stokes parameters
    out_jones =  np.ndarray(shape=(4,4,stokesI.shape[0]),dtype=complex)
    bjones =  np.ndarray(shape=(4,4,stokesI.shape[0]),dtype=complex)
    sinverse =  np.ndarray(shape=(4,4,stokesI.shape[0]),dtype=complex)
    sjones =  np.ndarray(shape=(4,4,stokesI.shape[0]),dtype=complex)
    b_m = np.ndarray(shape=(4,stokesI.shape[0]),dtype=complex) # array the stores the corrected Stokes parameters for each image pixel
    for i in range(0,stokesI.shape[0]):
      out_jones[:,:,i] = np.kron(jones[:,:,i],jones[:,:,i].conj()) # JxJ
      bjones[:,:,i] = np.dot(out_jones[:,:,i],stokes) # JxJS
      sinverse[:,:,i] = np.dot(np.array(np.linalg.inv(stokes)),out_jones[:,:,i]) #S^-1JxJ
      sjones[:,:,i] = np.dot(sinverse[:,:,i],stokes) #S^-1JxJS (Mueller matrix)
      b_m[:,i] = np.dot(bjones[:,:,i],b_s[:,i]) # JxJSb_m

    np.savez(open('jones-f%.4g_j%.5f.npz'%(nu*1e-6,np.float(jd)),'wb'), mueller=sjones)
    np.savez(open('%s-f%.4g_j%.5f.npz'%(filename,nu*1e-6,np.float(jd)),'wb'),intstokes=b_s,stokes=b_m)

#=====================================================================
def genCatStokes(data,nu,nu_0,alpha,RM,pfrac,pol=False):
    """
    Evaluates Stokes Q and U for a source catalgue given Stokes I, polarization fraction and Rotation Measure
   
    -data    : Total intensity; type:numpy.ndarray
    -nu      : Observation frequencies in Hz; type:float 
    -nu0     : Reference frequency in Hz; type:float
    -alpha   : Spectral Index; type:numpy.ndarray
    -RM      : Rotation Measure; type:float
    -pfrac   : Poalrization fraction, type:numpy.ndarray
    -pol     : Enable polarization; type:boolean 
    """
    rad2deg = lambda val: val * 180./np.pi
    stokesI = data*((nu/float(nu_0))**alpha) # power law formula to generate Stokes I at different frequencies
    stokesQ = np.zeros((stokesI.shape),dtype=float)
    stokesU = np.zeros((stokesI.shape),dtype=float)
    stokesV = np.zeros((stokesI.shape),dtype=float)

    if pol:
         P_ems = pfrac * stokesI
         phi = RM *(c/nu)**2 # computing polarization angle
         stokesQ = P_ems * np.cos(2*phi)
         stokesU = P_ems * np.sin(2*phi)

               
    #setting stokes I to zero, hack for paper to estimate leakage
    #stokesI = np.zeros((stokesI.shape),dtype=float)

    return stokesI, stokesQ, stokesU, stokesV
#=====================================================================
def genHealStokes(data,nu,nu0,RM,pol):
   """
   Evaluates Stokes Q and U from polarized emission

   -data    : Polarized; type:numpy.ndarray
   -nu      : Observation frequencies in Hz; type:flaot
   -nu0     : Reference frequency in Hz; type:float
   -RM      : Rotation Measure; type:float
   -pol     : Enable polarized outputs
   """
   phi = RM *(c/nu)**2 # computing polarization angle
   stokesQ = data * np.cos(2*phi)
   stokesU = data * np.sin(2*phi)
   stokesV = np.zeros((stokesI.shape),dtype=complex)

   return stokesI, stokesQ, stokesU, stokesV
#=====================================================================
def call_genForegrounds(data,nu,nu_0,ra,dec,xpol,ypol,jd,RM,filename,polmap=None,pfrac=None,nside=None,pol=None):
  """
    Calling genForegrounds to evaluate the sky as seen the telescope (I_xx,I_xy,I_yx,I_yy)

    -data     : Stokes I/total intensity data; type: numpy.ndarray
    -nu       : Range of frequencies in Hz; type: float
    -nu_0     : Reference frequency in Hz; type:float
    -ra       : Right ascension of sources/healpix pixels; type:numpy.ndarray
    -dec      : Declination of sources/healpix pixels; type:numpy.ndarray
    -xpol     : Input beam model for xx polarization; type:string
    -ypol     : Input beam model for yy polarization; type:string 
    -jd       : Julian date; dtype:float
    -RM       : Rotation Measure; type:float
    -filename : Output filename; type:string
    -polmap   : Polarized emission; type:numpy.ndarray
    -pfrac    : Polarization Fraction; type: numpy.ndarray 
    -nside    : Nside/Resolution of Healpix map; type:float
    -pol      : Enable polarized outputs
    -healpix  : True if input data is a healpix map

  """
  print "Generating Stokes parameters"
  if polmap==None:
     sI,sQ,sU,sV = genCatStokes(data,nu,nu_0,alpha,RM,pfrac,pol=pol)
  else:
     sI,sQ,sU,sV = genHealStokes(polmap,nu,nu0,RM)
  
  print "Generating beam"
  beam = bm.InterBeam(xpol,ypol,nu,jd,ra,dec)
  print "Generating Maps"
  genForegrounds(sI,sQ,sU,sV,beam,nu,jd,filename)

#=====================================================================
def loadCat(catalogue):
   """
   Extracting data from a catalogue
   
   - catalogue : Input catalogue; type:string
   """
   cat = np.loadtxt(catalogue,dtype='str') # loading catalogue
   if cat.ndim==1:
      ra =  deg2rad(np.array([float(cat[1])])) #converting RA to radians
      dec = deg2rad(np.array([float(cat[2])])) # converting DEC to radians
      data = np.array([float(cat[3])])
      alpha = np.array([float(cat[4])])
      RM = np.array([float(cat[5])])
      pfrac = np.array([float(cat[6])])
   else:
      ra =  deg2rad(np.array(map(float,cat[:,1]))) #converting RA to radians
      dec = deg2rad(np.array(map(float,cat[:,2]))) # converting DEC to radians
      data = np.array(map(float,cat[:,3]))
      alpha = np.array(map(float,cat[:,4]))
      RM = np.array(map(float,cat[:,5]))
      pfrac = np.array(map(float,cat[:,6]))
   
   return data,ra,dec,alpha,RM,pfrac
#=====================================================================
def convHealCoord(nside):
   """
   Convert healpix coordinates to ra and dec
   
   - nside : Nside of healpix map
   """
   coord = hp.pix2ang(nside,np.arange(12*nside*nside))
   theta = 90-(rad2deg(coord[0])) # galactic longitude
   phi = rad2deg(coord[1])  # galactic latitude
   ra = np.ndarray(shape=(12*nside*nside),dtype=float)
   dec = np.ndarray(shape=(12*nside*nside),dtype=float)
   for ii in range(len(theta)):
      eq_coord = ga2equ([phi[ii],theta[ii]]) # converting to equatorial coordinates
      ra[ii] = deg2rad(eq_coord[0])
      dec[ii] = deg2rad(eq_coord[1])

   return ra,dec
#=====================================================================
if __name__ == '__main__':
    from optparse import OptionParser
    o = OptionParser()
    o.set_usage('%prog [options] Catalogue of point sources or healpix map')
    o.set_description(__doc__)
    o.add_option('--healpix',dest='healpix',action='store_true',help='if input is a healpix map')
    o.add_option('--pol',dest='pol',action='store_true',help='set polarization mode: generates stokes Q and U')
    o.add_option('--jd',dest='jd',default=None,help='Npz file containing julian dates')
    o.add_option('--start',dest='start',type=float,default=1e9,help='Starting frequency in Hz')
    o.add_option('--stop',dest='stop',type=float,default=2e9,help='Stopping frequency in Hz')
    o.add_option('--chan',dest='chan',type=int,default=203,help='Number of frequency channels')
    o.add_option('--nside',dest='nside',type=int,default=128,help='Resolution of the map')
    o.add_option('--xpol',dest='xpol',default='PAPER_FF_X.ffe',help='FEKO model for the X linear polarization')
    o.add_option('--ypol',dest='ypol',default='PAPER_FF_Y.ffe',help='FEKO model for the Y linear polarization')
    o.add_option('--ref',dest='ref',type=float,default=1e9,help='Reference frequency ref')
    o.add_option('--filename',dest='filename',default='Stokes',help='Output npz file, DEFAULT:Stokes-f[FREQ]_j[JD].npz')
    o.add_option('--pmap',dest='pmap',default=None, help='Input polarized map')
    o.add_option('--cores',dest='cores',type=int,default=4,help='Number of jobs')
    opts, args = o.parse_args(sys.argv[1:])

    pol=True if opts.pol else False
    healpix=True if opts.healpix else False

    rad2deg = lambda val: val * 180./np.pi # function to convert from radians to degrees
    deg2rad = lambda val: val * np.pi/180 # function to convert from degrees to radians

    # array containing julian dates
    jds = np.load(opts.jd)['jd']
    # range of frequencies
    frequency = np.linspace(opts.start,opts.stop,opts.chan)
    np.savez(open('frequency.npz','wb'),frequency=frequency)
    
    if opts.healpix:
       ra,dec = convHealCoord(opts.nside)
       if opts.pol:
          if opts.pmap==None:
             # generates realizations of polarized map from angular power spectrum
             random.seed(1001)
             ll = np.arange(1,2700)
             A_700 =  0.64#0.05 # RM normalization factor
             P_cl = A_700*((ll/700.)**(-1.65)) # angular power spectrum
             _dat = hp.synfast(P_cl,nside) # generating healpix realization
             polmap = _dat.copy()
             pfrac=None;RM=None
          else:
             polmap = hp.read_map(opts.pmap)
             polmap = hp.ud_grade(polmap,opts.nside)
       else:
          data = hp.read_map(args[0])
          data = hp.ud_grade(data,opts.nside)
          pfrac=0.1;RM=12;alpha=-2.6 
          polmap=None
    else:
       data,ra,dec,alpha,RM,pfrac = loadCat(args[0])
       polmap=None
       print data,ra,dec,alpha

    np.savez(open('RADEC.npz','wb'),ra=ra,dec=dec)
    Parallel(n_jobs=opts.cores)(delayed(call_genForegrounds)(data,fq,opts.ref,ra,dec,opts.xpol,opts.ypol,jd,RM,opts.filename,polmap=polmap,pfrac=pfrac,pol=opts.pol) for jd in jds for fq in frequency)
