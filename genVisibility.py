#!/usr/bin/env python
"""
Generates correlator output/visibilities given the foregrounds
"""
import numpy as np
import aipy
import optparse
import os,sys

c = 3e8 # speed of light
rad2deg = lambda val: val * 180./np.pi # convert radians to degrees
deg2rad = lambda val: val * np.pi/180 # convert degrees to radians
DEC = -30.721527777777776 # latitude of PAPER
#======================================================
def raddec2lm(ra0,dec0,ra,dec):# ra and dec in radians
    """
    Converts ra/dec to direction cosines l/m

    ra0  : reference/phase right ascension; type: float
    dec0 : reference/phase declination; type:float
    ra   : right ascension in radians; type:float
    dec  : declination in radians; type:float
    """
    rad2deg = lambda val: val * 180./np.pi
    l = np.cos(dec)* np.sin(ra0 - ra)
    m = -1*(np.sin(dec)*np.cos(dec0) - np.cos(dec)*np.sin(dec0)*np.cos(ra-ra0))
    return l,m
#======================================================
if __name__=='__main__':
   o = optparse.OptionParser()
   o.set_usage('genVisibility.py [options]')
   o.set_description(__doc__)
   aipy.scripting.add_standard_options(o,cal=True,ant=True)
   o.add_option('--start',dest='start',type=float,default=1e8,help='Starting frequency in GHz')
   o.add_option('--stop',dest='stop',type=float,default=2e8,help='Stopping frequency in GHz')
   o.add_option('--chan',dest='chan',type=int,default=203,help='Number of frequency channels')
   o.add_option('--radec',dest='radec',default='RADEC.npz',help='Array with ra/dec coordiates')
   o.add_option('--jd',dest='jd',default='julian_dates.npz',help='Array with julian dates')
   o.add_option('--inname', dest='inname', default='stokes',help='Name of the array consisting the stokes parameter, default = stokes-f[FREQ]_j[JD].npz')
   o.add_option('--filename', dest='filename', default='test',help='Filename of created Miriad UV file(ex: test.uv).')
   o.add_option('--save', dest='save',action='store_true', help='If option set to save, visibilities will be stored in MIRIAD format')
   o.add_option('--stokes', dest='stokes',action='store_true', help='If option set you stokes, stokes visibilities will be stored in MIRIAD file')
   o.add_option('--healpix', dest='healpix',action='store_true', help='Enable healpis normalization')
   o.add_option('--nside', dest='nside',default=128, help='Nside to be use ifro normalization, should be as healpix map used to generate foregrounds')
   opts, args = o.parse_args(sys.argv[1:])

   freqs = np.linspace(opts.start,opts.stop,opts.chan)
   freqs = freqs * 1e-9 # frequency in GHz
   sfreq = freqs[0] 
   sdf = freqs[1]-freqs[0]
   nchan = len(freqs)

   jds = np.load(opts.jd)['jd']
   times = map(float,jds)
   inttime = (times[1]-times[0])*24*60*60 #integration time in seconds

   # u,v coordinates
   i,j = map(int,opts.ant.split('_'))
   aa = aipy.cal.get_aa(opts.cal, sfreq, sdf, nchan)
   bl = aa.get_baseline(i,j) # getting the bl coordinates

   shape = (len(times),len(freqs))
   flags = np.zeros(shape, dtype=np.int32)
   uvgridxx = np.zeros(shape, dtype=np.complex64)
   uvgridxy = np.zeros(shape, dtype=np.complex64)
   uvgridyx = np.zeros(shape, dtype=np.complex64)
   uvgridyy = np.zeros(shape, dtype=np.complex64)

   radec = np.load(opts.radec)
   ra = radec['ra'] # ra in radians
   dec = radec['dec'] # dec in radians

   if opts.healpix:
      dtheta = 4*np.pi/(12*opts.nside*opts.nside)
   else:
      dtheta = 1
   for ii, t in enumerate(times):
      lm = np.zeros((2,len(ra)),)
      aa.set_jultime(t + (inttime/2.)/aipy.const.s_per_day) # setting the phase reference to evaluate the direction cosines l and m
      lst = aa.sidereal_time()
      RA_r , DEC_r = lst,deg2rad(DEC)
      for kk in range(len(ra)):
         l, m = raddec2lm(RA_r,DEC_r, ra[kk],dec[kk]) # direction cosines
         lm[0,kk] = l; lm[1,kk] = m

         for jj, f in enumerate(freqs):
            fng = np.exp(-2j*np.pi*((bl[0]*lm[0] + bl[1]*lm[1])*f)) # evaluating the fringe
            stokes = np.load('%s-f%.4g_j%.5f.npz'%(opts.inname,f*1e3,t))['stokes'] #loading the stokes parameters

            # discrete Fourier transform
            visxx = np.sum(stokes[0,:]*fng*dtheta)
            visxy = np.sum(stokes[1,:]*fng*dtheta)
            visyx = np.sum(stokes[2,:]*fng*dtheta)
            visyy = np.sum(stokes[3,:]*fng*dtheta)

            uvgridxx[ii,jj] = visxx
            uvgridxy[ii,jj] = visxy
            uvgridyx[ii,jj] = visyx
            uvgridyy[ii,jj] = visyy

   if opts.save:
      uvfile = opts.filename + '.uv'
      print 'Writing to ', uvfile
      uv = aipy.miriad.UV(uvfile, status='new')

      uv._wrhd('obstype','mixed-auto-cross')
      uv._wrhd('history','MDLVIS: created file from scratch.\n')
      uv.add_var('telescop' ,'a'); uv['telescop'] = 'AIPY'
      uv.add_var('operator' ,'a'); uv['operator'] = 'AIPY'
      uv.add_var('version' ,'a'); uv['version'] = '0.0.1'
      uv.add_var('epoch' ,'r'); uv['epoch'] = 2000.
      uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
      uv.add_var('nchan' ,'i'); uv['nchan'] = nchan
      uv.add_var('sdf' ,'d'); uv['sdf'] = sdf
      uv.add_var('sfreq' ,'d'); uv['sfreq'] = sfreq
      uv.add_var('freq' ,'d'); uv['freq'] = sfreq
      uv.add_var('restfreq' ,'d'); uv['restfreq'] = sfreq
      uv.add_var('nschan' ,'i'); uv['nschan'] = uv['nchan']
      uv.add_var('inttime' ,'r'); uv['inttime'] = inttime
      uv.add_var('npol' ,'i'); uv['npol'] = 4
      uv.add_var('nspect' ,'i'); uv['nspect'] = 1
      uv.add_var('nants' ,'i'); uv['nants'] = 32

      uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
      uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
      uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
      uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.

      #variables to be updated

      uv.add_var('coord' ,'d')
      uv.add_var('time' ,'d')
      uv.add_var('lst' ,'d')
      uv.add_var('ra' ,'d')
      uv.add_var('obsra' ,'d')
      uv.add_var('baseline' ,'r')
      uv.add_var('pol' ,'i')

      #get antenna array
      uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
      uv.add_var('dec' ,'d'); uv['dec'] = aa.lat
      uv.add_var('obsdec' ,'d'); uv['obsdec'] = aa.lat
      uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
      uv.add_var('antpos' ,'d'); uv['antpos'] = (np.array([ant.pos for ant in aa], dtype = np.double)).transpose().flatten() #transpose is miriad convention

      for ii, t in enumerate(times):
         print '%d/%d' % (ii+1, len(times))+' done'
         aa.set_jultime(t)
         uv['time'] = t
         uv['lst'] = aa.sidereal_time()
         uv['ra'] = aa.sidereal_time()
         uv['obsra'] = aa.sidereal_time()

         preamble = (bl, t, (i,j))

         if opts.stokes:
            uv['pol'] = aipy.miriad.str2pol['I']
            uv.write(preamble, uvgridxx[ii]+uvgridyy[ii], flags[ii])
            uv['pol'] = aipy.miriad.str2pol['Q']
            uv.write(preamble, uvgridxx[ii]-uvgridyy[ii], flags[ii])
            uv['pol'] = aipy.miriad.str2pol['U']
            uv.write(preamble, uvgridxy[ii]+uvgridyx[ii],flags[ii])
            uv['pol'] = aipy.miriad.str2pol['V']
            uv.write(preamble, 1j*(uvgridyx[ii]-uvgridxy[ii]), flags[ii])

         else:
            uv['pol'] = aipy.miriad.str2pol['xx']
            uv.write(preamble, uvgridxx[ii], flags[ii])
            uv['pol'] = aipy.miriad.str2pol['xy']
            uv.write(preamble, uvgridxy[ii], flags[ii])
            uv['pol'] = aipy.miriad.str2pol['yx']
            uv.write(preamble, uvgridyx[ii],flags[ii])
            uv['pol'] = aipy.miriad.str2pol['yy']
            uv.write(preamble, uvgridyy[ii], flags[ii])
      del(uv)   
