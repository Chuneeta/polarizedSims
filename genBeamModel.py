#!/usr/bin/env python

import numpy as np
import healpy as hp
import pyfits
import ephem
import fmt
import math


def ga2equ(ga):
  """
  Convert Galactic to Equatorial coordinates (J2000.0)
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

  

def InterBeam(xpol,ypol,freq,jd,ra,dec):
  """
  Interpolate beam at given frequency and sky position (RA and DEC)

  - xpol : Input beam model for xx polarization; type:string
  - ypol : Input beam model for yy polarization; type:string
  - freq : Frequency in Hz; type:float
  - ra   : Right ascensions; type: numpy.ndarray
  - dec  : Declinations; type:numpy.ndarray
  """

  # Reading FEKO BEAMS ..."
  fekoX=fmt.FEKO(xpol)
  fekoY=fmt.FEKO(ypol)
  fields = fekoX.fields
  feko_xpol=fekoX.fields[0]
  feko_ypol=fekoY.fields[0]
  phi=feko_xpol.phi*np.pi/180. # azimuth
  theta=feko_xpol.theta*np.pi/180. # zenith
  theta = np.pi/2 - theta # pyephem wants the elevation rather than the zenith angle
  
  beamFreqs = np.zeros(11)
  beamFreqs[0] = 100e6
  for i in range(1,11):
     beamFreqs[i] =  beamFreqs[i-1] + 10e6
 
  # Computing BEAMS
  gxx=feko_xpol.etheta*np.conj(feko_xpol.etheta)+feko_xpol.ephi*np.conj(feko_xpol.ephi)
  beam = np.ndarray(shape=(gxx.shape[0],4,len(beamFreqs)),dtype=complex)
  for j in range(len(beamFreqs)):
    feko_xpol = fekoX.fields[j]
    feko_ypol = fekoY.fields[j]

    for k in range(len(phi)):
       R_phi = np.array([[np.cos(phi[k]),-1*np.sin(phi[k])],[np.sin(phi[k]),np.cos(phi[k])]]) # rotation matrix to convert for Ludwig-3 coordinates to normal alt-az coordinate system
       Ex = np.array([[feko_xpol.etheta[k]],[feko_xpol.ephi[k]]]).reshape(2)
       Ey = np.array([[feko_ypol.etheta[k]],[feko_ypol.ephi[k]]]).reshape(2)
       Ex_hv = np.dot(R_phi,Ex)
       Ey_hv = np.dot(R_phi,Ey)

       beam[k,0,j] = Ex_hv[0]
       beam[k,1,j] = Ex_hv[1]
       beam[k,2,j] = Ey_hv[0]
       beam[k,3,j] = Ey_hv[1]

    beam[:,0,j] = beam[:,0,j]/np.max(beam[:,0,j])
    beam[:,3,j] = beam[:,3,j]/np.max(beam[:,3,j])
  # Create OBSERVER
  paper = ephem.Observer()
  paper.lat, paper.long, paper.elevation = '-30:43:17', '21:25:40.08', 0.0
  j0 = ephem.julian_date(0)
  paper.date = float(jd) - j0 + 5./60/24

  beamRA = np.ndarray(shape=(phi.shape[0]),dtype=float)
  beamDEC = np.ndarray(shape=(phi.shape[0]),dtype=float)
  for k in range(beamRA.shape[0]):
    ra0,dec0 = paper.radec_of(phi[k],theta[k])
    beamRA[k] = ra0  # RA in radians
    beamDEC[k] = dec0    # DEC in radians

  print "FREQUENCY INTERPOLATION"
  InterBeamF = np.ndarray(shape=(beam.shape[0],beam.shape[1]),dtype=complex)
  dist = np.abs(freq - np.array(beamFreqs)) # calculating the frequency distance
  ind = np.argsort(dist)
  ind.flatten() 
  
  for p in range(4):
     if dist[ind[0]] ==0:
        InterBeamF[:,p] = beam[:,p,ind[0]]
     else:  
        # 2-closest point interpolation
        weights = dist[ind[0]]**(-2) + dist[ind[1]]**(-2)       
        InterBeamF[:,p] = (beam[:,p,ind[0]]*dist[ind[0]]**(-2) + beam[:,p,ind[1]]*dist[ind[1]]**(-2))/weights
 
  InterHealBeam = np.ndarray(shape=(ra.shape[0],beam.shape[1]),dtype=complex) 
  
  print "SPATIAL INTERPOLATION"
  for r in range(ra.shape[0]):
    ra[r] = ra[r] - 2*np.pi if ra[r]>(2*np.pi) else ra[r] # wrapping over ra
    dist = np.sqrt((ra[r] - beamRA)**2 + (dec[r] - beamDEC)**2) # calculating distance
    ind = np.argsort(dist)
    ind = ind.flatten()
    weights = dist[ind[0]]**(-1) + dist[ind[1]]**(-1) + dist[ind[2]]**(-1)
    
    for p in range(4):
       # 3-closest point interpolation
       InterHealBeam[r,p] = InterBeamF[ind[0],p]*dist[ind[0]]**(-1) + InterBeamF[ind[1],p]*dist[ind[1]]**(-1) + InterBeamF[ind[2],p]*dist[ind[2]]**(-1) 
       InterHealBeam[r,p] = InterHealBeam[r,p]/weights
  

  np.savez(open('Beam-f%.4g_j%.5f.npz'%(freq*1e-6,np.float(jd)),'wb'),beam=InterHealBeam)
  return InterHealBeam

