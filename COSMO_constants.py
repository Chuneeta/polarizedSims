#! /usr/bin/env python

import numpy as np
import astropy.cosmology as CS

C = 2.99e8 # SPEED OF LIGHT IN M/S
F21 = 1.42040575177 # FREQUENCY OF 21 CM HYDROGEN LINE IN GHZ
COSMO = CS.FlatLambdaCDM(H0=100.0, Om0=0.27)  # Using H0 = 100 km/s/Mpc

def fq2z(fq):
   '''
   Redshift corresponding the specified frequency

   Input(s)
      fq :  [scalar] frequency in Hz
   '''
   z = F21/fq-1
   return z

def z2fq(z):
   '''
   Frequency (in Hz) corresponding the specified redshift

   Input(s)
      z :  [scalar] redshift 
   '''
   fq = F21/(z+1)
   return fq

   
def transverse_comoving_distance(z):
   '''
   Transverse comoving distance at redshift z corresponding to an angular separation of 1 radian in Mpc/h 

   Input(s)
      z :  [scalar] redshift
 
   '''
   Dz =  COSMO.comoving_distance(z).value # Mpc/h 
   return Dz


def comoving_depth(B,z):
   '''
   Comoving line-of-sight depth corresponding to specified  redshift and bandwidth for redshifted
   21 cm line in Mpc/h
  
   Input(s)
      B :    [scalar] Observing bandwith in Hz

      z :    [scalar] redshift
   '''
   deltaD =  (C/1e3) * B * (1+z)**2 /F21/COSMO.H0.value/COSMO.efunc(z) # Mpc/h 
   return deltaD


def dkprll_deta(z):
   '''
   Constant to transform delays to line-of-sight wavenumbers corresponding to redshift and 21 CM HI line
   in h/Mpc
   
   Input(s)
      z :  [scalar] redshift
   '''
   #omega_m = 0.27
   #return 2*np.pi/((1.7 / 0.1) * ((1+z) / 10.)**.5 * (omega_m/0.15)**-0.5 * 1e3)
   return 2 * np.pi * COSMO.H0.value * F21 * COSMO.efunc(z) / C /  (1+z)**2 * 1e3

def k_parallel(delays, z):
   '''
   Compute line-of-sight wavenumbers corresponding to specified delays and redshift for redshifted 21 cm line in h/Mpc

   Input(s):
      z : [scalar] redshift
   '''
   k_pllel = dkprll_deta(z) * delays
   return k_pllel

def k_perp(z):
   '''
   Compute transverse wavenumbers corresponding for redshifted 21 cm line in h/Mpc

   Input(s)
      z              : [scalar] redshift
   '''
   kperp = 2 * np.pi / transverse_comoving_distance(z)
   return kperp

def horizon_limit(z):
   '''
   Compute the horizon limit in Mpc/h at a given redshift for a specified baseline length (from Thyagarajan 2013)

   Input(s):
      baseline_length : [scalar] baseline length in wavelengths

      z : [scalar] redshift
   '''
   hlimit = COSMO.H0.value * COSMO.efunc(z) * transverse_comoving_distance(z) / C / (1+z) * 1e3
   return hlimit

