"""
Wrapper to generate visibilities for all the baselines pairs of PAPER-32
"""
import os,sys
import optparse

if __name__=='__main__':
   o = optparse.OptionParser()
   o.set_usage('python wrapper.py [option]')
   o.set_description(__doc__)
   o.add_option('--start',dest='start',type=float,default=1e8,help='Starting frequency in GHz')
   o.add_option('--stop',dest='stop',type=float,default=2e8,help='Stopping frequency in GHz')
   o.add_option('--chan',dest='chan',type=int,default=203,help='Number of frequency channels')
   o.add_option('--radec',dest='radec',default='RADEC.npz',help='Array with ra/dec coordiates')
   o.add_option('--jd',dest='jd',default='julian_dates.npz',help='Array with julian dates')
   o.add_option('--inname', dest='inname', default='stokes',help='Name of the array consisting the stokes parameter, default = stokes-f[FREQ]_j[JD].npz')
   o.add_option('--stokes', dest='stokes',action='store_true', help='If option set you stokes, stokes visibilities will be stored in MIRIAD file')
   o.add_option('--healpix', dest='healpix',action='store_true', help='Enable healpis normalization')
   o.add_option('--nside', dest='nside',default=128, help='Nside to be use for normalization, should be as healpix map used to generate foregrounds')
   o.add_option('-C',dest='calfile',default='Calfile containing antenna positions ')
   o.add_option('-o','--output',dest='output',default='test',help='Name pf output uvfile')
   o.add_option('-n','--nants',dest='nants',type=int,default=32,help='Number of antennas')

   opts,args = o.parse_args(sys.argv[1:])

   for ii in range(opts.nants):
      for jj in range(ii+1,opts.nants):
         bl = '%s_%s'%(ii,jj)
         outfile = opts.output + '-' + bl
         if opts.stokes:
            os.system('python genVisibility.py --inname=%s --start=%s --stop=%s --chan=%s --save --stokes --jd=%s -C %s --radec=%s  -a %s --filename=%s'%(opts.inname,opts.start,opts.stop,opts.chan,opts.jd,opts.calfile,opts.radec,bl,outfile))
         else:
             os.system('python genVisibility.py --inname=%s --start=%s --stop=%s --chan=%s --save --jd=%s -C %s --radec=%s  -a %s --filename=%s'%(opts.inname,opts.start,opts.stop,opts.chan,opts.jd,opts.calfile,opts.radec,bl,opts.outfile))
