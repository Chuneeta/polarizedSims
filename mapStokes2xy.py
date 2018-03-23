"""
Mapping (I,Q,U V) --> (xx,xy,yx,yy)
"""
import aipy as a
import sys

from optparse import OptionParser
o = OptionParser()
o.set_usage('%prog [options] uv')
o.set_description(__doc__)
o.add_option('-p',dest='pol',default='xx',help='polarization e.g xx,xy,yx,yy')    

opts, args = o.parse_args(sys.argv[1:])
filename = sys.argv[1]
uvi = a.miriad.UV(filename)	

pol = opts.pol
if pol=='xx':
  print '%s -> %s'%(filename,filename+'.I2xx')
  uvo_I=a.miriad.UV(filename+'.I2xx',status='new')
  uvo_I.init_from_uv(uvi,override={'pol':-5})
  uvo_I.pipe(uvi,append2hist='MAP: I <--> xx \n')
  del(uvo_I)
  del(uvi)

if pol=='yy':
  print '%s -> %s'%(filename,filename+'.Q2yy')
  uvo_I=a.miriad.UV(filename+'.Q2yy',status='new')
  uvo_I.init_from_uv(uvi,override={'pol':-6})
  uvo_I.pipe(uvi,append2hist='MAP: Q <--> yy \n')
  del(uvo_I)
  del(uvi)


if pol=='xy':
  print '%s -> %s'%(filename,filename+'.U2xy')
  uvo_I=a.miriad.UV(filename+'.U2xy',status='new')
  uvo_I.init_from_uv(uvi,override={'pol':-7})
  uvo_I.pipe(uvi,append2hist='MAP: U <--> xy \n')
  del(uvo_I)
  del(uvi)


if pol=='yx':
  print '%s -> %s'%(filename,filename+'.V2yx')
  uvo_I=a.miriad.UV(filename+'.V2yx',status='new')
  uvo_I.init_from_uv(uvi,override={'pol':-8})
  uvo_I.pipe(uvi,append2hist='MAP: V <--> yx \n')
  del(uvo_I)
  del(uvi)


