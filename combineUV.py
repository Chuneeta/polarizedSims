#/usr/bin/env python

""" 
Combining multiple uv files into a single one file 
"""

import aipy
import os,sys
import optparse

o = optparse.OptionParser()
o.set_usage('combineUV.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-o', dest='outfile', default='combine.uv',help='Filename of created Miriad UV file (ex: combine.uv).')
opts, args = o.parse_args(sys.argv[1:])

uvo = aipy.miriad.UV(opts.outfile, status='new')

for uvfile in args:
   uvi = aipy.miriad.UV(uvfile)
   print uvfile,'->',opts.outfile

   uvo.init_from_uv(uvi)
   uvo.pipe(uvi)

del(uvo)

