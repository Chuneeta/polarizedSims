"""
Storing calibration solutions from CASA into a dictoionary
"""
from pyrap.tables import table
import numpy as np
import os, sys
import optparse
import pickle

o = optparse.OptionParser()
o.set_usage('python cal2npz.py [options] *.cal')
o.set_description(__doc__)
o.add_option('--smooth', dest='smooth', action='store_true', help='Smoothen the calibration solutions')
o.add_option('-o', dest='outfile', default='gain_solution', help='Name of output file. Default is gain_solution')
opts, args = o.parse_args(sys.argv[1:])

for cal in args:
    print ('Loading {}'.format(cal))
    tb = table(cal)
    gsoln = tb.getcol('CPARAM')
    flag = tb.getcol('FLAG')
    _sh =  gsoln.shape

    freqs = np.linspace(100,200,203)
    gain_dict = {}
    for ii in range(_sh[0]):
        if not ii in gain_dict.keys(): gain_dict[ii] = {}
        sgains = gsoln[ii, :, 0]
        flags = flag[ii, :, 0]
        gain_dict[ii]['flag'] = flags
        xdata = freqs[~flags]
        if len(xdata) > 0:
            gains = np.ones((len(sgains)), dtype=np.complex64)
            ydata = sgains[~flags]
            if opts.smooth:
                weights = np.polyfit(xdata, np.abs(ydata), 3)
                model = np.poly1d(weights)
                gains[~flags] = np.abs(model(xdata)) * np.exp(1j * np.angle(ydata)) # preserving the original phase 
                outfile = opts.outfile + '.smooth.pkl'
            else:
                gains[~flags] = ydata
                outfile = opts.outfile + '.pkl'
            gain_dict[ii]['gain'] = gains
        else:
            gain_dict[ii]['gain'] = np.ones((len(freqs)), dtype=np.complex64)

with open(outfile, 'wb') as fl:
    pickle.dump(gain_dict, fl)
