import numpy as np
import optparse
import aipy
import os, sys
from collections import OrderedDict
import pickle
import scipy

o = optparse.OptionParser()
o.set_usage('python corrupt_vis.py [options] uv')
aipy.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('-b', '--baseline', dest='baseline', help='Baseline or antenna pair e.g 0_1')
o.add_option('--pkl', dest='pkl', type='str', help='pickle file containing 2D array (time, freq) representing the complex gains for each antenna.')
o.add_option('--crange', dest='crange', type=str, default='40_160', help='Frequency range e.g 80_90')
o.add_option('--cen', dest='cen', type=int, default=60, help='Centre frequency channel.')
o.add_option('--outfile', dest='outfile', type='str', default='corrupted_vis.uv', help='Name of output uvfile name. Default is corrupted_vis.uv')
opts, args = o.parse_args(sys.argv[1:])

ants = map(int, opts.baseline.split('_'))
aa = aipy.cal.get_aa(opts.cal, np.array([0.15]))
bl = aa.get_baseline(ants[0], ants[1])
gain = pickle.load(open(opts.pkl))
keys = gain.keys()
new_dict = OrderedDict()
freqs = np.linspace(120, 180, 120)
chan = map(int, opts.crange.split('_'))
c0, c1 = chan[0], chan[1]
chans = np.arange(c0, c1)
cfreqs = np.linspace(100, 200, 203)
for uvfile in args:
    uv = aipy.miriad.UV(uvfile)
    aipy.scripting.uv_selector(uv, opts.ant, opts.pol)
    new_uvfile = opts.outfile
    uvo = aipy.miriad.UV(new_uvfile, status='new')
    k = 0
    for (uvw, t, (i, j)), _d, f in uv.all(raw=True):
	if not t in new_dict.keys(): new_dict[t]=OrderedDict()
        assert i != j
        gain0 = gain[i]['gain'][c0:c1]
        gain1 = gain[j]['gain'][c0:c1]
        flag0 = gain[i]['flag'][c0:c1]
        flag1 = gain[j]['flag'][c0:c1]
        g0_mx = gain0[opts.cen]
        g1_mx = gain1[opts.cen]
        inds0 = np.where(np.abs(gain0)!=1)
        inds1 = np.where(np.abs(gain1)!=1)
        if len(inds0[0]) == 0 and len(inds1[0]) == 0:
            new_dict[t]['flag'][:] = True 
        else:
            if len(inds0[0]) == 0:
                g0 = gain1[inds1[0]] / g1_mx
                g1 = gain1[inds1[0]] / g1_mx
                f0 = inds1[0][0] 
                f1 = inds1[0][-1]
                freqs0 = freqs[inds1]
                freqs1 = freqs[inds1]
            elif len(inds1[0]) == 0:
                g0 = gain0[inds0[0]] / g0_mx
                g1 = gain0[inds0[0]] / g0_mx
                f0 = inds0[0][0]
                f1 = inds0[0][-1]
                freqs0 = freqs[inds0]
                freqs1 = freqs[inds0]
            else:
                g0 = gain0[inds0[0]] / g0_mx
                g1 = gain1[inds1[0]] / g1_mx
                f0 = inds0[0][0]
                f1 = inds0[0][-1]
                freqs0 = freqs[inds0]
                freqs1 = freqs[inds1]
            newfreqs = cfreqs[f0:f1] 
            print f0, f1
            interp_gains0 = scipy.interp(newfreqs, freqs0, g0)
            interp_gains1 = scipy.interp(newfreqs, freqs1, g1)
            new_dict[t]['data'] = interp_gains0 * _d[f0:f1] * np.conj(interp_gains1)
            f[0:f0] = True; f[f1:] = True
            new_dict[t]['flag'] = f[f0:f1]
            k += 1
    uvo._wrhd('obstype','mixed-auto-cross')
    uvo._wrhd('history','MDLVIS: created file from scratch.\n')
    uvo.add_var('telescop' ,'a'); uvo['telescop'] = 'AIPY'
    uvo.add_var('operator' ,'a'); uvo['operator'] = 'AIPY'
    uvo.add_var('version' ,'a'); uvo['version'] = '0.0.1'
    uvo.add_var('epoch' ,'r'); uvo['epoch'] = 2000.
    uvo.add_var('source'  ,'a'); uvo['source'] = 'zenith'
    uvo.add_var('nchan' ,'i'); uvo['nchan'] = uv['nchan']
    uvo.add_var('sdf' ,'d'); uvo['sdf'] = uv['sdf']
    uvo.add_var('sfreq' ,'d'); uvo['sfreq'] = uv['sfreq']
    uvo.add_var('freq' ,'d'); uvo['freq'] = uv['freq']
    uvo.add_var('restfreq' ,'d'); uvo['restfreq'] = uv['restfreq']
    uvo.add_var('nschan' ,'i'); uvo['nschan'] = uvo['nchan']
    uvo.add_var('inttime' ,'r'); uvo['inttime'] = uv['inttime']
    uvo.add_var('npol' ,'i'); uvo['npol'] = 1
    uvo.add_var('nspect' ,'i'); uvo['nspect'] = 1
    uvo.add_var('nants' ,'i'); uvo['nants'] = 32

    uvo.add_var('vsource' ,'r'); uvo['vsource'] = 0.
    uvo.add_var('ischan'  ,'i'); uvo['ischan'] = 1
    uvo.add_var('tscale'  ,'r'); uvo['tscale'] = 0.
    uvo.add_var('veldop'  ,'r'); uvo['veldop'] = 0.

    #variables to be updated
    uvo.add_var('coord' ,'d')
    uvo.add_var('time' ,'d')
    uvo.add_var('lst' ,'d')
    uvo.add_var('ra' ,'d')
    uvo.add_var('obsra' ,'d')
    uvo.add_var('baseline' ,'r')
    uvo.add_var('pol' ,'i')

    times = new_dict.keys()
    for ii, tt in enumerate(times):
        print '%d/%d' % (ii+1, len(times))+' done'
        aa.set_jultime(tt)
        uvo['time'] = t
        uvo['lst'] = aa.sidereal_time()
        uvo['ra'] = aa.sidereal_time()
        uvo['obsra'] = aa.sidereal_time()

        preamble = (bl, tt, (i,j))
        uvo['pol'] = aipy.miriad.str2pol['xx']
        uvo.write(preamble, new_dict[tt]['data'], new_dict[tt]['flag'])
    
    del(uv)
