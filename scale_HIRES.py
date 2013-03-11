#!/usr/bin/env python
import numpy as np
import pylab as pl

import os
from glob import glob

from barak.convolve import convolve_psf
from barak.utilities import between
from barak.spline import fit_spline
from barak.plot import subplots

np.seterr(divide='ignore', invalid='ignore')

def scale_spec_single(wa, ratio, chunkwidth=10, plot=True):  
    """ find a smoothly varying function, such that when multiplied by
    each spectrum in spectra, they match the flux in spectrum s0.

    Essentially forces spectra to have the same large-scale power as
    s0.
    """
    # need 1 bin every chunkwidth Angstroms
    nbins = int((wa[-1] - wa[0]) / chunkwidth)

    if nbins < 4:
        return None
        #return np.ones_like(fl0) * np.median(ratio[c0])    
    try:
        return fit_spline(wa, ratio, bins=nbins)
    except RuntimeError:
        return None

if 1:
    sp = [np.load(n).view(np.recarray) for n in sorted(glob('./rebinned/*.npy'))]

    T = np.load('template.npy').view(np.recarray)

    print 'scaling to match template'
    log = True
    # multiply to match first spectrum.
    if log:
        if not os.path.lexists('./scaling_plots'):
            os.mkdir('./scaling_plots')
        ax1,ax2 = subplots(4,3,':,:2 :,2',height=4)
        ax1.set_autoscale_on(0)
        ax2.set_autoscale_on(0)
        
    mults = [np.ones_like(s.fl) for s in sp]
    i = 0
    # only use pixels with at least this (smoothed) SNR for measuring scaling.
    minSN = 1
    option = 's'
    while i < len(sp):
        #if not i % 100: print i
        wa = sp[i].wa
        fl = sp[i].fl
        er = sp[i].er

        gooder = er > 0
        snr = np.median(fl[gooder] / er[gooder])

        wa0 = wa[gooder]
        goodwa = between(wa, wa0[0], wa0[-1])

        fls = convolve_psf(fl, 10)
        ser = convolve_psf(er, 10)
        goodsn = (fls / ser) > minSN

        good =  goodsn & goodwa & gooder
        #print good.sum(), 'good points'

        if snr < 0.5:
            # really low S/N, don't do any scaling
            option = 'k'
        elif good.sum() < 200 or snr < 1.5:
            # don't scale using a spline fit
            option = 'm'

        if option == 'k':
            pass
        elif option == 'm':
            # take median of flux ratios
            ratio = T.fl / fl
            med = np.median(ratio[good])
            mults[i][:] = med
        elif option == 's':
            # fit a spline to flux ratios
            ratio = T.fl / fl
            spline = scale_spec_single(wa[good], ratio[good])
            if spline is None:
                print 'too few points to fit spline, using median'
                option = 'm'
                continue
            m = spline(wa[goodwa])
            mults[i][goodwa] = m 
            mults[i][wa < wa0[0]] = m[0]
            mults[i][wa > wa0[-1]] = m[-1]
        else:
            raise ValueError('unknown option', option)

        if log:
            # save plots
            ax1.cla()
            ax2.cla()
            ax1.plot(wa, convolve_psf(T.fl, 5), 'b')
            ax1.plot(wa, convolve_psf(fl, 5), 'g--', alpha=0.5)
            ax1.plot(wa, convolve_psf(fl*mults[i], 5), 'r', alpha=0.8)
            ymax = np.percentile(T.fl[between(wa,wa0[0],wa0[-1])], 90)
            ax1.set_ylim(-0.5*ymax, 2*ymax)
            ax1.set_xlim(wa0[0], wa0[-1])
            ax1.axhline(0, color='k', ls='--')

            if option in 'ms':
                ax2.plot(wa[goodwa], ratio[goodwa], 'b.', alpha=0.5)
                ax2.plot(wa[good], ratio[good], 'ko')
                ax2.plot(wa[goodwa], mults[i][goodwa], 'r')
                ymax = np.percentile(mults[i][goodwa], 90)
                ax2.axhline(1, color='k', ls='--')
                ax2.set_xlim(wa0[0], wa0[-1])
                ax2.set_ylim(-1*ymax, 3*ymax)

            text = '%i %s %.2f %s' % (i, option, snr, mults[i][:3])
            ax2.set_title(text)

            pl.savefig('./scaling_plots/%03i.png' % i, dpi=50)

            #print text 
            #option = raw_input('Accept, try [m]edian or [k]eep previous? (accept) ')

            #if option.strip() == '':
            #    i += 1
            #    option = 's'
        i += 1

    np.save('multipliers.npy', np.array(mults))
