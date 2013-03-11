#!/usr/bin/env python
import numpy as np

from math import sqrt
import os
from glob import glob

from barak.convolve import convolve_psf
from barak.utilities import between, stats
from barak.io import writetable

np.seterr(divide='ignore')

def make_sky(skyvals, ervals):
    """ Make a rough combination of the background for all orders.
    """
    sky2d = np.array(skyvals)
    er2d = np.array(ervals)    
    sky2d[~(er2d > 0)] = np.nan
    sky = np.nansum(sky2d, axis=0)
    sky2d[~np.isnan(sky2d)] = 1
    sky = sky / np.nansum(sky2d, axis=0)
    return sky

def combine(fl, er, ser, mask, NSIG=3):
    """ Combine spectra, weighting flux by the smoothed error and
    clipping pixels with high residuals, and masking the requested
    pixels.

    in:  fluxes (m x n), errors (m x n), smoothed errors (m x n), masks (m x n)
    out: flux (n), error (n), npts (n), clipped (m x n) 

    m is the number of spectra, n is the number of pixels in each spectrum.
    """
    median = np.median
    fl, er, ser, m  = map(np.asarray, (fl,er,ser,mask))
    good = (er > 0) & ~m

    Nspec, Npix = fl.shape

    clipped = np.zeros(fl.shape, bool)
    cfl = []
    cer = []
    npts = []
    for i in range(Npix):
        if not i % 10000: print i

        ftot = 0
        etot = 0
        wtot = 0
        n = 0

        c0 = good[:, i-1:i+2]
        vals = fl[:, i-1:i+2][c0]
 
        if len(vals) > 0:
            med = median(vals)
            for j in range(Nspec):
                if not er[j,i] > 0 or m[j,i]:
                    continue
                
                resid = (fl[j,i] - med)/er[j,i]
                if abs(resid) > NSIG:
                    clipped[j,i] = True 
                    continue
    
                w = ser[j,i]**-2
                ftot += fl[j,i] * w
                etot += (er[j,i] * w)**2
                wtot += w
                n += 1

        if wtot == 0:
            f = np.nan
            if len(vals) > 5:
                f = med
            e = 0
        else:
            f = ftot / wtot
            e = sqrt(etot) / wtot
        cfl.append(f)
        cer.append(e)
        npts.append(n)

    return np.array(cfl), np.array(cer), np.array(npts), clipped

if 1:
    sp0 = [np.load(n).view(np.recarray) for n in sorted(glob('./rebinned/*.npy'))]
    mults = np.load('multipliers.npy')

    sp = []
    for s in sp0:
        # smoothed error
        ser = np.zeros_like(s.fl)
        clipped = np.zeros(len(s.fl), bool)
        masked = np.zeros(len(s.fl), bool)
        sp.append(np.rec.fromarrays([s.wa, s.fl, s.er, ser, clipped, masked, s.bg],
                                    names='wa,fl,er,ers,clipped,masked,bg'))

    print 'Applying scaling and masking regions'

    masks = {
          #67: [(6320.8, 6321.4)],
             }
    wa = sp[0].wa
    for i,s in enumerate(sp):
        s.fl *= mults[i]
        s.er *= mults[i]
        s.ser = convolve_psf(s.er, 25)

        # mask the last 5 pixels of each order
        s.masked[-5:] = True

        # bad regions are mostly sky lines
        s.masked[between(wa, 4879.9, 4881.6)] = True
        s.masked[between(wa, 5429.7, 5431.4)] = True
        s.masked[between(wa, 5577.7, 5578.3)] = True
        s.masked[between(wa, 5599, 5601)] = True 
        s.masked[between(wa, 5890.4, 5890.9)] = True
        s.masked[between(wa, 5896.5, 5896.9)] = True
        s.masked[between(wa, 6798.4, 6798.6)] = True
        s.masked[between(wa, 6864.5, 6865.2)] = True
        if i in masks:
           for w0,w1 in masks[i]:
               s.masked |= between(wa, w0, w1)

if 1:
    print 'combining spectra'
    # combine spectra, weighting flux by the smoothed error

    # in: fluxes (m x n), errors (m x n), smoothed errors (m x n), masks (m x n)
    # out: flux (n), error (n), npts (n), clipped (m x n) 

    fl = [s.fl for s in sp]
    er = [s.er for s in sp]
    ser = [s.ser for s in sp]
    masks = [s.masked for s in sp]
    fl, er, npts, clipped = combine(fl, er, ser, masks)
    for i,s in enumerate(sp):
        s.clipped = clipped[i]
        np.save('./rebinned/%03i.npy' % i, s)
    
    # now combine the sky
    bg = make_sky([s.bg for s in sp], [s.er for s in sp])

    rec = np.rec.fromarrays([wa,fl,er,npts,bg],names='wa,fl,er,npts,bg')
    writetable('combined.fits', rec)
