#!/usr/bin/env python
""" Rebin individual HIRES exposures to a common wavelength scale, and
make a rough combination of all the spectra to act as a template for
use scaling the spectra.
"""

import numpy as np
try:
    from astropy.io import fits 
except ImportError:
    import pyfits as fits
from pylab import plot,show

import sys, os

from barak.utilities import between
from .core import read_HIRES_wascale, rebin


def make_template(wa, wavs, fluxs):
    """ This makes a quick, unweighted combination of all the
    orders. The individual orders will be scaled to match this
    template.
    """
    work = []
    for i in range(len(wavs)):
        f = np.interp(wa, wavs[i], flux[i])
        f[~between(wa, wavs[i][0], wavs[i][-1])] = np.nan
        work.append(f)

    work = np.array(work)
    fsum = np.nansum(work, axis=0)
    work[~np.isnan(work)] = 1
    
    return fsum / np.nansum(work, axis=0)


# read all the (extracted) spectra
if 1:
    fluxnames = sorted(sys.argv[1:])
    assert fluxnames
    if 'Flux-' in fluxnames[0]:
        errnames = [n.replace('Flux-', 'Err-') for n in fluxnames]
        skynames = [n.replace('Flux-', 'Sky-') for n in fluxnames]
    elif '_Flux.fits' in fluxnames[0]:
        errnames = [n.replace('_Flux.fits','_Err.fits') for n in fluxnames]
        skynames = [n.replace('_Flux.fits','_Sky.fits') for n in fluxnames]
    else:
        raise ValueError('Unknown nameing scheme!')

    Nexp = len(fluxnames)

    wavs = []
    flux = []
    error = []
    bgs = []
    IDs = []
    for i in range(Nexp):
        fluxname = fluxnames[i]
        errname = errnames[i]
        skyname = skynames[i]
        print 'reading', fluxname, '+ Err + Sky' 
        hd = fits.getheader(fluxname)
        wa = read_HIRES_wascale(hd)
        wavs.extend(wa)
        flux.extend(fits.getdata(fluxname))
        error.extend(fits.getdata(errname))
        bgs.extend(fits.getdata(skyname))
        IDs.extend([fluxname.rsplit('/',1)[0]]*len(wa))

    # for i in range(len(wavs)):
    #         plot(wavs[i], flux[i])
    #         plot(wavs[i], error[i])
    # show()

if 1:
    # rebin orders
    # sort in order of increasing first bin wavelength
    isort = np.argsort([wa[0] for wa in wavs])
    wavs = [wavs[i] for i in isort]
    flux = [flux[i] for i in isort]
    error = [error[i] for i in isort]
    bgs = [bgs[i] for i in isort]
    IDs = [IDs[i] for i in isort]
    wmin = min([wa[0] for wa in wavs])
    wmax = max([wa[-1] for w in wavs])

    # Create a wavelength scale with constant
    # that covers all the individual exposures
    delta = 2.2 / 3.e5
    print 'using %.2f km/s pixels' % (delta * 3e5)
    dlog10wa = np.log10(1 + delta)
    log10wa = np.arange(np.log10(wmin), np.log10(wmax), dlog10wa)
    wa = 10**log10wa

    print 'making template'
    template = make_template(wa, wavs, flux)
    np.save('template.npy',np.rec.fromarrays([wa, template], names='wa,fl'))

    print 'rebinning'
    if not os.path.lexists('./rebinned'):
        os.mkdir('./rebinned')
        
    for i in range(len(wavs)):
        if not i % 20: print i
        fl, er = rebin(wavs[i], flux[i], error[i], wa)
        bg = np.interp(wa, wavs[i], bgs[i])
        sp = np.rec.fromarrays([wa, fl, er, bg], names='wa,fl,er,bg')
        np.save('rebinned/%03i.npy' % i, sp)
