""" Code for dealing with HIRES spectra.
"""
import numpy as np
try:
    from astropy.io import fits
except ImportError:
    import pyfits as fits

from collections import OrderedDict
from barak.utilities import concat_recarrays, between
from barak.io import writetable
from barak.spec import find_bin_edges
from atpy import Table

from math import sqrt

import pylab as pl

def get_datetime(n):
    """ Read the date keyword from a raw HIRES file (old or new detector)
    """
    hd = fits.getheader(n)
    try:
        s = hd['DATE']
    except KeyError:
        s = hd['DATE-OBS'] + 'T' + hd['UT']
    return s

def read_HIRES_wascale(name):
    """
    Given a HIRES reduced fits filename, return the wavelength array
    for each order.

    From the makee web pages:

    The 2D spectra (output by makee) have a wavelength scale for each
    echelle order specified by a 6th order polynomial of the form:
    wavelength in Angstroms = coef(0) + coef(1)*pixel +
    coef(2)*pixel*pixel + ... , where pixel is the column number.

    For each order, two FITS header cards are added. The two card
    names have the format: WV_0_## and WV_4_## where '##' is the
    echelle order number on a scale where the bluest order of a given
    exposure is number 01. Each card contains four values:

    WV_0_## = ' coef(0) coef(1) coef(2) coef(3) ' 
    WV_4_## = ' coef(4) coef(5) coef(6) coef(7) ' 

    where coef(7) will always be set to 0. Example: 

    WV_0_05 = ' 3.798860003E+03 +2.866401731E-02 -7.524915946E-07
    -8.422777264E-11'

    WV_4_05 = ' 8.209366207E-14 -5.373050237E-17 +1.107868830E-20
    +0.000000000E+00'

    for echelle order number 5 (also row number 5 in the output spectrum.)
    """
    hd = fits.getheader(name, 0)
    if 'NAXIS2' not in hd:
        d = fits.getdata(name)
        Nord = len(d[0])
        Npix = len(d)
    else:
        Nord = hd['NAXIS2']
        Npix = hd['NAXIS1']
    wavs = []
    for n in range(Nord):
        order = n + 1
        key = 'WV_0_%02i' % order
        if key not in hd:
            #print 'warning, using CO instead of', key
            key = 'CO_0_%02i' % order
        coeff = map(float, hd[key].split())
        #print hd[key]
        key = 'WV_4_%02i' % order
        if key not in hd:
            #print 'warning, using CO instead of ', key
            key = 'CO_4_%02i' % order
        #print hd[key]
        coeff.extend(map(float, hd[key].split()))
        #print coeff

        # was + 1 instead of + CRVAL1
        pix = np.arange(Npix) + int(hd['CRVAL1'])
        wa = np.polyval(list(reversed(coeff)), pix)
        wavs.append(wa)

    return wavs


def join_koa(filenames, outname):
    """ Join two or more IPAC (e.g. sci_*.tbl and cal_*.tbl) files
    listing KOA entries into a single fits file called outname.
    """
    if not outname.endswith('.fits'):
        outname += '.fits'

    T = [Table(filename) for filename in filenames]
    Tout = concat_recarrays([t.data for t in T])

    writetable(outname, Tout)


def rebin(wa, fl, er, rwa, debug=False):
    """ Rebins spectrum to a new wavelength scale.

    wa, fl, er : Old spectrum
    rwa        : Wavelength scale to which we rebin

    Returns the rebinned flux and error.

    Note - will probably get the flux and errors for the first and
    last pixel of the rebinned spectrum wrong.

    General pointers about rebinning if you care about errors in the
    rebinned values:

    1. Don't rebin to a smaller bin size.
    2. Be aware when you rebin you introduce correlations between
       neighbouring points and between their errors.
    3. Rebin as few times as possible.
    """
    # Note: 0 suffix indicates the old spectrum, 1 the rebinned spectrum.
    colors= 'brgy'
    
    # rebinned spectrum
    rfl = np.zeros_like(rwa)
    rer = np.zeros_like(rwa)

    # find pixel edges, used when rebinning
    edges0 = find_bin_edges(wa)
    edges1 = find_bin_edges(rwa)
    if debug:
        pl.clf()
        x0,x1 = edges1[0:2]
        yh, = pl.bar(x0, 0, width=(x1-x0),color='gray',
                    linestyle='dotted',alpha=0.3)
    widths0 = edges0[1:] - edges0[:-1]
    npts0 = len(wa)
    npts1 = len(rwa)
    df = 0.
    de2 = 0.
    npix = 0    # number of old pixels contributing to rebinned pixel,
    j = 0                # index of rebinned array
    i = 0                # index of old array

    # sanity check
    if edges0[-1] < edges1[0] or edges1[-1] < edges0[0]:
        raise ValueError('Wavelength scales do not overlap!')
    
    # find the first contributing old pixel to the rebinned spectrum
    if edges0[i+1] < edges1[0]:
        # Old wa scale extends lower than the rebinned scale. Find the
        # first old pixel that overlaps with rebinned scale.
        while edges0[i+1] < edges1[0]:
            i += 1
        i -= 1
    elif edges0[0] > edges1[j+1]:
        # New rebinned wa scale extends lower than the old scale. Find
        # the first rebinned pixel that overlaps with the old spectrum
        while edges0[0] > edges1[j+1]:
            rfl[j] = np.nan
            rer[j] = np.nan
            j += 1
        j -= 1
    lo0 = edges0[i]      # low edge of contr. (sub-)pixel in old scale
    while True:
        hi0 = edges0[i+1]  # upper edge of contr. (sub-)pixel in old scale
        hi1 = edges1[j+1]  # upper edge of jth pixel in rebinned scale

        if hi0 < hi1:
            if er[i] > 0:
                dpix = (hi0 - lo0) / widths0[i]
                df += fl[i] * dpix
                # We don't square dpix below, since this causes an
                # artificial variation in the rebinned errors depending on
                # how the old wav bins are divided up into the rebinned
                # wav bins.
                #
                # i.e. 0.25**2 + 0.75**2 != 0.5**2 + 0.5**2 != 1**2
                de2 += er[i]**2 * dpix
                npix += dpix
            if debug:
                yh.set_height(df/npix)
                c0 = colors[i % len(colors)]
                pl.bar(lo0, fl[i], width=hi0-lo0, color=c0, alpha=0.3)
                pl.text(lo0, fl[i], 'lo0')
                pl.text(hi0, fl[i], 'hi0')
                pl.text(hi1, fl[i], 'hi1')
                raw_input('enter...')
            lo0 = hi0
            i += 1
            if i == npts0:  break
        else:
            # We have all old pixel flux values that contribute to the
            # new pixel; append the new flux value and move to the
            # next new pixel.
            if er[i] > 0:
                dpix = (hi1 - lo0) / widths0[i]
                df += fl[i] * dpix
                de2 += er[i]**2 * dpix
                npix += dpix
            if debug:
                yh.set_height(df/npix)
                c0 = colors[i % len(colors)]
                pl.bar(lo0,  fl[i], width=hi1-lo0, color=c0, alpha=0.3)
                pl.text(lo0, fl[i], 'lo0')
                pl.text(hi0, fl[i], 'hi0')
                pl.text(hi1, fl[i], 'hi1')
                raw_input('df, de2, npix: %s %s %s   enter...' %
                          (df, de2, npix))
            if npix > 0:
                # find total flux and error, then divide by number of
                # pixels (i.e. conserve flux density).
                rfl[j] = df / npix
                rer[j] = sqrt(de2) / npix
            else:
                rfl[j] = np.nan
                rer[j] = np.nan
            df = 0.
            de2 = 0.
            npix = 0.
            lo0 = hi1
            j += 1
            if j == npts1:  break
            if debug:
                x0,x1 = edges1[j:j+2]
                yh, = pl.bar(x0, 0, width=x1-x0, color='gray',
                       linestyle='dotted', alpha=0.3)
                raw_input('enter...')

    return rfl, rer

def find_trace(filename, tol=0.002, findtrace=True,
               KOAfilename='/home/nhmc/data/HIRES/koa/KOA_HIRES_20130311.fits',
               ):
    """ Given a HIRES raw filename, find a suitable exposure in the
    archive to use as a trace.
    """
    fh = fits.open(filename)
    hd = fh[0].header
    echangl = hd['ECHANGL']
    xdangl = hd['XDANGL']
    binning = hd['BINNING']
    dateobs = hd['DATE-OBS']
    if 'XDISPERS' in hd:
        xdispers = hd['XDISPERS']
     
    koa = fits.getdata(KOAfilename).view(np.recarray)

    cond = between(koa.echangl, echangl - tol, echangl + tol)
    cond &= between(koa.xdangl, xdangl - tol, xdangl + tol)
    cond &= koa.binning == binning
    if xdispers is not None:
        cond &= koa.xdispers == xdispers
    if findtrace:
        cond &= koa.elaptime < 600
        isflat = koa.imagetyp == 'flatlamp'
        # don't want a flat unless it has the pinhole decker
        cond &= ~isflat | (isflat & (koa.deckname == 'D5'))
        # don't want darks or arcs.
        isdark = (koa.imagetyp == 'dark') | (koa.imagetyp == 'dark_lamp_on')
        cond &= ~isdark
        cond &= (koa.imagetyp != 'arclamp')
         
    if not cond.sum():
        print 'None found!'
        sys.exit()
     
    koa1 = koa[cond]
    dates = np.array([int(d.replace('-', '')) for d in koa1.date_obs])
    datediff = np.abs(dates - int(dateobs.replace('-','')))
    isort = datediff.argsort()
    count = 1
    outstr = []
    retval = []
    for k in koa1[isort]:
        obj = k['object'].lower()
        if 'flat' in obj or 'dark' in obj:
            continue
        outstr.append('%-20s|%s| %s| %s' %(
            k['object'][:20], k['koaid'], k['ut'], k['deckname']))
        retval.append((k['object'].strip(), k['koaid'].strip()))
        if count == 10:
            break
        count +=1

    return retval
