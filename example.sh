# rebin the spectra and make a template. This creates template.npy
# and one file in rebinned/ for each order in each spectrum.
HIRES_rebin reduced/Flux*fits
# make a template, used to scale the orders. This also plots all the
# rebinned spectra.
HIRES_make_template
# scale the rebinned spectra to match the template. This reads
# template.npy, all the files rebinned/*.npy and creates plots in
# scaling_plots/ and multipliers.npy.
HIRES_scale
# apply the scaling and combine the rebinned spectra. This reads
# multipliers.npy and all the files rebinned/*.npy, and then
# overwrites rebinned/*.npy, adding information about which pixels
# were masked and clipped during the combination. The combined
# spectrum is written to combined.fits. You will probably want to edit
# this file and change which regions are masked. WARNING: this changes
# the rebinned spectra in rebinned/!
HIRES_combine
# plot the scaled, rebinned exposures and the combined spectrum
HIRES_plot_combined
# fit a continuum to combined.fits and write to final.fits. 3.63 is
# the QSO redshift.
HIRES_fitcontinuum 3.63 combined.fits final.fits
