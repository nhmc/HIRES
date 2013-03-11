# rebin the spectra and make a template. This creates template.npy
# and one file in rebinned/ for each order in each spectrum.
./rebin_HIRES.py reduced/Flux*fits
# scale the rebinned spectra to match the template. This reads
# template.npy, all the files rebinned/*.npy and creates plots in
# scaling_plots/ and multipliers.npy.
./scale_HIRES.py
# apply the scaling and combine the rebinned spectra. This reads
# multipliers.npy and all the files rebinned/*.npy, and then
# overwrites rebinned/*.npy, adding information about which pixels
# were masked and clipped during the combination. The combined
# spectrum is written to combined.fits. You will probably want to edit
# this file and change which regions are masked.
./combine_HIRES.py
# plot the individual rebinned exposures and the combined spectrum
./plot_combined_HIRES.py
# fit a continuum to combined.fits and write to final.fits. 3.63 is
# the QSO redshift.
./fit_continuum 3.63 combined.fits final.fits
