# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 19, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

# Input parameters
s = 1.5
mass_ratios = [1e-7]
res = int(50)
method =  ['SG12']
coordinates = ['geo_cent']
tolerance = 0.00006
cutoff = 1.5
region = 'custom'
region_lim = (-.3, .1, 1.0, 1.3)
param = []
plot = []

for solver in method:
	for origin in coordinates:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'res': res, 'origin': origin,
							'solver': solver, 'tolerance': tolerance}))
			plot.append(BL(**param[-1]))

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = False
if plot_on:
	for p in plot:
		p.plot_n_solns(errors_only=False, region=region, region_lim=region_lim,
					   save=False, print_errors=True, s=3)
		caustics = mm.Caustics(s=s, q=p.q)
		caustics.plot(s=1)
		plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = False
if plot_on:
	for p in plot:
		p.plot_magnification(outliers=False, region=region,
				region_lim=region_lim, log_colorbar=True, cutoff=cutoff,
				save=False)
		caustics = mm.Caustics(s=s, q=p.q)
		caustics.plot(s=1)
		plt.show()

plot_on = True
if plot_on:
	for p in plot:
		p.plot_magn_coeff(cutoff=cutoff, outliers=False, region=region,
						region_lim=region_lim, save=False)

# Writes data to fits table
plot_on = False
if plot_on:
	for p in plot:
		p.write_to_fits()
