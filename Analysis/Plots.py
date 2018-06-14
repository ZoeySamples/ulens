# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 13, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

# Input parameters
s = 1.9
mass_ratios = [5e-8]
res = int(100)
method =  ['SG12']
coordinates = ['geo_cent']
tolerance = 0.00007
region = 'caustic'
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
		p.plot_n_solns(region = region, save=False, print_errors=True)
		caustics = mm.Caustics(s=s, q=p.q)
		caustics.plot(s=2)
		plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = True
if plot_on:
	for p in plot:
		p.plot_magnification(outliers = True, region = region, log_colorbar = False,
									cutoff = None)
		caustics = mm.Caustics(s=s, q=p.q)
		caustics.plot(s=1)
		plt.show()

# Writes data to fits table
plot_on = False
if plot_on:
	for p in plot:
		p.write_to_fits()
