# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 04, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

# Input parameters
s = 1.9
q = 5.e-8
res = int(1500)
origin = 'plan'
solver = 'zroots'

param = []
plot = []

loop_on = True
if loop_on:
	method =  ['numpy', 'SG12']
	coordinates = ['geo_cent', 'caustic', 'plan']
	mass_ratios = [1e-7, 5e-8, 2e-8, 1e-8]
	for (i, solver) in enumerate(method):
		for (j, origin) in enumerate(coordinates):
			for (k, q) in enumerate(mass_ratios):
				param.append({'s': s, 'q': q, 'res': res, 'origin': origin,
								'solver': solver, 'tolerance': 0.0001})
				plot.append(BL(**param[-1]))
else:
	param.append({'s': s, 'q': q, 'res': res, 'origin': origin, 'solver': solver, 
					'tolerance': 0.0001})
	plot.append(BL(**param[-1]))

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = False
if plot_on:
	for (i, p) in enumerate(plot):
		p.print_input()
		p.plot_n_solns(save=False, print_errors=False)
		plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = False
if plot_on:
	for (i, p) in enumerate(plot):
		p.plot_magnification()
		caustics = mm.Caustics(s=s, q=q)
		caustics.plot(s=1)
		plt.show()

# Writes data to fits table
plot_on = True
if plot_on:
	for p in plot:
		p.write_to_fits()
