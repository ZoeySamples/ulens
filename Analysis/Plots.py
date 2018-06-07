# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 04, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

# Input parameters
s = 1.9
q = 1.e-7
res = int(50)
origin = 'geo_cent'
solver = 'SG12'

param = []
plot = []

loop_on = True
if loop_on:
	method =  ['numpy', 'SG12']
	coordinates = ['geo_cent', 'plan']
	for (i, solver) in enumerate(method):
		for (j, origin) in enumerate(coordinates):
			param.append({'s': s, 'q': q, 'res': res, 'origin': origin,
							'solver': solver, 'tolerance': 0.0001})
			plot.append(BL(**param[-1]))
else:
	param.append({'s': s, 'q': q, 'res': res, 'origin': origin, 'solver': solver, 
					'tolerance': 0.0001})
	plot.append(BL(**param[-1]))

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = True
if plot_on:
	for (i, p) in enumerate(plot):
		p.plot_n_solns(save=False, print_errors=True)
		plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = True
if plot_on:
	for (i, p) in enumerate(plot):
		p.plot_magnification()
		caustics = mm.Caustics(s=s, q=q)
		caustics.plot(s=1)
		plt.show()

# Writes data to fits table
plot_on = False
if plot_on:
	for p in plot:
		p.write_to_fits()
