# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 04, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
import Functions as blf
from MakePlots import Plots
import MulensModel as mm

# Input parameters
s = 1.9		# Separation between bodies. Assume s>1. Type = string
q = 1.e-7	# mass ratio between bodies in units of larger body's mass. Type = float
pts = int(2000)	# Number of data points on each side of the grid. Type = int
origin = 'geo_cent'		#Coordinate frame to carry out calculations. Type = string
solver = 'Skowron_and_Gould_12'	# Determines which method to solve polynomial with.

"""
origin options:

	'geo_cent' - the geometric center frame [default if not specified]
	'star' - the star's (or larger body's) frame
	'plan' - the planet's (or smaller body's) frame
	'com' - the center-of-mass frame


solver options:

	'numpy'					- Uses np.roots method
	'Skowron_and_Gould_12'	- Uses Skowron & Gould 2012 method
	'zroots'				- Uses zroots laguerre method (not precise enough
								to be useful on this computer

"""

params = []
plot = []

loop_on = True
if loop_on:
	method =  ['numpy', 'Skowron_and_Gould_12']
	coordinates = ['geo_cent', 'plan']
	for (i, solver) in enumerate(method):
		for (j, origin) in enumerate(coordinates):
			params.append({'s': s, 'q': q, 'pts': pts, 'origin': origin,
							'solver': solver})
			plot.append(Plots(**params[-1]))
else:
	params.append({'s': s, 'q': q, 'pts': pts, 'origin': origin, 'solver': solver})
	plot.append(Plots(**params[-1]))

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = False
if plot_on:
	for (i, p) in enumerate(plot):
		print((params[i]))
		p.plot_n_solns()
		plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = True
if plot_on:
	for (i, p) in enumerate(plot):
		print((params[i]))
		p.plot_magnification()
		caustics = mm.Caustics(s=s, q=q)
		caustics.plot(s=1)
		#plt.show()

# Writes data to fits table
plot_on = True
if plot_on:
	for p in plot:
		p.write_to_fits()
