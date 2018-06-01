# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: May 31, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm
import BinaryLensFunctions as blf
import BinaryLensMakePlots as blp

# Input parameters
s = 1.9		# Separation between bodies. Assume s>1. Type = string
q = 1.e-7	# mass ratio between bodies in units of larger body's mass. Type = float
pts = 20	# Number of data points on each side of the grid. Type = int
origin = 'plan'		#Coordinate frame to carry out calculations. Type = string

"""
Options:

	'geo_cent' - the geometric center frame [default if not specified]
	'star' - the star's (or larger body's) frame
	'plan' - the planet's (or smaller body's) frame
	'com' - the center-of-mass frame

"""
#method = ['numpy', 'Skowron_and_Gould_12']
solver = 'Skowron_and_Gould_12'
#solver = 'numpy'
#solver = 'zroots'		# Not precise enough to be useful on this computer

"""
Determines which method to solve polynomial with.
Options:

	'numpy'					- Uses np.roots method
	'Skowron_and_Gould_12'	- Uses Skowron & Gould 2012 method
	'zroots'				- Uses zroots laguerre method

"""

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = False
if plot_on:
	for solver in method:
		blp.plot_n_solns(s, q, origin, solver, pts)
		plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = False
if plot_on:
	blp.plot_magnification(s=s, q=q, origin=origin, solver=solver, pts=pts)
	plt.show()

blp.write_to_fits(s, q, origin = 'geo_cent', solver='numpy', pts=500)
