# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: May 25, 2018; 7:10PM

import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm
import BinaryLensFunctions as blf
import BinaryLensMakePlots as blp

def size_caustic(s, q):
	w = 4.*np.sqrt(q)*(1. + 1./(2.*(s**2))) / (s**2)
	h = 4.*np.sqrt(q)*(1. - 1./(2.*(s**2))) / (s**2)
	x = 0.5*s - 1.0/s
	return w, h, x

# Input parameters
s = 1.9		# Separation between bodies. Assume s>1. Type = string
q = 1e-7	# mass ratio between bodies in units of larger body's mass. Type = float
pts = 150	# Number of data points on each side of the grid. Type = int
origin = 'geo_cent'		#Coordinate frame to carry out calculations. Type = string

"""Options are:

	'geo_cent' - the geometric center frame [default if not specified]
	'star' - the star's (or larger body's) frame
	'plan' - the planet's (or smaller body's) frame
	'com' - the center-of-mass frame

"""

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = True
if plot_on:	
	blp.plot_n_solns(s, q, origin, pts)
	plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = True
if plot_on:
	blp.plot_magnification(s, q, origin, pts)
	plt.show()
		
