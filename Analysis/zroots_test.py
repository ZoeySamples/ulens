# Zoey Samples
# Created: Jun 08, 2018
# BinaryLensPlots.py
# Last Updated: Jun 15, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

# Input parameters
s = 1.8
q = 1.e-8
res = int(100)
origin = 'plan'
solver = 'SG12'
tolerance = 0.00007

"""

Note: Issues resolved when set polish to 0 and re-compiled.


Test runs:
	Works properly for me when I use:

		#Test 1
		s = 1.2
		q = 1.e-5
		res = int(30)
		origin = 'plan'
		solver = 'zroots'
		tolerance = 0.0001

		#Test 2
		s = 1.2
		q = 1.e-3
		res = int(100)
		origin = 'geo_cent'
		solver = 'zroots'
		tolerance = 0.01

		#Test 3
		s = 1.2
		q = 1.e-7
		res = int(15)
		origin = 'plan'
		solver = 'zroots'
		tolerance = 0.0001


	Runs but gives horrible result when I use:

		origin = 'com', 'star', or 'geo_cent' AND
		tolerance < 0.005 OR tolerance > 0.05;
		OR	when I use q < 1e-6


	Gives run-time error when I use:

		#Test 1*
		s = 1.9
		q = 1.e-2
		res = int(100)
		origin = 'plan'
		solver = 'zroots'

		#Test 2*
		s = 1.9
		q = 1.e-7
		res = int(50)
		origin = 'plan'
		solver = 'zroots'

		#Test 3*
		s = 1.2
		q = 1.e-3
		res = int(20)
		origin = 'plan'
		solver = 'zroots'

			*tolerance = (any value)
"""

param = {'s': s, 'q': q, 'res': res, 'origin': origin, 'solver': solver,
				'tolerance': tolerance}
plot = BL(**param)

# Plots the number of solutions in a grid of points centered on the caustic
plot_on = True
if plot_on:
	plot.plot_num_images(save=False, print_errors=True)
	plt.show()

# Plots magnification in a grid of points centered on the caustic
plot_on = False
if plot_on:
	plot.plot_magnification(log_colorbar=True)
	caustics = mm.Caustics(s=s, q=q)
	caustics.plot(s=2)
	plt.show()
