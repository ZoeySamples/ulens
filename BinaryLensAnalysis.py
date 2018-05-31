# Zoey Samples
# Created: May 22, 2018
# BinaryLensAnalysis.py
# Last Updated: May 31, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm
import BinaryLensFunctions as blf

#solver = 'Skowron_and_Gould_12'
#solver = 'numpy'
solver = 'zroots'

"""
Determines whether to solve polynomial with np.roots or Skowron & Gould 2012.
Options:

	'numpy'					- Uses np.roots method
	'Skowron_and_Gould_12'	- Uses Skowron & Gould 2012 method
	'zroots'				- Uses zroots laguerre method

"""

origin = 'geo_cent'

"""
Coordinate frame to carry out calculations. Options are:

	'geo_cent' - the geometric center frame [default if not specified]
	'star' - the star's (or larger body's) frame
	'plan' - the planet's (or smaller body's) frame
	'com' - the center-of-mass frame

"""

(origin, print_str) = blf.print_frame(origin)

tests = [
	[0., 0., 1., 1.], 
	[1.3219, -0.0771, 1.35, .00578],
	[1.0799, 0.0985, 1.1, 0.99],
	[1.2489, 0.0209, 0.9357, 0.99]		
		]
"""
Input parameters for 4 trials; defaults to geometric center frame,
while calculations are carried out in frame determined by 'origin'
"""

for test in tests:
	"""Prints input paramters, image locations, and magnification"""
	print("Input", tests.index(test) + 1, ":\nx = {:}\ny = {:}\ns = {:}\nq = {:}\n".format(*test))
	solutions = blf.solution(*test, origin, solver)
	print("Image locations:")
	(dm, m, zeta, z1, z2) = blf.assign(*test, origin)
	for z in solutions:
		if blf.check_solution(dm, m, zeta, z1, z2, z, origin) == True:
			print("{:.5f}".format(z))
	magn = blf.magnification(*test, origin, solutions)
	print("\nThe magnification is: {:.5f}".format(magn))
	print("_"*20 + "\n")

plot_on = False
for test in tests:
	"""Prints each caustic"""
	if plot_on == False:
		break
	caustic = mm.Caustics(q=test[3],s=test[2])
	caustic.plot()
	plt.show()

plot_on = False
# Trajectory Plots: Not very useful
for test in tests:
	"""Makes a plot of magnification vs position for each scenario"""
	if not plot_on:
		break
	a = 2.0
	d = 0.5*test[2] - 1.0/test[2]
	x = np.linspace(d-0.3, d+0.3, 7000)
	y = a*(x-d)+0.02
	magn_array = x*0
	for j in range(len(x)):
		magn_array[j] = blf.magnification(x[j], y[j], test[2], test[3], origin)
	plt.plot(x, magn_array, '.b')
	plt.xlabel("X")
	plt.ylabel("Magnification")
	plt.show()

print('Calculations carried out in the', print_str,'\n')
