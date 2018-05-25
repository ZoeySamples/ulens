# Zoey Samples
# Created: May 22, 2018
# BinaryLensAnalysis.py
# Last Updated: May 24, 2018; 6:42PM

import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm
import BinaryLensFunctions as blf

tests = [
	[0, 0, 1, 1], 
	[1.3219, -0.0771, 1.35, .00578],
	[1.0799, 0.0985, 1.1, 0.99],
	[1.2489, 0.0209, 0.9357, 0.99]
		]	# Input parameters for 4 trials

for test in tests:
	"""Prints input paramters, image locations, and magnification"""
	print("Input", tests.index(test) + 1, ":\nx = {:}\ny = {:}\ns = {:}\nq = {:}\n".format(*test))
	solutions = blf.solution(*test)
	print("Image locations:")
	dm, m, zeta, z1 = blf.assign(*test)
	for z in solutions:
		if blf.check_solution(dm, m, zeta, z1, z) == True:
			print("{:.5f}".format(z))
	magn = blf.magnification(*test)
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
	if plot_on == False:
		break
	a = 2.0
	d = 0.5*test[2] - 1.0/test[2]
	x = np.linspace(d-0.3, d+0.3, 7000)
	y = a*(x-d)+0.02
	magn_array = x*0
	for j in range(len(x)):
		magn_array[j] = blf.magnification(x[j], y[j], test[2], test[3])
	plt.plot(x, magn_array, '.b')
	plt.xlabel("X")
	plt.ylabel("Magnification")
	plt.show()
