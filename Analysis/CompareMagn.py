# Zoey Samples
# Created: Jun 12, 2018
# CompareMagn.py
# Last Updated: Jun 12, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np

"""
Notes:
This code plots the magnification of points on a grid using the zroots solver,
and the difference between magnification, comparing it to MulensModel. There is
no discernable caustic shape for the zroots method when q < 4e-4.

There are runtime errors when s approaches 1.0 or 2.0, when the coefficients in
the BinaryLens object are not reduced for the given origin.
"""


# Input parameters
s = 1.1
q = 1e-3
res = int(100)
coordinates = ['geo_cent']
method = ['zroots']
tolerance = 0.005		# The optimal tolerance is around 0.005 for zroots

param = []
plot = []

for solver in method:
	for origin in coordinates:
		param.append({'s': s, 'q': q, 'res': res, 'origin': origin,
						'solver': solver, 'tolerance': tolerance})
		plot.append(BL(**param[-1]))
		plot[-1].grid_plots()

x_array = plot[0].x_array
y_array = plot[0].y_array
x_com = x_array + s*(1.-q)/(2*(1.+q))
y_com = y_array
mag_BL = []
for p in plot:
	mag_BL.append(p.mag_1d)

blens = []
mag_MM = []
blens = (mm.BinaryLens(mass_1 = 1./(1. + q), mass_2 = q/(1. + q),
									separation = s))

for i in range(len(x_com)):
	mag_MM.append(blens.point_source_magnification(x_com[i], y_com[i]))

"""
plt.scatter(x_array, y_array, c = mag_BL[0])
plt.show()
plt.scatter(x_com, y_com, c = mag_MM[0])
plt.show()
"""

# This will get rid of extremely high errors so that we can focus on the common data.
# To look at all data including extremely high values, set to False.
remove_outliers = True
if remove_outliers:
	for (i, magn) in enumerate(mag_BL):
		for (j, m) in enumerate(magn):
			if m > max(mag_MM):
				mag_BL[i][j] = max(mag_MM)


for (i, p) in enumerate(plot):
	# This makes a plot of the difference in magnification from the parameters
	# at the top of the code vs. the MulensModel calculation
	print('The maximum value of magnification from MulensModel is', max(mag_MM))
	plt.scatter(x_array, y_array, c = (mag_MM - mag_BL[i])/mag_MM, s=15, marker = 'o',
								cmap='plasma', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.ax.tick_params(labelsize=10)
	mag_plot.set_label('Magnification (MulensModel - zroots)')
	plt.xlabel('X-position of source', fontsize = 12)
	plt.ylabel('Y-position of source', fontsize = 12)
	plt.xlim(min(x_array), max(x_array))
	plt.ylim(min(y_array), max(y_array))
	plt.title('Difference in Magnification:\n{} Frame, {} Solver, q={}'.format(
							p.origin_title, p.solver_title, p.q))
	plt.gcf().set_size_inches(8, 6)
	plt.show()

	# This makes a plot of the magnification from the parameters at the top of the code
	plt.scatter(x_array, y_array, c = mag_BL[i], s=15, marker = 'o',
								cmap='plasma', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.ax.tick_params(labelsize=10)
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source', fontsize = 12)
	plt.ylabel('Y-position of source', fontsize = 12)
	plt.xlim(min(x_array), max(x_array))
	plt.ylim(min(y_array), max(y_array))
	plt.title('Magnification:\n{} Frame, {} Solver, q={}'.format(
							p.origin_title, p.solver_title, p.q))
	plt.gcf().set_size_inches(8, 6)
	plt.show()

