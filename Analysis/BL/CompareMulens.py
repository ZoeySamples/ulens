# Zoey Samples
# Created: Jun 12, 2018
# CompareMagn.py
# Last Updated: Jun 21, 2018

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np

"""
Notes:
This code plots the magnification of points on a grid using the given
input parameters and compares it that determined by MulensModel.
"""

# Input parameters
s = 3
q = 1e-3
res = int(150)
origins = ['plan', 'geo_cent']
solvers = ['SG12']
"""
region = 'custom'
region_lim = [-400, 120, -200, 200]
"""
region='caustic'
region_lim=None
param = []
plot = []

for (j, solver) in enumerate(solvers):
	for (k, origin) in enumerate(origins):
		param.append({'s': s, 'q': q, 'res': res, 'origin': origin,
						'solver': solver})
		plot.append(BL(**param[-1]))
		plot[-1].get_position_arrays(region=region, region_lim=region_lim)
		plot[-1].get_magnification_array()

x_array = plot[0].x_array
y_array = plot[0].y_array
mag_BL = [None]*len(plot)
for (i,p) in enumerate(plot):
	mag_BL[i] = p.magn_array

blens = []
mag_MM = []


# Assuming the total system's mass is 1 unit:
blens = (mm.BinaryLens(mass_1 = 1./(1. + q), mass_2 = q/(1. + q),
									separation = s))
"""
# Assuming the star's mass is 1 unit:
blens = (mm.BinaryLens(mass_1=1, mass_2=q, separation=s))
"""

for i in range(len(x_array)):
	mag_MM.append(blens._point_source_WM95(source_x=-x_array[i], source_y=y_array[i]))

"""
plt.scatter(x_array, y_array, c = mag_BL[0])
plt.show()
plt.scatter(x_com, y_com, c = mag_MM[0])
plt.show()
"""

# This will get rid of extremely high errors so that we can focus on the 
# common data. To look at all data including extremely high values, set
# to False.
remove_outliers = False
if remove_outliers:
	for (i, magn) in enumerate(mag_BL):
		for (j, m) in enumerate(magn):
			if m > max(mag_MM):
				mag_BL[i][j] = max(mag_MM)

# Plot the relative magnification
if False:
	for (i, p) in enumerate(plot):
		# This makes a plot of the difference in magnification from the parameters
		# at the top of the code vs. the MulensModel calculation
		print('The maximum value of magnification from MulensModel is',
			  max(mag_MM))
		plt.scatter(x_array, y_array, c = (mag_MM / mag_BL[i]),
					s=(400/res)**2, marker = 'o', cmap='plasma', lw=None)
		mag_plot = plt.colorbar()
		mag_plot.ax.tick_params(labelsize=10)
		mag_plot.set_label('Relative Magnification (MulensModel / BinaryLens)')
		plt.xlabel('X-position of source', fontsize = 12)
		plt.ylabel('Y-position of source', fontsize = 12)
		plt.xlim(min(x_array), max(x_array))
		plt.ylim(min(y_array), max(y_array))
		plt.title('Ratio of Magnification:\n{} Frame, {} Solver, q={}'.format(
								p.origin_title, p.solver_title, p.q))
		plt.gcf().set_size_inches(8, 6)
		plt.show()

# Plot the individual magnifications:
if True:
	for (i, p) in enumerate(plot):
		# This makes a plot of the difference in magnification from the parameters
		# at the top of the code vs. the MulensModel calculation
		print('The maximum value of magnification from MulensModel is',
			  max(mag_MM))
		p.plot_magnification(region=region, region_lim=region_lim, log_colorbar=False)
		caustic = mm.Caustics(s=s, q=p.q)
		caustic.plot(s=1)
		plt.show()

#	kwarg['norm'] = colors.LogNorm()
	plt.scatter(x_array, y_array, c=mag_MM, s=(400/res)**2,
				cmap='plasma', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.ax.tick_params(labelsize=10)
	mag_plot.set_label('Magnification via MulensModel')
	plt.xlabel('X-position of source', fontsize = 12)
	plt.ylabel('Y-position of source', fontsize = 12)
	plt.xlim(min(x_array), max(x_array))
	plt.ylim(min(y_array), max(y_array))
	plt.title('Ratio of Magnification:\n{} Frame, {} Solver, q={}'.format(
							p.origin_title, p.solver_title, p.q))
	plt.gcf().set_size_inches(8, 6)
	caustic = mm.Caustics(s=s, q=p.q)
	caustic.plot(s=1)
	plt.show()


