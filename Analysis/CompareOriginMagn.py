# Zoey Samples
# Created: Jun 15, 2018
# CompareOriginMagn.py
# Last Updated: Jun 15, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np

# Input parameters
s = 1.5
mass_ratios = [1e-7]
res = int(150)
coordinates = ['geo_cent', 'caustic', 'plan']
method = ['SG12']
tolerance = 0.005
region = 'offax_cusp'
cutoff = 4

param = [[[None] * len(coordinates) for i in range(len(method))] 
										for j in range(len(mass_ratios))]
plot = [[[None] * len(coordinates) for i in range(len(method))] 
										for j in range(len(mass_ratios))]

print('Getting data...')
for (i, q) in enumerate(mass_ratios):
	for (j, solver) in enumerate(method):
		for (k, origin) in enumerate(coordinates):
			param[i][j][k] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
							'solver': solver, 'tolerance': tolerance})
			plot[i][j][k] = (BL(**param[i][j][k]))
			print('Data retrieved for: {}; {}.'.format(plot[i][j][k].solver_phrase,
													 plot[i][j][k].origin_phrase))

for (i, q) in enumerate(mass_ratios):
	for (j, solver) in enumerate(method):
		for (k, origin) in enumerate(coordinates):
			for l in range(k+1, len(coordinates)):
				plot[i][j][k].plot_outlier_coeff(region = region, other_BL =
												plot[i][j][l], cutoff = cutoff)
				import sys
				sys.exit()
#				plot[i][j][k].plot_rel_magnification(other_BL = plot[i][j][l], 
#							region = region, outliers = True, hl_out = False, cutoff =
#							cutoff, log_colorbar = True, save=False)
#				plt.show()
