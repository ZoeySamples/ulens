# Zoey Samples
# Created: Jun 15, 2018
# CompareOriginMagn.py
# Last Updated: Jun 19, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np

# Input parameters
s = 1.5
mass_ratios = [1e-7]
res = int(100)
coordinates = ['plan', 'geo_cent']
method = ['SG12']
tolerance = 0.0007
region = 'custom'
region_lim = (-.3, .3, .5, 1.5)
ratio_cutoff = 4

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
			print('Data retrieved for: {}; {}.'.format(
					plot[i][j][k].solver_phrase, plot[i][j][k].origin_phrase))

# Plot relative magnification value vs coefficient value.
plot_on = True
if plot_on:
	for (i, q) in enumerate(mass_ratios):
		for (j, solver) in enumerate(method):
			for (k, origin) in enumerate(coordinates):
				for l in range(k+1, len(coordinates)):
					plot[i][j][k].plot_rel_magn_coeff(outliers=True, 
							region=region, region_lim=region_lim,
							other_BL=plot[i][j][l], ratio_cutoff=ratio_cutoff,
							save=False)

#Plot relative magnification vs position.
plot_on = False
if plot_on:
	for (i, q) in enumerate(mass_ratios):
		for (j, solver) in enumerate(method):
			for (k, origin) in enumerate(coordinates):
				for l in range(k+1, len(coordinates)):
					plot[i][j][k].plot_rel_magnification(
							other_BL=plot[i][j][l], region=region,
							region_lim=region_lim, outliers=True,
							ratio_cutoff=ratio_cutoff, log_colorbar=True,
							save=False)
					plt.show()
