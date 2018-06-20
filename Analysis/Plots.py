# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 19, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

# Input parameters
s = 1.5
mass_ratios = [0.5*1e-7]
res = int(100)
method =  ['SG12']
coordinates = ['caustic']
tolerance = 0.00007
cutoff = 1.5
region = 'custom'
region_lim = (-.15, .1, 0.9, 1.3)
coeff_multiplier = None
param = []
plot = []

print('coeff_multiplier =', coeff_multiplier)

def make_plot():

	for p in plot:
		for plot_type in plot_types:

			caustic = mm.Caustics(s=s, q=p.q)

			if plot_type == 'num_images':
				p.plot_num_images(errors_only=False, region=region,
						region_lim=region_lim, save=False, print_errors=True)
				caustic.plot(s=1)
				plt.show()

			if plot_type == 'magn':
				p.plot_magnification(outliers=True, region=region,
						region_lim=region_lim, log_colorbar=True, cutoff=cutoff,
						save=False)
				caustic.plot(s=1)
				plt.show()

			if plot_type == 'num_iamges_coeff':
				p.plot_num_images_coeff(color_magn=True, log_colorbar=True,
						region='caustic', region_lim=None, save=False)

			if plot_type == 'magn_coeff':
				p.plot_magn_coeff(color_num=True, cutoff=cutoff,
						outliers=False, region=region, region_lim=region_lim,
						save=False)

			if plot_type == 'coeff':
				p.plot_coefficients(cutoff=cutoff, log_colorbar=False,
						outliers=False, region=region, region_lim=region_lim,
						save=False)

			if plot_type == 'fits':
				p.write_to_fits()


for solver in method:
	for origin in coordinates:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'res': res, 'origin': origin,
							'solver': solver, 'tolerance': tolerance,
							'coeff_multiplier': coeff_multiplier}))
			plot.append(BL(**param[-1]))

plot_types = []
#plot_types.append('num_images')
#plot_types.append('magn')
#plot_types.append('num_iamges_coeff')
plot_types.append('magn_coeff')
#plot_types.append('coeff')
#plot_types.append('fits')

make_plot()
print('\n')
"""
p = plot[0]
p.get_position_arrays(region=region, region_lim=region_lim)
p.get_coeff_array()
coeff = p.coeff_array

import numpy as np
mean_coeff = [[]]*len(coeff)
for (i, c) in enumerate(coeff):
	mean_coeff[i] = np.mean(c)

print(mean_coeff)

"""












