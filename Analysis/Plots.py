# Zoey Samples
# Created: May 22, 2018
# BinaryLensPlots.py
# Last Updated: Jun 26, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np

# Input parameters
s = 1.5
mass_ratios = [1e-15, 5e-16]
tolerance = 2e-8
solvers =  ['SG12']
origins = ['plan']

res = int(200)
region = 'caustic'
region_lim = (-.15, .15, 0.8, 1.3)

cutoff = 1.5
coeff_multiplier = None
specific_frame_derivation = True
#tolerance = ((specific_frame_derivation)*2e-8 +
#			(not specific_frame_derivation)*5e-5)

param = []
plot = []

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
				p.plot_magnification(outliers=False, region=region,
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


for solver in solvers:
	for origin in origins:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'res': res, 'origin': origin,
					'solver': solver, 'tolerance': tolerance,
					'specific_frame_derivation': specific_frame_derivation}))
			plot.append(BL(**param[-1]))


plot_types = []
plot_types.append('num_images')
#plot_types.append('magn')
#plot_types.append('num_iamges_coeff')
#plot_types.append('magn_coeff')
#plot_types.append('coeff')
#plot_types.append('fits')

make_plot()
print('\n')

