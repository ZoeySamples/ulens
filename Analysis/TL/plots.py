# Zoey Samples
# Created: May 22, 2018
# TLPlots.py
# Last Updated: Jun 26, 2018

import matplotlib.pyplot as plt
from TripleLens import TripleLens as TL
import MulensModel as mm
import numpy as np

# Input parameters
s1 = 1.5
s2 = 0.8
q1 = 1e-8
q2 = 5e-2
solvers =  ['SG12']
origins = ['body3']
phi = 130

res = int(60)
sample_res = 5
region = 'caustic'
region_lim = (0, 1.2, 0, 1.2)

cutoff = 1.5
specific_frame_derivation = True

param = []
plot = []

def make_plot():

	for p in plot:
		for plot_type in plot_types:

#			caustic = mm.Caustics(s=sPS, q=p.qPS)

			if plot_type == 'num_images':
				p.plot_num_images(errors_only=False, region=region,
						region_lim=region_lim, save=False, print_errors=True)
#				caustic.plot(s=1)
				plt.show()

			if plot_type == 'magn':
				p.plot_magnification(outliers=False, region=region,
						region_lim=region_lim, log_colorbar=False, cutoff=cutoff,
						save=False)
#				caustic.plot(s=1)
				plt.show()

			if plot_type == 'num_iamges_coeff':
				p.plot_num_images_coeff(color_magn=TPrue, log_colorbar=True,
						region='caustic', region_lim=None, save=False)

			if plot_type == 'magn_coeff':
				p.plot_magn_coeff(color_num=True, cutoff=cutoff,
						outliers=False, region=region, region_lim=region_lim,
						save=False)

			if plot_type == 'coeff':
				p.plot_coefficients(cutoff=cutoff, log_colorbar=False,
						outliers=False, region=region, region_lim=region_lim,
						save=False)

			if plot_type == 'coeff_tstat':
				p.plot_coeff_tstat(cutoff=None, outliers=False, region=region,
						region_lim=region_lim, sample_res=sample_res, save=False)

			if plot_type == 'position_tstat':
				p.plot_position_tstat(cutoff=None, outliers=False, region=region,
						region_lim=region_lim, sample_res=sample_res, save=False)

			if plot_type == 'fits':
				p.write_to_fits()


for solver in solvers:
	for origin in origins:
		param.append({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'solver': solver, 'specific_frame_derivation': 
					specific_frame_derivation})
		plot.append(TL(**param[-1]))

plot_types = []
plot_types.append('num_images')
#plot_types.append('magn')
#plot_types.append('num_iamges_coeff')
#plot_types.append('magn_coeff')
#plot_types.append('coeff')
#plot_types.append('fits')
#plot_types.append('coeff_tstat')
#plot_types.append('position_tstat')

make_plot()
print('\n')

