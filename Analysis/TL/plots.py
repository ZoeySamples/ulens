# Zoey Samples
# Created: May 22, 2018
# TLPlots.py
# Last Updated: Jul 09, 2018

import matplotlib.pyplot as plt
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
import MulensModel as mm
import numpy as np


# Input parameters
system = 'SPM'
plot_frame = 'geo_cent'
s1 = 0.8
s2 = 0.6
q1=1e-4
q2=1e-1
solvers =  ['SG12']
origins = ['body3']
phi = 0

res = int(120)
sample_res = 5
region = 'caustic2b'
region_lim = (5, 10, 0, 5)

cutoff = 1.5
SFD = True

param = []
plot = []

def make_plot():

	for p in plot:
		for plot_type in plot_types:

			caustic = caus(lens=p)

			if plot_type == 'num_images':
				p.plot_num_images(errors_only=False, save=False, print_errors=True)
				caustic.plot_caustic(s=0.2, color='yellow')
				plt.show()

			if plot_type == 'magn':
				p.plot_magnification(outliers=False, log_colorbar=True, cutoff=cutoff,
						save=False)
#				caustic.plot_caustic(s=0.5, color='red')
				plt.show()

			if plot_type == 'image_pos':
				p.plot_image_positions()
				caustic = caus(lens=plot[0])
				caustic.calculate()
				(x, y) = (caustic.critical_curve.x, caustic.critical_curve.y)
				plt.scatter(x, y, s=1, color='red')
				plt.show()

			if plot_type == 'num_images_coeff':
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
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame})
		plot.append(TL(**param[-1]))

plot_types = []
plot_types.append('num_images')
#plot_types.append('magn')
#plot_types.append('image_pos')
#plot_types.append('num_iamges_coeff')
#plot_types.append('magn_coeff')
#plot_types.append('coeff')
#plot_types.append('fits')
#plot_types.append('coeff_tstat')
#plot_types.append('position_tstat')

make_plot()
print('\n')
