# Zoey Samples
# Created: May 22, 2018
# TLPlots.py
# Last Updated: Jul 18, 2018

import matplotlib.pyplot as plt
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
import MulensModel as mm
import numpy as np


def make_plot():

	for p in plot:
		for plot_type in plot_types:

			if plot_type == 'num_images':
				p.plot_num_images(errors_only=False, save=False, print_errors=True)
				caustic = caus(lens=p, solver='SG12')
				caustic.plot_caustic(s=1, color='yellow', lw=0)
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

# Input parameters
test = 'test2'

if test == 'test1':
	q1 = 1e-2
	q2 = 1e-1
	s1 = 1.5
	s2 = 0.9
	phi = 40

elif test == 'test2':
	q1 = 1e-2
	q2 = 1e-1
	s1 = 1.0
	s2 = 1.3
	phi = 0

elif test == 'test3':
	q1 = 1e-1
	q2 = 1e-2
	s1 = 0.7
	s2 = 1.5
	phi = 0

system = 'SPM'
solvers =  ['SG12']
origins = ['body3']

plot_frame = 'caustic'

res = int(30)
sample_res = 5
cutoff = 1.5
region = 'custom_3a'
region_lim = (-10, 10, -10, 10)
region_lim = np.array(region_lim)*3
refine_region = False

SFD = True

param = []
plot = []

for solver in solvers:
	for origin in origins:
		param.append({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame, 'refine_region': refine_region})
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
