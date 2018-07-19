# Zoey Samples
# Created: May 22, 2018
# Plots.py
# Last Updated: Jul 18, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from Caustics import Caustics as caus
import MulensModel as mm
import numpy as np


# Input parameters
s = 1.5
mass_ratios = [1e-13]

solvers =  ['SG12']
origins = ['plan']

plot_frame = 'caustic'

res = int(60)
sample_res = 5
cutoff = 1.5
region = 'caustic_a'
region_lim = [-20, 20, -20, 20]
refine_region = True

SFD = True

param = []
plot = []

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
				caustic = caus(lens=p, solver='SG12')
				#caustic.plot_caustic(s=1, color='blue')
				plt.show()

			if plot_type == 'num_iamges_coeff':
				p.plot_num_images_coeff(color_magn=True, log_colorbar=True,
						save=False)

			if plot_type == 'magn_coeff':
				p.plot_magn_coeff(color_num=True, cutoff=cutoff,
						outliers=False, save=False)

			if plot_type == 'coeff':
				p.plot_coefficients(cutoff=cutoff, log_colorbar=False,
						outliers=False,	save=False)

			if plot_type == 'coeff_tstat':
				p.plot_coeff_tstat(cutoff=None, outliers=False,
						sample_res=sample_res, save=False)

			if plot_type == 'position_tstat':
				p.plot_position_tstat(cutoff=None, outliers=False,
						sample_res=sample_res, save=False)

			if plot_type == 'fits':
				p.write_to_fits()


for solver in solvers:
	for origin in origins:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'res': res, 'origin': origin,
					'solver': solver, 'plot_frame': plot_frame, 'region': region,
					'region_lim': region_lim, 'refine_region': refine_region, 'SFD': SFD}))
			plot.append(BL(**param[-1]))

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

"""
causticMM = mm.Caustics(s=s, q=plot[0].q)
causticUL = caus(lens=plot[0])
causticUL.plot_caustic(s=1, color='orange')
causticMM.plot(s=1, color='blue')
plt.show()
"""

print('\n')




