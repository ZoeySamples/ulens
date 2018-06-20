# Zoey Samples
# Created: Jun 15, 2018
# CompareOriginMagn.py
# Last Updated: Jun 19, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np

def initialize():

	for (i, q) in enumerate(mass_ratios):
		for (j, const) in enumerate(const_param):
			for (k, comp) in enumerate(comp_param):
				if comp_param == origins:
					origin = comp
					solver = const
				elif comp_param == solvers:
					origin = const
					solver = comp	
				param[i][j][k] = ({'s': s, 'q': q, 'res': res, 'origin':
								  origin, 'solver': solver, 'tolerance':
								  tolerance})
				lens[i][j][k] = (BL(**param[i][j][k]))

def prepare_plot(plot_type, **plot_args):

	for (i, q) in enumerate(mass_ratios):
		for (j, const) in enumerate(const_param):
			for (k, comp) in enumerate(comp_param):
				for l in range(k+1, len(comp_param)):
					lens1 = lens[i][j][k]
					lens2 = lens[i][j][l]
					make_plot(lens1, lens2, plot_type, **plot_args)

def make_plot(lens1, lens2, plot_type, **plot_args):

	if plot_type == 'rel_magn_coeff':
		lens1.plot_rel_magn_coeff(other_BL=lens2, **plot_args)
	elif plot_type == 'rel_magn':
		lens1.plot_rel_magnification(other_BL=lens2, **plot_args)
		plt.show()


# Input parameters
s = 1.5
mass_ratios = [1e-7]
res = int(40)
origins = ['plan', 'geo_cent', 'caustic']
solvers = ['SG12']
tolerance = 0.0007
region = 'custom'
region_lim = (-.3, .3, .5, 1.5)
ratio_cutoff = 4
compare = 'origins'

if compare == 'origins':
	const_param = solvers
	comp_param = origins
elif compare == 'solvers':
	comp_param = solvers
	const_param = origins
else:
	raise ValueError('Unknown value given for string variable, compare.')

param = [[[None] * len(comp_param) for i in range(len(const_param))]
		for j in range(len(mass_ratios))]
lens = [[[None] * len(comp_param) for i in range(len(const_param))]
		for j in range(len(mass_ratios))]

initialize()

if False:
	prepare_plot(plot_type = 'rel_magn', outliers=False,
			ratio_cutoff=ratio_cutoff, region=region, region_lim=region_lim,
			save=False, log_colorbar=True)

if False:
	prepare_plot(plot_type = 'rel_magn_coeff', outliers=True,
			ratio_cutoff=ratio_cutoff, region=region, region_lim=region_lim,
			save=False)

