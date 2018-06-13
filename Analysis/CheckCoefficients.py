# Zoey Samples
# Created: June 13, 2018
# CheckCoefficients.py
# Last Updated: June 13, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm

"""
coeff (list):
	coeff[i][j][k] represents a single coefficient value.
		i corresponds to a lens object given by any combination of the input
			parameters.
		j corresponds to a single point in a lens; each holds 6 coefficients.
		k corresponds to the coefficient of degree (5 - k) in the polynomial;
			i.e. k=0 corresponds to the coefficient for z**5.
	


"""
# Input parameters
s = 1.5
mass_ratios = [1e-7]
method =  ['SG12']
coordinates = ['geo_cent']
tolerance = 0.0001
sample_res = 10
region_res = 20
region = 'cusp'
param = []
lens = []
coefficients = []

for solver in method:
	for origin in coordinates:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'origin': origin,
							'solver': solver, 'tolerance': tolerance}))
			lens.append(BL(**param[-1]))

for lens in lens:
	lens.plot_coeff_tstat(region = region, region_res = region_res, sample_res = sample_res)





