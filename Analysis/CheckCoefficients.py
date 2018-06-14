# Zoey Samples
# Created: June 13, 2018
# CheckCoefficients.py
# Last Updated: June 13, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm


# Input parameters
s = 1.5
mass_ratios = [1e-7]
method =  ['SG12']
coordinates = ['geo_cent']
tolerance = 0.0001
sample_res = 4
res = 150
region = 'cusp'
param = []
lens = []
coefficients = []

for solver in method:
	for origin in coordinates:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'origin': origin, 'solver': solver,
								 'res': res, 'tolerance': tolerance}))
			lens.append(BL(**param[-1]))

for lens in lens:
	lens.plot_coeff_tstat(region = region, plot_position_tstat = True, cutoff = 150,
					plot_coeffs = False, sample_res = sample_res, outliers = False)
#	lens.plot_magnification(log_colorbar = True)
#	caustic = mm.Caustics(s=s, q=lens.q)
#	caustic.plot(s=1)
#	plt.show()

