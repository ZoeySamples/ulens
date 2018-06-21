# Zoey Samples
# Created: June 13, 2018
# PlotTstat.py
# Last Updated: June 21, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm


# Input parameters
s = 1.5
mass_ratios = [1e-7]
solvers =  ['zroots']
origins = ['plan']
tolerance = 0.00007
sample_res = 5
res = 50
region = 'caustic'
region_lim = (-1.2, 0, 0, 1.2)
param = []
lens = []
coefficients = []

for solver in solvers:
	for origin in origins:
		for q in mass_ratios:
			param.append(({'s': s, 'q': q, 'origin': origin, 'solver': solver,
								 'res': res, 'tolerance': tolerance}))
			lens.append(BL(**param[-1]))

for lens in lens:
	if False:
		lens.plot_coeff_tstat(cutoff=None, outliers=False, region=region,
				region_lim=region_lim, sample_res=sample_res, save=False)
	if True:
		lens.plot_position_tstat(cutoff=None, outliers=False, region=region,
				region_lim=region_lim, sample_res=sample_res, save=False)
	if True:
		lens.plot_magnification(log_colorbar=True)
		caustic = mm.Caustics(s=s, q=lens.q)
		caustic.plot(s=1)
		plt.show()

