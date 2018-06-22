# Zoey Samples
# Created: Jun 21, 2018
# SolverInfo.py
# Last Updated: Jun 22, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
from pathlib import Path

"""
This file highlights some observations I have made regarding the
coordinate frames.

Finding 1: The least error-prone frame is the one centered on the planet,
or the smallest body. When using the general form of the polynomial
equation, the planet frame performs better than the geometric center
frame and the caustic frame by a factor of a few. This is especially
noticeable when we approach mass ratios of 1e-7 and 1e-8, whereafter
we are given a run-time warning and very sloppy data.

However, when we use the specific form of the polynomial, derived for the
planet frame, we are allowed to model systems with mass ratios all the way
down to 1e-15. One strange observation is that, while using the specifically-
derived form of the polynomial frame is much better for the planet frame, it
is actually worse for the geometric center frame.

Here is a demonstration:
"""

def demonstration():

	param = [[[None] * len(origins) for i in range(len(solvers))]
					for j in range(len(mass_ratios))]
	plot = [[[None] * len(origins) for i in range(len(solvers))]
					for j in range(len(mass_ratios))]

	for (i, q) in enumerate(mass_ratios):
		for (j, solver) in enumerate(solvers):
			for (k, origin) in enumerate(origins):
				param[i][j][k] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
						'solver': solver, 'tolerance': tolerance[k],
						'specific_frame_derivation': specific_frame_derivation})
				plot[i][j][k] = BL(**param[i][j][k])
				plot[i][j][k].plot_num_images(errors_only=False, region=region,
						region_lim=region_lim, save=False, print_errors=True)
				caustic = mm.Caustics(s=s, q=q)
				caustic.plot(s=1)
				plt.suptitle('Number of Images; ' +
						'Derived for specific frame = {}'.format(
						specific_frame_derivation), x=0.435)
				plt.show()
				# 	FIXME:	Make this into 2 sets of 3x3 subplots.

	if False:
		saved = False
		for i in range(10):
			name = '../Tables/SolverError{}'.format(i)
			if Path(name).is_file():
				continue
			plt.savefig(name)
			print(name, 'has been saved')
			saved = True
			break
		if saved == False:
			print('Error: too many files of same name already exist. File not saved')

	plt.show()


# Input parameters
s = 1.5
mass_ratios = [1e-6, 1e-8, 1e-10]
res = int(50)
solvers =  ['SG12']
origins = ['plan', 'caustic', 'geo_cent']
tolerance = [0.00007, 0.00007, 0.00007]
region = 'caustic'
region_lim = None
specific_frame_derivation = False

demonstration()

s = 1.5
mass_ratios = [1e-6, 1e-10, 1e-14]
res = int(50)
solvers =  ['SG12']
origins = ['plan', 'caustic', 'geo_cent']
tolerance = [2e-8, 0.00007, 0.00007]
region = 'caustic'
region_lim = None
specific_frame_derivation = True

demonstration()






