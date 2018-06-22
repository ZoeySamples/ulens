# Zoey Samples
# Created: Jun 21, 2018
# SolverInfo.py
# Last Updated: Jun 22, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
from pathlib import Path

"""
This file highlights some observations I have made regarding the solvers.
"""

"""
Finding 1: The most accurate solver is consistently Skowron & Gould.

Finding 2: All solvers are equally good when working in the planet frame,
which produces less errors than all the other frames

These findings will be demonstrated by creating a grid of points in complex
space, centered on the planetary caustic, where each coordinate represents
the position of the source. For each point in this grid, the corresponding
number of images (for a source at that point) will be calculated. While we
expect there to be 3 images for source positions outside the caustic, and 5
images for locations inside the caustic, each of the solvers occassionally
make mistakes due to computational errors. This demonstration will compare
the number of errors resulting from calculations for each solver with varying
mass ratios and coordinate frames.
"""

# Input parameters
s = 1.5
mass_ratios = [1e-7, 5e-8, 2e-8]
res = int(80)
solvers =  ['SG12', 'numpy', 'zroots']
#origins = ['plan', 'caustic', 'geo_cent']
origins = ['plan']
tolerance = 0.00007
region = 'caustic'
region_lim = None
table = [None]*len(mass_ratios)
data = [[[None] * (1+len(solvers)) for i in range(len(origins))]
				for j in range(len(mass_ratios))]
param = [[[None] * len(solvers) for i in range(len(origins))]
				for j in range(len(mass_ratios))]
plot = [[[None] * len(solvers) for i in range(len(origins))]
				for j in range(len(mass_ratios))]
num_images_err = [[[None] * len(solvers) for i in range(len(origins))]
				for j in range(len(mass_ratios))]

for (i, q) in enumerate(mass_ratios):
	for (j, origin) in enumerate(origins):
		for (k, solver) in enumerate(solvers):
			param[i][j][k] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
							'solver': solver, 'tolerance': tolerance})
			plot[i][j][k] = BL(**param[i][j][k])
			plot[i][j][k].get_position_arrays(region=region)
			plot[i][j][k].get_num_images_array()
			num_images_err[i][j][k] = plot[i][j][k].get_num_images_errors()[2]
			data[i][j][k] = len(num_images_err[i][j][k])
		data[i][j][len(solvers)] = str('out of {}'.format(res**2))
		
fig, axs = plt.subplots(len(mass_ratios),1)
for (i, q) in enumerate(mass_ratios):
	collabel = []
	rowlabel = []
	temp_data = data[i]
	axs[i].axis('off')
	for (k, solver) in enumerate(solvers):
		collabel.append('{}'.format(plot[i][0][k].solver_title))
	for (j, origin) in enumerate(origins):
		rowlabel.append('{}'.format(plot[i][j][0].origin_title))
	plt.gcf().set_size_inches(8, 6)

	axs[i].set_title('Number of Errors\nMass Ratio = {}'.format(q), y=0.8,
					  fontsize=12)
	collabel.append('Number Points')
	table[i] = axs[i].table(cellText=data[i], colLabels=collabel,
							rowLabels=rowlabel,
			colWidths=[0.15 for x in range(len(collabel) - 1)]+[0.25],
						loc='center')
	table[i].scale(1, 1.2)

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




