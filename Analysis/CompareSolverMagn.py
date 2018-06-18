# Zoey Samples
# Created: Jun 12, 2018
# CompareSolverMagn.py
# Last Updated: Jun 12, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
import MulensModel as mm
import numpy as np


# Input parameters
s = 1.5
q = 1e-7
res = int(50)
coordinates = ['geo_cent', 'caustic', 'plan']
method = ['SG12', 'numpy']
tolerance = 0.0001

param = [[None] * len(coordinates) for i in range(len(method))]
plot = [[None] * len(coordinates) for i in range(len(method))]
mag_BL = [[None] * len(coordinates) for i in range(len(method))]

print('Getting data...')
for (i, solver) in enumerate(method):
	for (j, origin) in enumerate(coordinates):
		param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
						'solver': solver, 'tolerance': tolerance})
		plot[i][j] = (BL(**param[i][j]))
		plot[i][j].grid_plots()
		mag_BL[i][j] = (plot[i][j].mag_1d)
		print('Data retrieved for: {}; {}.'.format(plot[i][j].solver_phrase,
												 plot[i][j].origin_phrase))

x_array = plot[0][0].x_array
y_array = plot[0][0].y_array

print('Plotting...')
for (j, origin) in enumerate(coordinates):
	for (i, solver) in enumerate(method):
		for k in range(i+1, len(method)):
			plt.scatter(x_array, y_array, c = (mag_BL[i][j] / mag_BL[k][j]), 
									s=10, marker = 'o', cmap='plasma', lw=None)
			mag_plot = plt.colorbar()
			mag_plot.ax.tick_params(labelsize=10)
			mag_plot.set_label('Magnification Ratio')
			plt.xlabel('X-position of source', fontsize = 12)
			plt.ylabel('Y-position of source', fontsize = 12)
			plt.xlim(min(x_array), max(x_array))
			plt.ylim(min(y_array), max(y_array))
			plt.title('Ratio of Magnification:\n{} Frame; {} Solver / {} Solver'.format(
									plot[i][j].origin_title, plot[i][j].solver_title, 
									plot[k][j].solver_title))
			plt.gcf().set_size_inches(8, 6)
			plt.show()

print('Done.\n')
