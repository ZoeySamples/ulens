# Zoey Samples
# Created: Jun 21, 2018
# SolverInfo.py
# Last Updated: Jun 25, 2018

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

	param = [[None] * len(origins) for j in range(len(mass_ratios))]
	plot = [[None] * len(origins) for j in range(len(mass_ratios))]

	fig, ax = plt.subplots(len(mass_ratios), len(origins))

	for (i, q) in enumerate(mass_ratios):
		for (j, origin) in enumerate(origins):
			idx = 1 + j + len(origins)*i
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'solver': solver, 'specific_frame_derivation': 
					specific_frame_derivation})
			plot[i][j] = BL(**param[i][j])

			kwargs = plot[i][j].check_kwargs()
			kwargs['cmap'] = 'coolwarm'
			plot[i][j].get_position_arrays(region=region, region_lim=region_lim)
			plot[i][j].get_num_images_array()
			(x, y, num_images) = (plot[i][j].x_array, plot[i][j].y_array,
								  plot[i][j].num_images)

			ax[i][j] = plt.subplot(len(mass_ratios), len(origins), idx)
			ax[i][j].scatter(x, y, c=num_images, **kwargs)
			caustic = mm.Caustics(s=s, q=q)
			caustic.plot(s=1, color='yellow')
			(xmin, xmax) = (min(plot[i][j].x_array), max(plot[i][j].x_array))
			(ymin, ymax) = (min(plot[i][j].y_array), max(plot[i][j].y_array))
			dx = xmax - xmin
			ax[i][j].set_xlim(xmin, xmax)
			ax[i][j].set_ylim(ymin, ymax)
			ax[i][j].axes.get_xaxis().set_visible(False)
			ax[i][j].axes.get_yaxis().set_visible(False)
			ax[i][j].set_title('Mass ratio = {}; origin = {}'.format(q, origin),
								fontsize=8)
	plt.tight_layout()
	plt.suptitle('Number of Images; ' +	'Derived for specific frame = {}'.
				 format(specific_frame_derivation), x=0.435)
	plt.gcf().set_size_inches(3*len(mass_ratios), 3*len(origins))

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
mass_ratios = [1e-6, 1e-9, 1e-12]
origins = ['plan', 'caustic', 'geo_cent']
res = int(20)
solver =  'SG12'
region = 'caustic'
region_lim = None

specific_frame_derivation = False
demonstration()

specific_frame_derivation = True
demonstration()






