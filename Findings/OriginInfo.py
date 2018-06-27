# Zoey Samples
# Created: Jun 21, 2018
# SolverInfo.py
# Last Updated: Jun 26, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib
import numpy as np
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

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'solver': solver, 'specific_frame_derivation': 
					specific_frame_derivation})
			plot[i][j] = BL(**param[i][j])

			cmap = plt.cm.Blues
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			bounds = np.linspace(-0.5,5.5,7)
			norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
			ticks = np.linspace(0,5,6)

			# Get the data for the plots.
			kwargs = plot[i][j].check_kwargs()
			kwargs['cmap'] = cmap
			kwargs['norm'] = norm
			plot[i][j].get_position_arrays(region=region,
						region_lim=region_lim)
			plot[i][j].get_num_images_array()
			(x, y, num_images) = (plot[i][j].x_array, plot[i][j].y_array,
						plot[i][j].num_images)

			# Create and adjust the plots appropriately.
			ax[i][j] = plt.subplot(len(mass_ratios), len(origins), idx)
			sc = ax[i][j].scatter(x, y, c=num_images, vmin=0, vmax=5, **kwargs)
			caustic = mm.Caustics(s=s, q=q)
			caustic.plot(s=0.5, color='yellow')
			(xmin, xmax) = (min(plot[i][j].x_array), max(plot[i][j].x_array))
			(ymin, ymax) = (min(plot[i][j].y_array), max(plot[i][j].y_array))
			(dx, dy) = (xmax-xmin, ymax-ymin)
			plt.xlim(xmin, xmax)
			plt.ylim(ymin, ymax)
			plt.xticks(np.arange(xmin+0.2*dx, xmax, 0.6*dx))
			plt.yticks(np.arange(ymin+0.2*dy, ymax, 0.3*dy))
			ax[i][j].tick_params(axis='x', labelsize=11)
			ax[i][j].tick_params(axis='y', labelsize=11)
			ax[i][j].axes.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.1e'))
			ax[i][j].axes.xaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.3e'))
			plt.ylabel('Mass Ratio: {}'.format(plot[i][j].q), fontsize=14)
			if (j != 0):
				ax[i][j].axes.get_yaxis().set_visible(False)
			if (i == 0):
				ax[i][j].axes.set_title('{}\nFrame'.format(
						plot[i][j].origin_title), fontsize=14)

	# Get the string for the title.
	if specific_frame_derivation:
		sfd_str = ''
	else:
		sfd_str = 'Not '

	# Make final plot adjustments.
	cbar = fig.add_axes([0.88, 0.15, 0.03, 0.7])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks)
	num_color.set_label('Number of Images', fontsize=14, labelpad=10)
	cbar.axes.tick_params(labelsize=11) 
	plt.subplots_adjust(wspace=0.1, hspace=0.2, top=0.90, right=0.85)
#	plt.suptitle('Number of Images; ' +	'{}Using Calculation for Specific Frame'.
#				 format(sfd_str), x=0.5, y=0.97, fontsize=14)
	plt.gcf().set_size_inches(2.8*len(origins), 2.6*len(mass_ratios))

	# Save the plot as a .png file.
	if save_fig:
		saved = False
		for i in range(10):
			name = '../Tables/test_SFD_{}.png'.format(i)
			if Path(name).is_file():
				continue
			plt.savefig(name)
			print(name, 'has been saved')
			saved = True
			break
		if saved == False:
			print('Error: too many files of same name already exist. File not saved')

#	plt.show()

# Here are the input parameters for making the plots.
s = 1.5
mass_ratios = [1e-6, 1e-12]
origins = ['plan', 'caustic', 'geo_cent']
res = int(200)
solver =  'SG12'
region = 'caustic'
region_lim = None
save_fig = True

specific_frame_derivation = False
demonstration()

specific_frame_derivation = True
demonstration()



