# Zoey Samples
# Created: Jun 21, 2018
# OriginInfo.py
# Last Updated: Jul 2, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib
import pandas as pd
import numpy as np
from BinaryLens import BinaryLens as BL
from Caustics import Caustics as caus
import MulensModel as mm
from pathlib import Path

"""
This file shows plots for the planet frame, using the SG12 solver,for 
varying separations and mass ratios to show how well the modeling works.
"""

def num_images_demo():

	param = [[None] * len(separations) for j in range(len(mass_ratios))]
	plot = [[None] * len(separations) for j in range(len(mass_ratios))]
	fig, ax = plt.subplots(len(mass_ratios), len(separations))

	for (i, q) in enumerate(mass_ratios):
		for (j, origin) in enumerate(separations):
			idx = 1 + j + len(separations)*i

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD})
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
			plot[i][j].get_position_arrays()
			plot[i][j].get_num_images_array()
			(x, y, num_images) = (plot[i][j].x_array, plot[i][j].y_array,
						plot[i][j].num_images)

			# Create and adjust the plots appropriately.
			ax[i][j] = plt.subplot(len(mass_ratios), len(separations), idx)
			sc = ax[i][j].scatter(x, y, c=num_images, vmin=0, vmax=5, **kwargs)
			caustic = caus(lens=plot[i][j], solver='SG12')
			caustic.plot_caustic(s=1, color='yellow', points=5000)
			(xmin, xmax) = (min(plot[i][j].x_array), max(plot[i][j].x_array))
			(ymin, ymax) = (min(plot[i][j].y_array), max(plot[i][j].y_array))
			(dx, dy) = (xmax-xmin, ymax-ymin)
			plt.xlim(xmin, xmax)
			plt.ylim(ymin, ymax)
			plt.xticks(np.arange(xmin+0.2*dx, xmax, 0.6*dx))
			plt.yticks(np.arange(ymin+0.2*dy, ymax, 0.3*dy))
			ax[i][j].tick_params(axis='x', labelsize=8+len(separations))
			ax[i][j].tick_params(axis='y', labelsize=8+len(separations))
			ax[i][j].axes.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.1e'))
			ax[i][j].axes.xaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.2e'))
			if (i == 0):
				ax[i][j].axes.set_title('{} Frame\ns={}'.format(
						plot[i][j].origin_title, s), fontsize=12+len(separations))

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.20, 0.09, 0.60, 0.04])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks, orientation='horizontal')
	num_color.set_label('Number of Images', fontsize=12+len(separations), labelpad=1)
	cbar.axes.tick_params(labelsize=8+len(separations))

	for (i, q) in enumerate(mass_ratios):
		fig.text(0.93, 0.90 - .33/len(mass_ratios) - .79*i/len(mass_ratios), 'q={}'.format(q), ha='center', va='center', fontsize=12+len(separations))

	plt.subplots_adjust(wspace=0.5, hspace=0.15, top=0.90, bottom=0.19, left=0.12, right=0.88)
	plt.gcf().set_size_inches(3.0*len(separations)+1.0, 2.5*len(mass_ratios)+1.5)

	# Save the plot as a .png file.
	if save_fig:
		saved = False
		for i in range(10):
			name = '../Tables/images_SFD_{}.png'.format(i)
			if Path(name).is_file():
				continue
			plt.savefig(name)
			print(name, 'has been saved')
			saved = True
			break
		if saved == False:
			print('Error: too many files of same name already exist. File not saved')

	if show_fig:
		plt.show()

def magnification_demo():

	param = [[None] * len(separations) for j in range(len(mass_ratios))]
	plot = [[None] * len(separations) for j in range(len(mass_ratios))]
	fig, ax = plt.subplots(len(mass_ratios), len(separations))

	for (i, q) in enumerate(mass_ratios):
		for (j, origin) in enumerate(separations):
			idx = 1 + j + len(separations)*i

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim, 
					'solver': solver, 'SFD': SFD})
			plot[i][j] = BL(**param[i][j])

			# Get the data for the plots.

			cmap = plt.cm.YlOrRd
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			ticks = np.array([1,10,100])

			kwargs = plot[i][j].check_kwargs()
			kwargs['cmap'] = cmap
			kwargs['norm'] = colors.LogNorm()
			plot[i][j].get_position_arrays()
			plot[i][j].get_magnification_array()
			(x, y, magnification) = (plot[i][j].x_array, plot[i][j].y_array,
						plot[i][j].magn_array)

			# Create and adjust the plots appropriately.
			ax[i][j] = plt.subplot(len(mass_ratios), len(separations), idx)
			sc = ax[i][j].scatter(x, y, c=magnification, vmin=1, vmax=100, **kwargs)
#			caustic = caus(lens=plot[i][j], solver='SG12')
#			caustic.plot_caustic(s=1, color='blue', points=5000)
			(xmin, xmax) = (min(plot[i][j].x_array), max(plot[i][j].x_array))
			(ymin, ymax) = (min(plot[i][j].y_array), max(plot[i][j].y_array))
			(dx, dy) = (xmax-xmin, ymax-ymin)
			plt.xlim(xmin, xmax)
			plt.ylim(ymin, ymax)
			plt.xticks(np.arange(xmin+0.2*dx, xmax, 0.6*dx))
			plt.yticks(np.arange(ymin+0.2*dy, ymax, 0.3*dy))
			ax[i][j].tick_params(axis='x', labelsize=8+len(separations))
			ax[i][j].tick_params(axis='y', labelsize=8+len(separations))
			ax[i][j].axes.yaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.1e'))
			ax[i][j].axes.xaxis.set_major_formatter(
								mtick.FormatStrFormatter('%.3e'))
			if (i == 0):
				ax[i][j].axes.set_title('{} Frame\ns={}'.format(
						plot[i][j].origin_title, s), fontsize=12+len(separations))

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.20, 0.09, 0.60, 0.04])
	magn_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks, orientation='horizontal')
	magn_color.set_label('Magnification', fontsize=12+len(separations), labelpad=1)
	cbar.axes.tick_params(labelsize=8+len(separations))

	for (i, q) in enumerate(mass_ratios):
		fig.text(0.93, 0.90 - .33/len(mass_ratios) - .79*i/len(mass_ratios), 'q={}'.format(q), ha='center', va='center', fontsize=12+len(separations))

	plt.subplots_adjust(wspace=0.5, hspace=0.15, top=0.90, bottom=0.19, left=0.12, right=0.88)
	plt.gcf().set_size_inches(3.0*len(separations)+1.0, 2.5*len(mass_ratios)+1.5)

	# Save the plot as a .png file.
	if save_fig:
		saved = False
		for i in range(10):
			name = '../Tables/magn_SFD_{}.png'.format(i)
			if Path(name).is_file():
				continue
			plt.savefig(name)
			print(name, 'has been saved')
			saved = True
			break
		if saved == False:
			print('Error: too many files of same name already exist. File not saved')

	if show_fig:
		plt.show()


# Here are the input parameters for making the plots.
s = 1.5
mass_ratios = [1e-6, 1e-12]
origins = 'plan'
res = int(20)
solver =  'SG12'
region = 'caustic'
region_lim = [-.5, .5, 0.0, 2]
save_fig = False
show_fig = True

SFD = False
num_images_demo()
magnification_demo()


SFD = True
#num_images_demo()
#magnification_demo()



