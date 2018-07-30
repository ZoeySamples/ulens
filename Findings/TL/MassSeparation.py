# Zoey Samples
# Created: Jun 21, 2018
# OriginInfo.py
# Last Updated: Jul 2, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import numpy as np
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
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
		for (j, s) in enumerate(separations):
			idx = 1 + j + len(separations)*i

			# Initialize each binary lens system with the BinaryLens class.
			if idx == 2:
				region = regions[1]
			else:
				region = regions[0]

			param[i][j] = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame, 'refine_region': refine_region})
			plot[i][j] = TL(**param[i][j])

			cmap = plt.cm.Blues
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			bounds = np.linspace(-0.5,5.5,7)
			norm = colors.BoundaryNorm(bounds, cmap.N)
			ticks = np.linspace(0,5,6)

			# Get the data for the plots.
			kwargs = plot[i][j].check_kwargs()
			kwargs['cmap'] = cmap
			kwargs['norm'] = norm
			kwargs['s'] = 1
			kwargs['lw'] = 0
			plot[i][j].get_position_arrays()
			plot[i][j].get_num_images_array()
			(x, y, num_images) = (plot[i][j].x_array, plot[i][j].y_array,
						plot[i][j].num_images)

			# Create and adjust the plots appropriately.
			ax[i][j] = plt.subplot(len(mass_ratios), len(separations), idx)
			sc = ax[i][j].scatter(x, y, c=num_images, vmin=0, vmax=5, **kwargs)
		#	caustic = caus(lens=plot[i][j], solver='SG12')
		#	caustic.plot_caustic(s=1, color='yellow', points=2000, lw=0)
			get_plot_parameters(plot=plot[i][j], ax=ax[i][j], i=i, j=j)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.12, 0.90, 0.60, 0.035])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
				orientation='horizontal')
	num_color.set_label('Number of Images', fontsize=16, labelpad=-62)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL_images_plan_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def magnification_demo():

	param = [[None] * len(separations) for j in range(len(mass_ratios))]
	plot = [[None] * len(separations) for j in range(len(mass_ratios))]
	fig, ax = plt.subplots(len(mass_ratios), len(separations))

	for (i, q) in enumerate(mass_ratios):
		for (j, s) in enumerate(separations):
			idx = 1 + j + len(separations)*i

			# Initialize each binary lens system with the BinaryLens class.
			if idx == 2:
				region = regions[1]
			else:
				region = regions[0]

			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim, 
					'solver': solver, 'SFD': SFD, 'refine_region': refine_region})
			plot[i][j] = TL(**param[i][j])

			# Get the data for the plots.
			cmap = plt.cm.YlOrRd
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			ticks = np.array([1,10,100])

			kwargs = plot[i][j].check_kwargs()
			kwargs['cmap'] = cmap
			kwargs['norm'] = colors.LogNorm()
			kwargs['s'] = 1
			kwargs['lw'] = 0
			plot[i][j].get_position_arrays()
			plot[i][j].get_magnification_array()
			(x, y, magnification) = (plot[i][j].x_array, plot[i][j].y_array,
						plot[i][j].magn_array)

			# Create and adjust the plots appropriately.
			ax[i][j] = plt.subplot(len(mass_ratios), len(separations), idx)
			sc = ax[i][j].scatter(x, y, c=magnification, vmin=1, vmax=200, **kwargs)
			get_plot_parameters(plot=plot[i][j], ax=ax[i][j], i=i, j=j)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.12, 0.90, 0.60, 0.035])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
				orientation='horizontal')
	num_color.set_label('Magnification', fontsize=15, labelpad=-62)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL_magn_plan_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def get_plot_text(plot, fig):

	for (i, q) in enumerate(mass_ratios):
		fig.text(0.925, 0.81 - .33/len(mass_ratios) - .77*i/len(mass_ratios),
				'q={:.0e}'.format(q), ha='center', va='center', fontsize=16)
	fig.text(0.79, 0.905, '{} Frame\n{} Solver'.format(plot[0][0].origin_title,
				plot[0][0].solver_title), fontsize=16)
	plt.subplots_adjust(wspace=0.35, hspace=0.23, top=0.82, bottom=0.06,
				left=0.12, right=0.88)
	plt.gcf().set_size_inches(2.8*len(separations)+0.5, 1.8*len(mass_ratios)+1.0)

def get_plot_parameters(plot, ax, i, j):

	(xmin, xmax) = (min(plot.x_array), max(plot.x_array))
	(ymin, ymax) = (min(plot.y_array), max(plot.y_array))
	(dx, dy) = (xmax-xmin, ymax-ymin)
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	plt.xticks(np.arange(-0.3*dx, xmax, 0.6*dx))
	plt.yticks(np.arange(-0.25*dy, ymax, 0.5*dy))
	ax.tick_params(axis='x', labelsize=11)
	ax.tick_params(axis='y', labelsize=11)
	ax.axes.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
	ax.axes.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
	labels = ax.get_yticklabels()
	plt.setp(labels, rotation=60, verticalalignment='center', stretch='condensed',  x=0.05)
	labels = ax.get_xticklabels()
	plt.setp(labels, stretch='condensed')
	if (i == 0):
		ax.axes.set_title('s={}'.format(separations[j]), fontsize=16)

def save_png(file_name):

	for i in range(10):
		name = file_name[:-4] + '{}'.format(i) + file_name[-4:]
		if Path(name).is_file():
			continue
		plt.savefig(name, dpi=300)
		print(name, 'has been saved')
		return
	print('Error: too many files of same name already exist. File not saved')

# Here are the input parameters for making the plots.
separations = [0.6, 0.9, 1.1, 2.0, 5.0]
mass_ratios = [1e-1, 1e-3, 1e-7, 1e-12]
origin = 'plan'
res = int(20)
solver =  'SG12'
regions = ['caustic_a', 'custom_a']
region_lim = [-2, 2.8, -5.08, 2.2]
save_fig = False
show_fig = True

refine_region = True
SFD = True
num_images_demo()
magnification_demo()

