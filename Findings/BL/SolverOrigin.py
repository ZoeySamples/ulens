# Zoey Samples
# Created: Jun 21, 2018
# SolverInfo.py
# Last Updated: Jun 22, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap as LinSegCmap
import numpy as np
from BinaryLens import BinaryLens as BL
from Caustics import Caustics as caus
from pathlib import Path

"""
The most consistently accurate solver in frames other than the planet frame
is Skowron & Gould.
"""

def num_images_demo():

	param = [[None] * len(solvers) for j in range(len(origins))]
	plot = [[None] * len(solvers) for j in range(len(origins))]
	fig, ax = plt.subplots(len(origins), len(solvers))

	for (i, origin) in enumerate(origins):
		for (j, solver) in enumerate(solvers):
			idx = 1 + j + len(solvers)*i

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'refine_region': refine_region})
			plot[i][j] = BL(**param[i][j])

			# Use our two base colormaps
			orng = plt.cm.Oranges
			blues = plt.cm.Blues

			# Get the list values from each colormap
			ornglist = [orng(i) for i in range(orng.N)] # This list contains 256 colors.
			blueslist = [blues(i) for i in range(blues.N)] # This list contains 256 colors.

			# Select the regions of the colormaps we want, and slice them together.
			start = 0
			jump = 24
			clist = np.linspace(start, start+8*jump, 6)	# Slicing points for merged list.
			clist = [int(val) for val in clist]		# Convert the list into integers.

			# Create the new list with segments of the Oranges and Blues colormaps.
			colorlist = (ornglist[clist[0]:clist[3]] + blueslist[clist[2]:clist[3]] +
						 ornglist[clist[4]:clist[5]] + blueslist[clist[4]:clist[5]])

			# Create new colormap.
			cmap_images = LinSegCmap.from_list('Custom cmap', colorlist, 256)

			# Discretize the colormap.
			bounds = np.linspace(-0.5, 5.5, 7)	# This is the discretized boundary.
			norm = colors.BoundaryNorm(bounds, cmap_images.N) # This is the scale.
			ticks = np.linspace(0,5,6)	# These are the tickmark locations.

			kwargs = plot[i][j].check_kwargs()
			kwargs['cmap'] = cmap_images
			kwargs['norm'] = norm
			kwargs['s'] = 1
			kwargs['lw'] = 0
			plot[i][j].get_position_arrays()
			plot[i][j].get_num_images_array()
			if errors_only:
				(x, y, num_images) = plot[i][j].get_num_images_errors()
			else:
				(x, y, num_images) = (plot[i][j].x_array, plot[i][j].y_array,
						plot[i][j].num_images)

			# Create and adjust the plots appropriately.
			ax[i][j] = plt.subplot(len(origins), len(solvers), idx)
			sc = ax[i][j].scatter(x, y, c=num_images, vmin=0, vmax=5, **kwargs)
			caustic = caus(lens=plot[i][j], solver='SG12')
			caustic.plot_caustic(s=1, color='yellow', points=5000, lw=0)
			get_plot_parameters(plot=plot[i][j], ax=ax[i][j], i=i, j=j)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.25, 0.04, 0.55])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							 orientation='vertical')
	num_color.set_label('Number of Images', fontsize=15, labelpad=-70)
	cbar.axes.tick_params(labelsize=12)

	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/images_solver_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def magnification_demo():

	param = [[None] * len(solvers) for j in range(len(origins))]
	plot = [[None] * len(solvers) for j in range(len(origins))]
	fig, ax = plt.subplots(len(origins), len(solvers))

	for (i, origin) in enumerate(origins):
		for (j, solver) in enumerate(solvers):
			idx = 1 + j + len(solvers)*i

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim, 
					'solver': solver, 'SFD': SFD, 'refine_region': refine_region})
			plot[i][j] = BL(**param[i][j])

			# Get the data for the plots.
		#	cmap = plt.cm.YlOrRd
			cmap = plt.cm.gray
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmaplist = cmaplist[20:]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			ticks = np.array([1,10])

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
			ax[i][j] = plt.subplot(len(origins), len(solvers), idx)
			sc = ax[i][j].scatter(x, y, c=magnification, vmin=1, vmax=10, **kwargs)
			get_plot_parameters(plot=plot[i][j], ax=ax[i][j], i=i, j=j)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.25, 0.04, 0.55])
	magn_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							  orientation='vertical')
	magn_color.set_label('Magnification', fontsize=15, labelpad=-80)
	cbar.axes.tick_params(labelsize=12)

	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/magn_solver_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def get_plot_text(plot, fig):

	for (i, origin) in enumerate(origins):
		origin_title = plot[i][0].origin_title
		fontsize=16
		if format(origin_title) == 'Geometric Center':
			origin_title = 'Geometric\nCenter'
			fontsize -= 1
		if format(origin_title) == 'Center-of-Mass':
			origin_title = 'Center-\nof-Mass'
			fontsize -= 1
		fig.text(0.86, 0.94 - .48/len(origins) - .878*i/len(origins),
				'{}\nFrame'.format(origin_title), ha='left', va='center',
				fontsize=fontsize)
	fig.text(0.05, 0.12, 'q={}\ns={}'.format(q, s), fontsize=16)
	plt.subplots_adjust(wspace=0.12, hspace=0.12, top=0.93, bottom=0.06,
				left=0.28, right=0.84)
	plt.gcf().set_size_inches(3.0*len(solvers)+1.5, 2.0*len(origins)-0.5)

def get_plot_parameters(plot, ax, i, j):

	(xmin, xmax) = (min(plot.x_array), max(plot.x_array))
	(ymin, ymax) = (min(plot.y_array), max(plot.y_array))
	(dx, dy) = (xmax-xmin, ymax-ymin)
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	plt.xticks(np.arange(-0.3*dx, xmax, 0.6*dx))
	plt.yticks(np.arange(-0.3*dy, ymax, 0.3*dy))
	ax.tick_params(axis='x', labelsize=12)
	ax.tick_params(axis='y', labelsize=12)
	labels = ax.get_yticklabels()
	plt.setp(labels, verticalalignment='center', stretch='condensed')
	labels = ax.get_xticklabels()
	plt.setp(labels, stretch='condensed')
	ax.axes.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
	ax.axes.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
	if (i == 0):
		ax.axes.set_title('{} Solver'.format(
				plot.solver_title),
				fontsize=16)
	if (j!=0):
		ax.axes.get_yaxis().set_visible(False)
	if (i!=len(origins)-1):
		ax.axes.get_xaxis().set_visible(False)

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
s = 1.5
q = 4e-8
origins = ['plan', 'caustic', 'geo_cent', 'star', 'com']
res = int(400)
solvers =  ['SG12', 'zroots', 'numpy']
region = 'caustic'
region_lim = [0.6, 1.4, -0.20, 0.20]
save_fig = True
show_fig = True

refine_region = True
SFD = True
errors_only = False
#num_images_demo()
magnification_demo()


