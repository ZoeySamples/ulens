# Zoey Samples
# Created: Jun 21, 2018
# SolverEqual.py
# Last Updated: Jul 10, 2018

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
All solvers are equally successful when working in the planet frame.
"""

def num_images_demo():

	param = [[None] * len(solvers) for j in range(len(mass_ratios))]
	plot = [[None] * len(solvers) for j in range(len(mass_ratios))]
	fig, ax = plt.subplots(len(mass_ratios), len(solvers))

	for (i, q) in enumerate(mass_ratios):
		for (j, solver) in enumerate(solvers):
			idx = 1 + j + len(solvers)*i

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim, 'SFD': SFD,
					'refine_region': refine_region, 'solver': solver})
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
			ax[i][j] = plt.subplot(len(mass_ratios), len(solvers), idx)
			sc = ax[i][j].scatter(x, y, c=num_images, vmin=0, vmax=5, **kwargs)
	#		caustic = caus(lens=plot[i][j], solver='SG12')
	#		caustic.plot_caustic(s=1, color='yellow', points=5000, lw=0)
			get_plot_parameters(plot=plot[i][j], ax=ax[i][j], i=i, j=j)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.12, 0.90, 0.60, 0.04])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							 orientation='horizontal')
	num_color.set_label('Number of Images', fontsize=15, labelpad=-60)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/images_plan_solver_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def magnification_demo():

	param = [[None] * len(solvers) for j in range(len(mass_ratios))]
	plot = [[None] * len(solvers) for j in range(len(mass_ratios))]
	fig, ax = plt.subplots(len(mass_ratios), len(solvers))

	for (i, q) in enumerate(mass_ratios):
		for (j, solver) in enumerate(solvers):
			idx = 1 + j + len(solvers)*i

			# Initialize each binary lens system with the BinaryLens class.
			param[i][j] = ({'s': s, 'q': q, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim, 
					'solver': solver, 'SFD': SFD, 'refine_region': refine_region})
			plot[i][j] = BL(**param[i][j])

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
			ax[i][j] = plt.subplot(len(mass_ratios), len(solvers), idx)
			sc = ax[i][j].scatter(x, y, c=magnification, vmin=1, vmax=100, **kwargs)
			get_plot_parameters(plot=plot[i][j], ax=ax[i][j], i=i, j=j)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.12, 0.90, 0.60, 0.04])
	magn_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							  orientation='horizontal')
	magn_color.set_label('Magnification', fontsize=15, labelpad=-60)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/magn_plan_solver_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def get_plot_text(plot, fig):

	for (i, q) in enumerate(mass_ratios):
		fig.text(0.93, 0.80 - .33/len(mass_ratios) - .78*i/len(mass_ratios),
				'q={:.0e}'.format(q), ha='center', va='center', fontsize=16)
	fig.text(0.78, 0.89, '{} Frame\ns={}'.format(plot[0][0].origin_title, s),
				fontsize=16)
	plt.subplots_adjust(wspace=0.30, hspace=0.22, top=0.80, bottom=0.06,
				left=0.12, right=0.87)
	plt.gcf().set_size_inches(3.0*len(solvers)+1.0, 2.0*len(mass_ratios)+1.0)

def get_plot_parameters(plot, ax, i, j):

	(xmin, xmax) = (min(plot.x_array), max(plot.x_array))
	(ymin, ymax) = (min(plot.y_array), max(plot.y_array))
	(dx, dy) = (xmax-xmin, ymax-ymin)
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	plt.xticks(np.arange(-0.3*dx, xmax, 0.6*dx))
	plt.yticks(np.arange(-0.25*dy, ymax, 0.5*dy))
	ax.tick_params(axis='x', labelsize=12)
	ax.tick_params(axis='y', labelsize=12)
	labels = ax.get_yticklabels()
	plt.setp(labels, rotation=60, verticalalignment='center', x=0.05, stretch='condensed')
	labels = ax.get_xticklabels()
	plt.setp(labels, stretch='condensed')
	ax.axes.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
	ax.axes.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3e'))
	if (i == 0):
		ax.axes.set_title('{} Solver'.format(
				plot.solver_title),
				fontsize=16)

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
mass_ratios = [1e-14, 1e-15, 1e-16]
origin = 'plan'
res = int(250)
solvers =  ['SG12', 'zroots', 'numpy']
region = 'custom'
region_lim = [0.6, 1.4, -0.20, 0.20]
save_fig = True
show_fig = False

refine_region = False
SFD = True
errors_only = False
num_images_demo()
magnification_demo()


