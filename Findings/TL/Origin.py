# Zoey Samples
# Created: Jun 21, 2018
# OriginInfo.py
# Last Updated: Jul 17, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap as LinSegCmap
import matplotlib.gridspec as gridspec
from itertools import product
import numpy as np
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
from pathlib import Path


def num_images_demo():

	num_outer_plots = len(ps_mass_ratios)*len(mp_mass_ratios)
	num_inner_plots = len(origins)*len(angles)
	fig = plt.figure(figsize=(10, 8))
	outer = gridspec.GridSpec(len(ps_mass_ratios), len(mp_mass_ratios),
							  wspace=0.20, hspace=0.35)
	param = [[None]*num_inner_plots for j in range(num_outer_plots)]
	plot = [[None]*num_inner_plots for j in range(num_outer_plots)]

	for (i, q2), (j, q1) in product(enumerate(mp_mass_ratios), enumerate(ps_mass_ratios)):
		outer_idx = j + i*len(mp_mass_ratios)
		inner = gridspec.GridSpecFromSubplotSpec(len(angles), len(origins),
					subplot_spec=outer[outer_idx], wspace=0.1, hspace=0.25)

		for (k, phi), (l, origin) in product(enumerate(angles), enumerate(origins)):
			inner_idx = l + k*len(origins)
			(m, n) = (outer_idx, inner_idx)

			# Initialize each triple lens system with the TripleLens class.
			param[m][n] = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame, 'refine_region': refine_region})
			plot[m][n] = TL(**param[m][n])

			"""
			# Old color map -- all blue
			cmap = plt.cm.Blues
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			bounds = np.linspace(-0.5,10.5,12)
			norm = colors.BoundaryNorm(bounds, cmap.N)
			ticks = np.linspace(0,10,11)
			"""

			# Use our two base colormaps
			orng = plt.cm.Oranges
			blues = plt.cm.Blues

			# Get the list values from each colormap
			ornglist = [orng(i) for i in range(orng.N)] # This list contains 256 colors.
			blueslist = [blues(i) for i in range(blues.N)] # This list contains 256 colors.

			# Select the regions of the colormaps we want, and slice them together.
			start = 0
			jump = 24
			clist = np.linspace(start, start+10*jump, 11)	# Slicing points for merged list.
			clist = [int(val) for val in clist]		# Convert the list into integers.

			# Create the new list with segments of the Oranges and Blues colormaps.
			colorlist = (ornglist[clist[0]:clist[4]] + blueslist[clist[1]:clist[2]] +
						 ornglist[clist[4]:clist[5]] + blueslist[clist[3]:clist[4]] +
						 ornglist[clist[6]:clist[7]] + blueslist[clist[5]:clist[6]] +
						 ornglist[clist[8]:clist[9]] + blueslist[clist[7]:clist[8]])

			# Create new colormap.
			cmap_images = LinSegCmap.from_list('Custom cmap', colorlist, 256)

			# Discretize the colormap.
			bounds = np.linspace(-0.5,10.5,12)	# This is the discretized boundary.
			norm = colors.BoundaryNorm(bounds, cmap_images.N) # This is the scale.
			ticks = np.linspace(0,10,11)	# These are the tickmark locations.

			kwargs = plot[m][n].check_kwargs()
			kwargs['cmap'] = cmap_images
			kwargs['norm'] = norm
			kwargs['s'] = 1
			kwargs['lw'] = 0

			# Get the data for the plots.
			plot[m][n].get_position_arrays()
			plot[m][n].get_num_images_array()
			(x, y, num_images) = (plot[m][n].x_array, plot[m][n].y_array,
								  plot[m][n].num_images)

			# Create and adjust the plots appropriately.
			ax = plt.subplot(inner[n])

			# Set the limits of the colormap values in the scatter plot to 0 and 10.
			sc = ax.scatter(x, y, c=num_images, vmin=0, vmax=10, **kwargs)
		#	caustic = caus(lens=plot[m][n], solver='SG12')
		#	caustic.plot_caustic(s=1, color='yellow', points=5000, lw=0)
			get_plot_parameters(plot=plot[m][n], ax=ax, k=k, l=l)
			fig.add_subplot(ax)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.895, 0.50, 0.04])
	num_color = plt.colorbar(sc, cax=cbar, cmap=cmap_images, ticks=ticks,
							orientation='horizontal')

	num_color.set_label('Number of Images', fontsize=15, labelpad=-70)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL/num_orgn_{}_.png'.format(system)
		save_png(file_name)

	if show_fig:
		plt.show()

def magnification_demo():

	num_outer_plots = len(ps_mass_ratios)*len(mp_mass_ratios)
	num_inner_plots = len(origins)*len(angles)
	fig = plt.figure(figsize=(10, 8))
	outer = gridspec.GridSpec(len(ps_mass_ratios), len(mp_mass_ratios),
							  wspace=0.20, hspace=0.35)
	param = [[None]*num_inner_plots for j in range(num_outer_plots)]
	plot = [[None]*num_inner_plots for j in range(num_outer_plots)]

	for (i, q2), (j, q1) in product(enumerate(mp_mass_ratios), enumerate(ps_mass_ratios)):
		outer_idx = j + i*len(mp_mass_ratios)
		inner = gridspec.GridSpecFromSubplotSpec(len(angles), len(origins),
					subplot_spec=outer[outer_idx], wspace=0.1, hspace=0.25)

		for (k, phi), (l, origin) in product(enumerate(angles), enumerate(origins)):
			inner_idx = l + k*len(origins)
			(m, n) = (outer_idx, inner_idx)

			# Initialize each triple lens system with the TripleLens class.
			param[m][n] = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame, 'refine_region': refine_region})
			plot[m][n] = TL(**param[m][n])

			cmap = plt.cm.YlOrRd
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			ticks = np.array([1,10,100])

			kwargs = plot[m][n].check_kwargs()
			kwargs['cmap'] = cmap
			kwargs['norm'] = colors.LogNorm()
			kwargs['s'] = 1
			kwargs['lw'] = 0

			# Get the data for the plots.
			plot[m][n].get_position_arrays()
			plot[m][n].get_magnification_array()
			(x, y, magnification) = (plot[m][n].x_array, plot[m][n].y_array,
						plot[m][n].magn_array)

			# Create and adjust the plots appropriately.
			ax = plt.subplot(inner[n])
			sc = ax.scatter(x, y, c=magnification, vmin=1, vmax=100, **kwargs)
	#		caustic = caus(lens=plot[m][n], solver='SG12')
	#		caustic.plot_caustic(s=1, color='yellow', points=5000, lw=0)
			get_plot_parameters(plot=plot[m][n], ax=ax, k=k, l=l)
			fig.add_subplot(ax)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.89, 0.50, 0.04])
	magn_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							  orientation='horizontal')
	magn_color.set_label('Magnification', fontsize=15, labelpad=-70)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL/mag_orgn_{}_.png'.format(system)
		save_png(file_name)

	if show_fig:
		plt.show()

def get_plot_text(plot, fig):

	for (i, q) in enumerate(mp_mass_ratios):
		fig.text(0.97, 0.76 - .33/len(mp_mass_ratios) - .76*i/len(mp_mass_ratios),
				'q2={:.0e}'.format(q), stretch='ultra-expanded', ha='center',
				va='center', fontsize=18, rotation=90)

	for (i, q) in enumerate(ps_mass_ratios):
		fig.text(.53/len(ps_mass_ratios) + .90*i/len(ps_mass_ratios), 0.83,
				'q1={:.0e}'.format(q), stretch='ultra-expanded', ha='center',
				va='center', fontsize=18)

	fig.text(0.66, 0.95, 's1={}, s2={}\n{} Solver'.format(s1, s2, plot[0][0].solver_title),
				fontsize=16)
	fig.text(0.66, 0.91, '{}'.format(plot[0][0].sys_string), fontsize=16)
	fig.text(0.66, 0.88, '{}'.format(plot[0][0].caustic_phrase), fontsize=16)
	plt.subplots_adjust(top=0.75, bottom=0.06, left=0.08, right=0.92)
	plt.gcf().set_size_inches(1.8*len(origins)*len(mp_mass_ratios)+1.5,
							  1.3*len(angles)*len(ps_mass_ratios)+1.5)

def get_plot_parameters(plot, ax, k, l):

	(xmin, xmax) = (min(plot.x_array), max(plot.x_array))
	(ymin, ymax) = (min(plot.y_array), max(plot.y_array))
	(dx, dy) = (xmax-xmin, ymax-ymin)
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	plt.xticks(np.arange(-0.3*dx, xmax, 0.6*dx))
	plt.yticks(np.arange(-0.3*dy, ymax, 0.6*dy))
	ax.tick_params(axis='x', labelsize=12)
	ax.tick_params(axis='y', labelsize=12)

	ax.axes.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
	ax.axes.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
	labels = ax.get_yticklabels()
	plt.setp(labels, x=0.07, rotation=60, verticalalignment='center', stretch='extra-condensed')
	labels = ax.get_xticklabels()
	plt.setp(labels, stretch='extra-condensed')

	if (k == 0):
		ax.axes.set_title('{}\nFrame'.format(plot.origin_title), fontsize=16)
	if (l != 0) or (k != len(angles)-1):
		ax.axes.get_yaxis().set_visible(False)
		ax.axes.get_xaxis().set_visible(False)
	if (l == len(origins)-1):
		ax.axes.text(1.2*xmax, 0.0, 'phi={}'.format(angles[k]), ha='center', va='center',
				fontsize=14, rotation=90)


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
s1 = 1.5
s2 = 1.0
ps_mass_ratios = [1e-3, 1e-6]
mp_mass_ratios = [1e-2, 1e-3]
angles = [90, 135, 180]
system = 'SPM'

origins = ['geo_cent', 'body2', 'body3']
res = int(200)
solver =  'SG12'
region = 'caustic_2'
region_lim = [-.5, .5, 0.0, 2]
save_fig = True
show_fig = False

refine_region = False
plot_frame = 'caustic'

SFD = True
num_images_demo()
magnification_demo()


