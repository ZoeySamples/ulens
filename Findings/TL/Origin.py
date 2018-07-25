# Zoey Samples
# Created: Jun 21, 2018
# OriginInfo.py
# Last Updated: Jul 17, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
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

	for (i, q1), (j, q2) in product(enumerate(ps_mass_ratios), enumerate(mp_mass_ratios)):
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

			cmap = plt.cm.Blues
			cmaplist = [cmap(i) for i in range(cmap.N)]
			cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
			bounds = np.linspace(-0.5,10.5,12)
			norm = colors.BoundaryNorm(bounds, cmap.N)
			ticks = np.linspace(0,10,11)

			kwargs = plot[m][n].check_kwargs()
			kwargs['cmap'] = cmap
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
			sc = ax.scatter(x, y, c=num_images, vmin=0, vmax=10, **kwargs)
		#	caustic = caus(lens=plot[m][n], solver='SG12')
		#	caustic.plot_caustic(s=1, color='yellow', points=5000, lw=0)
			get_plot_parameters(plot=plot[m][n], ax=ax, k=k, l=l)
			fig.add_subplot(ax)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.895, 0.60, 0.04])
	num_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks, orientation='horizontal')
	num_color.set_label('Number of Images', fontsize=15, labelpad=-66)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL_images_origin_.png'
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

	for (i, q1), (j, q2) in product(enumerate(ps_mass_ratios), enumerate(mp_mass_ratios)):
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
	cbar = fig.add_axes([0.08, 0.89, 0.60, 0.04])
	magn_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							  orientation='horizontal')
	magn_color.set_label('Magnification', fontsize=15, labelpad=-66)
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL_magn_origin_.png'
		save_png(file_name)

	if show_fig:
		plt.show()

def get_plot_text(plot, fig):

	for (i, q) in enumerate(ps_mass_ratios):
		fig.text(0.97, 0.76 - .33/len(ps_mass_ratios) - .76*i/len(ps_mass_ratios),
				'q1={:.0e}'.format(q), stretch='ultra-expanded', ha='center',
				va='center', fontsize=18, rotation=90)

	for (i, q) in enumerate(mp_mass_ratios):
		fig.text(.55/len(ps_mass_ratios) + .91*i/len(ps_mass_ratios), 0.83,
				'q2={:.0e}'.format(q), stretch='ultra-expanded', ha='center',
				va='center', fontsize=18)

	fig.text(0.74, 0.905, 's1={}, s2={}\n{} Solver'.format(s1, s2, plot[0][0].solver_title),
				fontsize=16)
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
res = int(250)
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


