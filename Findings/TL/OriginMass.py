# Zoey Samples
# Created: Jun 21, 2018
# OriginInfo.py
# Last Updated: Jul 17, 2018

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap as LinSegCmap
from itertools import product
import numpy as np
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
from pathlib import Path


def num_images_demo():

	num_outer_plots = len(ps_mass_ratios)*len(mp_mass_ratios)
	num_inner_plots = len(origins)*len(solvers)
	fig = plt.figure(figsize=(10, 8))
	outer = gridspec.GridSpec(len(ps_mass_ratios), len(mp_mass_ratios),
							  wspace=0.50, hspace=0.25)
	param = [[None]*num_inner_plots for j in range(num_outer_plots)]
	plot = [[None]*num_inner_plots for j in range(num_outer_plots)]

	for (i, q2), (j, q1) in product(enumerate(mp_mass_ratios), enumerate(ps_mass_ratios)):
		outer_idx = j + i*len(mp_mass_ratios)
		inner = gridspec.GridSpecFromSubplotSpec(len(origins), len(solvers),
					subplot_spec=outer[outer_idx], wspace=0.1, hspace=0.15)

		for (k, origin), (l, solver) in product(enumerate(origins), enumerate(solvers)):
			inner_idx = l + k*len(solvers)
			(m, n) = (outer_idx, inner_idx)

			# Initialize each triple lens system with the TripleLens class.
			param[m][n] = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame, 'refine_region': refine_region})
			plot[m][n] = TL(**param[m][n])

			"""
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
			sc = ax.scatter(x, y, c=num_images, vmin=0, vmax=10, **kwargs)
	#		caustic = caus(lens=plot[m][n], solver='SG12')
	#		caustic.plot_caustic(s=1, color='yellow', points=5000, lw=0)
			get_inner_plot_parameters(plot=plot[m][n], ax=ax, k=k, l=l)
			fig.add_subplot(ax)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.10, 0.04, 0.70])
	num_color = plt.colorbar(sc, cax=cbar, cmap=cmap_images, ticks=ticks,
							 orientation='vertical')
	num_color.set_label('Number of Images', fontsize=15, labelpad=-40-15*len(ps_mass_ratios))
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL/num_qs_{}_.png'.format(system)
		save_png(file_name)

	if show_fig:
		plt.show()

def magnification_demo():

	num_outer_plots = len(ps_mass_ratios)*len(mp_mass_ratios)
	num_inner_plots = len(origins)*len(solvers)
	fig = plt.figure(figsize=(10, 8))
	outer = gridspec.GridSpec(len(ps_mass_ratios), len(mp_mass_ratios),
							  wspace=0.50, hspace=0.25)
	param = [[None]*num_inner_plots for j in range(num_outer_plots)]
	plot = [[None]*num_inner_plots for j in range(num_outer_plots)]

	for (i, q2), (j, q1) in product(enumerate(mp_mass_ratios), enumerate(ps_mass_ratios)):
		outer_idx = j + i*len(mp_mass_ratios)
		inner = gridspec.GridSpecFromSubplotSpec(len(origins), len(solvers),
					subplot_spec=outer[outer_idx], wspace=0.1, hspace=0.15)

		for (k, origin), (l, solver) in product(enumerate(origins), enumerate(solvers)):
			inner_idx = l + k*len(solvers)
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
			get_inner_plot_parameters(plot=plot[m][n], ax=ax, k=k, l=l)
			fig.add_subplot(ax)

	# Add an axis for the color bar.
	cbar = fig.add_axes([0.08, 0.10, 0.04, 0.70])
	magn_color = plt.colorbar(sc, cax=cbar, cmap=kwargs['cmap'], ticks=ticks,
							  orientation='vertical')
	magn_color.set_label('Magnification', fontsize=15, labelpad=-50-15*len(ps_mass_ratios))
	cbar.axes.tick_params(labelsize=12)
	get_plot_text(plot, fig)

	# Save the plot as a .png file.
	if save_fig:
		file_name = '../../Tables/TL/mag_qs_{}_.png'.format(system)
		save_png(file_name)

	if show_fig:
		plt.show()

def plot_images():

	num_outer_plots = len(ps_mass_ratios)*len(mp_mass_ratios)
	num_inner_plots = len(origins)*len(solvers)
	fig = plt.figure(figsize=(10, 8))
	outer = gridspec.GridSpec(len(ps_mass_ratios), len(mp_mass_ratios),
							  wspace=0.25, hspace=0.38)
	param = [[None]*num_inner_plots for j in range(num_outer_plots)]
	plot = [[None]*num_inner_plots for j in range(num_outer_plots)]

	for (i, q2), (j, q1) in product(enumerate(mp_mass_ratios), enumerate(ps_mass_ratios)):
		outer_idx = j + i*len(mp_mass_ratios)
		inner = gridspec.GridSpecFromSubplotSpec(len(origins), len(solvers),
					subplot_spec=outer[outer_idx], wspace=0.1, hspace=0.25)

		for (k, origin), (l, solver) in product(enumerate(origins), enumerate(solvers)):
			inner_idx = l + k*len(solvers)
			(m, n) = (outer_idx, inner_idx)

			# Initialize each triple lens system with the TripleLens class.
			param[m][n] = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
					'q1': q1, 'res': res, 'origin': origin,
					'region': region, 'region_lim': region_lim,
					'solver': solver, 'SFD': SFD, 'system': system,
					'plot_frame': plot_frame, 'refine_region': refine_region})
			plot[m][n] = TL(**param[m][n])

			roots = plot[m][n].get_roots(x=0., y=0.)
			accepted_images = plot[m][n].get_accepted_solutions(x=0., y=0.)
			rejected_images = []

			for root in roots:
				if root in accepted_images:
					continue
				else:
					rejected_images.append(root)

			accepted_images = np.array(accepted_images)
			rejected_images = np.array(rejected_images)

			caustic = caus(lens=plot[m][n], solver='SG12')
			z1 = caustic.z1
			z2 = caustic.z2
			z3 = caustic.z3

			ax = plt.subplot(inner[n])
			fig.add_subplot(ax)

			s = 15
		#	sc = ax.scatter(z1.real, z1.imag, marker='*', s=10*s, color='blue', lw=0)
		#	sc = ax.scatter(z2.real, z2.imag, marker='o', s=4*s, color='blue', lw=0)
		#	sc = ax.scatter(z3.real, z3.imag, marker='o', s=s, color='blue', lw=0)

			sc = ax.scatter(accepted_images.real, accepted_images.imag,
						s=s, color='black', lw=0)

			sc = ax.scatter(rejected_images.real, rejected_images.imag,
						s=s, color='red', lw=0)

		#	caustic.plot_caustic(s=1, color='orange', points=5000, lw=0)

			plt.xlim(-5, 5)
			plt.ylim(-5, 5)


			ax.axes.get_xaxis().set_visible(False)
			ax.axes.get_yaxis().set_visible(False)

		#	ax.axes.text(4.5, 4.5, 'num accepted: {}'.format(len(accepted_images)), ha='center',
		#			va='center', fontsize=8)

	get_plot_text(plot, fig)

	if save_fig:
		file_name = '../../Tables/TL/soln_slvr_q_{}_.png'.format(system)
		save_png(file_name)

	if show_fig:
		plt.show()


def get_plot_text(plot, fig):

	for (i, q) in enumerate(mp_mass_ratios):
		fig.text(0.95, 0.93 - .36/len(mp_mass_ratios) - .96*i/len(mp_mass_ratios),
				'q2={:.0e}'.format(q), stretch='ultra-expanded', ha='center',
				va='center', fontsize=18, rotation=90)

	for (i, q) in enumerate(ps_mass_ratios):
		fig.text(0.24 + .24/len(ps_mass_ratios) + .73*i/len(ps_mass_ratios), 0.94,
				'q1={:.0e}'.format(q), stretch='ultra-expanded', ha='center',
				va='center', fontsize=18)

	fig.text(0.02, 0.94, 's1={}, s2={}, phi={}'.format(s1, s2, phi), fontsize=15, stretch='ultra-condensed')
	fig.text(0.02, 0.91, '{}'.format(plot[0][0].sys_string), fontsize=15, stretch='ultra-condensed')
	fig.text(0.02, 0.88, '{}'.format(plot[0][0].caustic_phrase), fontsize=15, stretch='ultra-condensed')
	if len(ps_mass_ratios) == 2:
		plt.subplots_adjust(top=0.92, bottom=0.06, left=0.26, right=0.85)
	elif len(ps_mass_ratios) == 3:
		plt.subplots_adjust(top=0.92, bottom=0.06, left=0.23, right=0.88)

	plt.gcf().set_size_inches(2.5*len(ps_mass_ratios)+5.5,
							  1.8*len(origins)*len(mp_mass_ratios))

def get_inner_plot_parameters(plot, ax, k, l):


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
	plt.setp(labels, x=0.05, rotation=60, verticalalignment='center',
			 stretch='extra-condensed')
	labels = ax.get_xticklabels()
	plt.setp(labels, stretch='extra-condensed')

#	if (k == 0):
#		ax.axes.set_title('{}\nSolver'.format(plot.solver_title), fontsize=16)
	if (l != 0) or (k != len(origins)-1):
		ax.axes.get_yaxis().set_visible(False)
		ax.axes.get_xaxis().set_visible(False)
	if (l == len(solvers)-1):
		title = plot.origin_title
		xratio = 1.25
		if title == 'Geometric Center':
			title = 'Geometric\nCenter'
			xratio = 1.50
		ax.axes.text(xratio*xmax, 0.0, '{}\nFrame'.format(title), ha='center',
				va='center', fontsize=14, rotation=90)

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
s2 = 0.8
ps_mass_ratios = [1e-5, 1e-6]
mp_mass_ratios = [1e-2, 1e-3]
phi = 135
system = 'SPM'

origins = ['body2', 'body3']
res = int(250)
solvers =  ['SG12']
region = 'caustic_2'
region_lim = [-.5, .5, 0.0, 2]
refine_region = True
plot_frame = 'caustic'
save_fig = True
show_fig = False

refine_region = True
plot_frame = 'caustic'

SFD = True
num_images_demo()
magnification_demo()

s1 = 1.5
s2 = 1.5
ps_mass_ratios = [1e-2, 1e-4, 1e-6]
mp_mass_ratios = [1e-2, 1e-4, 1e-6]
system = 'SPP'

num_images_demo()
magnification_demo()
#plot_images()



