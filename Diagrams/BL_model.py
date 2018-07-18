# Zoey Samples
# Created: May 22, 2018
# diagrams.py
# Last Updated: Jul 10, 2018


import matplotlib.pyplot as plt
import matplotlib.axes as axes
from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
import numpy as np

def get_binary_lens_plot(save_fig=False, show_fig=True, arrow_on=True, caustic_on=False):

	param = []
	plot = []

	param = (({'s': s, 'q': q, 'origin': origin, 'solver': solver,
			'plot_frame': plot_frame, 'refine_region': refine_region,
			 'SFD': SFD}))
	plot = (BL(**param))


	causticUL = caus(lens=plot)
	z1 = causticUL.z1
	z2 = causticUL.z2

	plt.scatter(z1.real, z1.imag, s=1000, color='blue')
	plt.scatter(z2.real, z2.imag, s=2000, color='red', marker='*')

	if arrow_on:
		ax = plt.axes()
		ax.arrow(z2.real+0.05, z2.imag, z1.real - z2.real - 0.1, z1.imag - z2.imag, width=1e-6,
				 length_includes_head=True, head_width=2e-2, head_length=5e-2, fc='k', ec='k')
		ax.arrow(z1.real-0.05, z1.imag, z2.real - z1.real + 0.1, z1.imag - z2.imag, width=1e-6,
				 length_includes_head=True, head_width=2e-2, head_length=5e-2, fc='k', ec='k')

	if caustic_on:
		caustic = caus(lens=plot)
		caustic.plot_caustic(s=1, lw=0, color='orange')

	title = 'Binary Lens System'
	got_plot_adjustments(title=title)

	if save_fig:
		name='BL_model'
		_save_fig(name, arrow_on, caustic_on)

	if show_fig:
		plt.show()
	else:
		plt.clf()

def get_triple_lens_plot(save_fig=False, show_fig=True, arrow_on=False, caustic_on=False):

	param = []
	plot = []

	param = (({'s1': s1, 's2': s2, 'q1': q1, 'q2': q2, 'system': system,
			'origin': origin, 'solver': solver, 'phi': phi,
			'plot_frame': plot_frame, 'SFD': SFD}))
	plot = (TL(**param))


	causticUL = caus(lens=plot)
	z1 = causticUL.z1
	z2 = causticUL.z2
	z3 = causticUL.z3

	plt.scatter(z1.real, z1.imag, s=2000, color='red', marker='*')
	plt.scatter(z2.real, z2.imag, s=1000, color='blue')
	plt.scatter(z3.real, z3.imag, s=400, color='purple')

	if arrow_on:
		ax = plt.axes()
		ax.arrow(z2.real-0.05, z2.imag, z1.real - z2.real + 0.1, z1.imag - z2.imag, width=1e-6,
				 length_includes_head=True, head_width=2e-2, head_length=5e-2, fc='k', ec='k')
		ax.arrow(z1.real+0.05, z1.imag, z2.real - z1.real - 0.1, z1.imag - z2.imag, width=1e-6,
				 length_includes_head=True, head_width=2e-2, head_length=5e-2, fc='k', ec='k')

		dx = (z3.real - z2.real)
		dy = (z3.imag - z2.imag)
		ax.arrow(z2.real+0.10*dx, z2.imag+0.10*dy, 0.8*(dx), 0.8*(dy), width=1e-6,
				 length_includes_head=True, head_width=1.5e-2, head_length=3e-2, fc='k', ec='k')
		ax.arrow(z3.real-0.10*dx, z3.imag-0.10*dy, -0.8*(dx), -0.8*(dy), width=1e-6,
				 length_includes_head=True, head_width=1.5e-2, head_length=3e-2, fc='k', ec='k')

	if caustic_on:
		caustic = caus(lens=plot)
		caustic.plot_caustic(s=1, lw=0, color='orange')

	title = 'Triple Lens System'
	got_plot_adjustments(title=title)

	if save_fig:
		name='TL_model'
		_save_fig(name, arrow_on, caustic_on)

	if show_fig:
		plt.show()
	else:
		plt.clf()

def got_plot_adjustments(title):

	plt.xlim(-1.2, 1.2)
	plt.ylim(-0.4, 0.4)
	plt.xlabel('X-Position', fontsize=20, labelpad=10)
	plt.ylabel('Y-Position', fontsize=20, rotation=90)
	plt.xticks(np.arange(-1.0, 1.1, 0.5))
	plt.yticks(np.arange(-0.4, 0.41, 0.1))
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)
	plt.title('{}'.format(title), fontsize=32, y=1.04)
	plt.gcf().set_size_inches(12,8)

def _save_fig(name, arrow_on, caustic_on):

	if arrow_on:
		name = name + '_arr'
	if caustic_on:
		name = name +  '_caus'
	name = name + '.png'
	plt.savefig(name, dpi=200)
	print('file', name, 'saved')


# Input parameters for binary lens
s = 1.5
q = 1e-2
solver = 'SG12'
origin = 'plan'

plot_frame = 'caustic'
refine_region = True
SFD = True

save_fig = True
show_fig = False

get_binary_lens_plot(save_fig, show_fig, arrow_on=False, caustic_on=False)
get_binary_lens_plot(save_fig, show_fig, arrow_on=False, caustic_on=True)
get_binary_lens_plot(save_fig, show_fig, arrow_on=True, caustic_on=False)
get_binary_lens_plot(save_fig, show_fig, arrow_on=True, caustic_on=True)

# Input parameters for triple lens
s1 = 1.5
s2 = 2.0
q1 = 1e-2
q2 = 4e-2
phi = 155
solver = 'SG12'
origin = 'body2'
system = 'SPM'

plot_frame = 'caustic'
SFD = True

save_fig = True
show_fig = False

get_triple_lens_plot(save_fig, show_fig, arrow_on=False, caustic_on=False)
get_triple_lens_plot(save_fig, show_fig, arrow_on=False, caustic_on=True)
get_triple_lens_plot(save_fig, show_fig, arrow_on=True, caustic_on=False)
get_triple_lens_plot(save_fig, show_fig, arrow_on=True, caustic_on=True)



