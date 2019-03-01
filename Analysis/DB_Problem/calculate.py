#Zoey Samples
#Created: Sept 6, 2018
#calculate.py
# Last Updated: Oct 23, 2018

import sys
import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
import MulensModel as mm
import numpy as np
import math
import configparser


"""This is a script to convert the parameters from David Bennett's
Triple Lens calculations and convert them into the parameter space
defined in the ulens class, TripleLens. 
"""

def get_dictionaries(sect):

	parameters = dict()
	for key in config[sect]:
		if key == 'param_str':
			parameters[key] = config[sect][key]
		else:
			parameters[key] = config.getfloat(sect, key)

	source_dict = dict()
	TL_param = dict()
	for key, value in parameters.items():
		if key in ['param_str', 'v_e', 't0', 'u0', 'theta1cm']:
			source_dict[key] = parameters[key]
		else:
			TL_param[key] = parameters[key]

	key_change = [['ang23', 'eps1', 'eps2', 'sep1cm', 'sep23'],
				  ['phi', 'q1', 'q2', 's1', 's2']]

	for (i, old_key) in enumerate(key_change[0]):
		new_key = key_change[1][i]
		TL_param[new_key] = TL_param[old_key]
		del TL_param[old_key]

	return source_dict, TL_param

def get_Rhie_object(phi, s1, s2, q1, q2):
	
	param = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
				'q1': q1, 'res': res, 'origin': origin,
				'region': region, 'region_lim': region_lim,
				'solver': solver, 'system': system,
				'plot_frame': plot_frame, 'refine_region': False})
	plot = TL(**param)
	return plot

def get_new_object(plot, phi, s1, s2, q1, q2):
	q1_new = q1
	q2_new = q2
	s1_new = np.abs(plot.z2 - plot.z1)
	s2_new = np.abs(plot.z3 - plot.z1)
	s3_temp = np.abs(plot.z2 - plot.z3)
	phi_new = (180./math.pi)*math.acos((s1_new**2 + s2_new**2 - s3_temp**2) / 
									   (2.*s1_new*s2_new))
	refine_region = True
	SFD = True

	param_new = ({'s2': s2_new, 's1': s1_new, 'phi': phi_new, 'q2': q2_new,
				'q1': q1_new, 'res': res, 'origin': origin_new,
				'region': region, 'region_lim': region_lim,
				'solver': solver, 'SFD': SFD, 'system': system_new,
				'plot_frame': plot_frame, 'refine_region': refine_region})
	plot_new = (TL(**param_new))
	return plot_new

def plot_trajectory(plot, param_str, v_e, t0, u0, theta1cm, t_start, t_stop):

	z0 = u0*( math.cos(math.pi/2. - theta1cm) +
					1.j*math.sin(math.pi/2. - theta1cm) )
	t_traj = np.linspace(t_start, t_stop, 1000)
	z_traj = z0 + (t_traj - t0)*v_e*( math.cos(theta1cm) + 1.j*math.sin(theta1cm) )
	x_traj = z_traj.real
	y_traj = z_traj.imag
	caustic = caus(lens=plot, solver='SG12')
	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
	plt.scatter(x_traj, y_traj, color='blue', s=1, lw=0)
	plt.title('{} Trajectory with Caustics'.format(param_str))
	plt.xlabel('Position - Real')
	plt.ylabel('Position - Imaginary')
	plt.show()

def plot_light_curve(plot, param_str, v_e, t0, u0, theta1cm, t_start, t_stop,
					 t_pts=10000):
	z0 = u0*( math.cos(math.pi/2. - theta1cm) +
					1.j*math.sin(math.pi/2. - theta1cm) )
	t_traj = np.linspace(t_start, t_stop, t_pts)
	z_traj = z0 + (t_traj - t0)*v_e*( math.cos(theta1cm) + 1.j*math.sin(theta1cm) )
	x_traj = z_traj.real
	y_traj = z_traj.imag
	magn_traj = np.zeros(len(t_traj))
	for (i,t) in enumerate(t_traj):
		magn_traj[i] = plot.get_magnification(x_traj[i], y_traj[i])
	plt.scatter(t_traj, magn_traj, color='blue', s=1, lw=0)
	plt.ylim(1, 100)
	plt.title('{} Light Curve'.format(param_str))
	plt.ylabel('Magnification')
	plt.xlabel('Time (BJD\')')
	plt.yscale('log')
	plt.show()

def plot_num_images(plots, num_images_params):

		plots[0].plot_num_images(**num_images_params)
		caustic = caus(lens=plots[0], solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plots[1].plot_num_images(**num_images_params)
		caustic = caus(lens=plots[1], solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

def plot_magnification(plots, magn_params, show_caustic):

		plots[0].plot_magnification(**magn_params)
		caustic = caus(lens=plots[0], solver='SG12')
		if show_caustic:
			caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plots[1].plot_magnification(**magn_params)
		caustic = caus(lens=plots[1], solver='SG12')
		if show_caustic:
			caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

## Plotting parameters

# Calculation Parameters for Rhie (2002) TL object(s) (do not change)
origin = 'Rhie2002'
system = 'Rhie2002'
solver =  'SG12'

# Graphical parameters
res = int(100)
region = 'caustic_2a'
region_lim = (-8, 8, -8, 8)
plot_frame = 'caustic'

# Convert-to Parameters (in ulens/TripleLens.py parameter space)
origin_new = 'body2'
system_new = 'SSP'

t_start = 7470.
t_stop = 7510.

# Open the config file
if len(sys.argv) != 2:
    raise ValueError('Need to include cfg file')
config_file = sys.argv[1]
config = configparser.ConfigParser()
config.optionxform = str
config.read(config_file)

section = ['bsplwlb_1', 'bsplwmb_1', 'bsplwba_2']
source_dict = [dict()]*3
TL_param = [dict()]*3
for (i, sect) in enumerate(section):
	if sect not in config:
		raise KeyError('{} section not found in config file'.format(sect))
	source_dict[i], TL_param[i] = get_dictionaries(sect=sect)

## Create pairs of plot objects for each set of parameters.

plot_Rhie_1 = get_Rhie_object(**TL_param[0])
plot_SSP_1 = get_new_object(plot_Rhie_1, **TL_param[0])

plot_Rhie_2 = get_Rhie_object(**TL_param[1])
plot_SSP_2 = get_new_object(plot_Rhie_2, **TL_param[1])

plot_Rhie_3 = get_Rhie_object(**TL_param[2])
plot_SSP_3 = get_new_object(plot_Rhie_3, **TL_param[2])


## Show plots for each pair of plot objects

# Plot types
trajectory = True
light_curve = False
num_images = False
magn = False

# Plot pairs
plot1 = True
plot2 = True
plot3 = True

# Plot parameteres
save = False
size = 5


if trajectory:

	if plot1:
		plot_trajectory(plot=plot_Rhie_1, **source_dict[0], t_start=t_start, t_stop=t_stop)

	if plot2:
		plot_trajectory(plot=plot_Rhie_2, **source_dict[1], t_start=t_start, t_stop=t_stop)

	if plot3:
		plot_trajectory(plot=plot_Rhie_3, **source_dict[2], t_start=t_start, t_stop=t_stop)

if light_curve:

	if plot1:
		plot_light_curve(plot_Rhie_1, **source_dict[0], t_start=t_start, t_stop=t_stop)

	if plot2:
		plot_light_curve(plot_Rhie_2, **source_dict[1], t_start=t_start, t_stop=t_stop)

	if plot2:
		plot_light_curve(plot_Rhie_3, **source_dict[2], t_start=t_start, t_stop=t_stop)

if num_images:

	num_images_params = {'errors_only': False, 'print_errors': False,
						  'lw': 0, 'save': save, 's': size}

	if plot1:
		plots = [plot_Rhie_1, plot_SSP_1]
		plot_num_images(plots, num_images_params)

	if plot2:
		plots = [plot_Rhie_2, plot_SSP_2]
		plot_num_images(plots, num_images_params)

	if plot3:
		plots = [plot_Rhie_3, plot_SSP_3]
		plot_num_images(plots, num_images_params)

if magn:

	magn_params = {'outliers': False, 'log_colorbar': True, 'lw': 0,
				   'save': save, 's': size}
	magn_kwargs = {'show_caustic': False}

	if plot1:
		plots = [plot_Rhie_1, plot_SSP_1]
		plot_magnification(plots, magn_params, **magn_kwargs)

	if plot2:
		plots = [plot_Rhie_2, plot_SSP_2]
		plot_magnification(plots, magn_params, **magn_kwargs)

	if plot3:
		plots = [plot_Rhie_3, plot_SSP_3]
		plot_magnification(plots, magn_params, **magn_kwargs)


