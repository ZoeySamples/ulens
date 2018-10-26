#Zoey Samples
#Created: Sept 6, 2018
#calculate.py
# Last Updated: Oct 23, 2018

import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL
from Caustics import Caustics as caus
import MulensModel as mm
import numpy as np
import math


"""This is a script to convert the parameters from David Bennett's
Triple Lens calculations and convert them into the parameter space
defined in the ulens class, TripleLens. 
"""

def get_Rhie_object():
	
	param = ({'s2': s2, 's1': s1, 'phi': phi, 'q2': q2,
				'q1': q1, 'res': res, 'origin': origin,
				'region': region, 'region_lim': region_lim,
				'solver': solver, 'system': system,
				'plot_frame': plot_frame})
	plot = TL(**param)
	return plot

def get_new_object(plot):
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


# Calculation Parameters for Rhie (2002) TL object(s) (do not change)
origin = 'Rhie2002'
system = 'Rhie2002'
solver =  'SG12'

# Graphical parameters
res = int(150)
region = 'caustic_2a'
region_lim = (-8, 8, -8, 8)
plot_frame = 'caustic'

# Convert-to Parameters (in ulens/TripleLens.py parameter space)
origin_new = 'body2'
system_new = 'SSP'

## Rhie (2002) Lens parameters
#bsplwlb_1
q1 = 0.182580613773819
q2 = 5.771314214367329E-004
s1 = 2.08319873661748
s2 = 0.907858239420497 
phi = -0.808129441734630

# Make a Rhie (2002) TL object and a matching SSP pair
plot_Rhie_1 = get_Rhie_object()
plot_SSP_1 = get_new_object(plot_Rhie_1)

## Rhie (2002) Lens parameters
#bsplwmb_1
q1 = 0.184974402810256
q2 = 2.510336552345038E-002
s1 = 2.08381807183272
s2 = 1.51361923300292
phi = -0.824251125263340

# Make a Rhie (2002) TL object and a matching SSP pair
plot_Rhie_2 = get_Rhie_object()
plot_SSP_2 = get_new_object(plot_Rhie_2)

## Rhie (2002) Lens parameters
#bsplwba_2
q1 = 0.180596618168268
q1 = 1.077013484344390E-003
s1 = 2.10173374019687
s2 = 0.983945204539880
phi = -0.844107672043463

# Make a Rhie (2002) TL object and a matching SSP pair
plot_Rhie_3 = get_Rhie_object()
plot_SSP_3 = get_new_object(plot_Rhie_3)


## Make plots for each pair
num_images = True
magn = True
plot1 = True
plot2 = True
plot3 = True

if num_images:

	num_images_params = {'errors_only': False, 'save': False,
						 'print_errors': False, 's': 3}

	if plot1:

		plot_Rhie_1.plot_num_images(**num_images_params)
		caustic = caus(lens=plot_Rhie_1, solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plot_SSP_1.plot_num_images(**num_images_params)
		caustic = caus(lens=plot_SSP_1, solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

	if plot2:

		plot_Rhie_2.plot_num_images(**num_images_params)
		caustic = caus(lens=plot_Rhie_2, solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plot_SSP_2.plot_num_images(**num_images_params)
		caustic = caus(lens=plot_SSP_2, solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

	if plot3:

		plot_Rhie_3.plot_num_images(**num_images_params)
		caustic = caus(lens=plot_Rhie_3, solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plot_SSP_3.plot_num_images(**num_images_params)
		caustic = caus(lens=plot_SSP_3, solver='SG12')
		caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

if magn:

	magn_params = {'outliers': False, 'log_colorbar': True, 'save': False, 's': 4}

	if plot1:

		plot_Rhie_1.plot_magnification(**magn_params)
		caustic = caus(lens=plot_Rhie_1, solver='SG12')
	#	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plot_SSP_1.plot_magnification(**magn_params)
		caustic = caus(lens=plot_SSP_1, solver='SG12')
	#	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

	if plot2:

		plot_Rhie_2.plot_magnification(**magn_params)
		caustic = caus(lens=plot_Rhie_2, solver='SG12')
	#	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plot_SSP_2.plot_magnification(**magn_params)
		caustic = caus(lens=plot_SSP_2, solver='SG12')
	#	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

	if plot3:

		plot_Rhie_3.plot_magnification(**magn_params)
		caustic = caus(lens=plot_Rhie_3, solver='SG12')
	#	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()

		plot_SSP_3.plot_magnification(**magn_params)
		caustic = caus(lens=plot_SSP_3, solver='SG12')
	#	caustic.plot_caustic(points=10000, s=1, color='red', lw=0)
		plt.show()








