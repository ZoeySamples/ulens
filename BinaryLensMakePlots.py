# Zoey Samples
# Created: May 30, 2018
# BinaryLensMakePlots.py
# Last Updated: May 30, 2018; 11:35AM

import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm
import BinaryLensFunctions as blf

def size_caustic(s, q, origin = 'geo_cent'):
	w = 4.*np.sqrt(q)*(1. + 1./(2.*(s**2))) / (s**2)
	h = 4.*np.sqrt(q)*(1. - 1./(2.*(s**2))) / (s**2)
	x = 0.5*s - 1.0/s
	return w, h, x

	"""
	if origin == 'geo_cent':
		x = 0.5*s - 1.0/s
	elif origin == 'star':
		x = s - 1.0/s
	elif origin == 'plan':
		x = -1.0/s
	elif origin == 'com':
		x = 2.0*s/(1.0 + q) - 1.0/s
	else:
		print('Not a valid coordinate system. Converting to geometric center frame')
		x = 0.5*s - 1.0/s
	"""

def plot_n_solns(s, q, origin = 'geo_cent', solver='numpy', pts=150):
	"""
	Plot showing number of images on a grid (x,y) around planetary caustic

	"""
	(w_caustic, h_caustic, x_cent) = size_caustic(s, q, origin)
	x_factor = 15.
	y_factor = 15.
	x_grid = np.linspace(x_cent - x_factor*w_caustic, x_cent + x_factor*w_caustic, pts)
	y_grid = np.linspace(-y_factor*h_caustic, y_factor*h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	im_num = np.zeros(pts**2, dtype=int)
	color = np.zeros(pts**2)
	num_two = 0
	num_four = 0
	num_other = 0
	print('parameters:\ns={:}\nq={:}'.format(s, q))
	for i, xx in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			x_1d[idx] = xx
			y_1d[idx] = yy
			(dm, m, zeta, z1, z2) = blf.assign(xx, yy, s, q, origin)
			solutions = blf.solution(xx, yy, s, q, origin, solver)
			for z in solutions:
				if blf.check_solution(dm, m, zeta, z1, z2, z, origin):
					im_num[idx] += 1
			if im_num[idx] == 5:
				color[idx] = 255
			elif im_num[idx] == 3:
				color[idx] = 1
			elif im_num[idx] == 4:
				color[idx] = 120
				num_four += 1
			elif im_num[idx] == 2:
				color[idx] = 120
				num_two += 1
			else:
				num_other += 1
				print('Concern: number of images=', im_num[idx])
				print('x_source={:}\ny_source={:}'.format(xx, yy))
	print('Number of points where the number of images is not 3 or 5 is',
		num_four + num_two + num_other,'out of', pts**2)
	plt.scatter(x_1d, y_1d, c=im_num, s=8, cmap='jet')
	im_plot = plt.colorbar()
	im_plot.set_label('Num Images')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x_cent - x_factor*w_caustic, x_cent + x_factor*w_caustic)
	plt.ylim(-y_factor*h_caustic, y_factor*h_caustic)
	plt.title('Num Images using "{}" frame'.format(origin))

def plot_magnification(s, q, origin = 'geo_cent', solver='numpy', pts=150):
	"""Make square grid of points that shows the magnification at each point"""
	(w_caustic, h_caustic, x_cent) = size_caustic(s, q, origin)
	x_factor = 15.
	y_factor = 15.
	x_grid = np.linspace(x_cent - x_factor*w_caustic, x_cent + x_factor*w_caustic, pts)
	y_grid = np.linspace(-y_factor*h_caustic, y_factor*h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	mag_1d = np.zeros(pts**2)
	for i, xx in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			mag_1d[idx] = blf.magnification(x=xx, y=yy, s=s, q=q, origin=origin, solver=solver)
			x_1d[idx] = xx
			y_1d[idx] = yy
	plt.scatter(x_1d, y_1d, c=mag_1d, s=8, cmap='jet')
	mag_plot = plt.colorbar()
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x_cent - x_factor*w_caustic, x_cent + x_factor*w_caustic)
	plt.ylim(-y_factor*h_caustic, y_factor*h_caustic)
	plt.title('Magnification using "{}" frame'.format(origin))
