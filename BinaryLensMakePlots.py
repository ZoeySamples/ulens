# Zoey Samples
# Created: May 30, 2018
# BinaryLensMakePlots.py
# Last Updated: May 31, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm
import BinaryLensFunctions as blf
from astropy.io import fits

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
	x_grid = np.linspace(x_cent - w_caustic, x_cent + w_caustic, pts)
	y_grid = np.linspace(-h_caustic, h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	im_num = np.zeros(pts**2, dtype=int)
	color = np.zeros(pts**2)
	num_one = 0
	num_two = 0
	num_four = 0
	num_tot = 0
	print('parameters:\ns={:}\nq={:}\n'.format(s, q))
	for (i, xx) in enumerate(x_grid):
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
			elif im_num[idx] == 1:
				num_one += 1
			else:
				print('Concern: found {:} solutions'.format(im_num[idx]))
	num_tot = num_one + num_two + num_four
	print('Concern: number of points where the number of images is',
		'\n1: {:}\n2: {:}\n4: {:}\ntotal: {:} out of {:} points'
		.format(num_one, num_two, num_four, num_tot, pts**2))
	plt.scatter(x_1d, y_1d, c=im_num, s=((500./pts)**2), marker = 'o', cmap='jet', lw=None)
	im_plot = plt.colorbar()
	im_plot.set_label('Num Images')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x_cent - w_caustic, x_cent + w_caustic)
	plt.ylim(-h_caustic, h_caustic)
	plt.title('Num Images using "{}" frame'.format(origin))

def plot_magnification(s, q, origin = 'geo_cent', solver='numpy', pts=150):
	"""Make square grid of points that shows the magnification at each point"""
	(w_caustic, h_caustic, x_cent) = size_caustic(s, q, origin)
	x_grid = np.linspace(x_cent - w_caustic, x_cent + w_caustic, pts)
	y_grid = np.linspace(-h_caustic, h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	mag_1d = np.zeros(pts**2)
	for i, xx in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			mag_1d[idx] = blf.magnification(x=xx, y=yy, s=s, q=q, origin=origin, solver=solver)
			x_1d[idx] = xx
			y_1d[idx] = yy
	plt.scatter(x_1d, y_1d, c=mag_1d, s=((500./pts)**2), marker = 'o', cmap='jet', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x_cent - w_caustic, x_cent + w_caustic)
	plt.ylim(-h_caustic, h_caustic)
	plt.title('Magnification using "{}" frame'.format(origin))

def write_to_fits(s, q, origin = 'geo_cent', solver='numpy', pts=500):
	"""
	Writes information about grid to a .fits table for comparison of magnification
	and number of images between different coordinate systems and solving methods.
	"""
	(w_caustic, h_caustic, x_cent) = size_caustic(s, q, origin)
	x_grid = np.linspace(x_cent - w_caustic, x_cent + w_caustic, pts)
	y_grid = np.linspace(-h_caustic, h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	mag_1d = np.zeros(pts**2)
	im_num = np.zeros(pts**2, dtype=int)
	print('parameters:\ns={:}\nq={:}\n'.format(s, q))
	for (i, xx) in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			mag_1d[idx] = blf.magnification(x=xx, y=yy, s=s, q=q, origin=origin,
											solver=solver)
			x_1d[idx] = xx
			y_1d[idx] = yy
			(dm, m, zeta, z1, z2) = blf.assign(xx, yy, s, q, origin)
			solutions = blf.solution(xx, yy, s, q, origin, solver)
			for z in solutions:
				if blf.check_solution(dm, m, zeta, z1, z2, z, origin):
					im_num[idx] += 1
	col = []
	col.append(fits.Column(name='x', array=x_1d, format='D'))
	col.append(fits.Column(name='y', array=y_1d, format='D'))
	col.append(fits.Column(name='Magnification', array=mag_1d, format='D'))
	col.append(fits.Column(name='Number Images', array=im_num, format='I'))
	hdu1 = fits.BinTableHDU.from_columns(col)
	hdr = fits.Header()
	hdr['SEPARAT'] = '{:f}'.format(s)
	hdr['M_RATIO'] = '{:f}'.format(q)
	hdr['ORIGIN'] = origin
	hdr['SOLVER'] = solver
	hdu0 = fits.PrimaryHDU(header = hdr)
	hdus = fits.HDUList([hdu0, hdu1])
	hdus.writeto('BinaryLens2.fits')

