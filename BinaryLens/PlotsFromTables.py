# Zoey Samples
# Created: Jun 04, 2018
# BinaryLensMakePlots.py
# Last Updated: Jun 04, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
import Functions as blf
from astropy.io import fits

fits_table_filename = '../Tables/BL_plan_SG.fits'
hdul = fits.open(fits_table_filename)

def get_variables(hdul):
	x = hdul[1].data['x']
	y = hdul[1].data['y']
	num_images = hdul[1].data['Number Images']
	magn = hdul[1].data['Magnification']
	res = float(hdul[0].header['RES'])
	origin = hdul[0].header['ORIGIN']
	return x, y, num_images, magn, res, origin

def plot_magnification(x, y, num_images, magn, res, origin):
	plt.scatter(x, y, c=magn, s=((800./res)**2), marker = 'o', cmap='jet', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x[0], x[-1])
	plt.ylim(y[0], y[-1])
	plt.title('Magnification using "{}" frame'.format(origin))

def plot_n_solutions(x, y, num_images, magn, res, origin):
	plt.scatter(x, y, c=num_images, s=((800./res)**2), marker = 'o', cmap='jet', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.set_label('Num Images')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x[0], x[-1])
	plt.ylim(y[0], y[-1])
	plt.title('Magnification using "{}" frame'.format(origin))

plot_list = get_variables(hdul)
plot_n_solutions(*plot_list)
plt.show()
plot_magnification(*plot_list)
plt.show()
