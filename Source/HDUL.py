# Zoey Samples
# Created: Jun 11, 2018
# HDUL.py
# Last Updated: Jun 11, 2018

import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
import MulensModel as mm


class HDUL(object):

	"""
	Reads data from fits files, makes plots, assesses data for "outliers,"
	compares multiple sets of data, etc.

	Attributes:

		Required:
			file_name (string):
				The name of the fits file.

		Optional:
			cutoff (string):
				The minimum magnification value accepted when assessing "outliers."
	"""

	def __init__(self, file_name, cutoff=None):
		self.cutoff = cutoff
		self.file_name = file_name
		self.get_variables()
		self.strings()
		self.plot_limits()

	def get_variables(self):
		"""Assigns data to arrays."""

		self.hdul = fits.open(self.file_name)
		self.x = self.hdul[1].data['x']
		self.y = self.hdul[1].data['y']
		self.s = float(self.hdul[0].header['SEPARAT'])
		self.q = float(self.hdul[0].header['M_RATIO'])
		self.num_images = self.hdul[1].data['Number Images']
		self.magn = self.hdul[1].data['Magnification']
		self.res = float(self.hdul[0].header['RES'])
		self.origin = self.hdul[0].header['ORIGIN']
		self.solver = self.hdul[0].header['SOLVER']

	def plot_limits(self):
		"""Assigns the plot limits for the x- and y-axis."""

		self.xmin = self.x[0]
		self.xmax = self.x[-1]
		self.ymin = self.y[0]
		self.ymax = self.y[-1]

	def check_compatible(self, other_hdul):
		"""
		Determines if tables have the same size and are plotted over the same grid.
		"""

		(x1, y1, res1) = (self.x, self.y, self.res)
		(x2, y2, res2) = (other_hdul.x, other_hdul.y, other_hdul.res)

		# Here, we check if the resolution for each grid is the same, and print
		# an error if they aren't.
		print('Checking if hdu lists are compatible...')
		if int(res1) != int(res2):
			sys.exit('Error: Resolutions are not the same.')
		z1 = x1 + y1*1.j
		z2 = x2 + y2*1.j
		for i in range(len(x1)):
			if np.abs(z1[i] - z2[i]) > 1e-10:
				sys.exit('Error: Grids are not formed over the same region')
		print('Test passed for all points.')

	def plot_rel_magn(self, other_hdul, save = False, outliers = False):
		"""
		Plots the fractional difference in magnification between two sets of data
		"""

		self.check_compatible(other_hdul = other_hdul)
		print('Plotting the relative magnification...')

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.
		if outliers:
			self.get_magn_outliers()
			other_hdul.get_magn_outliers()
			(x, y, magn1, magn2) = (self.x_outliers, self.y_outliers,
							self.magn_outliers, other_hdul.magn_outliers)

		else:
			(x, y, magn1, magn2) = (self.x, self.y, self.magn, other_hdul.magn)

		rel_magn = (magn1 / magn2)
		plt.scatter(x, y, c=rel_magn, s=((800./self.res)**2), marker = 'o',
						cmap='plasma', lw=None)
		rel_plot = plt.colorbar()
		rel_plot.set_label('Fractional Difference')
		plt.xlabel('X-position of source', fontsize = 14)
		plt.ylabel('Y-position of source', fontsize = 14)
		plt.xlim(self.xmin, self.xmax)
		plt.ylim(self.ymin, self.ymax)
		plt.title('Relative Magnification')

	def plot_magnification(self, save = False, outliers = False):
		"""
		Makes a plot of the magnification vs. position.
		"""

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.
		if outliers:
			(x, y, magn) = (self.x_outliers, self.y_outliers, self.magn_outliers)
			print('Plotting the magnification of outliers...')
		else:
			(x, y, magn) = (self.x, self.y, self.magn)
			print('Plotting the magnification...')

		plt.scatter(x, y, c=magn, s=7, marker = '^', cmap='plasma', lw=None)
		mag_plot = plt.colorbar()
		mag_plot.ax.tick_params(labelsize=10)
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of source', fontsize = 12)
		plt.ylabel('Y-position of source', fontsize = 12)
		plt.xlim(self.xmin, self.xmax)
		plt.ylim(self.ymin, self.ymax)
		if not outliers:
			plt.title('Magnification\n{} Frame; {} Solver'.format(
								self.origin_title, self.solver_title))

	def plot_magn_outliers(self, save = False):
		"""
		Plots the magnification at each position where the magnification is
		greater than the specified cutoff value.
		"""

		self.get_magn_outliers()
		self.plot_magnification(save=save, outliers=True)

		caustics = mm.Caustics(s=(self.s), q=(self.q))
		caustics.plot(s=5, c='red')

		plt.gcf().set_size_inches(8, 6)
		plt.title('High Magnification with Caustic\nFrame: {}; Solver = {}; M > {}'.
						format(self.origin_title, self.solver_title, self.cutoff))
		if save:
			plt.savefig('../Tables/high_mag_{}_{}'.format(self.solver_file,
									self.origin_file))

	def get_magn_outliers(self):
		"""
		Creates new arrays of (x, y, magn) only for magnification values that
		are above the cutoff value.
		"""

		self.x_outliers = []
		self.y_outliers = []
		self.magn_outliers = []

		# If the cutoff value is not specified, default to the 90th percentile
		# of magnification.
		if self.cutoff==None:
			magn_sorted = sorted(magn_outliers)
			self.cutoff = magn_sorted[0.9 * self.res]
			print('No cutoff value specified; defaulting to 90th percentile.')

		print('Finding the magnification outliers...')

		for (i, magn) in enumerate(self.magn):
			if magn > self.cutoff:
				self.x_outliers.append(self.x[i])
				self.y_outliers.append(self.y[i])
				self.magn_outliers.append(magn)

	def plot_n_solutions(self, save = False, outliers = False):
		"""
		Makes a plot of the number of solutions vs. position.
		"""

		print('Plotting the number of images at each point...')

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.
		if outliers:
			(x, y, num_images) = (self.x_errors, self.y_errors,
											self.num_images_errors)
		else:
			(x, y, num_images) = (self.x, self.y, self.num_images)

		plt.scatter(x, y, c=num_images, s=1,
							marker = 'o', cmap='copper', lw=None)
		soln_plot = plt.colorbar()
		soln_plot.set_label('Num Images')

		plt.xlabel('X-position of source', fontsize = 12)
		plt.ylabel('Y-position of source', fontsize = 12)
		plt.xlim(self.xmin, self.xmax)
		plt.ylim(self.ymin, self.ymax)
		if not outliers:
			plt.title('Num Images\n{} Frame; {} Solver'.format(self.origin_title,
								 self.solver_title))

	def plot_n_solutions_errors(self, save = False):
		"""
		Plots the number of images at each position where the number of images is
		not equal to 3 or 5.
		"""

		self.get_num_images_errors()
		self.plot_n_solutions(save = save, outliers = True)

		caustics = mm.Caustics(s=(self.s), q=(self.q))
		caustics.plot(s=5, c='red')

		plt.gcf().set_size_inches(8, 6)
		plt.title('Erroneous Num Images with Caustic\nFrame: {}; Solver = {}; M > {}'.
					format(self.origin_title, self.solver_title, self.cutoff))
		if save:
			plt.savefig('../Tables/high_mag_{}_{}'.format(self.solver_file,
									self.origin_file))

	def get_num_images_errors(self):
		"""
		Plots where the number of images is not equal to what we expect; i.e.
		where num_images != (3 or 5)
		"""

		print('Finding erroneous numbers of images...')

		self.x_errors = []
		self.y_errors = []
		self.num_images_errors = []

		for (i, num) in enumerate(self.num_images):
			if (num != 3) and (num != 5):
				self.x_errors.append(self.x[i])
				self.y_errors.append(self.y[i])
				self.num_images_errors.append(num)

	def write_to_fits(self, outliers=False):
		"""
		Writes new table to fits file.
		"""

		if outliers:
			self.get_magn_outliers()
			x = self.x_outliers
			y = self.y_outliers
			magn = self.magn_outliers
		else:
			x = self.x
			y = self.y
			z = self.z

		print('Writing data to fits file...')

		col = []
		col.append(fits.Column(name='x', array = x, format='D'))
		col.append(fits.Column(name='y', array = y, format='D'))
		col.append(fits.Column(name='Magnification', array = magn, format='D'))
		hdu1 = fits.BinTableHDU.from_columns(col)
		hdr = fits.Header()
		hdr['SEPARAT'] = self.s
		hdr['M_RATIO'] = self.q
		hdr['ORIGIN'] = self.origin
		hdr['SOLVER'] = self.solver
		hdr['RES'] = '{:d}'.format(self.res)
		hdu0 = fits.PrimaryHDU(header = hdr)
		hdus = fits.HDUList([hdu0, hdu1])
		hdus.writeto('../Tables/BL_extreme_{}_{:}.fits'.format(
							self.origin_file, self.solver_file))

	def strings(self):
		"""
		Assign appropriate strings for file-naming, printing, and plot titles.
		"""

		if self.origin == 'geo_cent':
			self.origin_file = 'gcent'
			self.origin_title = 'Geo Center'
			self.origin_phrase = 'geometric center frame'
		elif self.origin == 'plan':
			self.origin_file = self.origin
			self.origin_title = 'Planet'
			self.origin_phrase = 'planet frame'
		elif self.origin == 'com':
			self.origin_file = self.origin
			self.origin_title = 'Center Mass'
			self.origin_phrase = 'center-of-mass frame'
		elif self.origin == 'star':
			self.origin_file = self.origin
			self.origin_title = 'Star'
			self.origin_phrase = 'star frame'
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(self.origin))

		if self.solver == 'numpy':
			self.solver_file = 'np'
			self.solver_title = 'Numpy'
			self.solver_phrase = 'numpy root finder'
	#	elif self.solver == 'SG12':
		elif self.solver ==	'Skowron_and_Gould_12':
			self.solver_file = 'SG'
			self.solver_title = 'SG 2012'
			self.solver_phrase = 'Skowron and Gould 2012 root finder'
		elif self.solver == 'zroots':
			self.solver_file = 'zr'
			self.solver_title = 'Zroots'
			self.solver_phrase = 'Numerical Recipes: zroots root finder'
		else:
			raise ValueError('Unknown solver: {:}'.format(self.solver))
