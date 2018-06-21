# Zoey Samples
# Created: Jun 11, 2018
# HDUL.py
# Last Updated: Jun 11, 2018

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
	"""

	def __init__(self, file_name):
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
		
	def check_kwargs(self, log_colorbar=False, **kwargs):
		"""
		Checks if the user specified keyword arguments for plots, and
		assigns default values for those not specified.

		Parameters:
			log_colorbar (bool):
				If True, kwargs will be assigned value for norm that displays
				the colormap in a logarithmic fashion.

			**kwargs (dictionary):
				Keyword arguments specified by the user.

		Returns:
			kwargs (dictionary):
				Keyword arguments after being assigned default values for
				those not specified by the user.
		"""

		if 's' not in kwargs:
			kwargs['s'] = (400. / self.res)**2
		if 'cmap' not in kwargs:
			kwargs['cmap'] = 'plasma'
		if 'linewidths' not in kwargs and 'lw' not in kwargs:
			kwargs['lw'] = None
		if 'marker' not in kwargs:
			kwargs['marker'] = 'o'
		if log_colorbar:
			kwargs['norm'] = colors.LogNorm()
		return kwargs
		
	def plot_magnification(self, cutoff=None, log_colorbar=False, 
						   outliers=False, save = False, **kwargs):
		"""
		Makes a plot of the magnification vs. position.
		"""

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.
		kwargs = self.check_kwargs(log_colorbar, **kwargs)

		if outliers:
			self.get_magnification_outliers(cutoff=cutoff)
			(x, y, magn) = (self.x_outliers, self.y_outliers, self.magn_outliers)
			print('Plotting the magnification of outliers...')
		else:
			(x, y, magn) = (self.x, self.y, self.magn)
			print('Plotting the magnification...')

		plt.scatter(x, y, c=magn, **kwargs)
		(xmin, xmax) = (min(self.x), max(self.x))
		(ymin, ymax) = (min(self.y), max(self.y))
		dx = xmax - xmin
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
		mag_plot = plt.colorbar()
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of Source', fontsize=12)
		plt.ylabel('Y-position of Source', fontsize=12)
		plt.gcf().set_size_inches(8, 6)

		if outliers:
			if cutoff == None:
				cutoff = int(min(magn))
			plt.suptitle('High Magnification', x=0.435)
			title = ('Frame: {}; Solver: {}\n'.format(
					self.origin_title, self.solver_title) + 
					's={}, q={}, M>{:.0f}'.format(self.s, self.q, 
					cutoff))
			plt.title(title, fontsize=11)
			file_name = ('../Tables/HighMagn_{}_{}.png'.format(
					self.solver_file, self.origin_file))
		else:
			plt.suptitle('Magnification', x=0.435)
			title = ('Frame: {}; Solver: {}\n'.format(
					self.origin_title, self.solver_title) + 
					's={}, q={}'.format(self.s, self.q))
			plt.title(title, fontsize=11)
			file_name = ('../Tables/Magn_{}_{}.png'.format(self.solver_file,
					self.origin_file))

		if save:
			self.save_png(file_name=file_name)

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

	def get_magnification_outliers(self, cutoff=None):
		"""
		Creates new arrays of (x, y, magn) only for magnification values that
		are above the cutoff value.

		Parameters:
			cutoff (float):
				The minimum magnification a point can have for it to be
				accepted as an 'outlier'.
		"""

		self.x_outliers = []
		self.y_outliers = []
		self.magn_outliers = []

		# If the cutoff value is not specified, the arrays will default
		# to the 10 largest points.
		if cutoff==None:
			magn_sorted = sorted(self.magn)
			cutoff = magn_sorted[int(0.95*(self.res**2))]
			print('No cutoff value specified; selecting only upper 5% of points')

		print('Finding the magnification outliers...')

		for (i, magn) in enumerate(self.magn):
			if magn > cutoff:
				self.x_outliers.append(self.x[i])
				self.y_outliers.append(self.y[i])
				self.magn_outliers.append(magn)

	def plot_num_images(self, errors_only=False, print_errors=True,
						save=False, **kwargs):
		"""
		Makes a plot of the number of solutions vs. position.
		"""

		print('Plotting the number of images at each point...')

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.

		kwargs = self.check_kwargs(**kwargs)
		kwargs['cmap'] = 'coolwarm'

		if print_errors:
			self.print_num_images_errors()

		if errors_only:
			(x, y, num_images) = self.get_num_images_errors()
			file_name = '../Tables/NumImErr_{}_{}.png'.format(
					self.origin_file, self.solver_file)
		else:
			(x, y, num_images) = (self.x, self.y, self.num_images)
			file_name = '../Tables/NumIm_{}_{}.png'.format(
					self.origin_file, self.solver_file)

		plt.scatter(x, y, c=num_images, **kwargs)
		im_plot = plt.colorbar()
		im_plot.set_label('Num Images')
		plt.xlabel('X-position of Source', fontsize=12)
		plt.ylabel('Y-position of Source', fontsize=12)
		plt.gcf().set_size_inches(8, 6)
		(xmin, xmax) = (min(self.x), max(self.x))
		(ymin, ymax) = (min(self.y), max(self.y))
		dx = xmax - xmin
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
		plt.suptitle('Number of Images', x=0.435)
		title = ('Frame: {}; Solver: {}\ns={}, q={}'.format(
				self.origin_title, self.solver_title, self.s, self.q))
		plt.title(title, fontsize=11)

		if save:
			self.save_png(file_name=file_name)

	def get_num_images_errors(self):
		"""
		Returns only the points in the grid where the number of images is
		not equal to 3 or 5.

		Returns:

			x_err (array of floats):
				The x-position for each erroneous point.

			y_err (array of floats):
				The y-position for each erroneous point.

			num_images_err (array of floats):
				The number of images for each erroneous point.
		"""

		x_err = []
		y_err = []
		num_images_err = []

		for (i, num) in enumerate(self.num_images):
			if (num == 0) or (num == 1) or (num == 2) or (num == 4):
				x_err.append(self.x[i])
				y_err.append(self.y[i])
				num_images_err.append(num)

		return (x_err, y_err, num_images_err)

	def print_num_images_errors(self):
		"""
		Prints the number of points that have n images, where n is an integer
		ranging from 0 to 5.
		"""

		num = np.zeros(6, dtype = int)
		for num_im in self.num_images:
			for i in range(len(num)):
				if num_im == i:
					num[i] += 1
		print('Concern: number of points where the number of images is',
			'\n0: {:}\n1: {:}\n2: {:}\n3: {:}\n4: {:}\n5: {:}\nTotal: {:}'
			.format(*num, sum(num)))

	def save_png(self, file_name):
		"""
		Saves a plot to a .png file.

		Required parameters:
			file_name (string):
				The name to which the file will be saved (exluding an integer
				before the .png extension, to denote which number of times
				the same file name has been saved). The file_name is
				specified within each plotting method, and the user does
				not need to type any information to obtain it.
		"""

		for i in range(10):
			name = file_name[:-4] + '{}'.format(i) + file_name[-4:]
			if Path(name).is_file():
				continue
			plt.savefig(name)
			print(name, 'has been saved')
			return
		print('Error: too many files of same name already exist. File not saved')

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
		elif self.origin == 'caustic':
			self.origin_file = 'caus'
			self.origin_title = 'Caustic'
			self.origin_phrase = 'caustic frame'
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(self.origin))

		if self.solver == 'numpy':
			self.solver_file = 'np'
			self.solver_title = 'Numpy'
			self.solver_phrase = 'numpy root finder'
		elif self.solver ==	'Skowron_and_Gould_12' or self.solver == 'SG12':
			self.solver_file = 'SG'
			self.solver_title = 'SG 2012'
			self.solver_phrase = 'Skowron and Gould 2012 root finder'
		elif self.solver == 'zroots':
			self.solver_file = 'zr'
			self.solver_title = 'Zroots'
			self.solver_phrase = 'Numerical Recipes: zroots root finder'
		else:
			raise ValueError('Unknown solver: {:}'.format(self.solver))
