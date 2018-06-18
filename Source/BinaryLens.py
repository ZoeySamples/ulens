# Zoey Samples
# Created: June 06, 2018
# BinaryLens.py
# Last Updated: June 13, 2018

import sys
import os
from pathlib import Path
import ctypes
import numpy as np
import cmath
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import MulensModel as mm
from astropy.io import fits

MODULE_PATH = os.path.abspath(__file__)
for i in range(3):
	MODULE_PATH = os.path.dirname(MODULE_PATH)
PATH = os.path.join(MODULE_PATH, 'MulensModel-master', 'source', 'VBBL',
		"VBBinaryLensingLibrary_wrapper.so")

# Here we attempt to access the Skowron & Gould 2012 root finder
try:
	vbbl = ctypes.cdll.LoadLibrary(PATH)
except OSError as error:
	msg = "Something went wrong with VBBL wrapping ({:})\n\n" + repr(error)
	print(msg.format(PATH))
else:
	vbbl.VBBL_SG12_5.argtypes = 12 * [ctypes.c_double]
	vbbl.VBBL_SG12_5.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, 
		shape=(10,))
	_vbbl_SG12_5 = vbbl.VBBL_SG12_5

MODULE_PATH = os.path.abspath(__file__)
for i in range(2):
	MODULE_PATH = os.path.dirname(MODULE_PATH)
PATH = os.path.join(MODULE_PATH, 'NumericalRecipes', "zrootsBinaryLens_wrapper.so")

# Here we attempt to access the zroots root finder
try:
	zroots = ctypes.cdll.LoadLibrary(PATH)
except OSError as error:
	msg = "Something went wrong with zroots wrapping ({:})\n\n" + repr(error)
	print(msg.format(PATH))
else:
	zroots.zroots_5.argtypes = 12 * [ctypes.c_double]
	zroots.zroots_5.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, 
		shape=(10,))
	_zroots_5 = zroots.zroots_5


class BinaryLens(object):

	"""
	Using the binary lens equation: Calculates solutions, magnification, makes
	plots, and writes to files. Works for point test runs or filling a grid of 
	points centered on the planetary caustic.

	Attributes:

		Required:
			s (float):
				The separation between the two bodies.

			q (float):
				The mass ratio of the smaller body to the larger body.

			origin (string):
				The coordinate frame in which calculations are carried out.
				Options are:
					'geo_cent'	- the geometric center frame
					'star'		- the star's (or larger body's) frame
					'plan'		- the planet's (or smaller body's) frame
					'com'		- the center-of-mass frame
					'caustic'	- the planetary caustic frame

			solver (string):
				Determines which root finder to use to solve the polynomial.
				Options are:
					'numpy'		- Uses np.roots method
					'SG12'		- Uses Skowron & Gould 2012 method
					'zroots'	- Uses zroots laguerre method

		Optional:
			x (float):
				The x-position of the source. Required if user is not plotting
				a grid with BinaryLens class.

			y (float):
				The y-position of the source. Required if user is not plotting
				a grid with BinaryLens class.

			res (float):
				The resolution of the grid (i.e. number of points on each side). 
				Required if user is plotting a grid with BinaryLens class.

			tolerance (float):
				The maximum distance away a calculated image position must be from
				an actual solution when substitued back into the binary lens
				equation, in order not to be rejected.

		Methodology:
			To find solutions of test parameters, follow these steps:
				1. Write a dictionary to call BinaryLens with, including 
					parameters:	{x, y, s, q, origin, solver, tolerance (opt)}.
				2. Call any of the following functions:
					image_positions()
					print_image_position(print_input)
					print_magnification(print_input)

			To make a plot with a grid of points or to write a table to a .fits
			file, follow these steps:
				1. Write a dictionary to call BinaryLens with, including 
					parameters:	{s, q, origin, solver, res, tolerance (opt)}.
				2. Call any of the following functions:
					plot_n_solns(save, print_errors)
					plot_magnification(save)
					write_to_fits():
	"""

	def __init__(self, s, q, origin, solver, x=None, y=None, res=None,
				 tolerance=0.0001):
		self.s = s
		self.q = q
		self.origin = origin
		self.solver = solver
		self.res = res
		self.x = x
		self.y = y
		self.tolerance = tolerance
		self.strings()

### The following functions assign values to variables pertaining to the class.

	def get_variables(self, x, y):
		"""
		Assign all variables to be used in polynomial. Assigns:
	
			dm (float):
				Half the difference of the masses (positive).

			m (float):
				Half the sum of the masses.

			zeta (complex):
				The source's position in complex space.

			zeta_conj (complex):
				The complex conjugate of the source's position.

			z1 (float):
				The planet's position (on the real axis by definition).

			z2 (float):
				The star's position (on the real axis by definition).
		"""
	
		self.mass_ratio()
		self.source_position(x=x, y=y)
		self.lensing_body_positions()

	def mass_ratio(self):
		"""Define m and dm."""

		self.dm = (1. - self.q) / 2.
		self.m = (1. + self.q) / 2.

	def source_position(self, x, y):
		"""Define zeta and zeta_conj."""

		if self.origin == 'geo_cent':
			self.zeta = x + y*1.j
		elif self.origin == 'star':
			self.zeta = (x + self.s/2.) + y*1.j
		elif self.origin == 'plan':
			self.zeta = (x - self.s/2.) + y*1.j
		elif self.origin == 'com':
			self.zeta = (x + self.s*self.dm/(2.*self.m)) + y*1.j
		elif self.origin == 'caustic':
			self.zeta = (x + 1./(self.s) - self.s/2.) + y*1.j
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))
		self.zeta_conj = self.zeta.conjugate()

	def lensing_body_positions(self):
		"""Define z1 and z2."""

		if self.origin == 'geo_cent':
			self.z1 = 0.5*self.s
			self.z2 = -0.5*self.s
		elif self.origin == 'star':
			self.z1 = self.s
			self.z2 = 0.
		elif self.origin == 'plan':
			self.z1 = 0.
			self.z2 = -self.s
		elif self.origin == 'com':
			self.z1 = self.s*(self.m + self.dm) / (2.*self.m)
			self.z2 = -self.s*(self.m - self.dm) / (2.*self.m)
		elif self.origin == 'caustic':
			self.z1 = 1./self.s
			self.z2 = -self.s + 1./self.s
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))

		self.z1_conj = self.z1.conjugate()
		self.z2_conj = self.z2.conjugate()

### The following functions calculate physical values for further analysis.

	def get_coefficients(self, x, y):
		"""Returns the coefficients for the polynomial equation."""

		self.get_variables(x=x, y=y)

		# Assign the values of the coefficients of the polynomial.
		# Be aware that these expressions are very long and messy. They trail off
		# the screen without text wrapping and look non-indented with text wrapping.

		"""
		coeff5 = (self.z1**2 - self.zeta_conj**2)

		coeff4 = (-(self.z1*(2*self.dm + self.z1*self.zeta)) - 2*self.m*self.zeta_conj + self.zeta*self.zeta_conj**2)

		coeff3 = (-2*self.z1**4 + 4*(self.dm*self.z1 + self.m*self.zeta)*self.zeta_conj + 2*self.z1**2*self.zeta_conj**2)

		coeff2 = (4*self.dm*self.z1*(self.m + self.z1**2) + 2*(2*self.m**2 + self.z1**4)*self.zeta - 4*self.dm*self.z1*self.zeta*self.zeta_conj - 2*self.z1**2*self.zeta*self.zeta_conj**2)

		coeff1 = self.z1*(-4*self.dm**2*self.z1 - 4*self.m**2*self.z1 + self.z1**5 - 8*self.dm*self.m*self.zeta - 4*self.z1*(self.dm*self.z1 + self.m*self.zeta)*self.zeta_conj - self.z1**3*self.zeta_conj**2)

		coeff0 = self.z1**2*(4*self.dm*self.m*self.z1 - 2*self.dm*self.z1**3 + 4*self.dm**2*self.zeta - self.z1**4*self.zeta + 2*self.z1*(self.m*self.z1 + 2*self.dm*self.zeta)*self.zeta_conj + self.z1**2*self.zeta*self.zeta_conj**2)
		"""

		coeff5 = (-self.zeta_conj + self.z1)*(self.zeta_conj- self.z2)

		coeff4 = (self.m*self.z1 + self.m*self.z2 + 2.*(self.z1**2)*self.z2 + 2.*self.z1*(self.z2**2) + self.dm*(-self.z1 + self.z2) + self.z1*self.z2*self.zeta + (self.zeta_conj**2)*(2.*self.z1 + 2.*self.z2 + self.zeta) - self.zeta_conj*(2.*self.m + (self.z1 + self.z2)*(2.*self.z1 + 2.*self.z2 + self.zeta)))

		coeff3 = (self.dm*(self.z1**2) - self.m*(self.z1**2) - 2.*self.m*self.z1*self.z2 - (self.z1**3)*self.z2 - self.dm*(self.z2**2) - self.m*(self.z2**2) - 4.*(self.z1**2)*(self.z2**2) - self.z1*(self.z2**3) - 2.*self.m*self.z1*self.zeta - 2.*self.m*self.z2*self.zeta - 2.*(self.z1**2)*self.z2*self.zeta - 2.*self.z1*(self.z2**2)*self.zeta - (self.zeta_conj**2)*((self.z1**2) + 2.*self.z1*(2.*self.z2 + self.zeta) + self.z2*(self.z2 + 2.*self.zeta)) + self.zeta_conj*(2.*self.dm*(self.z1 - self.z2) + 2.*self.m*(self.z1 + self.z2 + 2.*self.zeta) + (self.z1 + self.z2)*((self.z1**2) + 4.*self.z1*self.z2 + (self.z2**2) + 2.*self.z1*self.zeta + 2.*self.z2*self.zeta)))

		coeff2 = (-2.*(self.m**2)*(self.z1 + self.z2 - 2.*self.zeta) - 3.*self.m*(2.*self.zeta_conj - self.z1 - self.z2)*(self.z1 + self.z2)*self.zeta + self.dm*(self.z1 - self.z2)*(2.*self.m - 2.*self.z1*self.z2 + self.z1*self.zeta + self.z2*self.zeta - 2.*self.zeta_conj*(self.z1 + self.z2 + self.zeta)) + (self.zeta_conj - self.z1)*(self.zeta_conj - self.z2)*((self.z2**2)*self.zeta + (self.z1**2)*(2.*self.z2 + self.zeta) + 2.*self.z1*self.z2*(self.z2 + 2.*self.zeta)))

		coeff1 = ((-self.dm**2)*((self.z1 - self.z2)**2) + (self.m**2)*((self.z1**2) + 6.*self.z1*self.z2 + (self.z2**2) - 4.*self.z1*self.zeta - 4.*self.z2*self.zeta) - self.m*(2.*self.zeta_conj - self.z1 - self.z2)*(self.z1*self.z2*(self.z2 - 4.*self.zeta) + (self.z1**2)*(self.z2 - self.zeta) - (self.z2**2)*self.zeta) - (self.zeta_conj - self.z1)*self.z1*(self.zeta_conj - self.z2)*self.z2*(2.*self.z2*self.zeta + self.z1*(self.z2 + 2.*self.zeta)) + self.dm*(self.z1 - self.z2)*(self.z1*self.z2*(self.z2 - 2.*self.zeta) + (self.z1**2)*(self.z2 - self.zeta) - (4.*self.m + (self.z2**2))*self.zeta + 2.*self.zeta_conj*(self.z2*self.zeta + self.z1*(self.z2 + self.zeta))))

		coeff0 = (-2.*(self.m**2)*(self.z1**2)*self.z2 - 2.*(self.m**2)*self.z1*(self.z2**2) - self.m*(self.z1**3)*(self.z2**2) - self.m*(self.z1**2)*(self.z2**3) + (self.m**2)*(self.z1**2)*self.zeta + (self.dm**2)*((self.z1 - self.z2)**2)*self.zeta + 2.*(self.m**2)*self.z1*self.z2*self.zeta + self.m*(self.z1**3)*self.z2*self.zeta + (self.m**2)*(self.z2**2)*self.zeta + (self.zeta_conj**2)*(self.z1**2)*(self.z2**2)*self.zeta + 2.*self.m*(self.z1**2)*(self.z2**2)*self.zeta + self.m*self.z1*(self.z2**3)*self.zeta + (self.z1**3)*(self.z2**3)*self.zeta - self.dm*(self.z1 - self.z2)*(2.*self.m + self.z1*self.z2)*(self.z1*(self.z2 - self.zeta) - self.z2*self.zeta) - self.zeta_conj*self.z1*self.z2*((2.*self.dm*(self.z1 - self.z2) + self.z1*self.z2*(self.z1 + self.z2))*self.zeta + self.m*(-2.*self.z1*self.z2 + 2.*self.z1*self.zeta + 2.*self.z2*self.zeta)))

		coeff_list = np.array([coeff5, coeff4, coeff3, coeff2, coeff1, coeff0])
		return coeff_list

	def get_solutions(self, x, y):
		"""Return solutions of polynomial."""

		coeff_list = self.get_coefficients(x=x, y=y)

		# Return the roots of the polynomial via the given root finder
		if self.solver == 'SG12':
			rev_list = coeff_list[::-1]
			out = _vbbl_SG12_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return roots
		elif self.solver == 'zroots':
			rev_list = coeff_list[::-1]
			out = _zroots_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return roots
		elif self.solver == 'numpy':
			roots = np.roots(coeff_list).tolist()
			return roots

	def check_solution(self, solution):
		"""
		Check if the determined solution is consistent with the binary 
		lens equation.
		"""

		z = solution

		# This is the binary lens equation.
		zeta_actual = (z + (self.m - self.dm) / (self.z1_conj - z.conjugate()) +
					  (self.m + self.dm) / (self.z2_conj - z.conjugate()))
		if np.abs(self.zeta - zeta_actual) > self.tolerance:
			return False
		else:
			return True

	def image_positions(self):
		"""
		Calculates the image positions (i.e. checks which solutions pass the
		check). Returns a list of the positions.
		"""

		roots = self.get_solutions(x=self.x, y=self.y)
		image_positions = []
		for solution in roots:
			if self.check_solution(solution=solution):
				image_positions.append(solution)
		return image_positions

	def get_magnification(self, x, y):
		"""Returns the magnification for each configuration."""

		roots = self.get_solutions(x=x, y=y)
		magn = list(range(5))
		for (i, z) in enumerate(roots):
			detJ = (1. - ((self.m - self.dm) / ((z - self.z1)**2) + (self.m +
					self.dm) / ((z - self.z2)**2)) * ((self.m - self.dm) /
					((z.conjugate() - self.z1)**2) + (self.m + self.dm) / 
					((z.conjugate() - self.z2)**2)))
			if self.check_solution(solution=z):
				magn[i] = np.abs(1./detJ)
			else:
				magn[i] = 0.

		# This is the sum of the calculated magnitude after removing 
		# non-physical results.
		return sum(magn)

	def get_size_caustic(self):
		"""
		Determines the width, height, and position of the center of the
		caustic.
		"""

		self.width_caustic = 4.*np.sqrt(self.q)*(1. + 1./(2.*(self.s**2))) / (self.s**2)
		self.height_caustic = 4.*np.sqrt(self.q)*(1. - 1./(2.*(self.s**2))) / (self.s**2)
		self.xcenter_caustic = 0.5*self.s - 1.0/self.s

### The following functions are used for assigning data into grids.

	def get_position_arrays(self, region):
		"""
		Fills arrays for the x- and y-position to prepare grid plots.

		Parameters:
			region (string):
				The area on which the grid will be calculated.
				Accepted values:

				'caustic' - the zoomed-out view showing the full
						planetary caustic.
				'onax-cusp' - the zoomed-in view to the right of the cusp
						on the horizontal axis.
				'offax-cusp' - the zoomed-in view above the cusp on the
						vertical axis.
		"""

		self.get_size_caustic()

		if region == 'caustic':
			region_xmin = self.xcenter_caustic - 0.8*self.width_caustic
			region_xmax = self.xcenter_caustic + 0.8*self.width_caustic
			region_ymin = -0.8*self.height_caustic
			region_ymax = 0.8*self.height_caustic
		if region == 'onax_cusp':
			region_xmin = self.xcenter_caustic + 0.55*self.width_caustic
			region_xmax = self.xcenter_caustic + 0.8*self.width_caustic
			region_ymin = -0.10*self.height_caustic
			region_ymax = 0.10*self.height_caustic
		if region == 'offax_cusp':
			region_xmin = self.xcenter_caustic - 0.10*self.width_caustic
			region_xmax = self.xcenter_caustic + 0.10*self.width_caustic
			region_ymin = 0.55*self.height_caustic
			region_ymax = 0.8*self.height_caustic

		x_grid = np.linspace(region_xmin, region_xmax, self.res)
		y_grid = np.linspace(region_ymin, region_ymax, self.res)
		self.x_array = np.zeros(self.res**2)
		self.y_array = np.zeros(self.res**2)

		for (i, xx) in enumerate(x_grid):
			for (j, yy) in enumerate(y_grid):
				idx = self.res*i + j
				self.x_array[idx] = xx
				self.y_array[idx] = yy

	def get_magnification_array(self):
		"""Fills an array for the magnification through the grid."""

		self.magn_array = np.zeros(self.res**2, dtype=float)
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			self.magn_array[idx] = self.get_magnification(x=x, y=y)

	def get_num_images_array(self):
		"""Fills an array for the number of images through the grid."""

		self.num_images = np.zeros(self.res**2, dtype=int)
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			roots = self.get_solutions(x=x, y=y)
			for solution in roots:
				if self.check_solution(solution=solution):
					self.num_images[idx] += 1

	def get_coeff_array(self):
		"""
		Fills an array for the values of the coefficients through the grid.
		"""

		self.coeff_array = [[]*self.res**2 for i in range(12)]
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			roots = self.get_solutions(x=x, y=y)
			coeffs = self.get_coeff_list(x=x, y=y)
			for k in range(12):
				self.coeff_array[k].append(coeffs[k])
		self.get_coeff_strings()

	def get_tstat_array(self, sample_res):
		"""
		Fills an array for the t-stat values through the grid.

		Parameters:
			sample_res (float):
				The resolution of the sample area around each point
				used to calculate the t-stat.
		"""

		self.tstat_array = np.zeros(self.res**2)
		self.sample_magn = [[]*sample_res**2 for i in range(self.res**2)]
		self.get_magnification_array()
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			self.sample_magn[idx] = self.return_sample_magnification(
					x=x, y=y, sample_res = sample_res)
			self.tstat_array[idx] = self.get_tstat(
					point_magn = self.magn_array[idx],
					sample_magn_array = self.sample_magn[idx])

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
			magn_sorted = sorted(self.magn_array)
			cutoff = magn_sorted[(self.res**2) - 11]
			print('No cutoff value specified; selecting only upper 10 points')

		print('Finding the magnification outliers...')

		for (i, magn) in enumerate(self.magn_array):
			if magn > cutoff:
				self.x_outliers.append(self.x_array[i])
				self.y_outliers.append(self.y_array[i])
				self.magn_outliers.append(magn)

	def get_tstat_outliers(self):
		"""
		Fills an array with the values of t-stat for points determined
		to be outliers.
		"""

		self.tstat_outliers = []

		try:
			for i in range(len(self.magn_outliers)):
				self.tstat_outliers.append(self.tstat_array[i])
		except:
			raise ValueError('You must call get_magnification_outliers()',
									'before calling get_tstat_outliers()')

### The following functions return a single element to a grid variable.

	def return_sample_magnification(self, x, y, sample_res):
		"""
		Returns an array of magnification values for points in a small
		area around each sample point.

		Parameters:
			x (float):
				The x-position around which the sample area will be centered.

			y (float):
				The y-position around which the sample area will be centered.

			sample_res (float):
				The resolution (number of points on each side) of the
				sample area.

		Returns:
			sample_magn (array):
				An array containing local magnification values for each
				point centered around sample point.
		"""

		sample_xgrid = np.linspace(x - self.width_caustic/500.,
				x + self.width_caustic/500., sample_res)
		sample_ygrid = np.linspace(y - self.height_caustic/500.,
				y + self.height_caustic/500., sample_res)
		sample_magn = np.zeros(sample_res**2)

		for (i, sample_x) in enumerate(sample_xgrid):
			for (j, sample_y) in enumerate(sample_ygrid):
				sample_idx = sample_res*i + j
				sample_magn[sample_idx] = self.get_magnification(
						x=sample_x, y=sample_y)
		return sample_magn

	def get_tstat(self, point_magn, sample_magn_array):
		"""
		Returns the t-stat value for a given sample point.

		Parameters:
			point_magn (float):
				The magnification at the sample point.

			sample_magn_array (array):
				The array of magnification values centered around the
				sample point.

		Returns:
			tstat (float):
				The t-stat value determined for the sample point.
		"""

		mean_magn = sum(sample_magn) / len(sample_magn)
		stderr_magn = np.std(sample_magn) / np.sqrt(len(sample_magn))
		tstat = abs(magn - mean_magn) / stderr_magn
		return tstat

	def get_coeff_list(self, x, y):
		"""
		Returns a list of values corresponding to each coefficient in
		the binary lens polynomial.

		Parameters:
			x (float):
				The x-position of the source (i.e. the current value of x
				in the grid)

			y (float):
				The y-position of the source (i.e. the current value of y
				in the grid)

		Returns:
			c (list of floats):
				The value of all 6 real followed by all 6 imaginary
				coefficients for the given source position (x, y)
		"""

		coeffs = self.get_coefficients(x=x, y=y)
		c_real = []
		c_imag = []

		# This fills the c_real and c_imag lists in descending order. Specifically,
		# c_real[n] is the real component of the nth-degree term, etc.
		for i in range(len(coeffs)):
			c_r = float(np.real(coeffs[i]))
			c_i = float(np.imag(coeffs[i]))
			c_real.append(c_r)
			c_imag.append(c_i)

		# This reverses the order of the coefficients and places the real
		# coefficients in front of the imaginary. For example, c[4] is the real
		# component of the 4th degree term, and c[9] is the imaginary component
		# of the (9-6)=3rd degree term.
		c = c_real[::-1] + c_imag[::-1]

		return c

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

### The following functions make plots for grid data.

	def plot_n_solns(self, print_errors=True, region='caustic',
					 save=False, **kwargs):
		"""
		Creates a plot showing number of images on a grid (x,y) around
		the planetary caustic.

		Parameters:
			print_errors (bool):
				Prints the number of instances where the number of accepted
				solutions was equal to each integer value from 0 to 5
				inclusive.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		kwargs = self.check_kwargs(**kwargs)
		self.get_position_arrays(region = region)
		self.get_num_images_array()

		if print_errors:
			self.print_num_images_errors()

		plt.scatter(self.x_array, self.y_array, c=self.num_images, **kwargs)
		im_plot = plt.colorbar()
		im_plot.set_label('Num Images')
		plt.xlabel('X-position of Source', fontsize = 12)
		plt.ylabel('Y-position of Source', fontsize = 12)
		plt.gcf().set_size_inches(8, 6)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmax, (xmax - xmin) / 6))
		plt.suptitle('Number of Images', x=0.44)
		title = ('Frame: {}; Solver: {}; Region: {}\ns = {}, q = {}'.format(
				self.origin_title, self.solver_title, region, self.s, self.q))
		plt.title(title, fontsize = 11)

		if save:
			file_name = '../Tables/NumIm_{}_{}.png'.format(
					self.origin_str, self.solver_str)
			self.save_png(file_name = file_name)

	def plot_magnification(self, region = 'caustic', outliers = False, cutoff = None,
							log_colorbar = False, save = False, **kwargs):
		"""
		Make square grid of points that shows the magnification at each point.
		Parameters:
			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			save (bool):
				If True, file will be saved by convention determined below. If
				False, it will do nothing.
		"""

		kwargs = self.check_kwargs(log_colorbar, **kwargs)
		self.get_position_arrays(region = region)
		self.get_magnification_array()

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.
		if outliers:
			self.get_magnification_outliers(cutoff)
			(x, y, magn) = (self.x_outliers, self.y_outliers, self.magn_outliers)
			print('Plotting the magnification of outliers...')
		else:
			(x, y, magn) = (self.x_array, self.y_array, self.magn_array)
			print('Plotting the magnification...')

		plot = plt.scatter(x, y, c = magn, **kwargs)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmax, (xmax - xmin) / 6))
		mag_plot = plt.colorbar()
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of Source', fontsize = 12)
		plt.ylabel('Y-position of Source', fontsize = 12)
		plt.gcf().set_size_inches(8, 6)

		if outliers:
			plt.suptitle('High Magnification')
			title = ('Frame: {}; Solver: {}; Region: {}\n'.format(
					self.origin_title, self.solver_title, region) + 
					's = {}, q = {}, M > {:.0f}'.format(self.s, self.q, 
					min(magn)))
			plt.title(title, fontsize = 11)
			file_name = ('../Tables/HighMagn_{}_{}.png'.format(
					self.solver_file, self.origin_file))
		else:
			plt.suptitle('Magnification', x=0.44)
			title = ('Frame: {}; Solver: {}; Region: {}\n'.format(
					self.origin_title, self.solver_title, region) + 
					's = {}, q = {}'.format(self.s, self.q))
			plt.title(title, fontsize = 11)
			file_name = ('../Tables/Magn_{}_{}.png'.format(self.solver_file,
					self.origin_file))

		if save:
			self.save_png(file_name = file_name)

	def plot_rel_magnification(self, other_BL, region = 'caustic', outliers = False,
				ratio_cutoff = None, log_colorbar = False, save = False, **kwargs):
		"""
		Plots the fractional difference in magnification between two sets of data.

		Parameters:
			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			save (bool):
				If True, file will be saved by convention determined below. If
				False, it will do nothing.
		"""

		kwargs = self.check_kwargs(log_colorbar, **kwargs)
		self.get_position_arrays(region = region)
		self.get_magnification_array()
		other_BL.get_position_arrays(region = region)
		other_BL.get_magnification_array()

		print('Plotting the relative magnification...')

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.

		(x, y, magn1, magn2) = (self.x_array, self.y_array, self.magn_array,
													other_BL.magn_array)
		rel_magn = (magn1 / magn2)
		file_name = ('../Tables/RelMagn_{}_{}.png'.format(self.solver_file,
															self.origin_file))

		if outliers:
			rel_magn_outliers = []
			x_outliers = []
			y_outliers = []
			if ratio_cutoff == None:
				ratio_cutoff = 2.0
			for (i, magn) in enumerate(rel_magn):
				if (magn > ratio_cutoff) or (1./magn > ratio_cutoff):
					x_outliers.append(x[i])
					y_outliers.append(y[i])
					rel_magn_outliers.append(magn)
			x = x_outliers
			y = y_outliers
			rel_magn = rel_magn_outliers
			file_name = ('../Tables/RelMagnOut_{}_{}.png'.format(self.solver_file,
															self.origin_file))
			if len(x) == 0:
				print('No outliers found. Continuing with next plot...')
				return

		plot = plt.scatter(x, y, c = rel_magn, **kwargs)
		rel_plot = plt.colorbar()
		rel_plot.set_label('Magnification Ratio')
		plt.xlabel('X-position of Source', fontsize = 12)
		plt.ylabel('Y-position of Source', fontsize = 12)
		plt.gcf().set_size_inches(8, 6)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmax, (xmax - xmin) / 6))
		plt.title('Relative Magnification\n({}, {} Frame) / ({}, {} Frame)'.
							format(self.solver_title, self.origin_title,
								other_BL.solver_title, other_BL.origin_title))
		if save:
			self.save_png(file_name = file_name)

	def plot_outlier_coeff(self, other_BL, region = 'caustic', ratio_cutoff = None,
														save = False, **kwargs):

		if 's' not in kwargs:
			kwargs['s'] = 8
		kwargs = self.check_kwargs(**kwargs)
		self.get_position_arrays(region = region)
		self.get_magnification_array()
		self.get_coeff_array()
		other_BL.get_position_arrays(region = region)
		other_BL.get_magnification_array()
		other_BL.get_coeff_array()

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.

		(x, y, magn1, magn2, coeff) = (self.x_array, self.y_array, self.magn_array,
											other_BL.magn_array, self.coeff_array)
		rel_magn = (magn1 / magn2)

		rel_magn_outliers = []
		x_outliers = []
		y_outliers = []
		coeff_outliers = [[] for i in range(12)]
		if ratio_cutoff == None:
			ratio_cutoff = 2.0
		for (i, magn) in enumerate(rel_magn):
			if (magn > ratio_cutoff) or (1./magn > ratio_cutoff):
				x_outliers.append(x[i])
				y_outliers.append(y[i])
				rel_magn_outliers.append(magn)
				for k in range(12):
					coeff_outliers[k].append(coeff[k][i])

		if len(x_outliers) == 0:
			print('No outliers found. Continuing with next plot...')
			return

		print('Plotting...')
		for (i, coeff) in enumerate(self.coeff_array):
			file_name = ('../Tables/BL_{}_{}_{}.png'.format(self.coeff_file[i],
										 self.solver_file, self.origin_file))
			ax = plt.gca()
			kwargs['color'] = 'black'
			ax.scatter(coeff, rel_magn, **kwargs)
			kwargs['color'] = 'red'
			ax.scatter(coeff_outliers[i], rel_magn_outliers, **kwargs)
			ax.set_yscale('log')
			xmin = min(coeff)
			xmax = max(coeff)
			dx = xmax - xmin
			plt.xlabel(self.coeff_string[i])
			plt.ylabel('Relative Magnification')
			plt.xlim(xmin - 0.01*dx, xmax + 0.01*dx)
			plt.ylim(0, max(rel_magn_outliers))
			plt.xticks(np.arange(xmin, xmax, dx/6))
			plt.title('Outliers for {}'.format(self.coeff_string[i]))
			if save:
				self.save_png(file_name = file_name)
			plt.show()

	def plot_coeff_tstat(self, cutoff = None, region = 'caustic',
			plot_position_tstat = False, sample_res = 5, save = False,
			outliers = False, plot_coeffs = True, **kwargs):

		if 's' not in kwargs:
			kwargs['s'] = 5
		kwargs = self.check_kwargs(**kwargs)
		self.get_position_arrays(region = region)
		self.get_coeff_array()
		self.get_tstat_array(sample_res = sample_res)

		if outliers:
			self.get_magnification_outliers(cutoff, data = 'tstat')
			(x, y, magn, tstat) = (self.x_outliers, self.y_outliers,
								self.magn_outliers, self.tstat_outliers)
		else:
			(x, y, magn, tstat) = (self.x_array, self.y_array, self.magn_array,
																self.tstat_array)

		if plot_position_tstat:
			print('Plotting...')
			self.plot_position_tstat(x=x, y=y, tstat=tstat, outliers = outliers)
		if plot_coeffs:
			print('Plotting...')
			for (i, coeff) in enumerate(self.coeff_array):
				file_name = ('../Tables/BL_tstat_{}_{}_{}.png'.format(
						self.coeff_file[i], self.solver_file, self.origin_file))
				plt.scatter(coeff, self.tstat_array, color='black', **kwargs)
				xmin = min(coeff)
				xmax = max(coeff)
				dx = xmax - xmin
				ymax = max(self.tstat_array)
				plt.xlabel(self.coeff_string[i])
				plt.ylabel('t-Test Result')
				plt.xlim(xmin - 0.01*dx, xmax + 0.01*dx)
				plt.ylim(0, 1.01*ymax)
				plt.xticks(np.arange(xmin, xmax, dx/6))
				plt.title('t-Test Result vs. {}'.format(self.coeff_string[i]))
				if save:
					self.save_png(file_name = file_name)
				plt.show()

	def plot_position_tstat(self, x, y, tstat, outliers, cutoff = None, **kwargs):

		kwargs = self.check_kwargs(**kwargs)
		plt.scatter(x, y, c = tstat, **kwargs)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmax, (xmax - xmin) / 6))
		plt.xlabel('X-position')
		plt.ylabel('Y-position')
		region_plot = plt.colorbar()
		region_plot.set_label('t-value')
		if outliers:
			plt.title(('Magnification outliers t-test Score\n{} Frame; {} Solver; M > {:.0f}, q={}'.
				format(self.origin_title, self.solver_title, cutoff, self.q)))
		else:
			plt.title('t-test Score vs Position\nFrame: {}; Solver: {}'.format(
				self.origin_title, self.solver_title))
		caustic = mm.Caustics(s=self.s, q=self.q)
		caustic.plot(s=1)
		plt.show()

### The following functions save the relevant data to a .png or .fits file

	def write_to_fits(self, region = 'caustic'):
		"""
		Writes information about grid to a .fits table for comparison of magnification
		and number of images between different coordinate systems and solving methods.
		"""

		self.get_position_arrays(region = region)
		self.get_magnification_array()
		self.get_num_images_array()
		col = []
		col.append(fits.Column(name='x', array=self.x_array, format='D'))
		col.append(fits.Column(name='y', array=self.y_array, format='D'))
		col.append(fits.Column(name='Magnification', array=self.magn_array, format='D'))
		col.append(fits.Column(name='Number Images', array=self.num_images, format='I'))
		hdu1 = fits.BinTableHDU.from_columns(col)
		hdr = fits.Header()
		hdr['SEPARAT'] = '{:f}'.format(self.s)
		hdr['M_RATIO'] = '{:f}e-8'.format(1.e8*self.q)
		hdr['ORIGIN'] = self.origin
		hdr['SOLVER'] = self.solver
		hdr['RES'] = '{:d}'.format(self.res)
		hdu0 = fits.PrimaryHDU(header = hdr)
		hdus = fits.HDUList([hdu0, hdu1])

		for i in range(10):
			try:
				file_name = '../Tables/BL_{}_{}_{}.fits'.format(
									self.origin_file, self.solver_file, i)
				hdus.writeto(file_name)
				print(file_name, 'has been saved')
			except:
				continue
			break

	def save_png(self, file_name):
		for i in range(10):
			name = file_name[:-4] + '{}'.format(i) + file_name[-4:]
			if Path(name).is_file():
				continue
			plt.savefig(name)
			print(name, 'has been saved')
			return
		print('Error: too many files of same name already exist. File not saved')

### The following functions print optional information.

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

	def print_input(self):
		print('\nInput:\nx = {:}\ny = {:}\ns = {:}\nq = {:}\n'
			.format(self.x, self.y, self.s, self.q))
		print('Calculated in the {} using {}\n'.format(self.origin_phrase,
					self.solver_phrase))

	def print_image_position(self, print_input=True):
		"""
		Prints the image positions with the option to display the input parameters. 
		"""

		if print_input:
			self.print_input()
		self.get_solutions(x=self.x, y=self.y)
		image_pos = self.image_positions()
		print('Image locations:')
		for (i, pos) in enumerate(image_pos):
			print('{:.5f}'.format(pos))

	def print_magnification(self, print_input=True):
		"""
		Prints the total magnification with the option to display the input
		parameters.

		Parameters:
			print_input (bool):
				Prints the input parameteers specified by the user if True.
		"""

		if print_input:
			self.print_input()
		magn = self.get_magnification(x=self.x, y=self.y)
		print('Magnification: {:.5f}'.format(magn))

### The following functions gather string data that are used within the source code.

	def get_coeff_strings(self):
		self.coeff_string = [None]*12
		self.coeff_file = [None]*12
		for i in range(6):
			self.coeff_string[i] = ('Re(coeff{})'.format(i))
			self.coeff_string[i+6] = ('Im(coeff{})'.format(i))
			self.coeff_file[i] = ('Re_c{}'.format(i))
			self.coeff_file[i+6] = ('Im_c{}'.format(i))

	def strings(self):
		"""
		Assign appropriate strings for file-naming, printing, and plot titles.
		Return an error if an input string parameter is not recognized.
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
		elif self.solver == 'SG12':
			self.solver_file = 'SG'
			self.solver_title = 'SG 2012'
			self.solver_phrase = 'Skowron and Gould 2012 root finder'
		elif self.solver == 'zroots':
			self.solver_file = 'zr'
			self.solver_title = 'Zroots'
			self.solver_phrase = 'Numerical Recipes: zroots root finder'
		else:
			raise ValueError('Unknown solver: {:}'.format(self.solver))

