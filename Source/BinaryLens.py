# Zoey Samples
# Created: June 06, 2018
# BinaryLens.py
# Last Updated: June 13, 2018

import sys
import os
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

	def __init__(self, s, q, origin, solver, x=None, y=None, res=None, tolerance=0.0001):
		self.s = s
		self.q = q
		self.origin = origin
		self.solver = solver
		self.res = res
		self.x = x
		self.y = y
		self.tolerance = tolerance
		self.roots = None
		self.strings()

	def get_variables(self):
		"""
		Assign all variables to be used in polynomial. Assigns:
	
			dm	 - defined as half the difference of the masses (positive)
			m	 - defined as half the sum of the masses
			zeta - the source's position
			z1	 - the planet's position
			z2	 - the star's position
	
		"""
	
		self.mass_ratio()
		self.source_position()
		self.lensing_body_positions()

	def mass_ratio(self):
		"""Define m and dm for use in polynomial."""

		self.dm = (1. - self.q) / 2.
		self.m = (1. + self.q) / 2.

	def source_position(self):
		"""Assign the position of the source in complex space."""

		if self.origin == 'geo_cent':
			self.zeta = self.x + self.y*1.j
		elif self.origin == 'star':
			self.zeta = (self.x + self.s/2.) + self.y*1.j
		elif self.origin == 'plan':
			self.zeta = (self.x - self.s/2.) + self.y*1.j
		elif self.origin == 'com':
			self.zeta = (self.x + self.s*self.dm/(2.*self.m)) + self.y*1.j
		elif self.origin == 'caustic':
			self.zeta = (self.x + 1./(self.s) - self.s/2.) + self.y*1.j
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))
		self.zeta_conj = self.zeta.conjugate()

	def lensing_body_positions(self):
		"""
		Assign the positions of the lensing bodies (assumed to be on the 
		real axis).
		"""

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

	def get_coefficients(self):
		"""
		Returns the coefficients for the polynomial equation.
		"""

		self.get_variables()

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

	def solutions(self):
		"""Return solutions of polynomial -- written in general form."""

		# Assign the values of the coefficients of the polynomial.
		# Be aware that these expressions are very long and messy. They trail off
		# the screen without text wrapping and look non-indented with text wrapping.

		coeff_list = self.get_coefficients()

		# Return the roots of the polynomial via the given root finder
		if self.solver == 'SG12':
			rev_list = coeff_list[::-1]
			out = _vbbl_SG12_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			self.roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return self.roots
		elif self.solver == 'zroots':
			rev_list = coeff_list[::-1]
			#print(rev_list)
			#print(self.s, self.q, self.x, self.y)
			out = _zroots_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			self.roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return self.roots
		elif self.solver == 'numpy':
			self.roots = np.roots(coeff_list).tolist()
			return self.roots

	def check_solution(self, solution):
		"""
		Check if the determined solution is consistent with the binary 
		lens equation.
		"""

		# Make sure the roots have been determined before checking them
		if self.roots == None:
			raise ValueError('Solutions have not been found yet')

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

		self.solutions()
		image_pos = []
		for solution in self.roots:
			if self.check_solution(solution=solution):
				image_pos.append(solution)
		return image_pos

	def magnification(self):
		"""Returns the magnification for each configuration"""

		# Make sure the roots have been determined before calculating magnification
		if self.roots == None:
			raise ValueError('Solutions have not been found yet')

		magn = list(range(5))
		for (i, z) in enumerate(self.roots):
			detJ = (1. - ((self.m - self.dm) / ((z - self.z1)**2) + (self.m + 
				self.dm) / ((z - self.z2)**2)) * ((self.m - self.dm) / ((z.conjugate() - 
				self.z1)**2) + (self.m + self.dm) / ((z.conjugate() - self.z2)**2)))
			if self.check_solution(solution = z):
				magn[i] = np.abs(1./detJ)
			else:
				magn[i] = 0.

		# This is the sum of the calculated magnitude after removing 
		# non-physical results.
		self.tot_magn = sum(magn)

	def print_image_position(self, print_input=True):
		"""
		Prints the image positions with the option to display the input parameters. 
		"""

		if print_input:
			self.print_input()
		self.solutions()
		image_pos = self.image_positions()
		print('Image locations:')
		for (i, pos) in enumerate(image_pos):
			print('{:.5f}'.format(pos))

	def print_magnification(self, print_input=True):
		"""
		Prints the total magnification with the option to display the input 
		parameters.
		"""

		if print_input:
			self.print_input()
		self.solutions()
		self.magnification()
		print('Magnification: {:.5f}'.format(self.tot_magn))

	def size_caustic(self):
		"""
		Determines the width, height, and position of the center of the caustic.
		"""

		w = 4.*np.sqrt(self.q)*(1. + 1./(2.*(self.s**2))) / (self.s**2)
		h = 4.*np.sqrt(self.q)*(1. - 1./(2.*(self.s**2))) / (self.s**2)
		x = 0.5*self.s - 1.0/self.s
		self.delta = 2.*w / 200.
		self.epsilon = 2.*h / 200.
		return w, h, x

	def fill_grid_arrays(self, region):
		"""Fills arrays for the x- and y-position to prepare grid plots."""

		(w_caustic, h_caustic, x_center) = self.size_caustic()

		if region == 'caustic':
			region_xmin = x_center - w_caustic
			region_xmax = x_center + w_caustic
			region_ymin = -h_caustic
			region_ymax = h_caustic
		if region == 'onax_cusp':
			region_xmin = x_center + 0.3*w_caustic
			region_xmax = x_center + 0.9*w_caustic
			region_ymin = -0.15*h_caustic
			region_ymax = 0.15*h_caustic
		if region == 'offax_cusp':
			region_xmin = x_center - 0.15*w_caustic
			region_xmax = x_center + 0.15*w_caustic
			region_ymin = 0.3*h_caustic
			region_ymax = 0.9*h_caustic
		x_grid = np.linspace(region_xmin, region_xmax, self.res)
		y_grid = np.linspace(region_ymin, region_ymax, self.res)

		self.x_array = np.zeros(self.res**2)
		self.y_array = np.zeros(self.res**2)

		for (i, xx) in enumerate(x_grid):
			for (j, yy) in enumerate(y_grid):
				idx = self.res*i + j
				self.x_array[idx] = xx
				self.y_array[idx] = yy
				self.x = xx
				self.y = yy

	def grid_plots(self, data, region, sample_res = None):
		"""
		Fills arrays for the total magnification and number of images at each
		point in the grid arrays.

			data (string):
				'n_solns' - returns information for plot_n_solns
				'magn'	  - returns information for plot_magnification
				'tstat'	  - returns information for plot_coeff_tstat
		"""

		self.fill_grid_arrays(region = region)
		self.num_images = np.zeros(self.res**2, dtype=int)
		self.magn_array = np.zeros(self.res**2, dtype=float)
		self.tstat_array = np.zeros(self.res**2)
		self.coeff_array = [[]*self.res**2 for i in range(12)]
		if data == 'tstat':
			self.coeff_array = [[]*self.res**2 for i in range(12)]
			self.sample_magn = [[]*sample_res**2 for i in range(self.res**2)]
		self.get_coeff_strings()

		print('Getting data...')
		for idx in range(self.res**2):

			if idx == (self.res**2)/int(4):
				print('25%...')
			if idx == (self.res**2)/int(2):
				print('50%...')
			if idx == int(3)*(self.res**2)/int(4):
				print('75%...')

			self.x = self.x_array[idx]
			self.y = self.y_array[idx]
			roots = self.solutions()
			coeffs = self.get_coeff_list()
			for k in range(12):
				self.coeff_array[k].append(coeffs[k])
			if data == 'n_solns':
				for solution in self.roots:
					if self.check_solution(solution=solution):
						self.num_images[idx] += 1
			else:
				self.magnification()
				self.magn_array[idx] = self.tot_magn
				if data == 'tstat':
					coeffs = self.get_coeff_list()
					for k in range(12):
						self.coeff_array[k].append(coeffs[k])
					self.sample_magn[idx] = self.return_sample_magn(x = self.x,
											y = self.x, sample_res = sample_res)
					self.tstat_array[idx] = self.get_tstat(magn = self.magn_array[idx],
										sample_magn = self.sample_magn[idx])

	def return_sample_magn(self, x, y, sample_res):
		sample_xgrid = np.linspace(x - 0.2*self.delta, x + 0.2*self.delta, sample_res)
		sample_ygrid = np.linspace(y - 0.2*self.epsilon, y + 0.2*self.epsilon, sample_res)
		sample_magn = np.zeros(sample_res**2)
		for (i, xx) in enumerate(sample_xgrid):
			for (j, yy) in enumerate(sample_ygrid):
				sample_idx = sample_res*i + j
				self.x = xx
				self.y = yy
				roots = self.solutions()
				self.magnification()
				sample_magn[sample_idx] = self.tot_magn
		return sample_magn

	def get_tstat(self, magn, sample_magn):
		mean_magn = sum(sample_magn) / len(sample_magn)
		stderr_magn = np.std(sample_magn) / np.sqrt(len(sample_magn))
		tstat = abs(magn - mean_magn) / stderr_magn
		return tstat

	def get_coeff_list(self):
		coeffs = self.get_coefficients()
		c_real = []
		c_imag = []
		# This fills the c_real and c_imag lists in descending order. Specifically,
		# c_real[n] is the real component of the nth-degree term, etc.
		for i in range(len(coeffs)):
			c_r = float(np.real(coeffs[i]))
			c_i = float(np.imag(coeffs[i]))
			c_real.append(c_r)
			c_imag.append(c_i)

		"""
		This reverses the order of the coefficients and places the real
		coefficients in front of the imaginary. For example, c[4] is the real
		component of the 4th degree term, and c[9] is the imaginary component of
		the (9-6)=3rd degree term.
		"""
		c = c_real[::-1] + c_imag[::-1]
		return c

	def print_errors(self):
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

	def plot_n_solns(self, save=False, print_errors = True, region = 'caustic'):
		"""
		Plot showing number of images on a grid (x,y) around planetary caustic
		"""

		self.grid_plots(region = region, data = 'n_solns')
		if print_errors:
			self.print_errors()

		plt.scatter(self.x_array, self.y_array, c=self.num_images, s=((600./self.res)**2), 
						marker = 'o', cmap='plasma', lw=None)
		im_plot = plt.colorbar()
		im_plot.set_label('Num Images')
		plt.xlabel('X-position of Source', fontsize = 12)
		plt.ylabel('Y-position of Source', fontsize = 12)
		plt.gcf().set_size_inches(8, 6)
		plt.xlim(min(self.x_array), max(self.x_array))
		plt.ylim(min(self.y_array), max(self.y_array))
		plt.title('Number Images\nFrame: {}; Solver: {}'.format(
			self.origin_title, self.solver_title))

		if save:
			file_name = '../Tables/NumIm_{}_{}.png'.format(origin_str, solver_str)
			plt.savefig(file_name)
			print(file_name, 'has been saved')

	def plot_magnification(self, region = 'caustic', outliers = False, cutoff = None,
							log_colorbar = False, save=False):
		"""
		Make square grid of points that shows the magnification at each point.
		Attributes:
			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			save (bool):
				If True, file will be saved by convention determined below. If
				False, it will do nothing.
		"""

		self.grid_plots(region = region, data = 'magn')

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.
		if outliers:
			self.cutoff = cutoff
			self.get_magn_outliers()
			(x, y, magn) = (self.x_outliers, self.y_outliers, self.magn_outliers)
			size = (1000/self.res)**2
			print('Plotting the magnification of outliers...')
		else:
			(x, y, magn) = (self.x_array, self.y_array, self.magn_array)
			size = (500/self.res)**2
			print('Plotting the magnification...')

		if log_colorbar:
			norm = colors.LogNorm()
		else:
			norm = None
		plt.scatter(x, y, c = magn, norm=norm, s=size, marker = 'o',
														 cmap='plasma', lw=None)
		plt.xlim(min(self.x_array), max(self.x_array))
		plt.ylim(min(self.y_array), max(self.y_array))
		mag_plot = plt.colorbar()
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of Source', fontsize = 12)
		plt.ylabel('Y-position of Source', fontsize = 12)
		plt.gcf().set_size_inches(8, 6)
		if outliers:
			plt.title('High Magnification with Caustic:\n{} Frame; {} Solver; M > {:.0f}, q={}'.
						format(self.origin_title, self.solver_title, self.cutoff, self.q))
		else:
			plt.title('Magnification\nFrame: {}; Solver: {}'.format(
				self.origin_title, self.solver_title))

		if save:
			for i in range(10):
				try:
					if outliers:
						file_name = ('../Tables/HighMagn_{}_{}{}'.format(
									self.solver_file, self.origin_file, i))
						plt.savefig(file_name)
						print(file_name, 'has been saved')
					else:
						file_name = '../Tables/Magn_{}_{}{}.png'.format(origin_str,
																	solver_str, i)
						plt.savefig(file_name)
						print(file_name, 'has been saved')
				except:
					continue
				break

	def plot_rel_magnification(self, other_BL, region = 'caustic', outliers = False,
				cutoff = None, log_colorbar = False, save = False, hl_out = False):
		"""
		Plots the fractional difference in magnification between two sets of data.

		Attributes:
			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			save (bool):
				If True, file will be saved by convention determined below. If
				False, it will do nothing.
		"""

		self.grid_plots(region = region, data = 'magn')
		other_BL.grid_plots(region = region, data = 'magn')

		print('Plotting the relative magnification...')

		if log_colorbar:
			norm = colors.LogNorm()
		else:
			norm = None

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.

		(x, y, magn1, magn2) = (self.x_array, self.y_array, self.magn_array,
													other_BL.magn_array)
		rel_magn = (magn1 / magn2)

		size = (400/self.res)**2
		#(x, y, magn) = (self.x_array, self.y_array, self.magn_array)


		if outliers:
			rel_magn_outliers = []
			x_outliers = []
			y_outliers = []
			if cutoff == None:
				cutoff = 1.1
			for (i, magn) in enumerate(rel_magn):
				if (magn > cutoff) or (1./magn > cutoff):
					x_outliers.append(x[i])
					y_outliers.append(y[i])
					rel_magn_outliers.append(magn)
			x = x_outliers
			y = y_outliers
			rel_magn = rel_magn_outliers
			if len(x) == 0:
				print('No outliers found. Continuing with next plot...')
				return

		if (hl_out and outliers):
			plt.scatter(self.x_array, self.y_array, c='cyan', lw=None, s=size)
			plt.scatter(x, y, c='red', lw=None, s=size)
		else:
			plt.scatter(x, y, c = rel_magn, s=size, marker = 'o',
											cmap='plasma', lw=None, norm=norm)
			rel_plot = plt.colorbar()
			rel_plot.set_label('Fractional Difference')
		plt.xlabel('X-position of Source', fontsize = 12)
		plt.ylabel('Y-position of Source', fontsize = 12)
		plt.gcf().set_size_inches(8, 6)
		plt.xlim(min(self.x_array), max(self.x_array))
		plt.ylim(min(self.y_array), max(self.y_array))
		plt.title('Relative Magnification\n({}, {} Frame) / ({}, {} Frame)'.
							format(self.solver_title, self.origin_title,
								other_BL.solver_title, other_BL.origin_title))

	def plot_outlier_coeff(self, other_BL, region = 'caustic', cutoff = None):

		self.grid_plots(region = region, data = 'magn')
		other_BL.grid_plots(region = region, data = 'magn')

		print('Plotting the relative magnification...')


		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.

		(x, y, magn1, magn2, coeff) = (self.x_array, self.y_array, self.magn_array,
											other_BL.magn_array, self.coeff_array)
		rel_magn = (magn1 / magn2)
		size = (300/self.res)**2

		rel_magn_outliers = []
		x_outliers = []
		y_outliers = []
		coeff_outliers = [[] for i in range(12)]
		if cutoff == None:
			cutoff = 1.1
		for (i, magn) in enumerate(rel_magn):
			print(i)
			if (magn > cutoff) or (1./magn > cutoff):
				print('True')
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
			ax = plt.gca()
			ax.scatter(coeff, rel_magn, s=size, color='black', lw=None)
			ax.scatter(coeff_outliers[i], rel_magn_outliers, s=size, color='red', lw=None)
			ax.set_yscale('log')
			xmin = min(coeff)
			xmax = max(coeff)
			dx = xmax - xmin
			plt.xlabel(self.coeff_string[i])
			plt.ylabel('Relative Magnification')
			plt.xlim(xmin - 0.01*dx, xmax + 0.01*dx)
			plt.ylim(0, max(rel_magn_outliers))
			plt.title('Outliers for {}'.format(self.coeff_string[i]))
			plt.show()

	def get_magn_outliers(self, data = None):
		"""
		Creates new arrays of (x, y, magn) only for magnification values that
		are above the cutoff value.
		"""

		self.x_outliers = []
		self.y_outliers = []
		self.magn_outliers = []
		self.tstat_outliers = []

		# If the cutoff value is not specified, default to the 90th percentile
		# of magnification.
		if self.cutoff==None:
			magn_sorted = sorted(self.magn_array)
			self.cutoff = magn_sorted[(self.res**2) - 11]
			print('No cutoff value specified; selecting only upper 10 points')

		print('Finding the magnification outliers...')

		for (i, magn) in enumerate(self.magn_array):
			if magn > self.cutoff:
				self.x_outliers.append(self.x_array[i])
				self.y_outliers.append(self.y_array[i])
				self.magn_outliers.append(magn)
				if data == 'tstat':
					self.tstat_outliers.append(self.tstat_array[i])

	def print_input(self):
		print('\nInput:\nx = {:}\ny = {:}\ns = {:}\nq = {:}\n'
			.format(self.x, self.y, self.s, self.q))
		print('Calculated in the {} using {}\n'.format(self.origin_phrase,
					self.solver_phrase))

	#Step 1
	def plot_coeff_tstat(self, cutoff = None, region = 'caustic',
			plot_position_tstat = False, sample_res = 5, save = False,
			outliers = False, plot_coeffs = True):

#		self.tstat_grid_plot(region = region, region_res = region_res,
#								sample_res = sample_res)
		self.grid_plots(region = region, data = 'tstat', sample_res = sample_res)

		if outliers:
			self.cutoff = cutoff
			self.get_magn_outliers(data = 'tstat')
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
				plt.scatter(coeff, self.tstat_array, s=5, color='black', lw=None)
				xmin = min(coeff)
				xmax = max(coeff)
				dx = xmax - xmin
				ymax = max(self.tstat_array)
				plt.xlabel(self.coeff_string[i])
				plt.ylabel('t-Test Result')
				plt.xlim(xmin - 0.01*dx, xmax + 0.01*dx)
				plt.ylim(0, 1.01*ymax)
				plt.title('t-Test Result vs. {}'.format(self.coeff_string[i]))
				plt.show()
				if save:
					continue	#FIXME: Include save condition

	#Step 3 (opt)
	def plot_position_tstat(self, x, y, tstat, outliers):
		plt.scatter(x, y, c = tstat, marker = 'o', s=(300/self.res)**2)
		plt.xlim(min(self.x_array), max(self.x_array))
		plt.ylim(min(self.y_array), max(self.y_array))
		plt.xlabel('X-position')
		plt.ylabel('Y-position')
		region_plot = plt.colorbar()
		region_plot.set_label('t-value')
		if outliers:
			plt.title(('Magnification outliers t-test Score\n{} Frame; {} Solver; M > {:.0f}, q={}'.
				format(self.origin_title, self.solver_title, self.cutoff, self.q)))
		else:
			plt.title('t-test Score vs Position\nFrame: {}; Solver: {}'.format(
				self.origin_title, self.solver_title))
		caustic = mm.Caustics(s=self.s, q=self.q)
		caustic.plot(s=(40/self.res)**2)
		plt.show()

	"""
	#Step 2
	def tstat_grid_plot(self, region, region_res, sample_res):
		self.fill_region_grid_arrays(region = region, region_res = region_res)
		self.get_region_arrays(region_res = region_res, sample_res = sample_res)
		print('Getting t-test values...')
		for (idx, magn) in enumerate(self.region_magn):
			if idx == len(self.region_magn)/int(4):
				print('25%...')
			elif idx == len(self.region_magn)/int(2):
				print('50%...')
			elif idx == 3*len(self.region_magn)/int(4):
				print('75%...')
			self.sample_magn[idx] = self.return_sample_magn(x = self.region_xarray[idx],
							y = self.region_yarray[idx], sample_res = sample_res)
			self.region_tstat[idx] = self.get_tstat(magn = magn, sample_magn =
														self.sample_magn[idx])
			self.x = self.region_xarray[idx]
			self.y = self.region_yarray[idx]
		self.get_coeff_strings()

	#Step 2.5
	def get_coeff_strings(self):
		self.coeff_string = [None]*12
		for i in range(6):
			self.coeff_string[i] = ('Re(coeff{})'.format(i))
			self.coeff_string[i+6] = ('Im(coeff{})'.format(i))

	#Step 2.4
	def get_tstat(self, magn, sample_magn):
		mean_magn = sum(sample_magn) / len(sample_magn)
		stderr_magn = np.std(sample_magn) / np.sqrt(len(sample_magn))
		tstat = abs(magn - mean_magn) / stderr_magn
		return tstat

	#Step 2.3
	def return_sample_magn(self, x, y, sample_res):
		sample_xgrid = np.linspace(x - 0.2*self.delta, x + 0.2*self.delta, sample_res)
		sample_ygrid = np.linspace(y - 0.2*self.epsilon, y + 0.2*self.epsilon, sample_res)
		sample_magn = np.zeros(sample_res**2)
		for (i, xx) in enumerate(sample_xgrid):
			for (j, yy) in enumerate(sample_ygrid):
				sample_idx = sample_res*i + j
				self.x = xx
				self.y = yy
				roots = self.solutions()
				self.magnification()
				sample_magn[sample_idx] = self.tot_magn
		return sample_magn

	#Step 2.2
	def get_region_arrays(self, region_res, sample_res):		
		plt.xlim(min(self.region_xarray), max(self.region_xarray))
		plt.ylim(min(self.region_yarray), max(self.region_yarray))
		self.region_xarray = np.zeros(region_res**2)
		self.region_yarray = np.zeros(region_res**2)
		self.region_magn = np.zeros(region_res**2)
		self.region_tstat = np.zeros(region_res**2)
		self.coeff_array = [[]*region_res**2 for i in range(12)]
		self.sample_magn = [[]*sample_res**2 for i in range(region_res**2)]
		for (i, xx) in enumerate(self.region_xgrid):
			for (j, yy) in enumerate(self.region_ygrid):
				idx = region_res*i + j
				self.region_xarray[idx] = xx
				self.region_yarray[idx] = yy
				self.x = xx
				self.y = yy
				roots = self.solutions()
				self.magnification()
				self.region_magn[idx] = self.tot_magn
				coeffs = self.get_coeff_list()
				for k in range(12):
					self.coeff_array[k].append(coeffs[k])

	#Step 2.2.1
	print('Getting coefficient values...')
	def get_coeff_list(self):
		coeffs = self.get_coefficients()
		c_real = []
		c_imag = []
		# This fills the c_real and c_imag lists in descending order. Specifically,
		# c_real[n] is the real component of the nth-degree term, etc.
		for i in range(len(coeffs)):
			c_r = float(np.real(coeffs[i]))
			c_i = float(np.imag(coeffs[i]))
			c_real.append(c_r)
			c_imag.append(c_i)

		""
		This reverses the order of the coefficients and places the real
		coefficients in front of the imaginary. For example, c[4] is the real
		component of the 4th degree term, and c[9] is the imaginary component of
		the (9-6)=3rd degree term.
		""
		c = c_real[::-1] + c_imag[::-1]
		return c

	
	#Step 2.1
	def fill_region_grid_arrays(self, region, region_res):
		(w_caustic, h_caustic, x_center) = self.size_caustic()
		if region == 'cusp':
			region_xmin = x_center + 0.4*w_caustic
			region_xmax = x_center + 0.7*w_caustic
			region_ymin = -0.1*h_caustic
			region_ymax = 0.1*h_caustic
		if region == 'caustic':
			region_xmin = x_center - w_caustic
			region_xmax = x_center + w_caustic
			region_ymin = -h_caustic
			region_ymax = h_caustic
		self.region_xgrid = np.linspace(region_xmin, region_xmax, region_res)
		self.region_ygrid = np.linspace(region_ymin, region_ymax, region_res)
	"""

	def write_to_fits(self, region = 'caustic'):
		"""
		Writes information about grid to a .fits table for comparison of magnification
		and number of images between different coordinate systems and solving methods.
		"""

		self.grid_plots(region = region, data = 'magn')
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

	def get_coeff_strings(self):
		self.coeff_string = [None]*12
		for i in range(6):
			self.coeff_string[i] = ('Re(coeff{})'.format(i))
			self.coeff_string[i+6] = ('Im(coeff{})'.format(i))

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
