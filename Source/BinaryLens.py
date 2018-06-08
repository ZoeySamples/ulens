# Zoey Samples
# Created: June 6, 2018
# BinaryLens.py
# Last Updated: June 6, 2018

import sys
import os
import ctypes
import numpy as np
import cmath
import matplotlib.pyplot as plt
import MulensModel as mm

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

PATH = os.path.join(MODULE_PATH, 'zroots_codes', "zrootsBinaryLens_wrapper.so")

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
	Using the binary lens equation:	Calculates solutions, magnification, makes
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
					'caustic'	- the planetary caustic frame *SOON*

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
		Assign all variables to be used in polynomial. Returns:
	
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
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))

		self.z1_conj = self.z1.conjugate()
		self.z2_conj = self.z2.conjugate()

	def get_coefficients(self):
		"""
		Returns the coefficients for the polynomial equation.
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

		self.get_variables()
		coeff_list = self.get_coefficients()

		# Return the roots of the polynomial via the given root finder
		if self.solver == 'SG12':
			rev_list = coeff_list[::-1]
			out = _vbbl_SG12_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			self.roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return self.roots
		elif self.solver == 'zroots':
			rev_list = coeff_list[::-1]
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

	def print_input(self):
		print('Input:\nx = {:}\ny = {:}\ns = {:}\nq = {:}\n'
			.format(self.x, self.y, self.s, self.q))
		print('Calculated in the {} using {}\n'.format(self.origin_phrase,
					self.solver_phrase))

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
		return w, h, x

	def get_grid_arrays(self):
		"""Fills arrays for the x- and y-position to prepare grid plots."""

		(w_caustic, h_caustic, x_cent) = self.size_caustic()
		x_grid = np.linspace(x_cent - w_caustic, 
								x_cent + w_caustic, self.res)
		y_grid = np.linspace(-h_caustic, h_caustic, self.res)
		self.x_array = np.zeros(self.res**2)
		self.y_array = np.zeros(self.res**2)

		for (i, xx) in enumerate(x_grid):
			for (j, yy) in enumerate(y_grid):
				idx = self.res*i + j
				self.x_array[idx] = xx
				self.y_array[idx] = yy
				self.x = xx
				self.y = yy

	def grid_plots(self):
		"""
		Fills arrays for the total magnification and number of images at each
		point in the grid arrays.
		"""

		self.get_grid_arrays()
		self.num_images = np.zeros(self.res**2, dtype=int)
		self.mag_1d = np.zeros(self.res**2, dtype=float)
		for idx in range(self.res**2):
			self.x = self.x_array[idx]
			self.y = self.y_array[idx]
			roots = self.solutions()
			self.magnification()
			self.mag_1d[idx] = self.tot_magn
			for solution in self.roots:
				if self.check_solution(solution=solution):
					self.num_images[idx] += 1

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

	def plot_n_solns(self, save=False, print_errors = True):
		"""
		Plot showing number of images on a grid (x,y) around planetary caustic
		"""

		self.grid_plots()
		if print_errors:
			self.print_errors()

		plt.scatter(self.x_array, self.y_array, c=self.num_images, s=((500./self.res)**2), 
						marker = 'o', cmap='jet', lw=None)
		im_plot = plt.colorbar()
		im_plot.set_label('Num Images')
		plt.xlabel('X-position of source')
		plt.ylabel('Y-position of source')
		(w_caustic, h_caustic, x_cent) = self.size_caustic()
		plt.xlim(x_cent - w_caustic, x_cent + w_caustic)
		plt.ylim(-h_caustic, h_caustic)
		plt.title('Number Images\nFrame: {}; Solver: {}'.format(
			self.origin_title, self.solver_title))

		if save:
			file_name = '../Tables/NumIm_{}_{}.png'.format(origin_str, solver_str)
			plt.savefig(file_name)
			print(file_name, 'has been saved')

	def plot_magnification(self, save=False):
		"""
		Make square grid of points that shows the magnification at each point
		"""

		self.grid_plots()
		plt.scatter(self.x_array, self.y_array, c=self.mag_1d, s=((500./self.res)**2),
								marker = 'o', cmap='jet', lw=None)
		mag_plot = plt.colorbar()
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of source')
		plt.ylabel('Y-position of source')
		(w_caustic, h_caustic, x_cent) = self.size_caustic()
		plt.xlim(x_cent - w_caustic, x_cent + w_caustic)
		plt.ylim(-h_caustic, h_caustic)
		plt.title('Magnification\nFrame: {}; Solver: {}'.format(
			self.origin_title, self.solver_title))

		if save:
			file_name = '../Tables/Magn_{}_{}.png'.format(origin_str, solver_str)
			plt.savefig(file_name)
			print(file_name, 'has been saved')

	def write_to_fits(self):
		"""
		Writes information about grid to a .fits table for comparison of magnification
		and number of images between different coordinate systems and solving methods.
		"""

		self.grid_plots()
		col = []
		col.append(fits.Column(name='x', array=self.x_array, format='D'))
		col.append(fits.Column(name='y', array=self.y_array, format='D'))
		col.append(fits.Column(name='Magnification', array=self.mag_1d, format='D'))
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

		file_name = '../Tables/BL_{}_{:2}.fits'.format(self.origin_file, 
													self.solver_file)
		hdus.writeto(file_name)
		print(file_name, 'has been saved')

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
