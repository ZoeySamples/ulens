# Zoey Samples
# Created: June 06, 2018
# BinaryLens.py
# Last Updated: June 21, 2018

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

# Here we attempt to access the Numerical Recipes zroots solver
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
				The maximum distance away a calculated image position must be 
				from an actual solution when substitued back into the binary
				lens equation, in order not to be rejected. If not provided,
				tolerance is given by default value.

			coeff_multiplier (float):
				If provided, all polynomial coefficients will be multiplied
				by this number before the root solver solves the polynomial.

		Global Variables:
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
					plot_num_images(save, print_errors)
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
		self.get_lensing_body_positions()
		self.get_mass()
		self.tolerance = tolerance

### The following functions assign values to variables pertaining to the class.

	def get_mass(self):
		"""Define m and dm."""

		self.dm = (1. - self.q) / 2.
		self.m = (1. + self.q) / 2.

	def get_lensing_body_positions(self):
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

	def get_source_position(self, x, y):
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

### The following functions calculate physical values for further analysis.

	def get_coefficients(self, x, y):
		"""Returns the coefficients for the polynomial equation."""

		self.get_source_position(x=x, y=y)

		# Assign the values of the coefficients of the polynomial.
		# Be aware that these expressions are very long and messy. They trail off
		# the screen without text wrapping and look non-indented with text wrapping.

		coeff5 = (-self.zeta_conj + self.z1)*(self.zeta_conj- self.z2)

		coeff4 = (self.m*self.z1 + self.m*self.z2 + 2.*(self.z1**2)*self.z2 + 2.*self.z1*(self.z2**2) + self.dm*(-self.z1 + self.z2) + self.z1*self.z2*self.zeta + (self.zeta_conj**2)*(2.*self.z1 + 2.*self.z2 + self.zeta) - self.zeta_conj*(2.*self.m + (self.z1 + self.z2)*(2.*self.z1 + 2.*self.z2 + self.zeta)))

		coeff3 = (self.dm*(self.z1**2) - self.m*(self.z1**2) - 2.*self.m*self.z1*self.z2 - (self.z1**3)*self.z2 - self.dm*(self.z2**2) - self.m*(self.z2**2) - 4.*(self.z1**2)*(self.z2**2) - self.z1*(self.z2**3) - 2.*self.m*self.z1*self.zeta - 2.*self.m*self.z2*self.zeta - 2.*(self.z1**2)*self.z2*self.zeta - 2.*self.z1*(self.z2**2)*self.zeta - (self.zeta_conj**2)*((self.z1**2) + 2.*self.z1*(2.*self.z2 + self.zeta) + self.z2*(self.z2 + 2.*self.zeta)) + self.zeta_conj*(2.*self.dm*(self.z1 - self.z2) + 2.*self.m*(self.z1 + self.z2 + 2.*self.zeta) + (self.z1 + self.z2)*((self.z1**2) + 4.*self.z1*self.z2 + (self.z2**2) + 2.*self.z1*self.zeta + 2.*self.z2*self.zeta)))

		coeff2 = (-2.*(self.m**2)*(self.z1 + self.z2 - 2.*self.zeta) - 3.*self.m*(2.*self.zeta_conj - self.z1 - self.z2)*(self.z1 + self.z2)*self.zeta + self.dm*(self.z1 - self.z2)*(2.*self.m - 2.*self.z1*self.z2 + self.z1*self.zeta + self.z2*self.zeta - 2.*self.zeta_conj*(self.z1 + self.z2 + self.zeta)) + (self.zeta_conj - self.z1)*(self.zeta_conj - self.z2)*((self.z2**2)*self.zeta + (self.z1**2)*(2.*self.z2 + self.zeta) + 2.*self.z1*self.z2*(self.z2 + 2.*self.zeta)))

		coeff1 = ((-self.dm**2)*((self.z1 - self.z2)**2) + (self.m**2)*((self.z1**2) + 6.*self.z1*self.z2 + (self.z2**2) - 4.*self.z1*self.zeta - 4.*self.z2*self.zeta) - self.m*(2.*self.zeta_conj - self.z1 - self.z2)*(self.z1*self.z2*(self.z2 - 4.*self.zeta) + (self.z1**2)*(self.z2 - self.zeta) - (self.z2**2)*self.zeta) - (self.zeta_conj - self.z1)*self.z1*(self.zeta_conj - self.z2)*self.z2*(2.*self.z2*self.zeta + self.z1*(self.z2 + 2.*self.zeta)) + self.dm*(self.z1 - self.z2)*(self.z1*self.z2*(self.z2 - 2.*self.zeta) + (self.z1**2)*(self.z2 - self.zeta) - (4.*self.m + (self.z2**2))*self.zeta + 2.*self.zeta_conj*(self.z2*self.zeta + self.z1*(self.z2 + self.zeta))))

		coeff0 = (-2.*(self.m**2)*(self.z1**2)*self.z2 - 2.*(self.m**2)*self.z1*(self.z2**2) - self.m*(self.z1**3)*(self.z2**2) - self.m*(self.z1**2)*(self.z2**3) + (self.m**2)*(self.z1**2)*self.zeta + (self.dm**2)*((self.z1 - self.z2)**2)*self.zeta + 2.*(self.m**2)*self.z1*self.z2*self.zeta + self.m*(self.z1**3)*self.z2*self.zeta + (self.m**2)*(self.z2**2)*self.zeta + (self.zeta_conj**2)*(self.z1**2)*(self.z2**2)*self.zeta + 2.*self.m*(self.z1**2)*(self.z2**2)*self.zeta + self.m*self.z1*(self.z2**3)*self.zeta + (self.z1**3)*(self.z2**3)*self.zeta - self.dm*(self.z1 - self.z2)*(2.*self.m + self.z1*self.z2)*(self.z1*(self.z2 - self.zeta) - self.z2*self.zeta) - self.zeta_conj*self.z1*self.z2*((2.*self.dm*(self.z1 - self.z2) + self.z1*self.z2*(self.z1 + self.z2))*self.zeta + self.m*(-2.*self.z1*self.z2 + 2.*self.z1*self.zeta + 2.*self.z2*self.zeta)))


		coefficients = np.array([coeff5, coeff4, coeff3, coeff2, coeff1, coeff0])
		return coefficients

	def get_roots(self, x, y):
		"""Return solutions of polynomial."""

		coefficients = self.get_coefficients(x=x, y=y)


		# Return the roots of the polynomial via the given root finder
		if self.solver == 'SG12':
			rev_list = coefficients[::-1]
			out = _vbbl_SG12_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return np.array(roots)
		elif self.solver == 'zroots':
			rev_list = coefficients[::-1]
			out = _zroots_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return np.array(roots)
		elif self.solver == 'numpy':
			roots = np.roots(coefficients).tolist()
			return np.array(roots)

	# Old method for calculating image positions; somewhat reliable
	# When adapting to new method, re-work every spot that calls this function
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

	# New method for calculating image positions; not working
	def get_accepted_solutions(self, x, y):
		
		roots = self.get_roots(x=x, y=y)
		lensing_body1 = (self.m + self.dm) / np.conjugate(roots - self.z1)
		lensing_body2 = (self.m - self.dm) / np.conjugate(roots - self.z2)
		solutions = self.zeta + lensing_body1 + lensing_body2

	
		len1 = [None]*5
		len2 = [None]*5
		r = [None]*5
		soln = [None]*5
		for k in range(5):
			len1[k] = ('{:.5f}'.format(abs(lensing_body1[k])))
			len2[k] = ('{:.5f}'.format(abs(lensing_body2[k])))
			r[k] = ('{:.4f}'.format(roots[k].real))
			soln[k] = ('{:.4f}'.format(solutions[k]))
		"""
		print('The real parts of the roots are:\n', r)
		print('The position of body 1 is:', self.z1, '\n')
		print('Here are the values for terms 2 and 3 in binary lens eqn:')
		print(len1, '\n', len2, '\n')
		"""

		accepted_solutions = []
		distances = []
		for (i, root) in enumerate(roots):
			distances_from_root = abs((solutions - root)**2)
			min_distance_arg = np.argmin(distances_from_root)


			dis = []
			for k in distances_from_root:
				dis.append('{:.2f}'.format(k))

			if i == min_distance_arg:
				accepted_solutions.append(root)
				distances.append(distances_from_root[min_distance_arg])
				print('Root {} Accepted'.format(i+1))
			else:
				print('Root {} Rejected'.format(i+1))

			print('The distances away are:', dis)
			print('The solved-for root is		', '{:.5f}'.format(root))
			print('The binary lens solution is	', soln[i], '\n')

		return accepted_solutions


	def image_positions(self):
		"""
		Calculates the image positions (i.e. checks which solutions pass the
		check). Returns a list of the positions.
		"""

		image_positions = []

		# Old method
		roots = self.get_roots(x=self.x, y=self.y)
		for solution in roots:
			if self.check_solution(solution=solution):
				image_positions.append(solution)
#		print('Old method:', image_positions)

		# New method; Comment next line to use old method.
		image_positions = self.get_accepted_solutions(x=self.x, y=self.y)
		return image_positions

	def get_magnification(self, x, y):
		"""Returns the magnification for each configuration."""

		roots = self.get_roots(x=x, y=y)
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

	def print_input(self):
		"""Prints the input parameters for a single test point."""

		print('\nInput:\nx = {:}\ny = {:}\ns = {:}\nq = {:}\n'
			.format(self.x, self.y, self.s, self.q))
		print('Calculated in the {} using {}\n'.format(self.origin_phrase,
					self.solver_phrase))

	def print_image_position(self, print_input=True):
		"""
		Prints the positions of those images which pass the check_solutions
		test, for a single test point.

		Optional parameters:

			print_input (bool):
				If True, prints the user-specified input parameters before
				displaying the image positions.
		"""

		if print_input:
			self.print_input()
		self.get_roots(x=self.x, y=self.y)
		image_pos = self.image_positions()
		print('Image locations:')
		for (i, pos) in enumerate(image_pos):
			print('{:.5f}'.format(pos))

	def print_magnification(self, print_input=True):
		"""
		Prints the total magnification for single test point.

		Parameters:

			print_input (bool):
				If True, prints the user-specified input parameters before
				displaying the image positions.
		"""

		if print_input:
			self.print_input()
		magn = self.get_magnification(x=self.x, y=self.y)
		print('Magnification: {:.5f}'.format(magn))
