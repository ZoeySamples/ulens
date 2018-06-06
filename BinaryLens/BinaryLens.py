# Zoey Samples
# Created: June 6, 2018
# BinaryLens.py
# Last Updated: June 6, 2018

import sys
import os
import ctypes
import numpy as np
import cmath


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


Class BinaryLens(object):
	"""
	Information about this class goes here.
	"""

	def __init__(self, s, q, x, y, origin, solver):
		self.s = s
		self.q = q
		self.x = x
		self.y = y
		self.origin = origin
		self.solver = solver
		self.get_variables()
		self.roots = None

	def print_parameters(self, origin, solver)
		"""Print information about solver and origin."""

		# Assign the string specifying the coordinate frame in which calculations
		# will be done.
		if origin == 'geo_cent':
			origin_str = 'geometric center frame'
		elif origin == 'star':
			origin_str = 'star frame'
		elif origin == 'plan':
			origin_str = 'planet frame'
		elif origin == 'com':
			origin_str = 'center-of-mass frame'
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))

		# Assign the string specifying the root finder that will be used
		if solver == 'numpy':
			solver_str = 'numpy root finder'
		elif solver == 'Skowron_and_Gould_12':
			solver_str = 'Skowron_and_Gould_12 root finder'
		elif solver == 'zroots':
			solver_str = 'zroots root finder'
		else:
			raise ValueError('Unknown solver: {:}'.format(solver))

		print('Calculations will be carried out in the{}\nusing the {}'.format(
				origin_str, solver_str))


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

		if origin == 'geo_cent':
			self.zeta = self.x + self.y*1.j
		elif origin == 'star':
			self.zeta = (self.x + self.s/2.) + self.y*1.j
		elif origin == 'plan':
			self.zeta = (self.x - self.s/2.) + self.y*1.j
		elif origin == 'com':
			self.zeta = (self.x + self.s*self.dm/(2.*self.m)) + self.y*1.j
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))
		self.zeta_conj = zeta.conjugate()

	def lensing_body_positions(self):
		"""
		Assign the positions of the lensing bodies (assumed to be on the 
		real axis).
		"""

		if origin == 'geo_cent':
			self.z1 = 0.5*self.s
			self.z2 = -0.5*self.s
		elif origin == 'star':
			self.z1 = self.s
			self.z2 = 0.
		elif origin == 'plan':
			self.z1 = 0.
			self.z2 = -self.s
		elif origin == 'com':
			self.z1 = self.s*(self.m + self.dm) / (2.*self.m)
			self.z2 = -self.s*(self.m - self.dm) / (2.*self.m)
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))

		self.z1_conj = self.z1.conjugate()
		self.z2_conj = self.z2.conjugate()

	def solutions(self):
		"""Return solutions of polynomial -- written in general form."""

		# Assign the values of the coefficients of the polynomial.
		# Be aware that these expressions are very long and messy. They trail off
		# the screen without text wrapping and look non-indented with text wrapping.

		coeff5 = (-self.zeta_conj + self.z1)*(self.zeta_conj- self.z2)

		coeff4 = (self.m*self.z1 + self.m*self.z2 + 2.*(self.z1**2)*self.z2 + 2.*self.z1*(self.z2**2) + self.dm*(-self.z1 + self.z2) + self.z1*self.z2*zeta + (self.zeta_conj**2)*(2.*self.z1 + 2.*self.z2 + zeta) - self.zeta_conj*(2.*self.m + (self.z1 + self.z2)*(2.*self.z1 + 2.*self.z2 + zeta)))

		coeff3 = (self.dm*(self.z1**2) - self.m*(self.z1**2) - 2.*self.m*self.z1*self.z2 - (self.z1**3)*self.z2 - self.dm*(self.z2**2) - self.m*(self.z2**2) - 4.*(self.z1**2)*(self.z2**2) - self.z1*(self.z2**3) - 2.*self.m*self.z1*zeta - 2.*self.m*self.z2*zeta - 2.*(self.z1**2)*self.z2*zeta - 2.*self.z1*(self.z2**2)*zeta - (self.zeta_conj**2)*((self.z1**2) + 2.*self.z1*(2.*self.z2 + zeta) + self.z2*(self.z2 + 2.*zeta)) + self.zeta_conj*(2.*self.dm*(self.z1 - self.z2) + 2.*self.m*(self.z1 + self.z2 + 2.*zeta) + (self.z1 + self.z2)*((self.z1**2) + 4.*self.z1*self.z2 + (self.z2**2) + 2.*self.z1*zeta + 2.*self.z2*zeta)))

		coeff2 = (-2.*(self.m**2)*(self.z1 + self.z2 - 2.*zeta) - 3.*self.m*(2.*self.zeta_conj - self.z1 - self.z2)*(self.z1 + self.z2)*zeta + self.dm*(self.z1 - self.z2)*(2.*self.m - 2.*self.z1*self.z2 + self.z1*zeta + self.z2*zeta - 2.*self.zeta_conj*(self.z1 + self.z2 + zeta)) + (self.zeta_conj - self.z1)*(self.zeta_conj - self.z2)*((self.z2**2)*zeta + (self.z1**2)*(2.*self.z2 + zeta) + 2.*self.z1*self.z2*(self.z2 + 2.*zeta)))

		coeff1 = ((-self.dm**2)*((self.z1 - self.z2)**2) + (self.m**2)*((self.z1**2) + 6.*self.z1*self.z2 + (self.z2**2) - 4.*self.z1*zeta - 4.*self.z2*zeta) - self.m*(2.*self.zeta_conj - self.z1 - self.z2)*(self.z1*self.z2*(self.z2 - 4.*zeta) + (self.z1**2)*(self.z2 - zeta) - (self.z2**2)*zeta) - (self.zeta_conj - self.z1)*self.z1*(self.zeta_conj - self.z2)*self.z2*(2.*self.z2*zeta + self.z1*(self.z2 + 2.*zeta)) + self.dm*(self.z1 - self.z2)*(self.z1*self.z2*(self.z2 - 2.*zeta) + (self.z1**2)*(self.z2 - zeta) - (4.*self.m + (self.z2**2))*zeta + 2.*self.zeta_conj*(self.z2*zeta + self.z1*(self.z2 + zeta))))

		coeff0 = (-2.*(self.m**2)*(self.z1**2)*self.z2 - 2.*(self.m**2)*self.z1*(self.z2**2) - self.m*(self.z1**3)*(self.z2**2) - self.m*(self.z1**2)*(self.z2**3) + (self.m**2)*(self.z1**2)*zeta + (self.dm**2)*((self.z1 - self.z2)**2)*zeta + 2.*(self.m**2)*self.z1*self.z2*zeta + self.m*(self.z1**3)*self.z2*zeta + (self.m**2)*(self.z2**2)*zeta + (self.zeta_conj**2)*(self.z1**2)*(self.z2**2)*zeta + 2.*self.m*(self.z1**2)*(self.z2**2)*zeta + self.m*self.z1*(self.z2**3)*zeta + (self.z1**3)*(self.z2**3)*zeta - self.dm*(self.z1 - self.z2)*(2.*self.m + self.z1*self.z2)*(self.z1*(self.z2 - zeta) - self.z2*zeta) - self.zeta_conj*self.z1*self.z2*((2.*self.dm*(self.z1 - self.z2) + self.z1*self.z2*(self.z1 + self.z2))*zeta + self.m*(-2.*self.z1*self.z2 + 2.*self.z1*zeta + 2.*self.z2*zeta)))

		coeff_list = np.array([coeff5, coeff4, coeff3, coeff2, coeff1, coeff0])

		# Return the roots of the polynomial via the given root finder
		if solver == 'Skowron_and_Gould_12':
			rev_list = coeff_list[::-1]
			out = _vbbl_SG12_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			self.roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return self.roots
		elif solver == 'zroots':
			rev_list = coeff_list[::-1]
			out = _zroots_5(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			self.roots = [out[i] + out[i+5] * 1.j for i in range(5)]
			return self.roots
		elif solver == 'numpy':
			self.roots = np.roots(coeff_list).tolist()
			return self.roots

	def check_solution(self, tolerance = 0.0001):
		"""
		Check if the determined solution is consistent with the binary 
		lens equation.
		"""

		# Make sure the roots have been determined before checking them
		if self.roots == None:
			self.solutions()

		for z in self.roots:
			# This is the binary lens equation.
			zeta_actual = (z + (self.m - self.dm) / (self.z1_conj - z.conjugate()) +
						(self.m + self.dm) / (self.z2_conj - self.z.conjugate()))
			if np.abs(self.zeta - zeta_actual) > tolerance:
				return False
			else:
				return True

	def magnification(self, tolerance = 0.0001):
		"""Returns the magnification for each configuration"""

		# Make sure the roots have been determined before calculating magnification
		if self.roots == None:
			self.solutions()

		magn = list(range(5))
		for (i, z) in enumerate(self.roots):
			detJ = (1. - ((self.m - self.dm) / ((z - self.z1)**2) + (self.m + 
				self.dm) / ((z - self.z2)**2)) * ((self.m - self.dm) / ((z.conjugate() - 
				self.z1)**2) + (self.m + self.dm) / ((z.conjugate() - self.z2)**2)))
			if self.check_solution(tolerance = tolerance):
				magn[i] = np.abs(1./detJ)
			else:
				magn[i] = 0.
		return sum(magn)
	
############



