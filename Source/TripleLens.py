# Zoey Samples
# Created: June 08, 2018
# TripleLens.py
# Last Updated: June 08, 2018

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


class TripleLens(object):

	"""
	Using the triple lens equation: Calculates solutions, magnification, makes
	plots, and writes to files. Works for point test runs or filling a grid of 
	points centered on the planetary caustic.

	Attributes:

		Required:

			eps1 (float):
				The mass ratio of planet 1 to the star.

			eps2 (float):
				The mass ratio of planet 2 to the star.

			m (float):
				The mass of the star

			s1 (float):
				The separation between ...

			s2 (float):
				The separation between ...

			alpha (float):
				The angle between ...

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
	"""


	def __init__(self, eps1, eps2, m, s1, s2, alpha, origin, solver, x=None, y=None, res=None, tolerance=0.0001):
		self.eps1 = eps1
		self.eps2 = eps2
		self.m = m
		self.s1 = s1
		self.s2 = s2
		self.alpha = alpha
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

			zeta - the source's position
			z0	 - the star's position
			z1	 - one planet's position
			z2	 - the other planet's position
		"""

		self.source_position()
		self.lensing_body_positions()

	def source_position(self):
		#

	def lensing_body_positions():
		#











