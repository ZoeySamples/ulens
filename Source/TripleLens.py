# Zoey Samples
# Created: June 08, 2018
# TripleLens.py
# Last Updated: June 28, 2018

import sys
import os
import ctypes
import numpy as np
import math
import matplotlib.pyplot as plt
import MulensModel as mm
import TLGetCoeff as getc


MODULE_PATH = os.path.abspath(__file__)
for i in range(2):
	MODULE_PATH = os.path.dirname(MODULE_PATH)
MODULE_PATH = os.path.join(MODULE_PATH, 'Solvers')
PATH = os.path.join(MODULE_PATH, 'VBBL',
		"VBBinaryLensingLibrary_wrapper.so")

# Here we attempt to access the Skowron & Gould 2012 root finder
try:
	vbbl = ctypes.cdll.LoadLibrary(PATH)
except OSError as error:
	msg = "Something went wrong with VBBL wrapping ({:})\n\n" + repr(error)
	print(msg.format(PATH))
else:
	vbbl.VBBL_SG12_10.argtypes = 22 * [ctypes.c_double]
	vbbl.VBBL_SG12_10.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, 
		shape=(20,))
	_vbbl_SG12_10 = vbbl.VBBL_SG12_10

PATH = os.path.join(MODULE_PATH, 'NumericalRecipes',
					'zrootsBinaryLens_wrapper.so')

# Here we attempt to access the Numerical Recipes zroots solver
try:

	zroots = ctypes.cdll.LoadLibrary(PATH)
except OSError as error:
	msg = "Something went wrong with zroots wrapping ({:})\n\n" + repr(error)
	print(msg.format(PATH))
else:
	zroots.zroots_10.argtypes = 22 * [ctypes.c_double]
	zroots.zroots_10.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, 
		shape=(20,))
	_zroots_10 = zroots.zroots_10

class TripleLens(object):

	"""
	Using the triple lens equation: Calculates solutions, magnification, makes
	plots, and writes to files. Works for point test runs or filling a grid of 
	points centered on the planetary caustic.

	Attributes:

		Required:

			system (string):
				The types of bodies in the triple system. Options are:

					'SPM' - A system with a star, a planet, and a moon.
							{z1,m1}: star; {z2,m2}: plan; {z3,m3}: moon
							phi is the angle formed by the points z1,z2,z3
					'SPP' - A system with a star and 2 planets
							{z1,m1}: star, {z2,m2}: plan1, {z3,m3}: plan2
					'SSP' - A system with 2 stars and a planet
							{z1,m1}: star1, {z2,m2}: star2, {z3,m3}: plan

				*Definitions of each m and z variable are given below.
				*The letters of the string should be ordered from highest
				 mass to lowest mass.

			q1 (float):
				If system is 'SPM': The mass ratio of the planet to the star.
				If system is 'SPP': The mass ratio of planet1 to the star.
				If system is 'SSP': The mass ratio of star2 to star1.

			q2 (float):
				If system is 'SPM': The mass ratio of the moon to the planet.
				If system is 'SPP': The mass ratio of planet2 to the star.
				If system is 'SSP': The mass ratio of the planet to star1.

			s1 (float):
				If system is 'SPM': The separation of the planet to the star
						in units of the whole system's Einstein radius.
				If system is 'SPP': The separation of planet1 to the star
						in units of the whole system's Einstein radius.
				If system is 'SSP': The separation of star1 to star2
						in units of the whole system's Einstein radius.

				*A positive separation means body2 is in the positive
				 x-direction from body1. A negative separation means body2
				 is in the negative x-direction from body1. Both bodies 
				 lie on the same horizontal axis.

			s2 (float):
				If system is 'SPM': The separation of the moon to the planet
						in units of the planet-moon Einstein radius.
				If system is 'SPP': The separation of planet2 to the star
						in units of the whole system's Einstein radius.
				If system is 'SSP': The separation of the planet to star1
						in units of the whole system's Einstein radius.

				*This separation should always be positive, as the angle,
				 phi, determines the direction of separation.

			phi (float):
				If system is 'SPM': The angle formed by star-planet-moon.
				If system is 'SPP': The angle planet2 makes with the
						horizontal axis passing through the star and planet1.
				If system is 'SSP': The angle the planet makes with the
						horizontal axis passing through star1 and star2.

				*phi is given in degrees in all cases.

			origin (string):
				The coordinate frame in which calculations are carried out.
				Options are:
					'geo_cent' - the central point between body1 and body2
					'body2'    - 2nd most massive body's frame
					'body3'    - the least massive body's frame

			solver (string):
				Determines which root finder to use to solve the polynomial.
				Options are:
					'numpy'  - Uses np.roots method
					'SG12'   - Uses Skowron & Gould 2012 method
					'zroots' - Uses zroots laguerre method

			*q1 and q2 should be less than 1, so that star2 is less massive
			 than star1, or planet2 is less massive than planet1.

			*For s1 and s2, the separation can be positive or negative,
			 depending on the scenario.


		Optional:
			x (float):
				The x-position of the source for a single test point in units
				of the total system's Einstein radius.

			y (float):
				The y-position of the source for a single test point in units
				of the total system's Einstein radius.

			res (float):
				The resolution of the grid (i.e. number of points on each side). 
				Required if user is plotting a grid with TripleLens class.

			specific_frame_derivation (bool):
				If True, the form of the polynomial given by the origin
				will be used to calculate solutions.

		Global Variables:
			m1 (float):
				The mass of the largest body (the star in SPM or SPP system,
				or star1 in SSP system)	as a fraction of the whole system's
				mass.

			m2 (float):
				The mass of 2nd largest body (the planet in SPM system,
				planet1 in SPP system, or star2 in SSP system),	as a fraction
				of the whole system's mass.

			m3 (float):
				The mass of the least massive body (the moon in SPM system,
				planet2 in SPP system, or the planet in SSP system) as a
				fraction of the whole system's mass.

			z1 (complex):
				The position of the most massive body (the star in SPM or
				SPP system, or star1 in SSP system) in units of the total
				system's Einstein radius.

			z2 (complex):
				The position of the 2nd most massive body (the planet in SPM
				system, planet1 in SPP system, or star2 in SSP system) in
				units of the total system's Einstein radius.

			z3 (complex):
				The position of the least massive body (the moon in SPM
				system, planet2 in SPP system, or the planet in SSP system)
				in units of the total system's Einstein radius.
	"""


	def __init__(self, q2, q1, s2, s1, phi, origin, solver, system=None, x=None,
				 y=None, res=None, specific_frame_derivation=True):

		self.q2 = q2
		self.q1 = q1
		self.s2 = s2
		self.s1 = s1
		self.phi = phi
		self.origin = origin
		self.solver = solver
		self.system = system
		self.res = res
		self.x = x
		self.y = y
		self.roots = None
		self.specific_frame_derivation = specific_frame_derivation
		self.get_mass()
		self.get_lensing_body_positions()
		self.strings()

		if (self.q1 > 1.) or (self.q2 > 1.):
			raise ValueError('q1 and q2 are not allowed to be greater than 1.')

		if self.s2 < 0.:
			raise ValueError('s2 is not allowed to be negative.')

#		if (self.system != 'SPM') and (self.system != 'SPP') and (self.system != 'SSP'):
#			raise ValueError('Unknown value for string variable, system.')

### The following functions assign values to variables pertaining to the class.

	def get_mass(self):

		denominator = 1. + self.q1 + self.q2*self.q1
		self.m1 = 1. / denominator
		self.m2 = self.q1 / denominator
		self.m3 = self.q2*self.q1 / denominator

		# Convert the separation between the moon and planet into units
		# of the total system's Einstein radius
		self.s2 *= np.sqrt((self.m3 + self.m2))

	def get_lensing_body_positions(self):

		self.phi *= math.pi / 180
		self.displacementMP = self.s2*(math.cos(math.pi - self.phi) +
									1.j*math.sin(math.pi - self.phi))

		if self.origin == 'geo_cent':
			self.z1 = -0.5*self.s1 + 0j
			self.z2 = 0.5*self.s1 + 0j
			self.z3 = self.z2 + self.displacementMP

		elif self.origin == 'body2':
			self.z1 = -self.s1 + 0j
			self.z2 = 0j
			self.z3 = self.displacementMP

		elif self.origin == 'body3':
			self.z1 = -self.s1 - self.displacementMP
			self.z2 = -self.displacementMP
			self.z3 = 0j

		else:
			raise ValueError('Unknown coordinate system: {:}'.format(self.origin))


	def get_source_position(self, x, y):

		if self.origin == 'geo_cent':
			zeta = x + y*1.j

		elif self.origin == 'body2':
			zeta = (x - self.s1/2.) + y*1.j

		elif self.origin == 'body3':
			zeta = (x - self.s1/2.) + y*1.j - self.displacementMP

		else:
			raise ValueError('Unknown coordinate system: {:}'.format(origin))

		return zeta


	def get_coefficients(self, x, y):
		"""Returns the coefficients for the polynomial equation."""

		zeta = self.get_source_position(x=x, y=y)
		if self.specific_frame_derivation:
			calc = self.origin
		else:
			calc = 'general'

		coeff = getc.get_coefficients(calc=calc, zeta=zeta, m3=self.m3,
				m2=self.m2, m1=self.m1, z3=self.z3,
				z2=self.z2, z1=self.z1)

		return coeff

	def get_roots(self, x, y):
		"""Return solutions of polynomial."""

		coefficients = self.get_coefficients(x=x, y=y)

			# Return the roots of the polynomial via the given root finder
		if self.solver == 'SG12':
			rev_list = coefficients[::-1]
			out = _vbbl_SG12_10(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			roots = [out[i] + out[i+10] * 1.j for i in range(10)]
			return np.array(roots)
		elif self.solver == 'zroots':
			rev_list = coefficients[::-1]
			out = _zroots_10(*(rev_list.real.tolist() + rev_list.imag.tolist()))
			roots = [out[i] + out[i+10] * 1.j for i in range(10)]
			return np.array(roots)
		if self.solver == 'numpy':
			roots = np.roots(coefficients).tolist()
			return np.array(roots)
		else:
			raise ValueError('Unknown Solver')

	def get_accepted_solutions(self, x, y):
		
		zeta = self.get_source_position(x=x, y=y)
		roots = self.get_roots(x=x, y=y)
		star_term = self.m1 / np.conj(roots - self.z1)
		plan_term = self.m2 / np.conj(roots - self.z2)
		moon_term = self.m3 / np.conj(roots - self.z3)
		solutions = zeta + star_term + plan_term + moon_term

		accepted_solutions = []
		for (i, root) in enumerate(roots):
			distances_from_root = abs((solutions - root)**2)
			min_distance_arg = np.argmin(distances_from_root)

			if i == min_distance_arg:
				accepted_solutions.append(root)

		return accepted_solutions

	def get_magnification(self, x, y):
		"""Returns the magnification for each configuration."""

		magn = 0
		image_positions = self.get_accepted_solutions(x=x, y=y)
		for z in image_positions:
			z_conj = np.conj(z)
			partial_conj = ((self.m1 / (self.z1 - z)**2) + 
							(self.m2 / (self.z2 - z)**2) +
							(self.m3 / (self.z3 - z)**2))
			
			partial = np.conj(partial_conj)
			detJ = (1. - partial*partial_conj)
			magn += np.abs(1./detJ)
		return magn

	def get_size_caustic(self):
		"""
		Determines the width, height, and position of the center of the
		planetary caustic. Estimated as a binary lens.
		"""

		self.width_caustic = 4.*np.sqrt(self.q1)*(
					1. + 1./(2.*(self.s1**2))) / (self.s1**2)
		self.height_caustic = 4.*np.sqrt(self.q1)*(
					1. - 1./(2.*(self.s1**2))) / (self.s1**2)
		self.xcenter_caustic = 0.5*self.s1 - 1.0/self.s1

### The following functions are used for assigning data into grids.

	def get_position_arrays(self, region, region_lim=None):
		"""
		Fills arrays for the x- and y-position to prepare grid plots.

		Parameters:
			region (string):
				The area on which the grid will be calculated.
				Accepted values:

				'caustic'
					The zoomed-out view showing the full planetary caustic.

				'onax-cusp'
					The zoomed-in view to the right of the cusp on the
					horizontal axis.

				'offax-cusp'
					The zoomed-in view above the cusp on the vertical axis.

				'custom'
					The user-specified region determined by the parameters:
					xmin, xmax, ymin, ymax

			region_lim (tuple):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax).

				xmin, xmax, ymin, ymax (floats):
					The limitget_accepted_solutionss of the grid area when region is 'custom.'
					The on-axis cus1 are given by: (x, y) = {(-1, 0), (1, 0)}
					The off-axis cus1 are given by: (x, y) = {(0, -1), (0, 1)}
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
		if region == 'both':
			region_xmin = -0.5*self.s1
			region_xmax = 0.5*self.s1
			region_ymin = -0.5*self.s1
			region_ymax = 0.5*self.s1
		if region == 'custom':
			(xmin, xmax, ymin, ymax) = (*region_lim,)
			region_xmin = self.xcenter_caustic + 0.5*xmin*self.width_caustic
			region_xmax = self.xcenter_caustic + 0.5*xmax*self.width_caustic
			region_ymin = 0.5*ymin*self.height_caustic
			region_ymax = 0.5*ymax*self.height_caustic

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

		"""
		self.num_images = np.zeros(self.res**2, dtype=int)
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			roots = self.get_roots(x=x, y=y)
			for solution in roots:
				if self.check_solution(solution=solution):
					self.num_images[idx] += 1
		"""

		self.num_images = np.zeros(self.res**2, dtype=int)
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			image_positions = self.get_accepted_solutions(x=x, y=y)
			self.num_images[idx] = len(image_positions)
		

	def get_coeff_array(self):
		"""
		Fills an array for the values of the coefficients through the grid.
		"""

		self.coeff_array = [[]*self.res**2 for i in range(12)]
		for idx in range(self.res**2):
			x = self.x_array[idx]
			y = self.y_array[idx]
			roots = self.get_roots(x=x, y=y)
			coeffs = self.get_coeff_parts(x=x, y=y)
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

	"""
	def get_num_images_errors(self):
		
		Returns only the points in the grid where the number of images is
		not equal to 3 or 5.

		Returns:

			x_err (array of floats):
				The x-position for each erroneous point.

			y_err (array of floats):
				The y-position for each erroneous point.

			num_images_err (array of floats):
				The number of images for each erroneous point.
		

		x_err = []
		y_err = []
		num_images_err = []

		for (i, num) in enumerate(self.num_images):
			if (num == 0) or (num == 1) or (num == 2) or (num == 4):
				x_err.append(self.x_array[i])
				y_err.append(self.y_array[i])
				num_images_err.append(num)

		return (x_err, y_err, num_images_err)
	"""

	def get_ratio_outliers(self, other_BL, ratio_cutoff=None):
		"""
		Returns the values of x, y, and the relative magnification where the
		relative magnification is greater than the specified cutoff limit.

		Required parameters:

			other_BL (object):
				Another instance of the class BinaryLens against which
				magnification values will be compared.


		Optional parameters:

			ratio_cutoff (float):
				The minimum relative magnification a point can have to be
				included in the array. If not specified, the value will default
				to 2.0.
		"""


		(x, y, magn1, magn2) = (self.x_array, self.y_array, self.magn_array,
								other_BL.magn_array)
		rel_magn = magn1 / magn2

		rel_magn_outliers = []
		x_outliers = []
		y_outliers = []

		# If ratio_cutoff has not been specified, default to 2.0
		if ratio_cutoff == None:
			ratio_cutoff = 2.0
		for (i, magn) in enumerate(rel_magn):
			if (magn > ratio_cutoff) or (1./magn > ratio_cutoff):
				x_outliers.append(x[i])
				y_outliers.append(y[i])
				rel_magn_outliers.append(magn)
		return x_outliers, y_outliers, rel_magn_outliers

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

	def get_tstat_plot_data(self, cutoff=None, outliers=False,
			region='caustic', region_lim=None, sample_res=5):
		"""
		Optional parameters:

			cutoff (float):
				The minimum magnification a point can have for the t-stat to
				be plotted. If not specified, the default is to include only
				the 10 greatest	values. Feature only implemented when 
				outliers == True.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			sample_res (int):
				The resolution (side length) of the sample grid created
				around each point in the region grid. These points are used
				to determine the t-stat value. If not specified, the default
				grid size is 5x5.
		"""

		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_coeff_array()
		self.get_tstat_array(sample_res = sample_res)

		if outliers:
			self.get_magnification_outliers(cutoff)
			self.get_tstat_outliers()
			(x, y, magn, tstat) = (self.x_outliers, self.y_outliers,
								   self.magn_outliers, self.tstat_outliers)
		else:
			(x, y, magn, tstat) = (self.x_array, self.y_array, self.magn_array,
								   self.tstat_array)
		return (x, y, magn, tstat)


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

		mean_magn = sum(sample_magn_array) / len(sample_magn_array)
		stderr_magn = np.std(sample_magn_array) / np.sqrt(len(sample_magn_array))
		tstat = abs(point_magn - mean_magn) / stderr_magn
		return tstat

	def get_coeff_parts(self, x, y):
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
### Note: All plots are made in the geometric center frame, regardless
### of which frame the calculations were done in.

	def plot_num_images(self, errors_only=False, print_errors=True,
			region='caustic', save=False, region_lim=None,
			default_settings=True, **kwargs):
		"""
		Creates a plot showing number of images on a grid over the specified
		region.

		Optional parameters:

			errors_only (bool):
				If True, only plots the points where the number of solutions
				is not equal to 3 or 5.

			print_errors (bool):
				If True, prints the number of instances where the number of
				accepted solutions was equal to each integer value from 0 to 5
				inclusive.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		kwargs = self.check_kwargs(**kwargs)
		kwargs['cmap'] = 'coolwarm'
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_num_images_array()

#		if print_errors:
#			self.print_num_images_errors()

		if errors_only:
			(x, y, num_images) = self.get_num_images_errors()
			file_name = '../Tables/NumImErr_{}_{}.png'.format(
					self.origin_file, self.solver_file)
		else:
			(x, y, num_images) = (self.x_array, self.y_array, self.num_images)
			file_name = '../Tables/NumIm_{}_{}.png'.format(
					self.origin_file, self.solver_file)

		plt.scatter(x, y, c=num_images, **kwargs)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		dx = xmax - xmin
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)

		if default_settings:
			im_plot = plt.colorbar()
			im_plot.set_label('Num Images')
			plt.xlabel('X-position of Source', fontsize=12)
			plt.ylabel('Y-position of Source', fontsize=12)
			plt.gcf().set_size_inches(9, 6)
			plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
			plt.suptitle('Number of Images', x=0.435)
			title = ('Frame: {}; Solver: {}; Region: {}'.format(
					self.origin_title, self.solver_title, region))
			plt.title(title, fontsize=11)

			if save:
				self.save_png(file_name=file_name)

	def plot_magnification(self, cutoff=None, log_colorbar=False, 
			outliers=False, region='caustic', region_lim=None, save=False,
			**kwargs):
		"""
		Creates a plot showing the magnification on a grid over the specified
		region.

		Optional parameters:

			cutoff (float):
				The minimum magnification a point can have and still be
				plotted. If not specified, the default is to include only
				the 10 greatest	values. Feature only implemented when 
				outliers == True.

			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		kwargs = self.check_kwargs(log_colorbar, **kwargs)
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_magnification_array()

		# Assign plotting variables data according to whether we include
		# all points or only the outliers.
		if outliers:
			self.get_magnification_outliers(cutoff)
			(x, y, magn) = (self.x_outliers, self.y_outliers, self.magn_outliers)
		else:
			(x, y, magn) = (self.x_array, self.y_array, self.magn_array)

		plt.scatter(x, y, c=magn, **kwargs)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		dx = xmax - xmin
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
		mag_plot = plt.colorbar()
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of Source', fontsize=12)
		plt.ylabel('Y-position of Source', fontsize=12)
		plt.gcf().set_size_inches(9, 6)

		if outliers:
			if cutoff == None:
				cutoff = int(min(magn))
			plt.suptitle('High Magnification', x=0.435)
			title = ('Frame: {}; Solver: {}; Region: {}\n'.format(
					self.origin_title, self.solver_title, region) + 
					'M>{:.0f}'.format(cutoff))
			plt.title(title, fontsize=11)
			file_name = ('../Tables/HighMagn_{}_{}.png'.format(
					self.solver_file, self.origin_file))
		else:
			plt.suptitle('Magnification', x=0.435)
			title = ('Frame: {}; Solver: {}; Region: {}\n'.format(
					self.origin_title, self.solver_title, region))
			plt.title(title, fontsize=11)
			file_name = ('../Tables/Magn_{}_{}.png'.format(self.solver_file,
					self.origin_file))

		if save:
			self.save_png(file_name=file_name)

	def plot_coefficients(self, cutoff=None, log_colorbar=False,
			outliers=False, region='caustic', region_lim=None, save=False,
			**kwargs):
		"""
		Creates a plot showing the magnification on a grid over the specified
		region.

		Optional parameters:

			cutoff (float):
				The minimum magnification a point can have and still be
				plotted. If not specified, the default is to include only
				the 10 greatest	values. Feature only implemented when 
				outliers == True.

			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		kwargs = self.check_kwargs(log_colorbar, **kwargs)
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_magnification_array()
		self.get_coeff_array()
		(x, y, magn, coeff) = (self.x_array, self.y_array, self.magn_array,
							   self.coeff_array)

		if outliers:
			self.get_magnification_outliers(cutoff)
			(x, y, magn) = (self.x_outliers, self.y_outliers,
					self.magn_outliers)
			for i in range(len(magn_out)):
				coeff_temp = self.get_coeff_parts(x=x_out[i], y=y_out[i])
				for k in range(12):
					coeff_out[k].append(coeff_temp[k])
			coeff = coeff_out

		for (i, coeff_val) in enumerate(coeff):
			plt.scatter(x, y, c=coeff_val, **kwargs)
			(xmin, xmax) = (min(self.x_array), max(self.x_array))
			(ymin, ymax) = (min(self.y_array), max(self.y_array))
			dx = xmax - xmin
			plt.xlim(xmin, xmax)
			plt.ylim(ymin, ymax)
			plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
			mag_plot = plt.colorbar()
			mag_plot.set_label('Coefficient value')
			plt.xlabel('X-position of Source', fontsize=12)
			plt.ylabel('Y-position of Source', fontsize=12)
			plt.gcf().set_size_inches(9, 6)

			if outliers:
				if cutoff == None:
					cutoff = int(min(magn))
				plt.suptitle('{} vs. Position (Only Outliers)'.format(self.coeff_string[i]),
							 x=0.435)
				title = ('Frame: {}; Solver: {}; Region: {}\n'.format(
						self.origin_title, self.solver_title, region) + 
						's={}, q={}, M>{:.0f}'.format(self.s, self.q, 
						cutoff))
				plt.title(title, fontsize=11)
				file_name = ('../Tables/HighMagn_{}_{}.png'.format(
						self.solver_file, self.origin_file))
			else:
				plt.suptitle('{} vs. Position'.format(self.coeff_string[i]),
							 x=0.435)
				title = ('Frame: {}; Solver: {}; Region: {}\n'.format(
						self.origin_title, self.solver_title, region) + 
						's={}, q={}'.format(self.s, self.q))
				plt.title(title, fontsize=11)
				file_name = ('../Tables/Magn_{}_{}.png'.format(self.solver_file,
						self.origin_file))
			if save:
				self.save_png(file_name=file_name)
			plt.show()

	def plot_magn_coeff(self, color_num=False, cutoff=None, outliers=False,
				region='caustic', region_lim=None, save=False, **kwargs):
		"""
		Creates a plot showing the magnification on a grid over the specified
		region.

		Optional parameters:

			color_num (bool):
				If True, the plot will have colored points corresponding to
				the number of images for the point.

			cutoff (float):
				The minimum magnification a point can have and still be
				plotted. If not specified, the default is to include only
				the 10 greatest	values. Feature only implemented when 
				outliers == True.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		if 's' not in kwargs:
			kwargs['s'] = 8
		kwargs = self.check_kwargs(**kwargs)
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_magnification_array()
		self.get_num_images_array()
		self.get_coeff_array()
		(x, y, magn, num_images) = (self.x_array, self.y_array, self.magn_array,
									self.num_images)

		if (outliers and not color_num):
			coeff_out = [[] for i in range(12)]
			self.get_magnification_outliers(cutoff)
			(x_out, y_out, magn_out) = (self.x_outliers, self.y_outliers,
							self.magn_outliers)
			for i in range(len(magn_out)):
				coeff_temp = self.get_coeff_parts(x=x_out[i], y=y_out[i])
				for k in range(12):
					coeff_out[k].append(coeff_temp[k])

		for (i, coeff) in enumerate(self.coeff_array):
			if color_num:
				kwargs['c'] = num_images
				kwargs['cmap'] = 'coolwarm'
				plt.scatter(coeff, magn, **kwargs)
				file_name = ('../Tables/mag_coeff_color_{}_{}_{}.png'.format(
						self.coeff_file[i], self.solver_file,
						self.origin_file))
			else:
				kwargs['color'] = 'black'
				plt.scatter(coeff, magn, **kwargs)
				if outliers:
					file_name = ('../Tables/mag_coeff_out_{}_{}_{}.png'.
							format(self.coeff_file[i], self.solver_file,
							self.origin_file))
					kwargs['color'] = 'red'
					plt.scatter(coeff_out[i], magn_out, **kwargs)
				else:
					file_name = ('../Tables/mag_coeff_{}_{}_{}.png'.format(
							self.coeff_file[i], self.solver_file,
							self.origin_file))
			xmin = min(coeff)
			xmax = max(coeff)
			dx = xmax - xmin
			plt.xlabel(self.coeff_string[i], fontsize=12)
			plt.ylabel('Magnification', fontsize=12)
			plt.gcf().set_size_inches(9, 6)
			if color_num:
				mag_plot = plt.colorbar()
				mag_plot.set_label('Number Images')
			plt.xlim(xmin - 0.05*dx, xmax + 0.05*dx)
			plt.ylim(0.90*min(magn), 1.05*max(magn))
			plt.yscale('log')
			plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
			plt.suptitle('Magnification for {}'.format(self.coeff_string[i]),
						 x=(0.515 - 0.08*color_num))
			title = ('{} Solver, {} Frame\n'.format(self.solver_title,
					self.origin_title)) + ('Region: {}, s={}, q={}'.format(
					region, self.s, self.q))
			plt.title(title, fontsize=11)
			if save:
				self.save_png(file_name=file_name)
			plt.show()

	def plot_num_images_coeff(self, color_magn=False, log_colorbar=False, 
				region='caustic', region_lim=None, save=False, **kwargs):
		"""
		Creates a plot showing the magnification on a grid over the specified
		region.

		Optional parameters:

			color_magn (bool):
				If True, the plot will have colored points corresponding to
				the magnification of the point.

			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear. Only
				implemented if color_magn==True.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		if 's' not in kwargs:
			kwargs['s'] = 8
		kwargs = self.check_kwargs(log_colorbar=log_colorbar, **kwargs)
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_magnification_array()
		self.get_num_images_array()
		self.get_coeff_array()

		(x, y, num_images) = (self.x_array, self.y_array, self.num_images)
		file_name = '../Tables/NumIm_{}_{}.png'.format(
				self.origin_file, self.solver_file)

		for (i, coeff) in enumerate(self.coeff_array):
			if color_magn:
				kwargs['c'] = self.magn_array
			else:
				kwargs['color'] = 'black'
			plt.scatter(coeff, num_images, **kwargs)
			xmin = min(coeff)
			xmax = max(coeff)
			dx = xmax - xmin
			plt.xlabel(self.coeff_string[i], fontsize=12)
			plt.ylabel('Num Images', fontsize=12)
			plt.gcf().set_size_inches(9, 6)
			if color_magn:
				region_plot = plt.colorbar()
				region_plot.set_label('Magnification')
			plt.xlim(xmin - 0.05*dx, xmax + 0.05*dx)
			plt.ylim(0, 5.05)
			plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
			plt.suptitle('Num Images for {}'.format(self.coeff_string[i]),
						 x = (0.515 - 0.08*color_magn))
			title = ('{} Solver, {} Frame\n'.format(self.solver_title,
					self.origin_title)) + ('Region: {}, s={}, q={}'.format(
					region, self.s, self.q))
			plt.title(title, fontsize=11)
			if save:
				self.save_png(file_name=file_name)
			plt.show()

	def plot_coeff_tstat(self, cutoff=None, outliers=False,	region='caustic',
						region_lim=None, sample_res=5, save=False, **kwargs):
		"""
		Creates a plot showing the t-stat vs. coefficient value for each of 
		the 12 coefficients.

		Optional parameters:

			cutoff (float):
				The minimum magnification a point can have for the t-stat to
				be plotted. If not specified, the default is to include only
				the 10 greatest	values. Feature only implemented when 
				outliers == True.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			sample_res (int):
				The resolution (side length) of the sample grid created
				around each point in the region grid. These points are used
				to determine the t-stat value. If not specified, the default
				grid size is 5x5.

			save (bool):
				If True, saves all 12 plots to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.

		Note: This feature has been determined not to be particularly useful
		or informative.
		"""

		if 's' not in kwargs:
			kwargs['s'] = 5
		kwargs = self.check_kwargs(**kwargs)
		(x, y, magn, tstat) = self.get_tstat_plot_data(cutoff=cutoff,
				outliers=outliers, region=region, region_lim=region_lim,
				sample_res=sample_res)
		if cutoff==None:
			cutoff = int(min(magn))

		print('Plotting...')
		for (i, coeff) in enumerate(self.coeff_array):
			file_name = ('../Tables/BL_tstat_{}_{}_{}.png'.format(
					self.coeff_file[i], self.solver_file, self.origin_file))
			plt.scatter(coeff, self.tstat_array, color='black', **kwargs)
			xmin = min(coeff)
			xmax = max(coeff)
			dx = xmax - xmin
			ymax = max(self.tstat_array)
			plt.xlabel(self.coeff_string[i], fontsize=12)
			plt.ylabel('t-Test Result', fontsize=12)
			plt.gcf().set_size_inches(9, 6)
			plt.xlim(xmin - 0.01*dx, xmax + 0.01*dx)
			plt.ylim(0, 1.01*ymax)
			plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
			plt.suptitle('t-Test Result vs. {}'.format(self.coeff_string[i]),
						  x=0.515)
			title = ('{} Solver, {} Frame\nRegion: {}, M>{}, s={}, q={}'.
					format(self.solver_title, self.origin_title, region,
					cutoff, self.s, self.q))
			plt.title(title, fontsize=11)
			if save:
				self.save_png(file_name=file_name)
			plt.show()

	def plot_position_tstat(self, cutoff=None, outliers=False,
			region='caustic', region_lim=None, sample_res=5, save=False,
			**kwargs):
		"""
		Creates a plot showing the t-stat vs. position in the specified region.

		Optional parameters:

			cutoff (float):
				The minimum magnification a point can have for the t-stat to
				be plotted. If not specified, the default is to include only
				the 10 greatest	values. Feature only implemented when 
				outliers == True.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			sample_res (int):
				The resolution (side length) of the sample grid created around
				each point in the region grid. These points are used to
				determine the t-stat value. If not specified, the default is 5.

			save (bool):
				If True, saves all 12 plots to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		kwargs = self.check_kwargs(**kwargs)
		(x, y, magn, tstat) = self.get_tstat_plot_data(cutoff=cutoff,
				outliers=outliers, region=region, region_lim=region_lim,
				sample_res=sample_res)
		if cutoff==None:
			cutoff = int(min(magn))

		plt.scatter(x, y, c=tstat, **kwargs)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		dx = xmax - xmin
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
		plt.xlabel('X-position', fontsize=12)
		plt.ylabel('Y-position', fontsize=12)
		plt.gcf().set_size_inches(9, 6)
		region_plot = plt.colorbar()
		region_plot.set_label('t-value')
		if outliers:
			plt.suptitle('t-test Score vs. Position (Outliers Only)', x=0.435)
			title('{} Solver; {} Frame\nRegion: {}, M>{:.0f}, s={}, q={}'.
					format(self.solver_title, self.origin_title, region,
					cutoff, self.s, self.q))
			plt.title(title, fontsize=11)
		else:
			plt.suptitle('t-test Score vs. Position', x=0.435)
			title = ('{} Solver; {} Frame\nRegion: {}, s={}, q={}'.
					format(self.solver_title, self.origin_title, region,
					self.s, self.q))
			plt.title(title, fontsize=11)
		caustic = mm.Caustics(s=self.s, q=self.q)
		caustic.plot(s=1)
		plt.show()

	# The folowing functions require 2 instances of the class to plot

	def plot_rel_magnification(self, other_BL, log_colorbar=False,
			outliers=False, ratio_cutoff=None, region='caustic',
			region_lim=None, save=False, **kwargs):
		"""
		Creates a plot showing the relative magnification of two instances
		of the class on a grid over the specified region.

		Required parameters:

			other_BL (object):
				Another instance of the class BinaryLens against which
				magnification values will be compared.

		Optional parameters:

			ratio_cutoff (float):
				The minimum relative magnification a point can have and still
				be plotted. If not specified, the value will default to that
				assigned in get_ratio_outliers. Feature only implemented when
				outliers == True.

			log_colorbar (bool):
				If True, the magnification colorbar will go on log
				scale. Otherwise, the scale will be linear.

			outliers (bool):
				If true, only the points with magnification above the cutoff
				will be plotted.

			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves the plot to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		kwargs = self.check_kwargs(log_colorbar, **kwargs)
		kwargs['s'] *= 0.8
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_magnification_array()
		other_BL.get_position_arrays(region=region, region_lim=region_lim)
		other_BL.get_magnification_array()

		# Assign the appropriate data, based on whether we want to include all
		# the data, or just the outliers.

		if outliers:
			(x, y, rel_magn) = self.get_ratio_outliers(other_BL, ratio_cutoff)
			file_name = ('../Tables/RelMagnOut_{}_{}.png'.format(
						 self.solver_file, self.origin_file))
		else:
			(x, y, magn1, magn2) = (self.x_array, self.y_array, self.magn_array,
									other_BL.magn_array)
			rel_magn = (magn1 / magn2)
			file_name = ('../Tables/RelMagn_{}_{}.png'.format(self.solver_file,
						 self.origin_file))

		# If no outliers were found, this aborts the function.
		if len(x) == 0:
			print('No outliers found. Continuing with next plot...')
			return

		plot = plt.scatter(x, y, c = rel_magn, **kwargs)
		rel_plot = plt.colorbar()
		rel_plot.set_label('Magnification Ratio')
		plt.xlabel('X-position of Source', fontsize=12)
		plt.ylabel('Y-position of Source', fontsize=12)
		plt.gcf().set_size_inches(9, 6)
		(xmin, xmax) = (min(self.x_array), max(self.x_array))
		(ymin, ymax) = (min(self.y_array), max(self.y_array))
		dx = xmax - xmin
		plt.xlim(xmin, xmax)
		plt.ylim(ymin, ymax)
		plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
		plt.suptitle('Relative Magnification', x=0.435)
		title = ('({} Solver, {} Frame) / ({} Solver, {} Frame)\n'.format(
				self.solver_title, self.origin_title, other_BL.solver_title,
				other_BL.origin_title)) + ('Region: {}, s={}, q={}'.format(
				region, self.s, self.q))
		plt.title(title, fontsize=11)
							
		if save:
			self.save_png(file_name=file_name)

	def plot_rel_magn_coeff(self, other_BL, outliers=True, ratio_cutoff=None,
				region='caustic', region_lim=None, save=False, **kwargs):
		"""
		Creates a plot showing the relative magnification of two instances
		of the class vs. the coefficient value on a grid over the specified
		region.

		Required parameters:

			other_BL (object):
				Another instance of the class BinaryLens against which
				magnification values will be compared.

		Optional parameters:

			outliers (bool):
				If True, the points with magnification outside the range
				determined by ratio_cutoff will be colored red, while the
				rest will be colored black.

			ratio_cutoff (float):
				The relative magnification cutoff. Points with relative
				magnification outside the range will be colored red, while 
				the rest will be colored black. If not specified, the value
				will default to that assigned in get_ratio_outliers. Feature
				does not work if outliers==False.

			region (string):
				The region in which the grid is to be filled.m See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.

			save (bool):
				If True, saves all 12 plots to a .png file.

			**kwargs (dictionary):
				Keyword arguments for pyplot module.
		"""

		# Get data for plotting
		if 's' not in kwargs:
			kwargs['s'] = 8
		kwargs = self.check_kwargs(**kwargs)
		self.get_position_arrays(region=region, region_lim=region_lim)
		self.get_magnification_array()
		self.get_coeff_array()
		other_BL.get_position_arrays(region=region, region_lim=region_lim)
		other_BL.get_magnification_array()

		rel_magn = self.magn_array / other_BL.magn_array

		if outliers:
			coeff_out = [[] for i in range(12)]
			(x_out, y_out, rel_magn_out) = self.get_ratio_outliers(
					other_BL, ratio_cutoff)

			for i in range(len(rel_magn_out)):
				coeff_temp = self.get_coeff_parts(x=x_out[i], y=y_out[i])
				for k in range(12):
					coeff_out[k].append(coeff_temp[k])

			if len(x_out) == 0:
				print('No outliers found. Continuing with next plot...')
				return

		for (i, coeff) in enumerate(self.coeff_array):
			kwargs['color'] = 'black'
			plt.scatter(coeff, rel_magn, **kwargs)
			if outliers:
				kwargs['color'] = 'red'
				plt.scatter(coeff_out[i], rel_magn_out, **kwargs)
				file_name = ('../Tables/relcoeff_out_{}_{}_{}.png'.format(
						self.coeff_file[i], self.solver_file,
						self.origin_file))

			else:
				file_name = ('../Tables/relcoeff_{}_{}_{}.png'.format(
						self.coeff_file[i], self.solver_file,
						self.origin_file))

			xmin = min(coeff)
			xmax = max(coeff)
			dx = xmax - xmin
			plt.xlabel(self.coeff_string[i], fontsize=12)
			plt.ylabel('Relative Magnification', fontsize=12)
			plt.gcf().set_size_inches(9, 6)
			plt.xlim(xmin - 0.05*dx, xmax + 0.05*dx)
			plt.ylim(0.90*min(rel_magn), 1.10*max(rel_magn))
			plt.yscale('log')
			plt.xticks(np.arange(xmin, xmin + 1.2*dx, dx / 4))
			plt.suptitle('Relative Magnification for {}'.format(
						 self.coeff_string[i]), x=0.515)
			title = ('({} Solver, {} Frame) / ({} Solver, {} Frame)\n'.format(
					self.solver_title, self.origin_title, other_BL.solver_title,
					other_BL.origin_title)) + ('Region: {}, s={}, q={}'.format(
					region, self.s, self.q))
			plt.title(title, fontsize=11)
			if save:
				self.save_png(file_name=file_name)
			plt.show()

### The following functions save the relevant data to a .png or .fits file

	def write_to_fits(self, region='caustic', region_lim=None):
		"""
		Writes grid data (x, y, magnification, number of images) at each
		point to a .fits table. Useful for comparing results between
		different solvers and coordinate frames.

		Optional parameters:
			region (string):
				The region in which the grid is to be filled. See the
				function, "get_position_arrays," to see accepted values.

			region_lim (tuple containing floats):
				The limits of the grid area when region is given by 'custom.'
				It is given by (xmin, xmax, ymin, ymax). See the method, 
				'get_position_arrays' to see what the limits correspond to.				
		"""

		self.get_position_arrays(region=region, region_lim=region_lim)
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

		# This checks if file of same name already exists. If so, it 
		# increases the integer value preceding the .fits extension.
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


	def print_image_position(self, print_input=True):
		"""
		Prints the positions of those images which pass the check_solutions
		test, for a single test point.

		Optional parameters:

			print_input (bool):
				If True, prints the user-specified input parameters before
				displaying the image positions.
		"""

#		if print_input:
#			self.print_input()
		self.get_roots(x=self.x, y=self.y)
		image_positions = self.get_accepted_solutions(x=self.x, y=self.y)
		print('Image locations:')
		for (i, pos) in enumerate(image_positions):
			print('{:.5f}'.format(pos))

	def print_magnification(self, print_input=True):
		"""
		Prints the total magnification for single test point.

		Parameters:

			print_input (bool):
				If True, prints the user-specified input parameters before
				displaying the image positions.
		"""

#		if print_input:
#			self.print_input()
		magn = self.get_magnification(x=self.x, y=self.y)
		print('Magnification: {:.5f}'.format(magn))

### The following functions gather string data that are used within the source code.

	def get_coeff_strings(self):
		"""
		Assigns strings to global lists regarding coefficient names.
		"""

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
			self.origin_title = 'Geometric Center'
			self.origin_phrase = 'geometric center frame'
		elif self.origin == 'body2':
			self.origin_file = 'b2'
			self.origin_title = 'Body2'
			self.origin_phrase = 'body2 frame'
		elif self.origin == 'body3':
			self.origin_file = 'b3'
			self.origin_title = 'Body3'
			self.origin_phrase = 'body3 frame'
		else:
			raise ValueError('Unknown coordinate system: {:}'.format(self.origin))

		if self.solver == 'numpy':
			self.solver_file = 'np'
			self.solver_title = 'Numpy'
			self.solver_phrase = 'Numpy root finder'
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


