# Zoey Samples
# Created: Jul 03, 2018
# Caustics.py
# Last Updated: Jul 09, 2018

import sys
from pathlib import Path
import ctypes
import os
import numpy as np
import matplotlib.pyplot as plt
import math
from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL


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
	vbbl.VBBL_SG12_4.argtypes = 10 * [ctypes.c_double]
	vbbl.VBBL_SG12_4.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, 
		shape=(8,))
	_vbbl_SG12_4 = vbbl.VBBL_SG12_4

	vbbl.VBBL_SG12_6.argtypes = 14 * [ctypes.c_double]
	vbbl.VBBL_SG12_6.restype = np.ctypeslib.ndpointer(dtype=ctypes.c_double, 
		shape=(12,))
	_vbbl_SG12_6 = vbbl.VBBL_SG12_6


class Caustics(object):
	"""
	Writes a set of points that correspond to the caustic in the source
	plane and the critical curve in the image plane.

	Attributes:

		Required:
			lens (object):
				An instance of the BinaryLens or TripleLens class from
				which all Caustics class variables will be calculated.

		Global Variables:
			lens_type (string):
				The type of lens ab instance of the class represents.
				Options are 'BL' or 'TL' for a binary lens and a triple
				lens, respectively.

			origin (string):
				The coordinate frame in which calculations are carried out.

			solver (string):
				Determines which root finder to use to solve the polynomial.

			If lens_type is 'BL':
				s (float):
					The separation between the two bodies.

				q (float):
					The mass ratio of the smaller body to the larger body.

				dm (float):
					Half the difference of the masses (positive).

				m (float):
					Half the sum of the masses.

				z1 (float):
					The planet's position (on the real axis by definition).

				z2 (float):
					The star's position (on the real axis by definition).

			If lens_type is 'TL':
				system (string):
					The types of bodies in the triple system. Options are:

						'SPM' - A system with a star, a planet, and a moon.
								{z1,m1}: star; {z2,m2}: plan; {z3,m3}: moon
						'SPP' - A system with a star and 2 planets
								{z1,m1}: star, {z2,m2}: plan1, {z3,m3}: plan2
						'SSP' - A system with 2 stars and a planet
								{z1,m1}: star1, {z2,m2}: star2, {z3,m3}: plan

				q1 (float):
					If system is 'SPM': The mass ratio of the planet to the star.
					If system is 'SPP': The mass ratio of planet1 to the star.
					If system is 'SSP': The mass ratio of star2 to star1.

				q2 (float):
					If system is 'SPM': The mass ratio of the moon to the planet.
					If system is 'SPP': The mass ratio of planet2 to the star.
					If system is 'SSP': The mass ratio of the planet to star1.

				s1 (float):
					The separation between body1 and body2 in units of the
					whole system's Einstein radius.

				s2 (float):
					The separation between body2 and body3 in units of the
					whole system's Einstein radius.

				phi (float):
					If system is 'SPM': The angle formed by star-planet-moon.
					If system is 'SPP': The angle formed by planet1-star-planet2.
					If system is 'SSP': The angle formed by star2-star1-planet.

					*phi is expressed in radians.

				m1 (float):
					The mass of the largest body as a fraction of the whole
					system's mass.

				m2 (float):
					The mass of 2nd largest body as a fraction of the whole
					system's mass.

				m3 (float):
					The mass of the least massive body as a fraction of
					the whole system's mass.

				z1 (complex):
					The position of the most massive body (the star in SPM or
					SPP system, or star1 in SSP system) in units of the total
					system's Einstein radius.

				z2 (complex):
					The position of the 2nd most massive body in units of
					the total system's Einstein radius.

				z3 (complex):
					The position of the least massive body in units of the
					total system's Einstein radius.

		**See the BinaryLens and TripleLens source codes for more detailed
		  explanations of these variables.
	"""


	def __init__(self, lens, solver='SG12'):

		# Determine if the variable type for lens is valid.
		if not isinstance(lens, BL) and not isinstance(lens, TL):
			raise TypeError('wrong type of lens: {:}'.format(type(lens)))
		self.solver = solver
		self.lens = lens
		self.plot_frame = self.lens.plot_frame
		self.lens.get_caustic_param(refine_region=self.lens.refine_region)

		if isinstance(self.lens, BL):
			self.lens_type = 'BL'
			self.s = self.lens.s
			self.q = self.lens.q
			(self.m, self.dm) = (self.lens.m, self.lens.dm)
			self.get_BL_lensing_body_positions()

		elif isinstance(self.lens, TL):
			self.lens_type = 'TL'
			self.s1 = self.lens.s1
			self.s2 = self.lens.s2_actual
			self.q1 = self.lens.q1
			self.q2 = self.lens.q2
			self.phi = self.lens.phi

			#FIXME: This initialization works for SPM systems, but the
			# one below, which should be equivalent, does not work.
			denominator = 1. + self.q1 + self.q2*self.q1
			self.m1 = 1. / denominator
			self.m2 = self.q1 / denominator
			self.m3 = self.q2*self.q1 / denominator
			#self.m3 = self.q2 / denominator
			#(self.m1, self.m2, self.m3) = (self.lens.m1, self.lens.m2, self.lens.m3)
			self.get_TL_lensing_body_positions()

	def get_BL_lensing_body_positions(self):

		xshift = self.lens.xshift + self.lens._custom_xshift
		yshift = self.lens.yshift + self.lens._custom_yshift
		#(xcenter, ycenter) = (self.lens.xcenter_caustic, self.lens.ycenter_caustic)

		# Come back and fix this assignment. Assign it according to geo_cent. If
		# plot_frame == 'caustic', apply shift.
		if self.plot_frame == 'caustic':
			self.z1 = 0.5*self.s - (xshift + 1j*(yshift))
			self.z2 = -0.5*self.s - (xshift + 1j*(yshift))
		elif self.plot_frame == 'geo_cent':
			self.z1 = 0.5*self.s
			self.z2 = -0.5*self.s

	def get_TL_lensing_body_positions(self):

		xshift = self.lens.xshift + self.lens._custom_xshift
		yshift = self.lens.yshift + self.lens._custom_yshift
		#(xcenter, ycenter) = (self.lens.xcenter_caustic, self.lens.ycenter_caustic)

		if self.lens.system == 'SPM':
			self.displacement23 = self.s2*(math.cos(math.pi - self.phi) +
									   1.j*math.sin(math.pi - self.phi))

			self.z1 = -0.5*self.s1 + 0j
			self.z2 = 0.5*self.s1 + 0j
			self.z3 = self.z2 + self.displacement23

			if self.plot_frame == 'caustic':
				self.z1 -= xshift + 1j*yshift
				self.z2 -= xshift + 1j*yshift
				self.z3 -= xshift + 1j*yshift

		elif (self.lens.system == 'SPP') or (self.lens.system == 'SSP'):
			self.displacement13 = self.s2*(math.cos(self.phi) +
									   1.j*math.sin(self.phi))

			self.z1 = -0.5*self.s1 + 0j
			self.z2 = 0.5*self.s1 + 0j
			self.z3 = self.z1 + self.displacement13

			if self.plot_frame == 'caustic':
				self.z1 -= xshift + 1j*yshift
				self.z2 -= xshift + 1j*yshift
				self.z3 -= xshift + 1j*yshift

	def plot_caustic(self, points=5000, **kwargs):

		self.calculate(points=points)
		plt.scatter(self.x, self.y, **kwargs)

	def plot_lens_bodies(self):

		s = 70
		plt.scatter(self.z1.real, self.z1.imag, marker='*', s=10*s, color='orange', lw=0)
		plt.scatter(self.z2.real, self.z2.imag, marker='o', s=4*s, color='green', lw=0)
		plt.scatter(self.z3.real, self.z3.imag, marker='o', s=s, color='yellow', lw=0)

	def plot_critical_curve(self, points=5000, **kwargs):

		self.calculate(points=points)
		plt.scatter(self.CriticalCurve.x, self.CriticalCurve.y, **kwargs)

	def calculate(self, points=5000):

		self.x = []
		self.y = []
		self.critical_curve = self.CriticalCurve()

		if self.lens_type == 'BL':

			(z1, z2, m, dm) = (self.z1, self.z2, self.m, self.dm)
			for alpha in np.linspace(0, 2*np.pi, points, endpoint=False):
				exp_ialpha = np.e**(alpha*1j)

				### Binary Lens, general case
				coeff4 = 1.
				coeff3 = (-2*z1 - 2*z2)
				coeff2 = (z1**2 + 4*z1*z2 + z2**2 - 2*m*exp_ialpha)
				coeff1 = (-2*z1**2*z2 - 2*z1*z2**2 + 2*dm*z1*exp_ialpha - 2*dm*z2*exp_ialpha + 2*m*(z1 + z2)*exp_ialpha)
				coeff0 = (z1**2*z2**2 - dm*z1**2*exp_ialpha - m*z1**2*exp_ialpha + dm*z2**2*exp_ialpha - m*z2**2*exp_ialpha)

				coefficients = np.array([coeff4, coeff3, coeff2, coeff1, coeff0])

				if self.solver == 'numpy':
					roots = np.roots(coefficients)
				elif self.solver == 'SG12':
					rev_list = coefficients[::-1]
					out = _vbbl_SG12_4(*(rev_list.real.tolist() + rev_list.imag.tolist()))
					roots = [out[i] + out[i+4] * 1.j for i in range(4)]
				else:
					raise ValueError('Unknown solver {:}'.format(solver))

				for root in roots:
					self.critical_curve.x.append(root.real)
					self.critical_curve.y.append(root.imag)

					caustic_point = self.solve_lens_equation(root)
					self.x.append(caustic_point.real)
					self.y.append(caustic_point.imag)

		elif self.lens_type == 'TL':

			(z1, z2, z3, m1, m2, m3) = (self.z1, self.z2, self.z3, self.m1,
										self.m2, self.m3)
			for alpha in np.linspace(0, 2*np.pi, points, endpoint=False):
				exp_ialpha = np.e**(alpha*1j)

				### Triple Lens, general case
				coeff6 = 1.
				coeff5 = (-2*z1 - 2*z2 - 2*z3)
				coeff4 = (z1**2 + 4*z1*z2 + z2**2 + 4*z1*z3 + 4*z2*z3 + z3**2 - exp_ialpha)
				coeff3 = (-2*z1**2*z2 - 2*z1*z2**2 - 2*z1**2*z3 - 8*z1*z2*z3 - 2*z2**2*z3 - 2*z1*z3**2 - 2*z2*z3**2 + 2*m2*z1*exp_ialpha + 2*m3*z1*exp_ialpha + 2*m1*z2*exp_ialpha + 2*m3*z2*exp_ialpha + 2*m1*z3*exp_ialpha + 2*m2*z3*exp_ialpha) 
				coeff2 = (z1**2*z2**2 + 4*z1**2*z2*z3 + 4*z1*z2**2*z3 + z1**2*z3**2 + 4*z1*z2*z3**2 + z2**2*z3**2 - m2*z1**2*exp_ialpha - m3*z1**2*exp_ialpha - 4*m3*z1*z2*exp_ialpha - m1*z2**2*exp_ialpha - m3*z2**2*exp_ialpha - 4*m2*z1*z3*exp_ialpha - 4*m1*z2*z3*exp_ialpha - m1*z3**2*exp_ialpha - m2*z3**2*exp_ialpha)
				coeff1 = (-2*z1**2*z2**2*z3 - 2*z1**2*z2*z3**2 - 2*z1*z2**2*z3**2 + 2*m3*z1**2*z2*exp_ialpha + 2*m3*z1*z2**2*exp_ialpha + 2*m2*z1**2*z3*exp_ialpha + 2*m1*z2**2*z3*exp_ialpha + 2*m2*z1*z3**2*exp_ialpha + 2*m1*z2*z3**2*exp_ialpha)
				coeff0 = (z1**2*z2**2*z3**2 - m3*z1**2*z2**2*exp_ialpha - m2*z1**2*z3**2*exp_ialpha - m1*z2**2*z3**2*exp_ialpha)

				coefficients = np.array([coeff6, coeff5, coeff4, coeff3, coeff2, coeff1, coeff0])

				if self.solver == 'numpy':
					roots = np.roots(coefficients)
				elif self.solver == 'SG12':
					rev_list = coefficients[::-1]
					out = _vbbl_SG12_6(*(rev_list.real.tolist() + rev_list.imag.tolist()))
					roots = np.array([out[i] + out[i+6] * 1.j for i in range(6)])
				else:
					raise ValueError('Unknown solver {:}'.format(solver))

				for root in roots:
					self.critical_curve.x.append(root.real)
					self.critical_curve.y.append(root.imag)

					caustic_point = self.solve_lens_equation(root)
					self.x.append(caustic_point.real)
					self.y.append(caustic_point.imag)

	"""
	### Binary Lens, planet
	coeff4 = 1
	coeff3 = 2*s
	coeff2 = (s**2 - exp_ialpha)
	coeff1 = -((2*q*s*exp_ialpha)/(1 + q))
	coeff0 = -((q*s**2*exp_ialpha)/(1 + q))

	### Triple Lens, planet
	coeff6 = 1.
	coeff5 = 2*(s1 - s2)
	coeff4 = (s1**2 - 4*s1*s2 + s2**2 - exp_ialpha) 
	coeff3 = 2*(s1*s2*(-s1 + s2) + (-(m3*s1) + m1*s2 + m2*(-s1 + s2))*exp_ialpha) 
	coeff2 = (s1**2*s2**2 - (m3*s1**2 + m1*s2**2 + m2*(s1**2 - 4*s1*s2 + s2**2))*exp_ialpha)
	coeff1 = 2*m2*s1*(s1 - s2)*s2*z*exp_ialpha
	coeff0 = -m2*s1**2*s2**2*exp_ialpha
	"""
	"""
	r_array = np.linspace(1e-5,10,5000)
	for r in r_array:
		d = self.s
		m1 = (self.m - self.dm)
		m2 = (self.m + self.dm)

		G = r**2 - d**2
		H = 2*d*r
		g = 1. - m1**2 / r**4
		h = 2*m1*m2 / r**2

		coeff2 = H**2*g - 2.*d**2*h
		coeff1 = H*(h-2*G*g)
		coeff0 = G*(G*g - h) + 2.*d**2*h - m2**2

		coefficients = np.array([coeff2, coeff1, coeff0])
		cosine_phi_solns = np.roots(coefficients)
		phi = []
		for cosine_phi in cosine_phi_solns:
			try:
				phi.append(math.acos(cosine_phi))
			except:
				break

		for angle in phi:
			root = r*math.cos(angle) + 1j*r*math.sin(angle)
			self.x.append(root.real)
			self.y.append(root.imag)

			self.critical_curve.x.append(root.real)
			self.critical_curve.y.append(root.imag)

			caustic_point = self.solve_lens_equation(root)
			self.x.append(caustic_point.real)
			self.y.append(caustic_point.imag)
		"""

	def solve_lens_equation(self, image_position):
		z = image_position
		z_conj = np.conj(image_position)

		if self.lens_type == 'BL':
			(m, dm, z1, z2) = (self.m, self.dm, self.z1, self.z2)
			(z1_conj, z2_conj) = (np.conj(z1), np.conj(z2))
			caustic_point = (z + (m - dm) / (z1_conj - z_conj) +
								 (m + dm) / (z2_conj - z_conj))
			return caustic_point
		elif self.lens_type == 'TL':
			(z1, z2, z3, m1, m2, m3) = (self.z1, self.z2, self.z3,
										self.m1, self.m2, self.m3)
			(z1_conj, z2_conj, z3_conj) = (np.conj(z1), np.conj(z2), np.conj(z3))
			caustic_point = (z + m1 / (z1_conj - z_conj) +
								 m2 / (z2_conj - z_conj) +
								 m3 / (z3_conj - z_conj))

			return caustic_point

	class CriticalCurve(object):

		def __init__(self):
			self.x = []
			self.y = []
