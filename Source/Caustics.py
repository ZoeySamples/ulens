# Zoey Samples
# Created: Jul 03, 2018
# Caustics.py
# Last Updated: Jul 03, 2018

import numpy as np
import matplotlib.pyplot as plt
import math
from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL

class Caustics(object):

	def __init__(self, lens, lens_type=None, plot_frame='caustic'):


		self.plot_frame = plot_frame
		self.lens = lens

		if isinstance(self.lens, BL):
			self.lens_type = 'BL'
			self.s = self.lens.s
			self.q = self.lens.q
			(self.m, self.dm) = (self.lens.m, self.lens.dm)
			if self.plot_frame == 'caustic':
				self.z1 = 1./self.s
				self.z2 = -self.s + 1./self.s
			elif self.plot_frame == 'geo_cent':
				self.z1 = 0.5*self.s
				self.z2 = -0.5*self.s

		elif isinstance(self.lens, TL):
			self.lens_type = 'TL'
			self.s1 = self.lens.s1
			self.s2 = self.lens.s2
			self.q1 = self.lens.q1
			self.q2 = self.lens.q2
			self.phi = self.lens.phi
			self.get_TL_mass()
			self.get_TL_lensing_body_positions()

	def get_TL_mass(self):

		denominator = 1. + self.q1 + self.q2*self.q1
		self.m1 = 1. / denominator
		self.m2 = self.q1 / denominator
		self.m3 = self.q2*self.q1 / denominator

		if self.lens.system == 'SPM':
			# Convert the separation between the moon and planet into units
			# of the total system's Einstein radius
			self.s2 *= np.sqrt((self.m3 + self.m2))

	def get_TL_lensing_body_positions(self):

		print(self.phi)

		if self.lens.system == 'SPM':
			self.displacementMP = self.s2*(math.cos(math.pi - self.phi) +
										1.j*math.sin(math.pi - self.phi))

			if self.plot_frame == 'geo_cent':
				self.z1 = -0.5*self.s1 + 0j
				self.z2 = 0.5*self.s1 + 0j
				self.z3 = self.z2 + self.displacementMP

			elif self.plot_frame == 'body2':
				self.z1 = -self.s1 + 0j
				self.z2 = 0j
				self.z3 = self.displacementMP

			elif self.plot_frame == 'body3':
				self.z1 = -self.s1 - self.displacementMP
				self.z2 = -self.displacementMP
				self.z3 = 0j

			elif self.plot_frame == 'caustic':
				self.z1 = -self.s1 + 1./self.s1
				self.z2 = 1./self.s1
				self.z3 = self.z2 + self.displacementMP

		"""
		self.lens_type = lens_type
		if lens_type == 'BL':
			param = ({'s': s, 'q': q, 'origin': plot_frame, 'solver': 'numpy'})
			self.lens = BL(**param)
			(z1, z2, m, dm) = (self.lens.z1, self.lens.z2, self.lens.m, self.lens.dm)
		if lens_type == 'TL':
			param = ({'s1': s1, 's2': s2, 'q1': q1,  'q2': q2, 'phi': phi, 'origin': originTL, 'solver': solver}) 
			test_TLlens = TL(**param)
		"""

	def plot_caustic(self, points=5000, **kwargs):

		self.calculate(points=points)
		plt.scatter(self.x, self.y, **kwargs)

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

				roots = np.roots([coeff4, coeff3, coeff2, coeff1, coeff0])
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
				coeff4 = (z1**2 + 4*z1*z2 + z2**2 + 4*z1*z3 + 4*z2*z3 + z3**2 - m1*exp_ialpha - m2*exp_ialpha - m3*exp_ialpha) 
				coeff3 = (-2*z1**2*z2 - 2*z1*z2**2 - 2*z1**2*z3 - 8*z1*z2*z3 - 2*z2**2*z3 - 2*z1*z3**2 - 2*z2*z3**2 + 2*m2*z1*exp_ialpha + 2*m3*z1*exp_ialpha + 2*m1*z2*exp_ialpha + 2*m3*z2*exp_ialpha + 2*m1*z3*exp_ialpha + 2*m2*z3*exp_ialpha) 
				coeff2 = (z1**2*z2**2 + 4*z1**2*z2*z3 + 4*z1*z2**2*z3 + z1**2*z3**2 + 4*z1*z2*z3**2 + z2**2*z3**2 - m2*z1**2*exp_ialpha - m3*z1**2*exp_ialpha - 4*m3*z1*z2*exp_ialpha - m1*z2**2*exp_ialpha - m3*z2**2*exp_ialpha - 4*m2*z1*z3*exp_ialpha - 4*m1*z2*z3*exp_ialpha - m1*z3**2*exp_ialpha - m2*z3**2*exp_ialpha)
				coeff1 = (-2*z1**2*z2**2*z3 - 2*z1**2*z2*z3**2 - 2*z1*z2**2*z3**2 + 2*m3*z1**2*z2*exp_ialpha + 2*m3*z1*z2**2*exp_ialpha + 2*m2*z1**2*z3*exp_ialpha + 2*m1*z2**2*z3*exp_ialpha + 2*m2*z1*z3**2*exp_ialpha + 2*m1*z2*z3**2*exp_ialpha)
				coeff0 = (z1**2*z2**2*z3**2 - m3*z1**2*z2**2*exp_ialpha - m2*z1**2*z3**2*exp_ialpha - m1*z2**2*z3**2*exp_ialpha)

				roots = np.roots([coeff6, coeff5, coeff4, coeff3, coeff2, coeff1, coeff0])
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

	def solve_lens_equation(self, image_position):
		z = image_position
		z_conj = np.conj(image_position)

		if self.lens_type == 'BL':
			(m, dm, z1, z2) = (self.m, self.dm, self.z1, self.z2)
			(z1_conj, z2_conj) = (np.conj(z1), np.conj(z2))
			caustic_point = z + (m - dm) / (z1_conj - z_conj) + (m + dm) / (z2_conj - z_conj)
			return caustic_point
		elif self.lens_type == 'TL':
			(z1, z2, z3, m1, m2, m3) = (self.z1, self.z2, self.z3,
										self.m1, self.m2, self.m3)
			(z1_conj, z2_conj, z3_conj) = (np.conj(z1), np.conj(z2), np.conj(z3))
			caustic_point = z + m1 / (z1_conj - z_conj) + m2 / (z2_conj - z_conj) + m3 / (z3_conj - z_conj)
			return caustic_point

	class CriticalCurve(object):

		def __init__(self):
			self.x = []
			self.y = []
