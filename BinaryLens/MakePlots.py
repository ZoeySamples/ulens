# Zoey Samples
# Created: May 30, 2018
# BinaryLensMakePlots.py
# Last Updated: Jun 04, 2018

import numpy as np
import cmath
import matplotlib.pyplot as plt
import Functions as blf
from astropy.io import fits

class Plots(object):
	
	def __init__(self, s=None, q=None, origin='plan', pts=200, solver = 'numpy'):
		self.s = s
		self.q = q
		self.pts = pts
		self.origin = origin
		self.solver = solver
	
	def get_variables(self):
		(w_caustic, h_caustic, x_cent) = self.size_caustic()
		self.x_grid = np.linspace(x_cent - w_caustic, 
								x_cent + w_caustic, self.pts)
		self.y_grid = np.linspace(-h_caustic, h_caustic, self.pts)
		self.x_1d = np.zeros(self.pts**2)
		self.y_1d = np.zeros(self.pts**2)
		solutions = []
		for (i, xx) in enumerate(self.x_grid):
			for (j, yy) in enumerate(self.y_grid):
				idx = self.pts*i + j
				solutions.append(blf.solution(x=xx, y=yy, s=self.s, q=self.q,
										origin=self.origin, solver=self.solver))
				self.x_1d[idx] = xx
				self.y_1d[idx] = yy
		self.solutions = solutions

	def size_caustic(self):
		w = 4.*np.sqrt(self.q)*(1. + 1./(2.*(self.s**2))) / (self.s**2)
		h = 4.*np.sqrt(self.q)*(1. - 1./(2.*(self.s**2))) / (self.s**2)
		x = 0.5*self.s - 1.0/self.s
		return w, h, x

	def construct_plots(self):
		self.get_variables()
		self.num_images = np.zeros(self.pts**2, dtype=int)
		self.mag_1d = np.zeros(self.pts**2, dtype=float)
		for idx in range(self.pts**2):
			(dm, m, zeta, z1, z2) = blf.assign(x=self.x_1d[idx], y=self.y_1d[idx], 
									s=self.s, q=self.q, origin=self.origin)
			self.mag_1d[idx] = blf.magnification(x=self.x_1d[idx], y=self.y_1d[idx],
									s=self.s, q=self.q, origin=self.origin)
			for z in self.solutions[idx]:
				if blf.check_solution(dm, m, zeta, z1, z2, z, self.origin):
					self.num_images[idx] += 1

	def print_errors(self):
		num_one = num_two = num_four = 0
		for i in self.num_images:
			if (i==3) or (i==5):
				continue
			elif i == 4:
				num_four += 1
			elif i == 2:
				num_two += 1
			elif i == 1:
				num_one += 1
			else:
				print('Concern: found {:} solutions'.format(i))
		num_tot = num_one + num_two + num_four
		print('Concern: number of points where the number of images is',
				'\n1: {:}\n2: {:}\n4: {:}\ntotal: {:} out of {:} points'
				.format(num_one, num_two, num_four, num_tot, (self.pts)**2))

	def plot_n_solns(self):
		"""
		Plot showing number of images on a grid (x,y) around planetary caustic

		"""

		self.get_variables()
		self.construct_plots()
		self.print_errors()
		plt.scatter(self.x_1d, self.y_1d, c=self.num_images, s=((500./self.pts)**2), 
						marker = 'o', cmap='jet', lw=None)
		im_plot = plt.colorbar()
		im_plot.set_label('Num Images')
		plt.xlabel('X-position of source')
		plt.ylabel('Y-position of source')
		(w_caustic, h_caustic, x_cent) = self.size_caustic()
		plt.xlim(x_cent - w_caustic, x_cent + w_caustic)
		plt.ylim(-h_caustic, h_caustic)
		plt.title('Num Images using "{}" frame'.format(self.origin))

	def plot_magnification(self):
		"""
		Make square grid of points that shows the magnification at each point

		"""

		self.get_variables()
		self.construct_plots()
		plt.scatter(self.x_1d, self.y_1d, c=self.mag_1d, s=((500./self.pts)**2),
								marker = 'o', cmap='jet', lw=None)
		mag_plot = plt.colorbar()
		mag_plot.set_label('Magnification')
		plt.xlabel('X-position of source')
		plt.ylabel('Y-position of source')
		(w_caustic, h_caustic, x_cent) = self.size_caustic()
		plt.xlim(x_cent - w_caustic, x_cent + w_caustic)
		plt.ylim(-h_caustic, h_caustic)
		plt.title('Magnification using "{}" frame'.format(self.origin))

	def write_to_fits(self):
		"""
		Writes information about grid to a .fits table for comparison of magnification
		and number of images between different coordinate systems and solving methods.
		"""

		self.get_variables()
		self.construct_plots()
		col = []
		col.append(fits.Column(name='x', array=self.x_1d, format='D'))
		col.append(fits.Column(name='y', array=self.y_1d, format='D'))
		col.append(fits.Column(name='Magnification', array=self.mag_1d, format='D'))
		col.append(fits.Column(name='Number Images', array=self.num_images, format='I'))
		hdu1 = fits.BinTableHDU.from_columns(col)
		hdr = fits.Header()
		hdr['SEPARAT'] = '{:f}'.format(self.s)
		hdr['M_RATIO'] = '{:f}'.format(self.q)
		hdr['ORIGIN'] = self.origin
		hdr['SOLVER'] = self.solver
		hdr['RES'] = '{:f}'.format(self.pts)
		hdu0 = fits.PrimaryHDU(header = hdr)
		hdus = fits.HDUList([hdu0, hdu1])

		# Abbreviate file name
		if self.origin == 'geo_cent':
			origin_str = 'gcent'
		else:
			origin_str = self.origin
		if self.solver == 'numpy':
			solver_str = 'np'
		elif self.solver == 'Skowron_and_Gould_12':
			solver_str = 'SG'
		else:
			 solver_str = solver
		hdus.writeto('../Tables/BL_{}_{:2}.fits'.format(origin_str, solver_str))
		
