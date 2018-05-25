# Zoey Samples
# Created: May 22, 2018
# BinaryLensAnalysis.py
# Last Updated: May 24, 2018; 6:59PM

import numpy as np
import cmath
import matplotlib.pyplot as plt
from matplotlib import cm
import MulensModel as mm
import BinaryLensFunctions as blf
from scipy.interpolate import griddata

def size_caustic(s, q):
	w = 4*np.sqrt(q)*(1 + 1/(2*(s**2))) / (s**2)
	h = 4*np.sqrt(q)*(1 - 1/(2*(s**2))) / (s**2)
	x = 0.5*s - 1.0/s
	return w, h, x

tests2 = [
	[1.35, 1e-7]
		]	# Input parameters for trials where s > 1

plot_on = True
for test in tests2:
	"""Plot showing number of images for each test, on a grid (x,y) around planetary caustic"""
	if not plot_on:
		break	
	pts = 150	# Number of data points on each side of the grid
	w_caustic, h_caustic, x_center = size_caustic(test[0], test[1])
	x_factor = 15.
	y_factor = 15.
	x_grid = np.linspace(x_center - x_factor*w_caustic, x_center + x_factor*w_caustic, pts)
	y_grid = np.linspace(-y_factor*h_caustic, y_factor*h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	im_num = np.zeros(pts**2, dtype=int)
	color = np.zeros(pts**2)
	print('parameters:\ns={:}\nq={:}'.format(*test))
	for i, xx in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			x_1d[idx] = xx
			y_1d[idx] = yy
			dm, m, zeta, z1 = blf.assign(xx, yy, test[0], test[1])
			solutions = blf.solution(xx, yy, test[0], test[1])
			for z in solutions:
				if blf.check_solution(dm, m, zeta, z1, z) == True:
					im_num[idx] += 1
			if im_num[idx] == 5:
				color[idx] = 255
			elif im_num[idx] == 3:
				color[idx] = 1
			else:
				color[idx] = 120
				print('Concern: number of images=', im_num[idx])
				print('x_source={:}\ny_source={:}'.format(xx, yy))
	plt.scatter(x_1d, y_1d, c=im_num, s=8, cmap='jet')
	im_plot = plt.colorbar()
	im_plot.set_label('Num Images')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.title('Number of Images vs. Position')
	plt.show()

# Scatter plot using 1D arrays only
plot_on = True
for test in tests2:
	"""Make square grid of points that shows the magnification at each point; assume s>1"""
	if not plot_on:
		break
	pts = 150	# Number of data points on each side of the grid
	w_caustic, h_caustic, x_center = size_caustic(test[0], test[1])
	x_factor = 1.
	y_factor = 1.
	x_grid = np.linspace(x_center - x_factor*w_caustic, x_center + x_factor*w_caustic, pts)
	y_grid = np.linspace(-y_factor*h_caustic, y_factor*h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	mag_1d = np.zeros(pts**2)
	for i, xx in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			mag_1d[idx] = blf.magnification(xx, yy, test[0], test[1])
			x_1d[idx] = xx
			y_1d[idx] = yy
	plt.scatter(x_1d, y_1d, c=mag_1d, s=8, cmap='jet')
	mag_plot = plt.colorbar()
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.title('Magnification of Image vs. Position')
	plt.show()

# Interpolated scatter plot; Not working & ugly
plot_on = False
for test in tests2:
	"""Make square grid of points that shows the magnification at each point; assume s>1"""
	if plot_on == False:
		break
	pts = 100	# Number of data points on each side of the grid
	w_caustic, h_caustic, x_center = size_caustic(test[0], test[1])
	x_factor = 1.
	y_factor = 1.
	x_grid = np.linspace(x_center - x_factor*w_caustic, x_center + x_factor*w_caustic, pts)
	y_grid = np.linspace(-y_factor*h_caustic, y_factor*h_caustic, pts)
	x_1d = np.zeros(pts**2)
	y_1d = np.zeros(pts**2)
	mag_1d = np.zeros(pts**2)
	mag_grid = np.zeros((pts, pts))
	for i, xx in enumerate(x_grid):
		for j, yy in enumerate(y_grid):
			idx = pts*i + j
			mag_grid[i,j] = blf.magnification(xx, yy, test[0], test[1])
			mag_1d[idx] = blf.magnification(xx, yy, test[0], test[1])
			x_1d[idx] = xx
			y_1d[idx] = yy
#	grid_x, grid_y = np.mgrid[0:1:pts, 0:1:pts]
#	grid_mag = griddata((x_1d, y_1d), mag_color, (grid_x, grid_y), method='cubic')
#	plt.imshow((grid_mag), cmap='viridis')
#	plt.scatter(x_1d, y_1d, c=mag_color, cmap='viridis', zorder=mag_color)
#	plt.show()
	xi = yi = np.arange(0,1.01,0.01)
	xi, yi = np.meshgrid(xi,yi)
	magi = griddata((x_1d, y_1d), mag_color, (xi, yi), method='linear')
#	plt.contour(xi, yi, magi, np.arange(0,1.01,0.01))
	im = plt.imshow(mag_grid, cmap='hot')
	plt.colorbar(im, orientation='horizontal')
	plt.scatter(x_1d, y_1d, c=mag_color, cmap='viridis', zorder=mag_color)
	plt.show()
#	points = np.vstack((x_1d, y_1d)).T
#	points = np.asarray(points)
#	values = mag_color
#	DEM = interpolate.griddata(points, values, (XI,YI), method='linear')
#	levels = np.arange(np.min(DEM),np.max(DEM),25)
#	plt.contour(DEM, levels,linewidths=0.2,colors='k')
#	plt.imshow(DEM,cmap ='RdYlGn_r',origin='lower')
#	plt.colorbar()
