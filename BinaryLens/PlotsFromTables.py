# Zoey Samples
# Created: Jun 04, 2018
# BinaryLensMakePlots.py
# Last Updated: Jun 04, 2018

import numpy as np
import cmath
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
import Functions as blf
import MulensModel as mm

# Assign file names to variables
if len(sys.argv) != 3:
    raise ValueError('Exactly two arguments needed')
table1 = '../Tables/' + sys.argv[1]
table2 = '../Tables/' + sys.argv[2]

# Open files to list variable
hdul = []
hdul.append(fits.open(table1))
hdul.append(fits.open(table2))

# Initialize
x=[]
y=[]
s=[]
q=[]
num_images = []
magn = []
res = []
origin = []

def get_variables(hdul):
	"""
	Assigns data to arrays.
	"""

	for i in range(len(hdul)):
		x.append(hdul[i][1].data['x'])
		y.append(hdul[i][1].data['y'])
		num_images.append(hdul[i][1].data['Number Images'])
		magn.append(hdul[i][1].data['Magnification'])
		res.append(float(hdul[i][0].header['RES']))
		origin.append(hdul[i][0].header['ORIGIN'])
		s.append(float(hdul[i][0].header['SEPARAT']))
		q.append(float(hdul[i][0].header['M_RATIO']))
	return x, y, num_images, magn, res, origin

def check_compatible(x, y, num_images, magn, res, origin):
	"""
	Determines if tables have the same size and are plotted over the same grid.
	"""

	if int(res[0]) != int(res[1]):
		sys.exit('Error: Resolutions are not the same. Exiting script')
	z = [0j]*(len(x[0]))
	z = [z, z]
	ok = 0
	for i in range(len(x[0])):
		for j in range(2):
			z[j][i] = x[j][i] + y[j][i]*1.j
		if np.abs(z[0][i] - z[1][i]) > 1e-10:
			print('Error: grids are not formed over the same region')
			print('At index {}, realtive position between points is {}'.format(i,
									np.abs(z[0][i]-z[1][i])))
		else:
			ok += 1
	print('Test passed for {} out of {} points'.format(ok, int(res[0]**2)))

def plot_magnification(x, y, num_images, magn, res, origin):
	"""
	Makes a plot of the magnification vs. position
	"""

	plt.scatter(x, y, c=magn, s=((800./res)**2), marker = 'o', cmap='jet', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x[0], x[-1])
	plt.ylim(y[0], y[-1])
	plt.title('Magnification using "{}" frame'.format(origin))

def plot_n_solutions(x, y, num_images, magn, res, origin):
	"""
	Makes a plot of the number of solutions vs. position
	"""

	plt.scatter(x, y, c=num_images, s=((800./res)**2), marker = 'o',
					cmap='jet', lw=None)
	soln_plot = plt.colorbar()
	soln_plot.set_label('Num Images')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x[0], x[-1])
	plt.ylim(y[0], y[-1])
	plt.title('Magnification using "{}" frame'.format(origin))

def plot_rel_magn(x, y, num_images, magn, res, origin):
	"""
	Plots the fractional difference in magnification between two sets of data
	"""
	print(repr(magn[0][0]), repr(magn[1][0]))
	rel_magn = (magn[0] / magn[1])
	print(repr(rel_magn[0]))
	plt.scatter(x[0], y[0], c=rel_magn, s=((800./res[0])**2), marker = 'o',
					cmap='jet', lw=None)
	rel_plot = plt.colorbar()
	rel_plot.set_label('Fractional Difference')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x[0][0], x[0][-1])
	plt.ylim(y[0][0], y[0][-1])
	plt.title('Relative Magnification')

plot_list = get_variables(hdul)
check_compatible(*plot_list)
caustics = mm.Caustics(s=s[0], q=q[0])
caustics.plot()
plot_rel_magn(*plot_list)
plt.show()


################################

#file_name_1 = sys.argv[1]
#file_name_2 = sys.argv[2]
#file_out = sys.argv[3]

#data_1 = BinaryLensSolutions(file_name=file_name_1)
#data_2 = BinaryLensSolutions(file_name=file_name_2)

#plot_magnification_ratio(data_1, data_2)
#plt.savefig(file_out)
