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
from astropy.io import fits
from itertools import repeat

# Definitions
def get_variables(hdul):
	"""
	Assigns data to arrays.
	"""

	x = hdul[1].data['x']
	y = hdul[1].data['y']
	num_images = hdul[1].data['Number Images']
	magn = hdul[1].data['Magnification']
	res = float(hdul[0].header['RES'])
	origin = hdul[0].header['ORIGIN']
	solver = hdul[0].header['SOLVER']
	s = float(hdul[0].header['SEPARAT'])
	q = float(hdul[0].header['M_RATIO'])
	return (x, y, s, q, num_images, magn, res, origin, solver)

def check_compatible(x1, y1, res1, x2, y2, res2):
	"""
	Determines if tables have the same size and are plotted over the same grid.
	"""

	if int(res1) != int(res2):
		sys.exit('Error: Resolutions are not the same. Exiting script')
	z = [0j]*(len(x1))
	z = [z, z]
	ok = 0
	z1 = x1 + y1*1.j
	z2 = x2 + y2*1.j
	for i in range(len(x1)):
		if np.abs(z1[i] - z2[i]) > 1e-10:
			print('Error: grids are not formed over the same region')
			print('At index {}, realtive position between points is {}'.format(i,
									np.abs(z1[i]-z2[i])))
		else:
			ok += 1
	print('Test passed for {} out of {} points'.format(ok, int(res1**2)))

def plot_magnification(x, y, magn, res, origin, solver, xmin=None, xmax=None,
									ymin=None, ymax=None):
	"""
	Makes a plot of the magnification vs. position
	"""

	if xmin==None:
		xmin = x[0]
	if xmax==None:
		xmax = x[-1]
	if ymin==None:
		ymin = y[0]
	if ymax==None:
		ymax = y[-1]

	plt.scatter(x, y, c=magn, s=((800./res)**2), marker = 'o', cmap='jet', lw=None)
	mag_plot = plt.colorbar()
	mag_plot.set_label('Magnification')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	(origin_str, solver_str) = abbreviate(origin=origin, solver=solver)
	plt.title('Magnification: {} Frame; {} Solver'.format(origin_str, solver_str))

def plot_n_solutions(x, y, num_images, res, origin, solver):
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
	(origin_str, solver_str) = abbreviate(origin=origin, solver=solver)
	plt.title('Num Images: {} Frame; {} Solver'.format(origin_str, solver_str))

def plot_rel_magn(x1, y1, magn1, x2, y2, magn2, res):
	"""
	Plots the fractional difference in magnification between two sets of data
	"""
	rel_magn = (magn1 / magn2)
	plt.scatter(x1, y1, c=rel_magn, s=((800./res)**2), marker = 'o',
					cmap='jet', lw=None)
	rel_plot = plt.colorbar()
	rel_plot.set_label('Fractional Difference')
	plt.xlabel('X-position of source')
	plt.ylabel('Y-position of source')
	plt.xlim(x1[0], x1[-1])
	plt.ylim(y1[0], y1[-1])
	plt.title('Relative Magnification')

def assess_outliers(x, y, num_images, magn, res, origin, cutoff):
	"""
	Creates new list of (x, y, magn) only for magn value that are above cutoff
	"""

	num_outliers = int(0)
	x_new = []
	y_new = []
	magn_new = []
	for (i, magn) in enumerate(magn):
		if magn > cutoff:
			x_new.append(x[i])
			y_new.append(y[i])
			magn_new.append(magn)
	return x_new, y_new, magn_new

def write_to_fits(x, y, magn, origin, solver, s, q, res):
	"""
	Writes new fits file only with extreme data
	"""

	col = []
	col.append(fits.Column(name='x', array=x, format='D'))
	col.append(fits.Column(name='y', array=y, format='D'))
	col.append(fits.Column(name='Magnification', array=magn, format='D'))
	hdu1 = fits.BinTableHDU.from_columns(col)
	hdr = fits.Header()
	hdr['SEPARAT'] = s
	hdr['M_RATIO'] = q
	hdr['ORIGIN'] = origin
	hdr['SOLVER'] = solver
	hdr['RES'] = '{:d}'.format(res)
	hdu0 = fits.PrimaryHDU(header = hdr)
	hdus = fits.HDUList([hdu0, hdu1])
	(origin_str, solver_str) = abbreviate(origin=origin, solver=solver)
	hdus.writeto('../Tables/BL_extreme_{}_{:2}.fits'.format(origin_str, solver_str))

def abbreviate(origin, solver):
	"""
	Abbreviate file name
	"""

	if origin == 'geo_cent':
		origin_str = 'gcent'
	else:
		origin_str = origin

	if solver == 'numpy':
		solver_str = 'np'
	elif solver == 'Skowron_and_Gould_12':
		solver_str = 'SG'
	else:
		 solver_str = solver
	return origin_str, solver_str


"""
Option 1: Input 2 tables to analyze & a fits file to save it to
"""
"""
# Assign file names to variables
if len(sys.argv) != (3 or 4):
    raise ValueError('Two arguments needed; Optional: One save file name')
table1 = '../Tables/' + sys.argv[1]
table2 = '../Tables/' + sys.argv[2]
if len(sys.argv) == 4:
	file_out = sys.argv[3]

# Open files to list variable
hdul = []
hdul.append(fits.open(table1))
hdul.append(fits.open(table2))
"""

"""
Option 2: Input any number of fits files; no output file
"""
table = []
hdul = []
for i in range(len(sys.argv) - 1):
	table.append('../Tables/' + sys.argv[i+1])
	hdul.append(fits.open(table[i]))

# Initialize
x = [[]] * len(hdul)
y = [[]] * len(hdul)
s = [[]] * len(hdul)
q = [[]] * len(hdul)
num_images = [[]] * len(hdul)
magn = [[]] * len(hdul)
res = [[]] * len(hdul)
origin = [[]] * len(hdul)
solver = [[]] * len(hdul)

# Main
cutoff = 1000
for (i, hdul) in enumerate(hdul):
	(x[i], y[i], s[i], q[i], num_images[i], magn[i], res[i], origin[i], 
				solver[i]) = get_variables(hdul = hdul)
	(x_new, y_new, magn_new) = assess_outliers(x=x[i], y=y[i], num_images=num_images[i],
					magn=magn[i], res=res[i], origin=origin[i], cutoff=cutoff)
	caustics = mm.Caustics(s=(s[i]), q=(q[i]))
	caustics.plot(s=5)
	xmin=x[0][0]
	xmax=x[0][-1]
	ymin=y[0][0]
	ymax=y[0][-1]
	plot_magnification(x=x_new, y=y_new, magn=magn_new, res=0.2*res[i],
				origin=origin[i], solver=solver[i], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
	plt.show()
	#write_to_fits(x=x_new, y=y_new, magn=magn_new, origin=origin[i],
	#			solver=solver[i], s=s[i], q=q[i], res=res[i])
	


	"""
	check_compatible(x[0], y[0], res[0], x[1], y[1], res[1])
	caustics = mm.Caustics(s=s[0], q=q[0])
	caustics.plot()
	plot_rel_magn(x[0], y[0], magn[0], x[1], y[1], magn[1], res[0])
	plt.show()
	#plt.savefig(file_out)
	"""


################################

#file_name_1 = sys.argv[1]
#file_name_2 = sys.argv[2]
#file_out = sys.argv[3]

#data_1 = BinaryLensSolutions(file_name=file_name_1)
#data_2 = BinaryLensSolutions(file_name=file_name_2)

#plot_magnification_ratio(data_1, data_2)
#plt.savefig(file_out)
