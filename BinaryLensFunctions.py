# Zoey Samples
# Created: May 22, 2018
# BinaryLensFunctions.py
# Last Updated: May 31, 2018

import sys
import os
import ctypes
import numpy as np
import cmath

MODULE_PATH = os.path.abspath(__file__)
MODULE_PATH = os.path.dirname(MODULE_PATH)
MODULE_PATH = os.path.dirname(MODULE_PATH)
PATH = os.path.join(MODULE_PATH, 'MulensModel-master', 'source', 'VBBL',
		"VBBinaryLensingLibrary_wrapper.so")

# Access SG12 code
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

# Access zroots code
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

def print_frame(origin):
	if origin == 'geo_cent':
		return origin, 'geometric center frame'
	elif origin == 'star':
		return origin, 'star frame'
	elif origin == 'plan':
		return origin, 'planet frame'
	elif origin == 'com':
		return origin, 'center-of-mass frame'
	else:
		origin = 'geo_cent'
		print('Error: unknown coordinate system given by string variable "origin"')
		print('Calculations will be carried out in the geometric center frame\n')
		return origin, 'geometric center frame'

def mass_ratio(q):
	"""Define m and dm for use in polynomial"""
	dm = (1. - q)/2.
	m = (q + 1.)/2.
	return dm, m

def source_position(x, y, s, m, dm, origin='geo_cent'):
	"""Assign the position of the source in complex space"""
	if origin == 'geo_cent':
		return x + y*1.j
	elif origin == 'star':
		return (x + s/2.) + y*1.j
	elif origin == 'plan':
		return (x - s/2.) + y*1.j
	elif origin == 'com':
		return (x + s*dm/(2.*m)) + y*1.j
	else:
		print('Error: unknown coordinate system given by string variable "origin"\n')
		print('Calculations will be carried out in the geometric center frame)')
		return x + y*1.j

def assign(x, y, s, q, origin):
	"""
	Assign all variables to be used in polynomial. Returns:

		dm	 - defined as half the difference of the masses (positive)
		m	 - defined as half the sum of the masses
		zeta - the source's position
		z1	 - the planet's position
		z2	 - the star's position

	"""

	(dm, m) = mass_ratio(q)   	 # Difference & sum of masses
	zeta = source_position(x, y, s, m, dm, origin)  # Source's position
	if origin == 'geo_cent':
		z1 = 0.5*s
		z2 = -0.5*s
	elif origin == 'star':
		z1 = s
		z2 = 0.
	elif origin == 'plan':
		z1 = 0.
		z2 = -s
	elif origin == 'com':
		z1 = s*(m + dm) / (2.*m)
		z2 = -s*(m - dm) / (2.*m)
	return dm, m, zeta, z1, z2

def solution(x, y, s, q, origin, solver='numpy'):
	"""
	Return solutions of polynomial.
	"""

	p = np.zeros(6, dtype=complex) # Polynomial to solve
	(dm, m, zeta, z1, z2) = assign(x, y, s, q, origin)

	"""
	The following is the general solution to the binary lens equation.
	"""
	
	p[0] = (-zeta.conjugate() + z1)*(zeta.conjugate()- z2)

	p[1] = (m*z1 + m*z2 + 2.*(z1**2)*z2 + 2.*z1*(z2**2) + dm*(-z1 + z2) + 
	z1*z2*zeta + (zeta.conjugate()**2)*(2.*z1 + 2.*z2 + zeta) -
	zeta.conjugate()*(2.*m + (z1 + z2)*(2.*z1 + 2.*z2 + zeta)))

	p[2] = (dm*(z1**2) - m*(z1**2) - 2.*m*z1*z2 - (z1**3)*z2 - dm*(z2**2) - 
	m*(z2**2) - 4.*(z1**2)*(z2**2) - z1*(z2**3) - 2.*m*z1*zeta - 2.*m*z2*zeta - 
	2.*(z1**2)*z2*zeta - 2.*z1*(z2**2)*zeta - (zeta.conjugate()**2)*((z1**2) + 
	2.*z1*(2.*z2 + zeta) + z2*(z2 + 2.*zeta)) + zeta.conjugate()*(2.*dm*(z1 - z2) + 
	2.*m*(z1 + z2 + 2.*zeta) + (z1 + z2)*((z1**2) + 4.*z1*z2 + (z2**2) + 
	2.*z1*zeta + 2.*z2*zeta)))

	p[3] = (-2.*(m**2)*(z1 + z2 - 2.*zeta) - 3.*m*(2.*zeta.conjugate() - z1 - z2)*
	(z1 + z2)*zeta + dm*(z1 - z2)*(2.*m - 2.*z1*z2 + z1*zeta + z2*zeta - 
	2.*zeta.conjugate()*(z1 + z2 + zeta)) + (zeta.conjugate() - z1)*(zeta.conjugate() -
	 z2)*((z2**2)*zeta + (z1**2)*(2.*z2 + zeta) + 2.*z1*z2*(z2 + 2.*zeta)))

	p[4] = ((-dm**2)*((z1 - z2)**2) + (m**2)*((z1**2) + 6.*z1*z2 + (z2**2) - 
	4.*z1*zeta - 4.*z2*zeta) - m*(2.*zeta.conjugate() - z1 - z2)*(z1*z2*(z2 - 
	4.*zeta) + (z1**2)*(z2 - zeta) - (z2**2)*zeta) - (zeta.conjugate() - z1)*z1*
	(zeta.conjugate() - z2)*z2*(2.*z2*zeta + z1*(z2 + 2.*zeta)) + dm*(z1 - z2)*
	(z1*z2*(z2 - 2.*zeta) + (z1**2)*(z2 - zeta) - (4.*m + (z2**2))*zeta + 
	2.*zeta.conjugate()*(z2*zeta + z1*(z2 + zeta))))

	p[5] = (-2.*(m**2)*(z1**2)*z2 - 2.*(m**2)*z1*(z2**2) - m*(z1**3)*(z2**2) - 
	m*(z1**2)*(z2**3) + (m**2)*(z1**2)*zeta + (dm**2)*((z1 - z2)**2)*zeta + 
	2.*(m**2)*z1*z2*zeta + m*(z1**3)*z2*zeta + (m**2)*(z2**2)*zeta + 
	(zeta.conjugate()**2)*(z1**2)*(z2**2)*zeta + 2.*m*(z1**2)*(z2**2)*zeta + 
	m*z1*(z2**3)*zeta + (z1**3)*(z2**3)*zeta - dm*(z1 - z2)*(2.*m + z1*z2)*
	(z1*(z2 - zeta) - z2*zeta) - zeta.conjugate()*z1*z2*((2.*dm*(z1 - z2) + 
	z1*z2*(z1 + z2))*zeta + m*(-2.*z1*z2 + 2.*z1*zeta + 2.*z2*zeta)))
	
	"""
	The following is the simplified solution to the binary lens equation in the
	geometric center frame.
	"""
	"""
	p[5] = (z1**2)*(4*(dm**2)*zeta + 4*m*dm*z1 + 4*dm*zeta*zeta.conjugate()*z1 + 2*m*zeta.conjugate()*(z1**2) + zeta*(zeta.conjugate()**2)*(z1**2) - 2*dm*(z1**3) - zeta*(z1**4))
	p[4] = -8*m*dm*zeta*z1 - 4*(dm**2)*(z1**2) - 4*(m**2)*(z1**2) - 4*m*zeta*zeta.conjugate()*(z1**2) - 4*dm*zeta.conjugate()*(z1**3) - (zeta.conjugate()**2)*(z1**4) + (z1**6)
	p[3] = 4*(m**2)*zeta + 4*m*dm*z1 - 4*dm*zeta*zeta.conjugate()*z1 - 2*zeta*(zeta.conjugate()**2)*(z1**2) + 4*dm*(z1**3) + 2*zeta*(z1**4)
	p[2] = 4*m*zeta*zeta.conjugate() + 4*dm*zeta.conjugate()*z1 + 2*(zeta.conjugate()**2)*(z1**2) - 2*(z1**4)
	p[1] = -2*m*zeta.conjugate() + zeta*(zeta.conjugate()**2) - 2*dm*z1 - zeta*(z1**2)
	p[0] = z1**2 - zeta.conjugate()**2
	"""

	if solver == 'Skowron_and_Gould_12':
		pp = p[::-1]
		out = _vbbl_SG12_5(*(pp.real.tolist() + pp.imag.tolist()))
		roots = [out[i] + out[i+5] * 1.j for i in range(5)]
		return roots
	elif solver == 'zroots':
		pp = p[::-1]
		out = _zroots_5(*(pp.real.tolist() + pp.imag.tolist()))
		roots = [out[i] + out[i+5] * 1.j for i in range(5)]
		return roots
	elif solver == 'numpy':
		r = np.roots(p)
		return r.tolist()

def check_solution(dm, m, zeta, z1, z2, z, origin):
	zeta_actual = (z + (m-dm)/(z1.conjugate()-z.conjugate()) +
		(m+dm)/(z2.conjugate()-z.conjugate()))
	if np.abs(zeta - zeta_actual) > 0.0001:	# Tolerance for solution check
		return False
	else:
		return True

def magnification(x, y, s, q, origin, solutions=None, solver='numpy'):
	"""Returns the magnification for each configuration"""
	if solutions is None:
		solutions = solution(x, y, s, q, origin, solver)
	(dm, m, zeta, z1, z2) = assign (x, y, s, q, origin)
	magn = list(range(5))
	for (i, z) in enumerate(solutions):
		detJ = (1. - ((m-dm)/((z-z1)**2) + (m+dm)/((z-z2)**2)) *
			((m-dm)/((z.conjugate()-z1)**2) + (m+dm)/((z.conjugate()-z2)**2)))
		if check_solution(dm, m, zeta, z1, z2, z, origin):
			magn[i] = np.abs(1./detJ)
		else:
			magn[i] = 0.
	return sum(magn)
