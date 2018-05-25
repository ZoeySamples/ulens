# Zoey Samples
# Created: May 22, 2018
# BinaryLensFunctions.py
# Last Updated: May 24, 2018; 6:41PM

import numpy as np
import cmath

def mass_ratio(q):
	"""Define m and dm for use in polynomial"""
	dm = (1. - q)/2.
	m = (q + 1.)/2.
	return dm, m

def source_position(x, y):
	"""Assign the position of the source in complex space"""
	return x + y*1j

def assign(x, y, s, q):
	""""Assign all variables to be used in polynomial"""
	dm, m = mass_ratio(q)   	 # Difference & sum of masses
	zeta = source_position(x,y)  # Source's position
	z1 = 0.5*s		    		 # Position of lensing body 1
	return dm, m, zeta, z1

def solution(*test):
	"""Return solutions of polynomial"""
	p = list(range(6)) # Polynomial to solve
	dm, m, zeta, z1 = assign (*test)
	p[5] = (z1**2)*(4*(dm**2)*zeta + 4*m*dm*z1 + 4*dm*zeta*zeta.conjugate()*z1 + 2*m*zeta.conjugate()*(z1**2) + zeta*(zeta.conjugate()**2)*(z1**2) - 2*dm*(z1**3) - zeta*(z1**4))
	p[4] = -8*m*dm*zeta*z1 - 4*(dm**2)*(z1**2) - 4*(m**2)*(z1**2) - 4*m*zeta*zeta.conjugate()*(z1**2) - 4*dm*zeta.conjugate()*(z1**3) - (zeta.conjugate()**2)*(z1**4) + (z1**6)
	p[3] = 4*(m**2)*zeta + 4*m*dm*z1 - 4*dm*zeta*zeta.conjugate()*z1 - 2*zeta*(zeta.conjugate()**2)*(z1**2) + 4*dm*(z1**3) + 2*zeta*(z1**4)
	p[2] = 4*m*zeta*zeta.conjugate() + 4*dm*zeta.conjugate()*z1 + 2*(zeta.conjugate()**2)*(z1**2) - 2*(z1**4)
	p[1] = -2*m*zeta.conjugate() + zeta*(zeta.conjugate()**2) - 2*dm*z1 - zeta*(z1**2)
	p[0] = z1**2 - zeta.conjugate()**2
	return np.roots(p)

def check_solution(dm, m, zeta, z1, z):
	zeta_actual = z + (m-dm)/(z1.conjugate()-z.conjugate()) + (m+dm)/(-z1.conjugate()-z.conjugate())
	if np.abs(zeta - zeta_actual) > 0.01:
		return False
	else:
		return True

def magnification(x, y, s, q):
	"""Returns the magnification for each configuration"""
	dm, m, zeta, z1 = assign (x, y, s, q)
	solutions = solution(x, y, s, q)
	magn = list(range(5))
	i=0
	for z in solutions:
		detJ = 1. - ((m-dm)/((z-z1)**2) + (m+dm)/((z+z1)**2)) * ((m-dm)/((z.conjugate()-z1)**2) + (m+dm)/((z.conjugate()+z1)**2))
		if check_solution(dm, m, zeta, z1, z) == True:
			magn[i] = np.abs(1./detJ)
		else:
			magn[i] = 0
		i+=1
	return sum(magn)
