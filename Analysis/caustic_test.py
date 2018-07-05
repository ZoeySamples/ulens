# Zoey Samples
# Created: Jun 7, 2018
# caustic_test.py
# Last Updated: Jun 25, 2018

from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL
from Caustics import Caustics as Caus
import MulensModel as mm
import matplotlib.pyplot as plt

solver = 'SG12'
originBL = 'star'
originTL = 'body3'

s=0.5
q=1e-3

s1=0.5
s2=0.5
q1=1e-3
q2=1e-7
phi=90

param = ({'s': s, 'q': q, 'origin': originBL, 'solver': solver}) 
test_BLlens = BL(**param)

param = ({'s1': s1, 's2': s2, 'q1': q1,  'q2': q2, 'phi': phi,
		'origin': originTL, 'solver': solver, 'system': 'SPP'}) 
test_TLlens = TL(**param)

causparam = ({'lens': test_BLlens})
test_BLcaustic = Caus(**causparam)

causparam = ({'lens': test_TLlens})
test_TLcaustic = Caus(**causparam)

#test_BLcaustic.plot_caustic(s=1)
#plt.show()

test_TLcaustic.plot_caustic(s=1)
plt.scatter(test_TLcaustic.z1.real, test_TLcaustic.z1.imag, color='red', s=50)
plt.scatter(test_TLcaustic.z2.real, test_TLcaustic.z2.imag, color='red', s=10)
plt.scatter(test_TLcaustic.z3.real, test_TLcaustic.z3.imag, color='red', s=3)
plt.show()
