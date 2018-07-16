# Zoey Samples
# Created: Jun 26, 2018
# ToleranceTest.py
# Last Updated: Jun 26, 2018

import sys
import os
import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from HDUL import HDUL

"""
The purpose of this script is to demonstrate whether it is worth using a
new algorithm to check solutions of the polynomial. The original method
involves using a tolerance to determine whether each solution satisfies
the polynomial within a certain range when substitued back in. The new
method instead uses an algorithm that does not rely on a tolerance.
According to these results, the new method is slightly more reliable.
Moreover, the new method allows the user to run their scripts without
specifying a tolerance. This is beneficial because there is no objective
way to determine the best tolerance, so I have been relying on trial and
error; that favors methods that are not replicable nor substantiable.
"""

file_name = (['../Tables/tolerance_on_0.fits', '../Tables/tolerance_off_0.fits', 
			 '../Tables/tolerance_on_1.fits', '../Tables/tolerance_off_1.fits'])
hdul = []
for i in range(len(file_name)):
	param = {'file_name': file_name[i]}
	hdul.append(HDUL(**param))

file_str = (['With Tolerance On; Mass Ratio 1e-15',
			'With Tolerance Off; Mass Ratio: 1e-15',
			'With Tolerance On; Mass Ratio: 5e-16',
			'With Tolerance Off; Mass Ratio: 5e-16'])

for i in range(len(hdul)):
	print('\n',file_str[i])
	hdul[i].plot_num_images(print_errors=True)
