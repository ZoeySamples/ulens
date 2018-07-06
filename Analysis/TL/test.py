# Zoey Samples
# Created: Jun 7, 2018
# TL_Test.py
# Last Updated: Jun 25, 2018

from TripleLens import TripleLens as TL
import MulensModel as mm
import matplotlib.pyplot as plt

solver = 'SG12'
origin = 'body3'
system = 'SPM'
SFD = True

x = [0.0, 1.3219, 1.0799, 1.2489]
y = [0.0, -0.0771, 0.0985, 0.0209]
s1 = [1.0, 1.35, 1.1, 0.9357]
s2 = [1.5, 1.5, 1.5, 1.5]
q1 = [1.0, 0.00578, 0.99, 0.99]
q2 = [1e-1, 1e-1, 1e-1, 1e-1]
phi = 0

param = []
test = []
for i in range(len(x)):
	param.append({'x': x[i], 'y': y[i], 's1': s1[i], 's2': s2[i],
				'phi': phi, 'q1': q1[i], 'q2': q2[i], 'origin': origin,
				'solver': solver, 'system': system, 'SFD': SFD}) 
	test.append(TL(**param[-1]))

for t in test:
#	t.print_image_position(print_input=False)
	t.print_magnification(print_input=False)
	print('_'*40)

