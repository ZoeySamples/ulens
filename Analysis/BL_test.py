# Zoey Samples
# Created: Jun 7, 2018
# BL_Test.py
# Last Updated: Jun 21, 2018

from BinaryLens import BinaryLens as BL

solver = 'SG12'
origin = 'caustic'

x = [0.0, 1.3219, 1.0799, 1.2489]
y = [0.0, -0.0771, 0.0985, 0.0209]
s = [1.0, 1.35, 1.1, 0.9357]
q = [1.0, 0.00578, 0.99, 0.99]


param = []
test = []
for i in range(len(x)):
	param.append({'x': x[i], 'y': y[i], 's': s[i], 'q': q[i], 'origin': origin,
					 'solver': solver}) 
	test.append(BL(**param[-1]))

for t in test:
	t.print_image_position(print_input=True)
	t.print_magnification(print_input=False)
	print('_'*40)
