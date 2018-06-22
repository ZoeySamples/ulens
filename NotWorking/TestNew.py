# Zoey Samples
# Created: Jun 7, 2018
# BL_Test.py
# Last Updated: Jun 21, 2018

from BL import BinaryLens as BL


"""
This file demonstrates that the new method for checking roots is not working
as expected.
"""

"""
# All 4 initial test points

x = [0.0, 1.3219, 1.0799, 1.2489]
y = [0.0, -0.0771, 0.0985, 0.0209]
s = [1.0, 1.35, 1.1, 0.9357]
q = [1.0, 0.00578, 0.99, 0.99]

"""

# Test point that isn't producing the correct results
x = [1.3219]
y = [-0.0771]
s = [1.35]
q = [0.00578]

solver = 'SG12'
origin = 'plan'

param = []
test = []
for i in range(len(x)):
	param.append({'x': x[i], 'y': y[i], 's': s[i], 'q': q[i], 'origin': origin,
					 'solver': solver}) 
	test.append(BL(**param[-1]))


region = 'custom'
region_lim = (-5, 15, -3, 2)

for t in test:
	t.print_image_position(print_input=False)
#	t.print_magnification(print_input=False)
	print('_'*40, '\n')


"""
The final "correct" image location, which is currently being rejected, is actually
supposed to come from root 4. However, root 4 has no chance of being accepted
because the "solution" is very far away from the solved-for root value. The
issue is that the solved-for root value is very close to the value of z1, the
position of the first lensing body. This makes the denominator of that term
very small, making the expression for the binary lens solution very large.
See the expression below.

# As seen in the source code, line 275 & after:

roots = self.get_roots(x=x, y=y)
lensing_body1 = (self.m + self.dm) / np.conjugate(roots - self.z1) *
lensing_body2 = (self.m - self.dm) / np.conjugate(roots - self.z2)
solutions = self.zeta + lensing_body1 + lensing_body2 **

* (root4 - self.z1) ~ 0  --->  lensing_body1(root4) ~ very large

** Because lensing_body1(root4) is very large, solutions(root4) is very large.
   However, the solved-for root value is approximately equal to z1, so this
   should not be the case.
"""






