# Zoey Samples
# Created: Jun 7, 2018
# caustic_test.py
# Last Updated: Jul 09, 2018

from BinaryLens import BinaryLens as BL
from TripleLens import TripleLens as TL
from Caustics import Caustics as Caus
import MulensModel as mm
import matplotlib.pyplot as plt

solver = 'SG12'
originBL = 'plan'
originTL = 'body3'

plot_frame = 'caustic'
region = 'caustic2a'

s=1.4
q=1e-14

s1=s
s2=1.2
q1=q
q2=1e-3
phi=0
system='SPM'

param = ({'s': s, 'q': q, 'origin': originBL, 'region': region, 'solver': solver,
		  'plot_frame': plot_frame})
testBL = BL(**param)
causticBL = Caus(lens=testBL)

param = ({'s1': s1, 's2': s2, 'q1': q1,  'q2': q2, 'phi': phi,
		'origin': originTL, 'region': region, 'solver': solver, 'system': system,
		'plot_frame': plot_frame})
testTL = TL(**param)
causticTL = Caus(lens=testTL, solver=solver)

# Plot the caustic alone
causticBL.plot_caustic(s=1, color='blue')
#causticTL.plot_caustic(s=1, color='red')
plt.show()

#Interesting parameters:
solver = 'SG12'
originBL = 'star'
originTL = 'body3'

plot_frame='caustic'

s=0.8
q=1e-5

s1=0.8
s2=0.8
q1=1e-5
q2=1e-1
phi=135
system='SPM'





