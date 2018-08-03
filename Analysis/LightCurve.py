#LightCurve.py

"""
This file shows how to define and plot single or binary lens and point
source models. The binary lens model has a planetary mass ratio and
one can clearly see that the planetary anomaly is a short perturbation
on top of smooth single point lens light curve. Also the lens source
geometry (i.e., the source trajectory and the caustic position are
plotted).

"""
import matplotlib.pyplot as pl
import numpy as np

import MulensModel
from MulensModel.model import Model


# Create a PSPL model
t_0 = 3583.
u_0 = 0.38
t_E = 12

pspl = Model({'t_0': t_0, 'u_0': u_0, 't_E': t_E})

# Create a planet model with same PSPL parameters
s = 1.5
q = 1e-2
alpha = 332

planet = Model(
	{'t_0': t_0, 'u_0': u_0, 't_E': t_E, 's': s, 'q': q, 'alpha': alpha})

# Plot planet and PSPL models and show difference in magnification at
# planet perturbation
pl.figure()
pspl.plot_magnification(
    color='blue', linestyle=':', zorder=1, label='Point Lens', lw=2)
planet.plot_magnification(
    color='red', linestyle='-', zorder=2, label='Planet', lw=1)
pl.ylabel('Magnification', fontsize=18)
pl.legend(loc='best', fontsize=18)

# Plot source trajectory and caustic
pl.figure(figsize=(6, 6))
planet.plot_trajectory(t_range=[t_0 - t_E, t_0], caustics=True, color='red', lw=2)
planet.plot_trajectory(t_range=[t_0, t_0 + t_E], caustics=True, color='blue', lw=2)
pl.xlim(-0.25, 1.0)
pl.ylim(-0.25, 1.0)
pl.title('Source Trajectory')

pl.show()
