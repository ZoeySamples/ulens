# Zoey Samples
# Created: Jun 04, 2018
# ReadTables.py
# Last Updated: Jun 11, 2018

import sys
import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from HDUL import HDUL


file_name = sys.argv[1:]
hdul = []
for i in range(len(file_name)):
	param = {'file_name': file_name[i], 'cutoff':1000}
	hdul.append(HDUL(**param))

hdul[0].plot_n_solutions_errors()
hdul[0].plot_magn_outliers()
plt.show()
