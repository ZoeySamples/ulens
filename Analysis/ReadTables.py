# Zoey Samples
# Created: Jun 04, 2018
# ReadTables.py
# Last Updated: Jun 21, 2018

import sys
import os
import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from HDUL import HDUL


file_name = sys.argv[1:]
hdul = []
for i in range(len(file_name)):
	param = {'file_name': file_name[i]}
	hdul.append(HDUL(**param))

for i in range(len(hdul)):

	hdul[i].plot_num_images(errors_only=False, print_errors=True)
	plt.gcf().set_size_inches(10, 7)
	plt.show()
#	hdul[i].plot_magnification(cutoff=None, log_colorbar=True,
#							   outliers=False, save = False)
#	plt.gcf().set_size_inches(10, 7)
#	plt.show()



