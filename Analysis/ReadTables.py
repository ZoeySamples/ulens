# Zoey Samples
# Created: Jun 04, 2018
# ReadTables.py
# Last Updated: Jun 11, 2018

import sys
import os
import matplotlib.pyplot as plt
from BinaryLens import BinaryLens as BL
from HDUL import HDUL


file_name = sys.argv[1:]
hdul = []
for i in range(len(file_name)):
	param = {'file_name': file_name[i], 'cutoff':1000}
	hdul.append(HDUL(**param))

for i in range(len(hdul)):
	hdul[i].plot_n_solutions_errors()
	hdul[i].plot_magn_outliers()
	plt.gcf().set_size_inches(10, 6)
	j = 0
	while True:
		j += 1
		file_name = ('../Tables/errors_{}_{}{}.png'.format(hdul[i].solver_file,
								hdul[i].origin_file, j))
		path = os.path.abspath(__file__)
		for k in range(2):
			path = os.path.dirname(path)
		file_path = os.path.join(path, file_name[3:])
		if os.path.exists(file_path):
		    continue
		plt.savefig(file_name)
		plt.show()
		print(file_name, 'has been saved')
		break


