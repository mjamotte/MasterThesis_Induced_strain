import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import cmath
import time
from mpl_toolkits import mplot3d

import nbimporter
import Fun_Honeycomb as FHoney

########################## Graphics for Honeycomb Xper ##########################

# Create the figure

fig_Xper_taus_Nys = pyplot.figure(figsize=(20,8))
ax1 = pyplot.subplot(111)
ax1.set_xlabel('$Numbers$')
ax1.set_ylabel('$E$')
pyplot.grid(axis='both')

# load the files to plot with corresponding parameters

tau_list = [0.005,0.015,0.02]
Nx_s = [2,121]
Ny_s = [2,99]
BC = 'Trigonal'

for Nx in Nx_s:
	for Ny in Ny_s:
		if (Nx==2 and Ny!=2) or (Ny==2 and Nx!=2) and Nx*Ny<10**4:
			for tau in tau_list:

				if tau<=2/(Ny-1): # if no Lifshitz transition

					temporary = '%s_Nx_%.1iNy_%.1i_tau_%.3f' % (BC,Nx,Ny,tau)
					data = np.load('Datas/'+temporary+'.npy',allow_pickle=True)

					# Definition of variables

					system = data[0]
					eigVal = data[1]

					numbers = np.arange(0,len(eigVal))
					
					# plot
					ax1.plot(numbers,eigVal.real,'.',label = 'Nx_%.1iNy_%.1i_tau_%.3f' % \
						(Nx,Ny, tau))

ax1.legend(loc=1);
pyplot.show()
fig_Xper_taus_Nys.savefig("Figures_Honey/fig_Honey_Xper_taus_Ns.pdf")