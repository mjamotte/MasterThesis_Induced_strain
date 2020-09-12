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

############ Comparison between different BC ############

Nx = 241
tau = 0.008

BC_toCompare = ['Yper','Ypert2t3','Yper_abs_strain']

fig_Yper_BCcompar = pyplot.figure(figsize=(8,8))
ax1 = pyplot.subplot(111)
ax1.set_xlabel('$k_y a$')
ax1.set_ylabel('$E$')

for i in range(len(BC_toCompare)):

	temporary = '%s_Nx_%.1i_tau_%.3f' % (BC_toCompare[i],Nx,tau) # put the correct BC to load the correct file
	data = np.load('Datas/'+temporary+'.npy',allow_pickle=True)

	system = data[0]
	kya_list = data[1]
	eigVal = data[2]
	eigVect = data[3]
	BC = system[0]
	Ny = system[2]
	gamma2 = system[5]
	gamma3 = system[6]

	# choose the 
	index_lvl_to_print = np.array([int(len(eigVal[0])/2)-2,int(len(eigVal[0])/2)-2,\
						int(len(eigVal[0])/2)-4])
	j = index_lvl_to_print[i]

	ax1.plot(kya_list,eigVal[:,j].real,'.',label= '%s' % (BC_toCompare[i]))
	ax1.grid(axis='both')

ax1.legend(loc=1);
pyplot.show()

temp = 'Yper_Comparison'
fig_Yper_BCcompar.savefig("Figures_Honey/fig_"+temp+".pdf")
