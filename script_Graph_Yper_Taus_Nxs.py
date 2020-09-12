import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import rc
from numpy import linalg
import cmath
import time
from mpl_toolkits import mplot3d

import os
import sys
path_above_Codes = (sys.path[0]).split('Codes')[0]
sys.path.append(path_above_Codes+'Codes/')
## Import the libraries in the parent folder "Codes"
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf
import lib_Optical_Honeycomb_Lattice as OHL
import lib_Honeycomb_Lattice as FHoney


rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

########################## Graphics for Honeycomb Yper ##########################

Save = True

# Create the figure

fig_Yper_taus_Nxs = pyplot.figure(figsize=(12,6))
ax1 = pyplot.subplot(111)
ax1.set_xlabel('$k_y a$',fontsize=20)
ax1.set_ylabel('$E/t$',fontsize=20)
pyplot.grid(axis='both')

# load the files to plot with corresponding parameters

tau_list = [0.005,0.01,0.02]
Nx_s = [121,241,361]
BC = 'Yper' # put the correct BC to load the correct file

for Nx in Nx_s:
	for tau in tau_list:

		if tau<=4/(Nx-1)*0.85: # if no Lifshitz transition

			temp = '%s_Nx_%.1i_tau_%.3f' % (BC,Nx,tau)
			data = np.load('Datas/'+temp+'.npy',allow_pickle=True)

			# Definition of variables

			system = data[0]
			kya_list = data[1]
			eigVal = data[2]
			eigVect = data[3]

			# plot

			select_eigVal = int(len(eigVal[0])/2)+1
			ax1.plot(kya_list,eigVal[:,select_eigVal].real,label = '$N_x$=%.1i, $\\tau$=%.3f' % \
				(Nx, tau))

ax1.set_xlim(0.7,1.7)
ax1.set_ylim(-0.05,0.45)
ax1.legend(loc='best',fontsize=17);
pyplot.show()

if Save==True:
	temp = '%s_Nx_%.1i_tau_%.3f' % (BC,Nx, tau)
	path_above_Maxime = (sys.path[0]).split('Maxime')[0]

	# path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
	#fig_Yper1.savefig(path_above_Maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

	fig_Yper_taus_Nxs.savefig(path_above_Maxime+"\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_"+BC+"_taus_Nxs.pdf") # Windows