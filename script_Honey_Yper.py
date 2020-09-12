import numpy as np
import scipy as sc
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
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

##################################### TEST HONEY Yper #####################################


tau_list = [0.01] #[0,0.005,0.015,0.03]
Nx_s = [241] #[121,241]
Ny = 2
gamma2 = 0.2
gamma3 = 0.2

BoundCond = ['Yper','Ypert2t3','Yper_abs_strain','Yper_opp_strain','Yper_Lifsh']
BC = BoundCond[0]

t = 1
a1 = np.array([1,0])
a2 = np.array([0,1])
My = 400
#kya_list = np.linspace(-4*np.pi/3/np.sqrt(3),4*np.pi/3/np.sqrt(3),My)
kya_list = np.linspace(2*np.pi/3/np.sqrt(3)-0.3,2*np.pi/3/np.sqrt(3)+0.3,My)
kxa = 0

select_eigVal =  Nx_s[0]*Ny# must be even
colors_all = np.zeros((My,select_eigVal))

start = time.time()
for Nx in Nx_s:
	for tau in tau_list:

		sites_dic = FHoney.square_lattice(Nx,Ny,a1,a2)
		x_s = OHL.extract_coord_from_dict(sites_dic)[0]
		system = [BC,Nx,Ny,t,tau,gamma2,gamma3,sites_dic]

		if tau<=4/(Nx-1): # if no Lifshitz transition

			eigVal = 1j*np.zeros((My,select_eigVal))
			eigVectors = 1j*np.zeros((My,Nx*Ny,select_eigVal))
			colors = np.zeros(len(eigVectors[0]))			

			i = 0
			for kya in kya_list:
				H_ky = FHoney.H_honey(system,kxa,kya,sites_dic)
				#eigVal[i],eigVectors[i] = sc.sparse.linalg.eigsh(H_ky,select_eigVal,which='SM')

				# select eigen values around 0 (that is put in the middle of eigVal_temp by eigh)

				eigVal[i],eigVectors[i] = linalg.eigh(H_ky)

				# eigVal[i] = eigVal_temp[int(len(eigVal_temp)/2)-int(select_eigVal/2):\
				# 			int(len(eigVal_temp)/2)+int(select_eigVal/2)] 
				# eigVectors[i] = eigVect_temp[:,int(len(eigVect_temp[0])/2)-int(select_eigVal/2):\
								#int(len(eigVect_temp[0])/2)+int(select_eigVal/2)]

				
				for m in range(len(eigVectors[0])):
					colors[m] = OHL.color_position(x_s,eigVectors[i,:,m],Nx)
				#print(colors)
				colors_all[i] = colors#[int(len(eigVal)/2)-int(select_eigVal/2):\
									#int(len(eigVal)/2)+int(select_eigVal/2)]

				i += 1

			# Save in a file
			data = [system, kya_list, eigVal, eigVectors,colors_all]
			dataID = '%s_Nx_%.1i_tau_%.3f' % (system[0],system[1], system[4])
			np.save('Datas/'+dataID+'.npy',data)

print("Execution time:",time.time()-start,"secondes")
			
# pyplot.pcolormesh(abs(H_ky))
# pyplot.colorbar()

# pyplot.show()