import numpy as np
import scipy as sc
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import rc

from numpy import linalg
import time

import os
import sys
path_above_Codes = (sys.path[0]).split('Codes')[0]
sys.path.append(path_above_Codes+'Codes/')
## Import the libraries in the parent folder "Codes"
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf
import lib_Optical_Honeycomb_Lattice as OHL

##################Script  OL 1D IT ##########################

print('"script_OL_Honey_bos_SC.py" is running')

#path = '/home/maxime/Desktop/Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_SC'
#list_args = ['N','V','U']
#split = ['N_','V_','U_']
params_interv = {'U': [1e-5],'V': [1e-6],\
				'N': [1e5]}
Nx_s = [201]
#files_names = Mf.select_files(path,extension,list_args,params_interv,split)
for Nx in Nx_s:
	k = 0
	for N in params_interv['N']:
		for V in params_interv['V']:
			for U in params_interv['U']:
				start = time.time()

				args_syst = {
				'syst' : 'Honey',	
				'J' : 1,
				'N' : N,
				'V' : V,
				'Nx' : Nx,
				'Ny' : 2,
				'U' : U,
				'Method' : 'SC',
				'Trap' : 'Harmonic',
				'Symm' : 'Isotropic', # could be "x_strain"
				'BC' : 'Yper'
				}
				args_syst.update({'sites_dic' : OHL.honey_lattice(args_syst)})

				## Kinetic + Trap part of Hamiltonian
				ka = 0
				H_KV = OHL.H_KV_honey_BEC(args_syst,ka)
				
				if k==0 or U==0:
				 	args_init = {
				 	'H_KV' : H_KV,
				 	'err_SC' : 1e-9
				 	}

				if k>0 and U!=0:
					args_init.update({'psi_init' : psi0})

				psi0 = GP.calc_psi0(args_syst,args_init)

				mu = linalg.eigh(H_KV+GP.H_int(psi0,args_syst))[0][0]

				E_funct_SC = GP.energy_functional(psi0,args_syst)

				#Trap = OHL.trap(args_syst)
				
				# U = args_syst['U']
				# N = args_syst['N']
				Nx = args_syst['Nx']
				# n_TF = (mu-Trap)/U/N
				# diff_TF_n0 = np.sum(np.abs(n_TF-abs(psi0)**2))

				data = [args_syst,args_init,mu,psi0,E_funct_SC,H_KV]#,diff_TF_n0]
				dataID = 'Honey_%s_%s_%s_Nx_%.1i_J_%.2f_N_%.4e_V_%.3e_U_%.3e' %\
						(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
						args_syst['Nx'],args_syst['J'],N,args_syst['V'],args_syst['U'])
						
				temp = "Datas_SC_Nx"+str(Nx)+"/"
				if not os.path.exists(temp):
					os.makedirs(temp)
				np.save("Datas_SC_Nx"+str(Nx)+"/"+dataID,data)

				print("For Nx = ", Nx, ", N = ", N, ", U = ", U, ", V = ", V,\
					"it took",time.time()-start,"secondes")

				k += 1
# if not os.path.exists('Hello'):
#     os.makedirs('Hello')

# pyplot.pcolormesh(np.abs(H_KV))
# pyplot.colorbar()
# pyplot.show()