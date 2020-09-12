import numpy as np
import scipy as sc
from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import time

import sys
sys.path.append('/home/maxime/Desktop/Codes/')
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf

##################Script  OL 1D IT ##########################
'''
#path = '/home/maxime/Desktop/Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC'
#list_args = ['N','V','U']
#split = ['N_','V_','U_']
params_interv = {'U': np.logspace(-6,-4,10),'V': [1e-5],\
				'N': [1e4-1,1e4,1e4+1]}
#extension = '.npy'
#files_names = Mf.select_files(path,extension,list_args,params_interv,split)
k = 0
for U in params_interv['U']:
	for V in params_interv['V']:
		Efun_SC = np.array([]) # to compute mu = dE_SC/dN
		for N in params_interv['N']:
			start = time.time()

			args_syst = {
			'syst' : '1D',
			'J' : 1,
			'N' : N,
			'V' : V,
			'Nx' : 201,
			'U' : U,
			'Method' : 'SC',
			'Trap' : 'Harmonic',
			'Symm' : 'Isotropic',
			}

			args_syst.update({'sites_dic' : GP.lattice_1D(args_syst)})

			## Kinetic + Trap part of Hamiltonian

			H_KV = GP.H_1D(args_syst)

			if k==0 or U==0:
			 	args_init = {
			 	'H_KV' : H_KV,
			 	'err_SC' : 1e-10
			 	}

			if k>0 and U!=0:
				args_init.update({'psi_init' : psi0})

			psi0 = GP.calc_psi0(args_syst,args_init)

			mu = linalg.eigh(H_KV+GP.H_int(psi0,args_syst))[0][0]

			E_funct_SC = GP.energy_functional(psi0,args_syst)
			Efun_SC = np.append(Efun_SC,E_funct_SC)
		
		dEfun_SCdN = GP.dEdN_O2(Efun_SC,1)

		dataID = 'U_'+str(U)
		data = [U,mu,dEfun_SCdN]
		
		np.save('Datas_SC_comp_mu_dEdN/'+dataID,data)

		print(time.time()-start,"secondes")

		k += 1
'''
path = '/home/maxime/Desktop/Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC_comp_mu_dEdN/'
list_args = ['U']
split = ['U_']
extension = '.npy'
params_interv = {'U': np.logspace(-6,-4,10)}
files_names = Mf.select_files(path,extension,list_args,params_interv,split)
split = ['U_']
list_args = ['U']
files_names = Mf.reorder_filesnames(files_names,list_args,split)

fig_OL1D = pyplot.figure(figsize=(8,8))
dEfun_SCdN_s = []
U_s = []
mu_s = []
for filename in files_names:

	data = np.load(filename+'.npy',allow_pickle=True)
	U_s.append(data[0])
	mu_s.append(data[1])
	dEfun_SCdN_s.append(data[2])

pyplot.semilogx(U_s,mu_s,'^g',fillstyle='none',label="$\mu$")
pyplot.semilogx(U_s,dEfun_SCdN_s,'sr',fillstyle='none',label="$dE/dN$")

args_syst = {
			'syst' : '1D',
			'J' : 1,
			'N' : 1e4,
			'V' : 1e-5,
			'Nx' : 201,
			'Method' : 'SC',
			'Trap' : 'Harmonic',
			'Symm' : 'Isotropic',
			}

pyplot.suptitle('Compar. $dE/dN$ and $\mu$ (SC) for %s, %s,\
	 Nx = %.1i, J = %.2f, V = %.3e' % \
	 (args_syst['Trap'],args_syst['Symm'],\
	args_syst['Nx'],args_syst['J'],args_syst['V']))

pyplot.xlabel('$U$')
pyplot.legend(loc=2);
pyplot.grid(axis='both')
pyplot.show()

temp = 'Comp_mu_dEdN'
fig_OL1D.savefig("Figures_OL_1D_BEC_SC/fig_"+temp+".pdf")
