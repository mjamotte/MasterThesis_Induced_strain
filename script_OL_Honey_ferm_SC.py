import numpy as np
#import scipy as sc
#from scipy.sparse.linalg import eigsh
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
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
print('"script_OL_Honey_ferm_SC.py" is running')

Nx = 201
alpha_s = [0.1]#np.linspace(0,0.2,5)
## Import datas
path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_SC_Nx'+str(Nx)
list_args = ['N','V','U']
split = ['N_','V_','U_']
params_interv = {'U': [1e-4],'V': [1e-5],\
				'N': [1e5]}
0.00
extension = '.npy'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)

for file_name in files_names:
	## Initialize datas
	data_bos =  np.load(file_name+extension,allow_pickle=True)
	args_syst = data_bos[0]
	N = args_syst['N']
	U = args_syst['U']
	J = args_syst['J']
	V = args_syst['V']
	Ny = args_syst['Ny']

	mu = data_bos[2]+3*J
	psi0 = data_bos[3]
	n0 = np.abs(psi0)**2

	x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]

	E_k = GP.E_kin(psi0,args_syst)
	E_U = GP.E_int(psi0,args_syst)
	print('E_kin/E_U = ',(E_k+3*N)/E_U)
	E_tr = GP.E_trap(psi0,args_syst)
	print('E_kin/E_tr = ',(E_k+3*N)/E_tr)

	n_TF = (mu-OHL.trap(args_syst))/U/N
	n_TF = OHL.remove_under0(n_TF) # negative values are put to 0
	R,Lx = GP.size_cloud(args_syst,mu,n0)
	print("The size of the cloud is", int(R/Lx*100)%100, "%","of the size of the system")
	OHL.compare2TF(args_syst,n0,n_TF) # returns a plot
	#pyplot.show()
	#exit()
	for alpha in alpha_s:

		start = time.time()

		args_syst.update({'alpha' : alpha})
		
		## Generate the spectrum/eigenstates
		DP_y = 2*np.pi/3/np.sqrt(3)
		My = 400
		kya_list = np.linspace(2*np.pi/3/np.sqrt(3)-1,2*np.pi/3/np.sqrt(3)+1,My)

		select_eigVal = Nx*Ny # must be even
		eigVal_all = 1j*np.zeros((My,select_eigVal))

		colors_all = np.zeros((My,select_eigVal)) # 3 = len(RGB)

		eigVal_TF_all = 1j*np.zeros((My,select_eigVal))

		i = 0
		for kya in kya_list:

			H_c = OHL.H_honey_c(n0*N,args_syst,kya)
			H_c_TF = OHL.H_honey_c(n_TF*N,args_syst,kya)
			eigVal,eigVectors = linalg.eigh(H_c)

			colors = np.zeros(len(eigVectors[0]))
			for m in range(len(eigVectors[0])):
				colors[m] = OHL.color_position(x_s,eigVectors[:,m],R)

			eigVal_TF = linalg.eigh(H_c_TF)[0]

			eigVal_all[i] = eigVal[int(len(eigVal)/2)-int(select_eigVal/2):\
									int(len(eigVal)/2)+int(select_eigVal/2)] 
			colors_all[i] = colors[int(len(eigVal)/2)-int(select_eigVal/2):\
									int(len(eigVal)/2)+int(select_eigVal/2)]

			eigVal_TF_all[i] = eigVal_TF[int(len(eigVal_TF)/2)-int(select_eigVal/2):\
									int(len(eigVal_TF)/2)+int(select_eigVal/2)] 

			i += 1


		data_ferm = [args_syst,eigVal_all,eigVal_TF_all,kya_list,mu,n0,colors_all]
		dataID = 'Ferm_Honey_%s_%s_%s_Nx_%.1i_J_%.2f_N_%.4e_V_%.3e_U_%.3e' %\
				(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
				args_syst['Nx'],args_syst['J'],args_syst['N'],\
				args_syst['V'],args_syst['U'])
		print('Coucou')
		temp = "Datas_Ferm_SC_Nx"+str(Nx)+"/"+"alpha"+str(alpha)+"/"
		if not os.path.exists(temp):
		    os.makedirs(temp)

		np.save(temp+dataID,data_ferm)


# print("For Nx = ", Nx, ", N = ", N, ", U = ", U, ", V = ", V,\
# 	"it took",time.time()-start,"secondes")

