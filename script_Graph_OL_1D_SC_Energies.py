import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from numpy import linalg
import cmath
import time
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

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

########################## Graphics for OPTICAL LATTICE ##########################
'''
This code is used to compare to analytical relations between the energies, 
the number of particles and the convergence of the self-consistent method 
in minimizing E-mu*N 
'''
plot_Vs = False
plot_Us = True


if plot_Us==True:

	fig_OL1D_Us = pyplot.figure(figsize=(8,8))
	ax2 = pyplot.subplot(111)

	E_funct_SC_s = np.array([])
	E_k_s = []
	E_U_s = []
	E_tr_s = []
	U_s = []
	mu_s = []

	## Select the files you want to plot
	path = path_above_Codes+'Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC'
	list_args = ['N','V','U']
	split = ['N_','V_','U_']
	params_interv = {'U': [1e-6,1e-3],'V': [0.5e-4],\
				'N': [1e4]}
	extension = '.npy'
	files_names = Mf.select_files(path,extension,list_args,params_interv,split)

	split = ['U_']
	list_args = ['U']
	files_names = Mf.reorder_filesnames(files_names,list_args,split)
	#print(files_names)
	
	for filename in files_names:

		data = np.load(filename+'.npy',allow_pickle=True)
		args_syst = data[0]
		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']

		args_init = data[1]
		#mu_all = data['mu_all']	
		mu = data[2]
		psi0 = data[3]
		E_funct_SC = data[4]

		## Compute contributions to energy

		E_k = GP.E_kin(psi0,args_syst)
		E_U = GP.E_int(psi0,args_syst)
		E_tr = GP.E_trap(psi0,args_syst)
		E_k_s.append(E_k)
		E_U_s.append(E_U)
		E_tr_s.append(E_tr)
		U_s.append(U)
		mu_s.append(mu)

	E_k_s = np.array(E_k_s)
	E_U_s = np.array(E_U_s)
	E_tr_s = np.array(E_tr_s)
	mu_s = np.array(mu_s)
	U_s = np.array(U_s)
	#print(U_s)

	## To see how E_kin and E_trap vary in fct of U
	
	ax2.semilogx(U_s,E_k_s/N+2*J,'-b',label="$E_k/N + 2J$",linewidth=2)#,fillstyle='none')
	ax2.semilogx(U_s,E_tr_s/N,'-g',label="$E_{trap}/N$",linewidth=2)#,fillstyle='none')
	ax2.semilogx(U_s,E_U_s/N,'-r',label="$E_U/N$",linewidth=2)#,fillstyle='none')
	#ax2.semilogx(U_s,(2*E_tr_s-2*E_k_s)/N-2*(2*J),'-k',label="Vir. theor.")#,fillstyle='none')
	#ax2.semilogx(U_s,2/5*(mu_s+2*J),'--r',label="$2(\mu+2J)/5$")
	#ax2.semilogx(U_s,1/5*(mu_s+2*J),'--g',label="$(\mu+2J)/5$")
	print("N = ", N, "U = ", U, "V0 = ", V,\
		": U*N*|phi_max|^2 = ",U*N*np.max(np.abs(psi0)**2))

	ax2.set_xlabel('$U$',fontsize=18)

	# axins = ax2.inset_axes([0.33, 0.6, 0.45, 0.33])
	# #axins.imshow(Z2, extent=extent, origin="lower")
	# # sub region of the original image
	# x1, x2, y1, y2 = 1e-5, 5e-5, 0.001, 0.007
	# axins.set_xlim(x1, x2)
	# axins.set_ylim(y1, y2)
	# axins.semilogx(U_s,E_k_s/N+2*J,'-b',fillstyle='none')#,label="$E_k/N + 2J$")
	# axins.semilogx(U_s,E_tr_s/N,'-g',fillstyle='none')#,label="$E_{trap}/N$")
	# axins.semilogx(U_s,E_U_s/N,'-r',fillstyle='none')#,label="$E_U/N$")

	# axins.semilogx(U_s,2/9*(mu_s+2*J),'--m',label="$2(\mu+2J)/9$",linewidth=1)
	# axins.semilogx(U_s,2/7*(mu_s+2*J),'--c',label="$2(\mu+2J)/7$",linewidth=1)
	# axins.set_yticklabels([])

	# axins.legend(fontsize=17);

	# elif what2plot=='mu_E':			
	# 	ax1.plot(np.log10(N),E_k/N+2*J,'*b',label="$E_k$")
	# 	ax1.plot(np.log10(N),E_tr/N,'*g',label="$E_t$")
	# 	ax1.plot(np.log10(N),E_U/N,'*r',label="$E_U$")
	# 	ax1.plot(np.log10(N),E_funct_IT.real/N,'o',label="$N = {}$".format(N))
	# 	for x in range(len(mu_all[0])):
	# 		ax2.plot(np.arange(len(mu_all)),mu_all[:,x].real) # mu.imag is zerO

	## Save the figure

	# pyplot.suptitle('%s, %s, %s,\
	# 	 Nx = %.1i, N = %.3e, J = %.2f, V = %.3e' % \
	# 	 (args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
	# 	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V']))

	ax2.grid(axis='both')
	ax2.legend(bbox_to_anchor=(0.1, 0.17),fontsize=17);
	pyplot.xticks(fontsize=16)
	pyplot.yticks(fontsize=16)
	pyplot.show()

	temp = '1D_Energies_Us_%s_%s_%s_Nx_%.1f_J_%.2f_V_%.1e' %\
			(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
			args_syst['Nx'],args_syst['J'],args_syst['V'])

	#path_above_Maxime = (sys.path[0]).split('Maxime')[0]
	#fig_OL1D_Us.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\Figures_OL_1D_BEC_SC'+temp+".pdf") # Windows
	#path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
	#fig_OL1D_Us.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

###################################
###################################

if plot_Vs==True:
	## Initialize the figure

	fig_OL1D_Vs = pyplot.figure(figsize=(8,8))
	ax1 = pyplot.subplot(111)

	E_funct_SC_s = np.array([])
	E_k_s = np.array([])
	E_U_s = np.array([])
	E_tr_s = np.array([])
	w0_s = np.array([])

	## Select the files you want to plot

	path = path_above_Codes+'Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC'
	list_args = ['N','V','U']
	split = ['N_','V_','U_']
	params_interv = {'U': [0],'V': [2e-6,1e-2],'N': [1e4]}
	extension = '.npy'
	files_names = Mf.select_files(path,extension,list_args,params_interv,split)
	
	split = ['V_']
	list_args = ['V']
	files_names = Mf.reorder_filesnames(files_names,list_args,split)

	for filename in files_names:

		data = np.load(filename+'.npy',allow_pickle=True)
		args_syst = data[0]
		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']
		m = 1/(2*J)

		args_init = data[1]
		#mu_all = data['mu_all']	
		mu = data[2]
		psi0 = data[3]
		E_funct_SC = data[4]

		## Compute contributions to energy

		E_k = GP.E_kin(psi0,args_syst)
		E_U = GP.E_int(psi0,args_syst)
		E_tr = GP.E_trap(psi0,args_syst)
		E_k_s = np.append(E_k_s,E_k)
		E_U_s = np.append(E_U_s,E_U)
		E_tr_s = np.append(E_tr_s,E_tr)
		w0_s = np.append(w0_s,np.sqrt(V/m))

	## To see how E_kin and E_trap vary in fct of V
	
	#print("w0/2 = ", w0/2)
	ax1.semilogx(w0_s,E_k_s/N+2*J,'sb',markersize=10,fillstyle='none',label="$E_{kin}/N+2J$")
	ax1.semilogx(w0_s,E_tr_s/N,'+g',markersize=10,fillstyle='none',label="$E_{trap}/N$")
	#ax1.plot(w0_s,E_U_s/N,'^r',fillstyle='none',label="$E_U$")
	ax1.semilogx(w0_s,w0_s/4,'--k',fillstyle='none',label="$\\hbar\\omega/4$",linewidth=1)

	ax1.set_xlabel('$\omega$',fontsize=20)

	# Save the figure

	# pyplot.suptitle('1D Ener. fct of V : %s, %s,\
	# 	 Nx = %.1i, N = %.3e, J = %.2f, U = %.3e' % \
	# 	 (args_syst['Trap'],args_syst['Symm'],\
	# 	args_syst['Nx'],args_syst['N'],args_syst['J'],\
	# 	args_syst['U']))

	ax1.grid(axis='both')
	ax1.legend(loc=2,fontsize=18);
	# ax1.set_xticklabels(['2 $\\times 10^{-3}$','$10^{-2}$ '])
	# ax1.set_xticks([0.002,0.01])
	pyplot.xticks(fontsize=16)
	pyplot.yticks(fontsize=16)
	pyplot.show()

	temp = '1D_Energies_Vs_%s_%s_%s_Nx_%.1f_J_%.2f_U_%.3e' %\
			(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
			args_syst['Nx'],args_syst['J'],args_syst['U'])
	#path_above_Maxime = (sys.path[0]).split('Maxime')[0]
	#fig_OL1D_Vs.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\Figures_OL_1D_BEC_SC'+temp+".pdf") # Windows
	path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
	fig_OL1D_Vs.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux