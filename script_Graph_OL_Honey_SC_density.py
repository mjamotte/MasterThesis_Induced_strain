import numpy as np
import matplotlib.colors as color
import matplotlib.axes as axes
import matplotlib.pyplot as pyplot
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import rc
from numpy import linalg
import cmath
import time
from mpl_toolkits import mplot3d
from random import randint

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

Nx = 401

# E_funct_IT_s = np.array([])
# diff_all = np.array([])
# args_diff_all = np.array([])

## Figure

fig_OL_Honey = pyplot.figure(figsize=(8,8))
ax2 = pyplot.subplot(111)
ax2.set_xlabel('$x/a$',fontsize=20)
#ax2.set_ylabel('$|\psi_0|^2$')#,color='r')

## Select the files you want to plot
path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_SC_Nx'+str(Nx)
list_args = ['N','V','U']
split = ['N_','V_','U_']
params_interv = {'U': [1e-3],'V': [3.5e-5],'N': [1e5]}
extension = '.npy'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)

if len(params_interv['U'])==1:

	for filename in files_names:

		data = np.load(filename+'.npy',allow_pickle=True)
		args_syst = data[0]
		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']
		sites_dic = args_syst['sites_dic']

		args_init = data[1]
		#mu_all = data['mu_all']	
		mu = data[2]+3*J
		psi0 = data[3]
		x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
		x_s,psi0 = OHL.reorder(x_s,psi0)

		E_funct_SC = data[4]
		
		n0 = np.abs(psi0)**2*N

		## Analytical gaussian 

		Trap = OHL.trap(args_syst)
		x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
		x_s,Trap = OHL.reorder(x_s,Trap)

		m = 1/(2*J)
		w0 = np.sqrt(V/m)
		R = np.sqrt(8*mu/V)

		ax2.plot(x_s,n0,'ok',label="$n_0, U = %.1e$" % (U),fillstyle='none')

		## Thomas-Fermi
		if U!=0:
			n_TF = (mu+3*J-Trap)/U# divide by N because here, normalised to 1
			n_TF = n_TF-np.max(n_TF)+np.max(n0)

			## Plot
			
			#Hermite_0 = (m*w0/np.pi)**0.25*np.exp(-m*w0*(positions-x0)**2/2)
			#ax2.plot(x_s,Hermite_0**2,'*g')
			n_TF = OHL.remove_under0(n_TF) # negative values are put to 0

			ax2.plot(x_s,n_TF,'-r',label="$n_{TF}$, $U = %.1e$" % (U),linewidth=2)

		axins = ax2.inset_axes([0.66, 0.66, 0.33, 0.33]) #([0.66,0.66,0.33,0.33]) #
		#axins.imshow(Z2, extent=extent, origin="lower")
		# sub region of the original image
		x1, x2, y1, y2 = np.max(x_s)/2+R/2-20, np.max(x_s)/2+R/2+10,-1, 80 #np.max(x_s)/2+R/2-15, np.max(x_s)/2+R/2+5,-1, 80 #
		axins.set_xlim(x1, x2)
		axins.set_ylim(y1, y2)
		#axins.set_yticklabels([])
		axins.plot(x_s,n0,'ok',fillstyle='none')
		axins.plot(x_s,n_TF,'-r')
	
	#axins.legend(loc=2,fontsize=13);

	# pyplot.suptitle('T-F approx. %s, %s,\
	# 	 Nx = %.1i, N = %.3e, J = %.2f, \n V = %.3e, U = %.3e' % \
	# 	 (args_syst['Trap'],args_syst['Symm'],\
	# 	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V']\
	# 	,args_syst['U']))

	#ax2.set_ylim(0,0.015)
	ax2.grid(axis='both')
	ax2.legend(loc=2,fontsize=16);
	pyplot.xticks(fontsize=16)
	pyplot.yticks(fontsize=16)
	pyplot.show()

	## Save the figure

	temp = 'Honey_density_%s_%s_%s_Nx_%.1i_J_%.3f_V_%.3e_U_%.3e' %\
		(args_syst['Method'],args_syst['Trap'],\
		args_syst['Symm'],args_syst['Nx'],args_syst['J'],\
		args_syst['V'],args_syst['U'])

	#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
	#fig_OL_Honey.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
	path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
	fig_OL_Honey.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

else:

	for filename in files_names[::7]:

		data = np.load(filename+'.npy',allow_pickle=True)
		args_syst = data[0]
		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']
		sites_dic = args_syst['sites_dic']

		args_init = data[1]
		#mu_all = data['mu_all']	
		mu = data[2]+3*J
		psi0 = data[3]
		x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
		x_s,psi0 = OHL.reorder(x_s,psi0)

		E_funct_SC = data[4]
		
		n0 = np.abs(psi0)**2*N

		## Analytical gaussian 

		Trap = OHL.trap(args_syst)
		x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
		x_s,Trap = OHL.reorder(x_s,Trap)

		m = 1/(2*J)
		w0 = np.sqrt(V/m)
		R = np.sqrt(8*mu/V)

		ax2.plot(x_s,n0,'o',label="$n_0, U = %.1e$" % (U),fillstyle='none')

		## Thomas-Fermi
		if U!=0:
			n_TF = (mu-Trap)/U# divide by N because here, normalised to 1
			n_TF = n_TF-np.max(n_TF)+np.max(n0)

			## Plot
			
			#Hermite_0 = (m*w0/np.pi)**0.25*np.exp(-m*w0*(positions-x0)**2/2)
			#ax2.plot(x_s,Hermite_0**2,'*g')
			n_TF = OHL.remove_under0(n_TF) # negative values are put to 0
			ax2.plot(x_s,n_TF,'-',label="$n_{TF}, U = %.1e$" % (U),linewidth=1)
	
	# pyplot.suptitle('T-F approx. %s, %s,\
	# 	 Nx = %.1i, N = %.3e, J = %.2f, \n V = %.3e, U = %.3e' % \
	# 	 (args_syst['Trap'],args_syst['Symm'],\
	# 	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V']\
	# 	,args_syst['U']))

	ax2.set_xlim(150)
	ax2.grid(axis='both')
	ax2.legend(loc=1,fontsize=14);
	pyplot.xticks(fontsize=16)
	pyplot.yticks(fontsize=16)
	pyplot.show()

	temp = 'Honey_densities_Us_%s_%s_%s_Nx_%.1i_J_%.3f_V_%.3e_U_%.3e' %\
		(args_syst['Method'],args_syst['Trap'],\
		args_syst['Symm'],args_syst['Nx'],args_syst['J'],\
		args_syst['V'],args_syst['U'])
	#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
	#fig_OL_Honey.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
	#path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
	#fig_OL_Honey.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux