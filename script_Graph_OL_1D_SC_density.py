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

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

########################## Graphics for OPTICAL LATTICE ##########################

E_funct_IT_s = np.array([])
diff_all = np.array([])
args_diff_all = np.array([])

## Figure

fig_OL1D = pyplot.figure(figsize=(8,8))
ax2 = pyplot.subplot(111)
ax2.set_xlabel('$x/a$',fontsize=20)
#ax2.set_ylabel('$|\psi_0|^2$',fontsize=20)#,color='r')

## Select the files you want to plot
path = path_above_Codes+'Codes/Optical_1D_Lattice_BEC/Self_consistent/Datas_SC'
list_args = ['U','V','N']
split = ['U_','V_','N_']
params_interv = {'U': [1e-5,2e-4],'V': [3.5e-5],'N': [1e5]}
extension = '.npy'
files_names = Mf.select_files(path,extension,list_args,params_interv,split)

for filename in files_names[::3]:

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

	n0 = np.abs(psi0)**2

	## Analytical gaussian 

	Trap = GP.trap_1D(args_syst)
	m = 1/(2*J)
	w0 = np.sqrt(V/m)
	xc = (Nx-1)/2
	sites_dic = args_syst['sites_dic']
	positions = OHL.extract_coord_from_dict(sites_dic)[0]
	positions,psi = OHL.reorder(positions,n0)

	## Thomas-Fermi

	#mu = mu_all[-1,int(len(mu_all[0])/2)]
	n_TF = (mu-Trap)/U/N# divide by because here, normalised to 1
	n_TF = n_TF-np.max(n_TF)+np.max(n0)

	## Plot

	#Hermite_0 = (m*w0/np.pi)**0.25*np.exp(-m*w0*(positions-x0)**2/2)
	#ax2.plot(positions,Hermite_0**2,'*g')
	ax2.plot(positions,n0*N,'-',label="$n_0$, $U = %.2e$" %(U),linewidth=3)#,label="$N = {}$".format(N))
	ax2.plot(positions,n_TF*N,'--',label="$n_{TF}$, $U = %.2e$" %(U))

## Save the figure
# pyplot.suptitle('T-F approx. %s, %s,\
# 	 Nx = %.1i, N = %.3e, J = %.2f, V = %.3e, U = %.3e' % \
# 	 (args_syst['Trap'],args_syst['Symm'],\
# 	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V']\
# 	,args_syst['U']))

ax2.set_ylim(0,2200)
ax2.set_xlim(xc)
pyplot.xticks(fontsize=16)
pyplot.yticks(fontsize=16)
ax2.grid(axis='both')
ax2.legend(fontsize=18)
pyplot.show()

temp = '1D_density_Ns_%s_%s_%s_Nx_%.1i_J_%.3f_V_%.3e_U_%.3e' %\
		(args_syst['Method'],args_syst['Trap'],\
		args_syst['Symm'],args_syst['Nx'],args_syst['J'],\
		args_syst['V'],args_syst['U'])

#fig_OL1D.savefig("Figures_OL_1D_BEC_SC\\comp_U"+".pdf")
#fig_OL1D.savefig("Figures_OL_1D_BEC_IT\\fig_"+temp+".pdf")

#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
#fig1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
fig_OL1D.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux