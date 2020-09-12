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

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def color2grey(colors_all,R,Lx):
	'''
	This function is specific to this script. It transforms the colours attributed
	to the mean positions of the eigenvectors. The variable colors_all contains 
	values from -x/(Lx/2) to +x/(Lx/2) and we want to transform them into 0 to x/R
	with different shades of grey. This allows to know if the eigenvectors are 
	centered around the edges of the cloud. 
	'''
	out = np.zeros((len(colors_all),len(colors_all[0])))
	for i in range(len(colors_all)):
		for j in range(len(colors_all[0])):
			out[i,j] = abs(1-2*(abs(abs(colors_all[i,j])-R/2/(Lx/2))))
	return out

########################## Graphics for OPTICAL LATTICE ##########################
print('"script_Graph_OL_Honey_spectYper.py" is running')
'''
This code is used to compare to analytical relations between the energies, 
the number of particles and the convergence of the self-consistent method 
in minimizing E-mu*N 
'''

only_strain = False
no_strain = False
in_out = True
is_edge = False
LL = False


Nx = 201
alpha = 0.1

## Select the files you want to plot

path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_SC_Nx'+str(Nx)
list_args = ['N','V','U']
split = ['N_','V_','U_']
params_interv = {'U': [1e-4],'V': [1e-5],'N': [1e5]}
extension = '.npy'
files_names_bos = Mf.select_files(path,extension,list_args,params_interv,split)

path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_Ferm_SC_Nx'+str(Nx)+'/alpha'+str(alpha)
files_names = Mf.select_files(path,extension,list_args,params_interv,split)

#split = ['U_']
#list_args = ['U']
files_names_ferm = Mf.reorder_filesnames(files_names,list_args,split)
#print(files_names)

fig1 = pyplot.figure(figsize=(10,8))
ax1 = pyplot.subplot(111)

for i,filename_ferm in enumerate(files_names_ferm):

	data_ferm = np.load(filename_ferm+'.npy',allow_pickle=True)
	args_syst = data_ferm[0]
	eigVal_all = data_ferm[1].real
	eigVal_TF_all = data_ferm[2].real
	kya_list = data_ferm[3]
	mu = data_ferm[4]
	n0 = data_ferm[5]
	colors_all = data_ferm[6]

	J = args_syst['J']
	Nx = args_syst['Nx']
	V = args_syst['V']
	U = args_syst['U']
	N = args_syst['N']
	alpha = args_syst['alpha']
	R = np.sqrt(mu*8/V)
	kxa = 2*np.pi/3#0
	sites_dic = args_syst['sites_dic']

	DP_y = 2*np.pi/3/np.sqrt(3)
	zoom_DPy = np.linspace(DP_y-0.1,DP_y+0.1,50)

	for j in range(len(eigVal_all[0])):
		if j==0:
			#ax1.plot(kya_list,eigVal_all[:,j],'-k')#,linewidth=3,label='$n_0$')
			ax1.plot(kya_list,eigVal_TF_all[:,j],color='green')#lighten_color('k', 0.4), lw=1,label='$n_{TF}$')
		else:
			#ax1.plot(kya_list,eigVal_all[:,j],'-k',linewidth=3)
			ax1.plot(kya_list,eigVal_TF_all[:,j],color='green')#lighten_color('k', 0.4), lw=1)

	## Lines with gradient of color

	if is_edge==True:
		colors_all = color2grey(colors_all,R,Nx*3/2)
		norm = pyplot.Normalize(colors_all.min(),colors_all.max())
		Cmap = 'coolwarm'
		for j in range(len(eigVal_all[0])):
			

			points = np.array([kya_list,eigVal_all[:,j]]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)
			lc = LineCollection(segments, cmap=Cmap, norm=norm)
			lc.set_array(colors_all[:,j])
			lc.set_linewidth(3)
			line = ax1.add_collection(lc)

		cbar = fig1.colorbar(line,ticks=[0.05, 0.95],ax=ax1)
		cbar.ax.set_yticklabels(['Not edge', 'Edge'])

	if in_out==True:
		norm = pyplot.Normalize(colors_all.min(),colors_all.max())
		Cmap = 'cool'

		for j in range(len(eigVal_all[0])):
			
			points = np.array([kya_list,eigVal_all[:,j]]).T.reshape(-1, 1, 2)
			segments = np.concatenate([points[:-1], points[1:]], axis=1)
			lc = LineCollection(segments, cmap=Cmap, norm=norm)
			lc.set_array(colors_all[:,j])
			lc.set_linewidth(3)
			line = ax1.add_collection(lc)
		cbar = fig1.colorbar(line,ticks=[-0.95, 0.95],ax=ax1)
		cbar.ax.set_yticklabels(['Left', 'Right'])

		
		# if j==0:
		# 	ax1.plot(kya_list,eigVal_TF_all[:,j], color=lighten_color('k', 0.4), lw=1,label='TF-regime')
		# else:
		# 	ax1.plot(kya_list,eigVal_TF_all[:,j], color=lighten_color('k', 0.4), lw=1)
	if LL==True:
		for j in range(len(eigVal_all[0])):
			frac = 3/2
			if j==0:
				ax1.plot(zoom_DPy,np.sqrt(j*alpha*V/U*frac)*np.sqrt(1-(zoom_DPy-DP_y)),\
					'--r',linewidth=1,label='$E \\propto \\sqrt{V_0 (1-q_ya)/U}$')
				ax1.plot(zoom_DPy,-np.sqrt(j*alpha*V/U*frac)*np.sqrt(1-(zoom_DPy-DP_y)),\
					'--r',linewidth=1)
			else:
				ax1.plot(zoom_DPy,np.sqrt(j*alpha*V/U*frac)*np.sqrt(1-(zoom_DPy-DP_y)),\
					'--r',linewidth=1)
				ax1.plot(zoom_DPy,-np.sqrt(j*alpha*V/U*frac)*np.sqrt(1-(zoom_DPy-DP_y)),\
					'--r',linewidth=1)
	
	ax1.xaxis.set_tick_params(labelsize=20)
	ax1.yaxis.set_tick_params(labelsize=20)
## Plot the eigenvector for specified values of ky and E.
	for filename_bos in files_names_bos:
		my_ky_s = [0.99,1.1503,1.21,1.3] # to pick a value near the one indicated
		my_eigVal_s = [0.1006,0.10650,0.0731,0.0052]
		data_bos = np.load(filename_bos+'.npy',allow_pickle=True)
		n0 = np.abs(data_bos[3])**2
		x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
		x_s,n0 = OHL.reorder(x_s,n0)

		fig_Yper2, axs = pyplot.subplots(int(len(my_ky_s)/2), int(len(my_ky_s)/2), figsize=(int(len(my_ky_s)/2)*6,(len(my_ky_s)/2)*6), sharex=True, sharey=True)
		
		m = 0
		for i in range(int(len(my_ky_s)/2)):
			for k in range(int(len(my_ky_s)/2)):
				kya0,index_kya0 = GP.nearest_value(kya_list,my_ky_s[m]) # select a value of ky to plot an edge state
				#print(eigVal_all[index_kya0])
				E,index_eigVal = GP.nearest_value(eigVal_all[index_kya0],my_eigVal_s[m])
				#print('E = ',E,'kya = ',kya0)
				my_eigVect,E,ky = OHL.get_eigVect(filename_ferm,\
									filename_bos,my_ky_s[m],my_eigVal_s[m])
				x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
				x_s,my_eigVect = OHL.reorder(x_s,my_eigVect)

				#print(OHL.color_outCloud(x_s,my_eigVect,0))
				#pyplot.plot(x_s,my_eigVect)

				abs_psi_A = np.abs(my_eigVect[0:-1:2])
				abs_psi_B = np.abs(my_eigVect[1:len(my_eigVect)+2:2])

				R,Lx = GP.size_cloud(args_syst,mu,n0)
				#print(OHL.get_color(my_eigVect,x_s,R))

				axs[i,k].plot(x_s[0:-1:2],OHL.relat2max(abs_psi_A**2),\
					'--r',label='$|\psi_A|^2$')
				axs[i,k].plot(x_s[1:len(x_s)+2:2],OHL.relat2max(abs_psi_B**2)\
					,'--b',label='$|\psi_B|^2$')

				## TF
				n_TF = 1/U*(mu-V/2*(x_s-np.max(x_s)/2)**2)
				shifted_n_TF = n_TF/N-np.max(n_TF)/N+np.max(n0)
				axs[i,k].plot(x_s,shifted_n_TF/np.max(shifted_n_TF),'m')

				
				axs[i,k].plot(x_s,n0/np.max(n0),'g',label='$n_0$')

				# if i<=2:
				# 	ax1.plot(ky,E,'^',fillstyle='none',markersize=15)
				# if i>2:
				# 	ax1.plot(ky,E,'s',fillstyle='none',markersize=15)

				l_B = 3/np.sqrt(2*alpha*V/U)
				x0 = x_s[int(len(x_s)/2)]#+l_B**2*(kya0-1.21)

				#print(len(eigVal_all[0]))
				n = index_eigVal-int(len(eigVal_all[0])/2)
				print('N=',N,'Index eig Val = ',index_eigVal)
				an_psi_A = abs(OHL.LL_eigVect(x_s[0:-1:2],n,x0,l_B))
				an_psi_B = abs(OHL.LL_eigVect(x_s[1:len(x_s)+2:2],n-1,x0,l_B))

				axs[i,k].xaxis.set_tick_params(labelsize=20)
				axs[i,k].yaxis.set_tick_params(labelsize=20)
				axs[i,k].legend(loc='upper right',fontsize=20)
				axs[i,k].set_ylim(0,1.05)
				axs[i,k].set_xlabel('$x/a$',fontsize=20)
				m += 1
				# axs[i].plot(x_s[0:-1:2],OHL.relat2max(an_psi_A**2)\
				# 	,'-r',label="Analyt. $|\psi_A|^2$")
				# axs[i].plot(x_s[1:len(x_s)+2:2],OHL.relat2max(an_psi_B**2)\
				# 	,'-b',label="Analyt. $|\psi_B|^2$")

	#pyplot.xlabel('$x$',fontsize='25')
	#pyplot.ylabel('$|\Psi|^2/\max(|\Psi|^2)$',fontsize='25')
	# pyplot.suptitle('Yper eigenstate OL Honey %s, %s, %s,Nx = %.1i, N = %.3e, \nJ = %.2f, V = %.3e, U = %.3e, alpha = %.3e' % \
	# 	 (args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
	# 	args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V'],\
	# 	args_syst['U'],args_syst['alpha']))
	#axs[i].grid(axis='both')
	

	ax1.set_xlim(0.7,1.7)
	ax1.set_ylim(-0.01,0.3)
	ax1.set_xlabel('$k_y a$',fontsize=20)
	ax1.set_ylabel('$E/t$',fontsize=20)
	
	pyplot.xticks(fontsize=18)
	pyplot.yticks(fontsize=18)

# if in_out==True:
# 	fig3 = pyplot.figure(figsize=(8,8))
# 	ax3 = pyplot.subplot(111)

# 	cmap = ListedColormap(['g', 'r','g'])
# 	val = ['-Lx/2','-R/2','R/2','Lx/2']
# 	Norm = [-Lx/2/Lx,-R/2/Lx, +R/2/Lx,Lx/2/Lx]
# 	norm = BoundaryNorm(Norm, cmap.N)
	

# 	for j in range(len(eigVal_all[0])):
		
# 		points = np.array([kya_list,eigVal_all[:,j]]).T.reshape(-1, 1, 2)
# 		segments = np.concatenate([points[:-1], points[1:]], axis=1)
# 		lc = LineCollection(segments, cmap=cmap, norm=norm)
# 		lc.set_array(colors_all[:,j])
# 		lc.set_linewidth(2)
# 		line = ax3.add_collection(lc)

# 	ax3.set_xlim(1,1.45)
# 	ax3.set_ylim(-0.01,0.15)
# 	ax3.set_xlabel('$k_y a$',fontsize=20)
# 	ax3.set_ylabel('$E/t$',fontsize=20)

# ## Superpose with only-strain part (the not-strained honeycomb is removed)
if only_strain==True:

	Nx = 201
	alpha = 0.1
	## Select the files you want to plot

	path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_SC_Nx'+str(Nx)
	list_args = ['N','V','U']
	split = ['N_','V_','U_']
	params_interv = {'U': [1e-3],'V': [3.5e-5],'N': [1e5]}
	extension = '.npy'
	files_names_bos = Mf.select_files(path,extension,list_args,params_interv,split)

	path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_Ferm_SC_Nx'+str(Nx)+'/alpha'+str(alpha)
	files_names = Mf.select_files(path,extension,list_args,params_interv,split)

	#split = ['U_']
	#list_args = ['U']
	files_names_ferm = Mf.reorder_filesnames(files_names,list_args,split)
	#print(files_names)

	for i,filename_ferm in enumerate(files_names_ferm):

		data_ferm = np.load(filename_ferm+'.npy',allow_pickle=True)
		args_syst = data_ferm[0]
		eigVal_all = data_ferm[1].real
		eigVal_TF_all = data_ferm[2].real
		kya_list = data_ferm[3]
		mu = data_ferm[4]
		n0 = data_ferm[5]
		colors_all = data_ferm[6]

		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']
		alpha = args_syst['alpha']
		kxa = 2*np.pi/3#0
		sites_dic = args_syst['sites_dic']

		for j in range(len(eigVal_all[0])):
			if j==0:
				ax1.plot(kya_list,eigVal_all[:,j],'--k',linewidth=1,label='Only strained')
			else:
				ax1.plot(kya_list,eigVal_all[:,j],'--k',linewidth=1)

## Superpose with no-strain part
if no_strain==True:

	Nx = 105
	alpha = 0

	## Select the files you want to plot

	path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_SC_Nx'+str(Nx)
	list_args = ['N','V','U']
	split = ['N_','V_','U_']
	params_interv = {'U': [1e-3],'V': [3.5e-5],'N': [1e5]}
	extension = '.npy'
	files_names_bos = Mf.select_files(path,extension,list_args,params_interv,split)

	path = path_above_Codes+'Codes/Optical_Honeycomb_Lattice/Self_consistent/Datas_Ferm_SC_Nx'+str(Nx)+'/alpha'+str(alpha)
	files_names = Mf.select_files(path,extension,list_args,params_interv,split)

	#split = ['U_']
	#list_args = ['U']
	files_names_ferm = Mf.reorder_filesnames(files_names,list_args,split)
	#print(files_names)

	for i,filename_ferm in enumerate(files_names_ferm):

		data_ferm = np.load(filename_ferm+'.npy',allow_pickle=True)
		args_syst = data_ferm[0]
		eigVal_all = data_ferm[1].real
		eigVal_TF_all = data_ferm[2].real
		kya_list = data_ferm[3]
		mu = data_ferm[4]
		n0 = data_ferm[5]
		colors_all = data_ferm[6]

		J = args_syst['J']
		Nx = args_syst['Nx']
		V = args_syst['V']
		U = args_syst['U']
		N = args_syst['N']
		alpha = args_syst['alpha']
		kxa = 2*np.pi/3#0
		sites_dic = args_syst['sites_dic']

		for j in range(len(eigVal_all[0])):
			if j==0:
				ax1.plot(kya_list,eigVal_all[:,j],'-.g',linewidth=1,label='Unstrained')
			else:
				ax1.plot(kya_list,eigVal_all[:,j],'-.g',linewidth=1)

		# ax1.xlabel('$k_y a$',fontsize='25')
		# ax1.ylabel('$E$',fontsize='25')
		# pyplot.xticks(fontsize=16)
		# pyplot.yticks(fontsize=16)
		# pyplot.grid(axis='both')
		#pyplot.legend(loc=2);
#ax1.legend(loc='lower left',fontsize=18)		
#ax1.legend(fontsize=16,bbox_to_anchor=(0.35, 0.22));



# if only_strain==False and no_strain==False:

# 	pyplot.show()

# 	## Save the figure
# 	temp = 'Ferm_spectrum_Yper_Honey_%s_%s_%s_Nx_%.1f_J_%.2f_V_%.1e_U_%.3e_a_%.2e' %\
# 			(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
# 			args_syst['Nx'],args_syst['J'],args_syst['V'],args_syst['U'],\
# 			args_syst['alpha'])
# 	#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
# 	#fig1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
# 	path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
# 	fig1.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

# if only_strain==True or no_strain==True:

# 	pyplot.show()

# 	## Save the figure
# 	temp = 'Ferm_spectrum_Honey_with_nostrain_%s_%s_%s_Nx_%.1f_J_%.2f_V_%.1e_U_%.3e_a_%.2e' %\
# 			(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
# 			args_syst['Nx'],args_syst['J'],args_syst['V'],args_syst['U'],\
# 			args_syst['alpha'])
# 	#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
# 	#fig1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
# 	path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
# 	fig1.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

if in_out==True:

	pyplot.show()

# 	## Save the figure
# 	temp = 'Ferm_spectrum_Honey_inoutCloud_%s_%s_%s_Nx_%.1f_J_%.2f_V_%.1e_U_%.3e_a_%.2e' %\
# 			(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
# 			args_syst['Nx'],args_syst['J'],args_syst['V'],args_syst['U'],\
# 			args_syst['alpha'])
# 	#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
# 	#fig1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
# 	path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
# 	fig1.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

# if is_edge==True:
# 	pyplot.show()

# 	## Save the figure
# 	temp = 'Ferm_spectrum_Honey_edgeCloud_%s_%s_%s_Nx_%.1f_J_%.2f_V_%.1e_U_%.3e_a_%.2e' %\
# 			(args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
# 			args_syst['Nx'],args_syst['J'],args_syst['V'],args_syst['U'],\
# 			args_syst['alpha'])
# 	#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
# 	#fig1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
# 	path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux
# 	fig1.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf") # Linux

# #fig_Yper2.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney_eigVect'+'.pdf') # Linux

#pyplot.show()