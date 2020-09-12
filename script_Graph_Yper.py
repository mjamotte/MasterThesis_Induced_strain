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
import lib_Honeycomb_Lattice as FHoney

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)


rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

########################## Graphics for Honeycomb Yper ##########################
Save = True
edge = False
edge_LL = False
localized = False

tau = 0.005
Nx = 241

BoundCond = ['Yper','Ypert2t3','Yper_abs_strain','Yper_opp_strain','Yper_Lifsh']
BC = BoundCond[0]

temporary = '%s_Nx_%.1i_tau_%.3f' % (BC,Nx,tau) # put the correct BC to load the correct file
data = np.load('Datas/'+temporary+'.npy',allow_pickle=True)

system = data[0]
kya_list = data[1]
eigVal = data[2]
eigVect = data[3]
colors_all = data[4]
BC = system[0]
Nx = system[1]
Ny = system[2]
gamma2 = system[5]
gamma3 = system[6]
sites_dic = FHoney.honey_lattice(Nx,Ny)

fig_Yper1 = pyplot.figure(figsize=(4,4))
ax1 = pyplot.subplot(111)
#axins = ax1.inset_axes([0.33, 0.55, 0.33, 0.33])


my_ky_s = [1.21] # to pick a value near the one indicated
my_eigVal_s = [0.001]
fig_Yper2, axs = pyplot.subplots(1, len(my_ky_s), figsize=(len(my_ky_s)*6,6), sharex=True, sharey=True)

for j in range(len(eigVal[0])): 

	if localized==True:
		norm = pyplot.Normalize(colors_all.min(),colors_all.max())
		Cmap = 'coolwarm'
		
		points = np.array([kya_list,eigVal[:,j]]).T.reshape(-1, 1, 2)
		segments = np.concatenate([points[:-1], points[1:]], axis=1)
		lc = LineCollection(segments, cmap=Cmap, norm=norm)
		lc.set_array(colors_all[:,j])
		lc.set_linewidth(3)

		# lc_axins = LineCollection(segments, cmap=Cmap, norm=norm)
		# lc_axins.set_array(colors_all[:,j])
		# lc_axins.set_linewidth(3)
		line = ax1.add_collection(lc)
		# line_axins = axins.add_collection(lc_axins)

	if localized==False:
		ax1.plot(kya_list,eigVal[:,j],'-k')
	
	for i in range(len(my_ky_s)):
		kya0,index_kya0 = FHoney.nearest_value(kya_list,my_ky_s[i]) # select a value of ky to plot an edge state
		E,index_eigVal = FHoney.nearest_value(eigVal[index_kya0],my_eigVal_s[i])

		if j==index_eigVal:# and edge_LL==False:
			x_s = FHoney.extract_coord_from_dict(sites_dic)[0]
			x_s,abs_psi = FHoney.reorder(x_s,abs(eigVect[index_kya0,:,j]))
			#ax1.plot(kya0,E,'x',markersize=10)
			#axins.plot(kya0,E,'x',markersize=10)

			abs_psi_A = abs_psi[0:-1:2]
			abs_psi_B = abs_psi[1:len(abs_psi)+2:2]

			l_B = 3/np.sqrt(2*tau)
			x0 = x_s[int(len(x_s)/2)]#+l_B**2*(kya0-1.21)
			print('B =',2/9*tau)

			if tau>0:
				N = index_eigVal-int(len(eigVal[0])/2)
				if N>=0:
					an_psi_A = abs(FHoney.LL_eigVect(x_s[0:-1:2],N,x0,l_B))
					an_psi_B = abs(FHoney.LL_eigVect(x_s[1:len(x_s)+2:2],N-1,x0,l_B))

				if N<0:
					an_psi_A = abs(FHoney.LL_eigVect(x_s[0:-1:2],N+1,x0,l_B))
					an_psi_B = abs(FHoney.LL_eigVect(x_s[1:len(x_s)+2:2],N+2,x0,l_B))


			if len(my_ky_s)>1:
				axs[i].plot(x_s[0:-1:2],abs_psi_A**2/np.max(abs_psi_A**2),\
					'--r',fillstyle='none',label='Num. $|\psi_A|^2$')
				axs[i].plot(x_s[1:len(abs_psi)+2:2],abs_psi_B**2/np.max(abs_psi_B**2)\
					,'--b',fillstyle='none',label='Num. $|\psi_B|^2$')

				if tau>0:
					N = index_eigVal-int(len(eigVal[0])/2)
					if N>=0:
						an_psi_A = abs(FHoney.LL_eigVect(x_s[0:-1:2],N,x0,l_B))
						an_psi_B = abs(FHoney.LL_eigVect(x_s[1:len(x_s)+2:2],N-1,x0,l_B))

					if N<0:
						an_psi_A = abs(FHoney.LL_eigVect(x_s[0:-1:2],N+1,x0,l_B))
						an_psi_B = abs(FHoney.LL_eigVect(x_s[1:len(x_s)+2:2],N+2,x0,l_B))
						axs[i].plot(x_s[0:-1:2],an_psi_A**2/np.max(an_psi_A**2)\
							,'-r',label="Analyt. $|\psi_A|^2$")
						axs[i].plot(x_s[1:len(x_s)+2:2],an_psi_B**2/np.max(an_psi_B**2)\
							,'-b',label="Analyt. $|\psi_B|^2$")


				axs[i].set_xlabel('$x/a$',fontsize='20')
				axs[i].grid(axis='both')
				axs[i].xaxis.set_tick_params(labelsize=20)
				axs[i].yaxis.set_tick_params(labelsize=20)


			if len(my_ky_s)==1:
				axs.plot(x_s[0:-1:2],abs_psi_A**2/np.max(abs_psi_A**2),\
					'--r',fillstyle='none',label='Num. $|\psi_A|^2$')
				axs.plot(x_s[1:len(abs_psi)+2:2],abs_psi_B**2/\
					np.max(abs_psi_B**2),'--b',fillstyle='none',label='Num. $|\psi_B|^2$')

				if tau>0:
					axs.plot(x_s[0:-1:2],an_psi_A**2/np.max(an_psi_A**2)\
						,'-r',label="Analyt. $|\psi_A|^2$")
					axs.plot(x_s[1:len(x_s)+2:2],an_psi_B**2/np.max(an_psi_B**2)\
						,'-b',label="Analyt. $|\psi_B|^2$")


				axs.set_xlabel('$x/a$',fontsize='20')
				axs.grid(axis='both')
				axs.xaxis.set_tick_params(labelsize=20)
				axs.yaxis.set_tick_params(labelsize=20)

			
			#ax1.plot(kya_list,eigVal[:,j].real,'k')

			# ax2.plot(x_s[0:-1:2],an_psi_A**2/np.max(an_psi_A**2),'-r',label="$|\psi|^2(x)$")
			# ax2.plot(x_s[1:len(abs_psi)+2:2],an_psi_B**2/np.max(an_psi_B**2),'-b',label="$|\psi|^2(x)$")	
if localized==True:
	cbar = fig_Yper1.colorbar(line,ticks=[colors_all.min(),colors_all.max()],ax=ax1)
	cbar.ax.set_yticklabels(['Center','Edge'],fontsize=15)

#axins.imshow(Z2, extent=extent, origin="lower")
# sub region of the original image

# x1, x2, y1, y2 = 1.21-0.1, 1.21+0.1, -0.1, 0.1
# axins.set_xlim(x1, x2)
# axins.set_ylim(y1, y2)
# axins.set_yticklabels([])
					
	

if tau>0 and edge_LL==False:
	flag = 0
	#plot the fits of relativistic Landau Levels
	DP1 = np.array([0,-4*np.pi/3/np.sqrt(3),+1])
	DP2 = np.array([0,-2*np.pi/3/np.sqrt(3),-1])
	DP3 = np.array([0,+2*np.pi/3/np.sqrt(3),+1])
	DP4 = np.array([0,+4*np.pi/3/np.sqrt(3),-1])
	DP_s = [DP1,DP2,DP3,DP4]
	#DP_s = [DP3]
	for DP in DP_s:
		kya_list_approx = np.linspace(DP[1]-0.1,DP[1]+0.1,100)
		kxa = 0
		ns = np.arange(0,10)
		E_p = np.zeros((len(ns),len(kya_list_approx)))
		E_m = np.zeros((len(ns),len(kya_list_approx)))
		E_p_tilt = np.zeros((len(ns),len(kya_list_approx)))
		E_m_tilt = np.zeros((len(ns),len(kya_list_approx)))
		for n in ns:
			i = 0
			for kya in kya_list_approx:
				E_p[n,i] = +np.sqrt(abs(n)*tau)
				E_m[n,i] = -E_p[n,i]
				E_p_tilt[n,i],E_m_tilt[n,i] = FHoney.E_relat_LL(n,system,kxa,kya,DP)
				i += 1

		# for n in ns:
		# 	if flag==0:
		# 		ax1.plot(kya_list_approx,E_p[n],'r-',linewidth=1,label="$E \\propto \\sqrt{N\\tau}$")
		# 		ax1.plot(kya_list_approx,E_m[n],'r-',linewidth=1)
		# 		ax1.plot(kya_list_approx,E_p_tilt[n],'g--',linewidth=1,label="$E \\propto \\sqrt{N\\tau(1-q_y a)}$")
		# 		ax1.plot(kya_list_approx,E_m_tilt[n],'g--',linewidth=1)
		# 		flag = 1
		# 	else:
		# 		ax1.plot(kya_list_approx,E_p[n],'r--',linewidth=1)
		# 		ax1.plot(kya_list_approx,E_m[n],'r--',linewidth=1)
		# 		ax1.plot(kya_list_approx,E_p_tilt[n],'g--',linewidth=1)
		# 		ax1.plot(kya_list_approx,E_m_tilt[n],'g--',linewidth=1)

# add more features to the figure

#if BC!='Ypert2t3':
#	pyplot.suptitle('%s, Nx = %.1i, $\\tau$ = %.3f' % (BC,Nx, tau))
#else:
#	pyplot.suptitle('%s, Nx = %.1i, $\\tau$ = %.3f,gam_2 = %.2f, gam_3 = %.2f' % (BC,Nx, tau,gamma2,gamma3))

ax1.grid(axis='both')
ax1.set_xlim([-1.8,1.8])
ax1.set_ylim([-0.2,0.2])
ax1.set_xlabel('$k_y a$',fontsize='25')
ax1.set_ylabel('$E/t$',fontsize='25')
ax1.tick_params(axis='both',labelsize=16)

# if edge_LL==False:
# 	ax1.legend(loc=1,fontsize='15')
# 	#ax1.legend(bbox_to_anchor=(0.5, 0.4),fontsize='18')
# if len(my_ky_s)>1:
# 	axs[-1].legend(bbox_to_anchor=(0.17, 0.5),fontsize='15')
# if len(my_ky_s)==1:
# 	axs.legend(loc=2,fontsize='15')

# axins3 = ax1.inset_axes([0.4, 0.4, 0.2, 0.33])

# kya0 = 1.21
# kya0,index_kya0 = FHoney.nearest_value(kya_list,kya0) # select a value of ky to plot an edge state
# axins3.plot(np.arange(len(eigVal[0])),eigVal[index_kya0],'o',markersize=3,label='$k_y a = %.2f$'%(kya0))
# axins3.set_ylim(-0.4,0.4)
# axins3.set_xlim(Nx*Ny/2-10,Nx*Ny/2+10)
# axins3.set_xlabel('\\# eigenvalues',fontsize='14')
# axins3.set_ylabel('$E/t$',fontsize='14')

# axins3.tick_params(axis='x',labelsize=16)
# axins3.tick_params(axis='y',labelsize=16)

# axins3.legend(fontsize=12,loc='upper center')


#pyplot.plot(np.arange(Nx*Ny),abs(eigVect[0])**2)
pyplot.show()

if Save==True:
	if edge_LL==True:
		temp = '%s_Nx_%.1i_tau_%.3f' % (BC,Nx, tau)
		#path_above_Maxime = (sys.path[0]).split('Maxime')[0]
		#fig_Yper1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
		path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux

		#fig_Yper1.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig_Honey_edgeLL'+temp+".pdf") # Linux
		#fig_Yper2.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_Graph'+temp+".pdf") # Windows
		#fig_Yper2.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig_HoneyGraph_edge_LL'+temp+".pdf") # Linux
		#fig_Yper2.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_Graph'+temp+".pdf") # Windows	else:
		# temp = '%s_Nx_%.1i_tau_%.3f' % (BC,Nx, tau)
		# #path_above_Maxime = (sys.path[0]).split('Maxime')[0]
		# #fig_Yper1.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
		# path_above_maxime = (sys.path[0]).split('maxime')[0] # Linux

		# fig_Yper1.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig_Honey'+temp+".pdf") # Linux
		# #fig_Yper2.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_Graph'+temp+".pdf") # Windows
		# fig_Yper2.savefig(path_above_maxime+'/maxime/ownCloud/MasterThesis/figures/fig_HoneyGraph'+temp+".pdf") # Linux
