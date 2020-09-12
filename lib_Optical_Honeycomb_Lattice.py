import numpy as np
import sympy as sy
from sympy import *
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.axes as axes
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib import rc
from numpy import linalg
import cmath
import time
from mpl_toolkits import mplot3d

import sys
sys.path.append('/home/maxime/Desktop/Codes/')
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf

##############################################################

def honey_lattice(args_syst): # boundaries = zigzag

	####################################################################################
	#																				   #	
	#	Constructs a square lattice with Bravais vectors a1 and a2, and composed of Nx # 
	# 	sites in the x direction and Ny sites in the y direction.					   #
	#																				   #	
	####################################################################################

	Nx = args_syst['Nx']
	Ny = args_syst['Ny']
	sites_dic = {}
	n = 0
	x = 0

	for j in range(Ny):

		## create the scheme to repeat along y axis
		x = 0
		y = j*np.sqrt(3)/2

		for i in range(Nx): # fill line by line, from left to right
			
			if i%2==0:
				sites_dic[n] = np.array([x,y])
				x = x+2

			else:
				sites_dic[n] = np.array([x,y])
				x = x+1
				
			n += 1	

		if j%2!=0:
			for k in range(Nx):
				L = len(sites_dic)
				if k%2==0:
					sites_dic[L-Nx+k][0] += 0.5
				else:
					sites_dic[L-Nx+k][0] -= 0.5
			
	return sites_dic


def draw_Honey(args_syst,drawAB,put_delta,put_t,params,Brillouin):

	fs = params['fontsize']
	lw = params['linewidth']

	pyplot.axes()
	ax = pyplot.gca()

	sites_dic = honey_lattice(args_syst)
	x_s,y_s = extract_coord_from_dict(sites_dic)

	Nx = args_syst['Nx']
	Ny = args_syst['Ny']

	

	if args_syst['BC']=='Xper':

		for i in range(0,len(sites_dic),2):

			if int(y_s[i]/(np.sqrt(3)/2))%2==0:

				if x_s[i]+1/2<np.max(x_s):
					line = pyplot.Line2D((x_s[i],x_s[i]-1/2),\
							(y_s[i],y_s[i]),lw=lw,color='green',ls='--')
					ax.add_line(line)
					line = pyplot.Line2D((x_s[i]+2,x_s[i]+2+1/2),\
							(y_s[i],y_s[i]),lw=lw,color='green',ls='--')
					ax.add_line(line)

				if y_s[i]+np.sqrt(3)/2<=np.max(y_s):
					line = pyplot.Line2D((x_s[i],x_s[i]+1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),color='black',lw=2+y_s[i])
					ax.add_line(line)
					line = pyplot.Line2D((x_s[i]+2,x_s[i]+1+1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),color='black',lw=2+y_s[i])
					ax.add_line(line)

				if drawAB==True:
					
					pyplot.plot(x_s[i],y_s[i],'or',markersize=12)
					pyplot.plot(x_s[i]+2,y_s[i],'ob',markersize=12)

			if int(y_s[i]/(np.sqrt(3)/2))%2!=0:

				line = pyplot.Line2D((x_s[i],x_s[i]+1),\
						(y_s[i],y_s[i]), lw=lw,color='green')
				ax.add_line(line)
				
				if y_s[i]+np.sqrt(3)/2<=np.max(y_s):
						line = pyplot.Line2D((x_s[i],x_s[i]-1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),color='black',lw=2+y_s[i])
						ax.add_line(line)
						line = pyplot.Line2D((x_s[i]+1,x_s[i]+1+1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),color='black',lw=2+y_s[i])
						ax.add_line(line)

				if drawAB=='True':
					
					pyplot.plot(x_s[i],y_s[i],'ob')
					pyplot.plot(x_s[i]+1,y_s[i],'or')


	if args_syst['BC']=='Yper':

		for i in range(len(x_s)):

			if i%2!=0: # draw aline to the right from a B site

				if x_s[i]+1<=np.max(x_s):

					line = pyplot.Line2D((x_s[i],x_s[i]+1),\
							(y_s[i],y_s[i]),lw=(1+x_s[i]),color='black')
					ax.add_line(line)

					if i<Nx:
						line = pyplot.Line2D((x_s[i],x_s[i]-1/2),\
						(y_s[i],y_s[i]-np.sqrt(3)/2), lw=lw,color='cyan',ls='--')
						ax.add_line(line)	

				if i>=Nx:		
					line = pyplot.Line2D((x_s[i],x_s[i]-1/2),\
					(y_s[i],y_s[i]+np.sqrt(3)/2),lw=lw,color='green',ls='--')
					ax.add_line(line)					

				if y_s[i]+np.sqrt(3)/2<=np.max(y_s):
					line = pyplot.Line2D((x_s[i],x_s[i]-1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),color='green',lw=lw)
					ax.add_line(line)
			else:
				if y_s[i]+np.sqrt(3)/2<=np.max(y_s):
					line = pyplot.Line2D((x_s[i],x_s[i]+1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),color='cyan',lw=lw)
					ax.add_line(line)

					if i<Nx:
						line = pyplot.Line2D((x_s[i],x_s[i]+1/2),\
						(y_s[i],y_s[i]-np.sqrt(3)/2),lw=lw,color='green',ls='--')
						ax.add_line(line)

				if i>=Nx:		
					line = pyplot.Line2D((x_s[i],x_s[i]+1/2),\
					(y_s[i],y_s[i]+np.sqrt(3)/2),lw=lw,color='cyan',ls='--')
					ax.add_line(line)	

			if drawAB==True:
				for i in range(len(x_s)):
					if i%2==0:
						pyplot.plot(x_s[i],y_s[i],'or',markersize=12)

					else:
						pyplot.plot(x_s[i],y_s[i],'ob',markersize=12)

	if args_syst['BC']=='Open':

		for i in range(len(x_s)):

			if i%2!=0: # draw aline to the right from a B site

				if x_s[i]+1<=np.max(x_s):

					line = pyplot.Line2D((x_s[i],x_s[i]+1),\
							(y_s[i],y_s[i]),lw=lw,color='black')
					ax.add_line(line)					

				if y_s[i]+np.sqrt(3)/2<=np.max(y_s):
					line = pyplot.Line2D((x_s[i],x_s[i]-1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),lw=lw,color='k')
					ax.add_line(line)
			else:
				if y_s[i]+np.sqrt(3)/2<=np.max(y_s):
					line = pyplot.Line2D((x_s[i],x_s[i]+1/2),\
							(y_s[i],y_s[i]+np.sqrt(3)/2),lw=lw,color='k')
					ax.add_line(line)
					
			if drawAB==True:
				for i in range(len(x_s)):
					if i%2==0:
						pyplot.plot(x_s[i],y_s[i],'ok',markersize=12)

					else:
						pyplot.plot(x_s[i],y_s[i],'ok',markersize=12)

			if put_delta==True:
				pyplot.text(0.95,0.95,'$\delta_1$',fontsize=fs)
				pyplot.text(1.8,0.5,'$\delta_3$',fontsize=fs)
				pyplot.text(1.8,1.1,'$\delta_2$',fontsize=fs)

				pyplot.arrow(1.5,np.sqrt(3)/2,-0.95,0,\
					length_includes_head=True,head_width=0.1,color='k')
				pyplot.arrow(1.5,np.sqrt(3)/2,0.5-0.02,np.sqrt(3)/2-0.03,\
					length_includes_head=True,head_width=0.1,color='k')
				pyplot.arrow(1.5,np.sqrt(3)/2,0.5-0.02,-np.sqrt(3)/2+0.03,\
					length_includes_head=True,head_width=0.1,color='k')

				pyplot.arrow(1.5,3/2*np.sqrt(3),3/2-0.05,np.sqrt(3)/2-0.05,\
					length_includes_head=True,head_width=0.1,color='k')
				pyplot.arrow(1.5,3/2*np.sqrt(3),3/2-0.05,-np.sqrt(3)/2+0.05,\
					length_includes_head=True,head_width=0.1,color='k')
				pyplot.text(2.25,2.85,'$\mathbf{a}_1$',fontsize=fs)
				pyplot.text(2.25,2.25,'$\mathbf{a}_2$',fontsize=fs)
				pyplot.text(0.98,2.4,'$a$',fontsize=fs)

			if put_t==True:
				pyplot.text(0.95,0.95,'$t_1(x)$',fontsize=fs)
				pyplot.text(1.8,0.5,'$t_3$',fontsize=fs)
				pyplot.text(1.8,1.1,'$t_2$',fontsize=fs)

	# if drawAB==True:
	# 	pyplot.text(1.4,0.6,'A',color='r',fontsize=fs)
	# 	pyplot.text(0.5,0.6,'B',color='b',fontsize=fs)

	if Brillouin==True:
		line = pyplot.Line2D((0,0),(0,1),lw=lw,color='k')
		ax.add_line(line)
		line = pyplot.Line2D((0,np.sqrt(3)/2),(1,3/2),lw=lw,color='k')
		ax.add_line(line)
		line = pyplot.Line2D((np.sqrt(3)/2,np.sqrt(3)),(3/2,1),lw=lw,color='k')
		ax.add_line(line)
		line = pyplot.Line2D((np.sqrt(3),np.sqrt(3)),(1,0),lw=lw,color='k')
		ax.add_line(line)
		line = pyplot.Line2D((np.sqrt(3),np.sqrt(3)/2),(0,-1/2),lw=lw,color='k')
		ax.add_line(line)
		line = pyplot.Line2D((np.sqrt(3)/2,0),(-1/2,0),lw=lw,color='k')
		ax.add_line(line)

		pyplot.text(1.8,0,'$K\'$',fontsize=fs)
		pyplot.text(1.8,1,'$K$',fontsize=fs)
		pyplot.text(0.8,-0.65,'$K$',fontsize=fs)
		pyplot.text(-0.2,0,'$K\'$',fontsize=fs)
		pyplot.text(-0.2,1,'$K$',fontsize=fs)
		pyplot.text(0.8,1.6,'$K\'$',fontsize=fs)



	#pyplot.axis('equal')
	pyplot.axis('scaled')

	return ax


def corr_hopping(i1,i2,density,args_syst,delta_link):

	####################################################################################
	# 																				   #
	#	Computes the Hamiltonian element H[i1,i2] where i1 and i2 are the index of the #
	#	sites in the sites dictionnary; 								   			   #
	#																				   #
	####################################################################################
	
	J = args_syst['J']
	alpha = args_syst['alpha']

	if args_syst['BC']=='Yper':

		if delta_link=='delta_1': # horizontal link, i1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = J*(1-alpha/3*(density[i1]-density[i2]))

		elif delta_link=='delta_2':
			hop = J*(1-alpha/3*(density[i1]-density[i2]))

		elif delta_link=='delta_3':
			hop = J*(1-alpha/3*(density[i1]-density[i2]))

	if args_syst['BC']=='Ypert2t3':

		Nx = system[1]
		Lx = 3/2*(Nx-1) # Nx must be odd
		xA = (sites_dic[i1][0])*3/2
		gamma2 = system[5]
		gamma3 = system[6]

		if delta_link=='delta_1': # horizontal link, i1 is always even
			#hop = 2*t # to observe Lifshitz transition
			hop = t*(1-tau/3*(xA-Lx/2))

		if delta_link=='delta_2':
			hop = t*(1+gamma2*tau/3*(xA-Lx/2))

		if delta_link=='delta_3':
			hop = t*(1+gamma3*tau/3*(xA-Lx/2))
	
	if args_syst['BC']=='Xper':

		yA = sites_dic[i1][1]*np.sqrt(3)/2
		Ny = system[2]
		Ly = Ny*np.sqrt(3)/2

		if delta_link=='delta_1':
			hop = t

		elif delta_link=='delta_2':
			hop = t*(1+tau*(yA-Ly/2)*np.sqrt(3)/3)

		elif delta_link=='delta_3':
			hop = t*(1-tau*(yA-Ly/2)*np.sqrt(3)/3)

	return hop


def hop_x(args_syst,density): # computes only on A sites!
	
	sites_dic = args_syst['sites_dic']
	x_s = extract_coord_from_dict(sites_dic)[0]
	x_s,density = reorder(x_s,density)
	hop_s = np.array([])

	for n in range(2,len(sites_dic),2): # x-strain
		hop_s = np.append(hop_s,corr_hopping(n,n-1,density,args_syst,'delta_1')) 

	return x_s,hop_s

def H_honey_c(density,args_syst,ka):

	d1 = np.array([-1,0])
	d2 = np.array([0.5,np.sqrt(3)/2])
	d3 = np.array([0.5,-np.sqrt(3)/2])
	deltas = [d1,d2,d3]
	
	J = args_syst['J']
	V = args_syst['V']
	U = args_syst['U']
	Nx = args_syst['Nx']
	Ny = args_syst['Ny']

	if args_syst['BC']=='Yper' or args_syst['BC']=='Ypert2t3':

		kya = ka
		H = 1j*np.zeros((Nx*Ny,Nx*Ny))
		
		for n in range(Nx*Ny):
			if n%Nx!=0 and n%2==0 and n>0: # n is a A site, this A isn't on the left boundry
				H[n,n-1] = corr_hopping(n,n-1,density,args_syst,'delta_1')
				H[n-1,n] = np.conj(H[n,n-1])

			if n<Nx*(Ny-1) and n%2==0:  # n is a A site, this A is not on the top boundary
				H[n,n+Nx] = corr_hopping(n,n+Nx,density,args_syst,'delta_2')*\
						np.exp(1j*kya*deltas[1][1])+corr_hopping(n,n+Nx,density,args_syst,'delta_3')*\
						np.exp(1j*kya*deltas[2][1])
				H[n+Nx,n] = np.conj(H[n,n+Nx])

			if n<Nx*(Ny-1) and n%2!=0:  # n is B site, this B is not on the top boundary
				H[n,n+Nx] = corr_hopping(n+Nx,n,density,args_syst,'delta_2')*\
						np.exp(-1j*kya*deltas[1][1])+corr_hopping(n+Nx,n,density,args_syst,'delta_3')*\
						np.exp(-1j*kya*deltas[2][1])
				H[n+Nx,n] = np.conj(H[n,n+Nx])


	if args_syst['BC']=='Xper': # ARMCHAIR
		H = 1j*np.zeros((Nx*Ny,Nx*Ny))
		kxa = ka

		for n in range(0,Nx*Ny,2):

			if n%(2*Nx)==0:
				# n = A site to n+1 = B site
				H[n,n+1] = honey_hopping(n,n+1,sites_dic,system,'delta_1')\
							*np.exp(+1j*kxa*3) # link t1 
				H[n+1,n] = np.conj(H[n,n+1])

				if n<Nx*(Ny-1):
					# n = A site to n+Nx = B site
					H[n,n+Nx] = honey_hopping(n,n+Nx,sites_dic,system,'delta_2') # link t2
					H[n+Nx,n] = np.conj(H[n,n+Nx])
					
					# n+1 = B site to n+1+Nx = A site
					H[n+1+Nx,n+1] = honey_hopping(n+1+Nx,n+1,sites_dic,system,'delta_3') # link t3
					H[n+1,n+1+Nx] = np.conj(H[n+1+Nx,n+1])


			elif n%(2*Nx)!=0:
				# n = B site to n+1 = A site
				H[n+1,n] = honey_hopping(n+1,n,sites_dic,system,'delta_1')
				H[n,n+1] = np.conj(H[n+1,n])

				if n<Nx*(Ny-1):
					# n = B site to n+Nx = A site
					H[n+Nx,n] = honey_hopping(n+Nx,n,sites_dic,system,'delta_3') # link t3
					H[n,n+Nx] = np.conj(H[n+Nx,n])

					# n+1 = A site to n+1+Nx = B site
					H[n+1,n+1+Nx] = honey_hopping(n+1,n+1+Nx,sites_dic,system,'delta_2') # link t2
					H[n+1+Nx,n+1] = np.conj(H[n+1,n+1+Nx])

	return H


def extract_coord_from_dict(dictio):

	L = len(dictio)
	x_s = np.zeros(L)
	y_s = np.zeros(L) # revlevant for 2D systems

	for n in range(L):

		if len(dictio[n])==1:
			x_s[n] = dictio[n]
		else:
			x_s[n],y_s[n] = dictio[n]

	return x_s,y_s

def trap(args_syst):

	V0 = args_syst['V']
	sites_dic = args_syst['sites_dic']
	Nx = args_syst['Nx']
	x_max = np.max(extract_coord_from_dict(sites_dic)[0])

	x_center = x_max/2
	V = np.zeros(len(sites_dic)) 

	for n in range(len(sites_dic)):
		x = sites_dic[n][0]
		V[n] = 1/2*V0*(x-x_center)**2

	return V

def remove_under0(n_TF):
	
	for i in range(len(n_TF)):
		if n_TF[i]<0:
			n_TF[i] = 0

	return n_TF

def H_KV_honey_BEC(args_syst,ka):

	d1 = np.array([-1,0])
	d2 = np.array([0.5,np.sqrt(3)/2])
	d3 = np.array([0.5,-np.sqrt(3)/2])
	deltas = [d1,d2,d3]
	J = args_syst['J']
	sites_dic = args_syst['sites_dic']
	Nx = args_syst['Nx']
	Ny = args_syst['Ny']	

	## without strain for bosons
	if args_syst['BC']=='Yper': #or system[0]=='Ypert2t3':

		kya = ka

		H = np.diag(trap(args_syst))*(1+0*1j)
		
		for n in range(len(sites_dic)):
			if n%Nx!=0 and n%2==0 and n>0: # this A isn't on the left boundry
				H[n,n-1] = -J #corr_hopping(n,n-1,sites_dic,system,'delta_1')
				H[n-1,n] = np.conj(H[n,n-1])

			if n<Nx*(Ny-1) and n%2==0:  # this A is not on the top boundary
				#H[n,n+Nx] = corr_hopping(n,n+Nx,sites_dic,system,'delta_2')*\
							#(1+np.exp(1j*kya*2*deltas[1][1]))	
				H[n,n+Nx] = -J*(1+np.exp(1j*kya*2*deltas[1][1]))	
				H[n+Nx,n] = np.conj(H[n,n+Nx])

			if n<Nx*(Ny-1) and n%2!=0:  # this B is not on the top boundary
				#H[n,n+Nx] = corr_hopping(n,n+Nx,sites_dic,system,'delta_3')*\
							#(1+np.exp(1j*kya*2*deltas[2][1]))
				H[n,n+Nx] = -J*(1+np.exp(1j*kya*2*deltas[2][1]))
				H[n+Nx,n] = np.conj(H[n,n+Nx])

	return H


def reorder(x_s,psi): 
	'''
	The "x" coord. are reordered from smallest to greatest and to avoid 
	to lose the correspondant psi, we conjointly keep the couples (x,psi(x))
	intact.
	'''
	
	Len = len(x_s)
	indices_values = []
	for i,value in enumerate(x_s):
		index_value = [i,value]
		indices_values.append(index_value)

	copy = np.array(indices_values.copy())
	temp = np.zeros((Len,2))
	out_x = x_s.copy()
	out_psi = psi.copy()
	for i in range(Len):
		index_min = int(np.argmin(copy[:,1]))
		temp[i] = indices_values[index_min]
		index_max = int(np.argmax(copy[:,1]))
		copy[index_min,1] = copy[index_max][1]+1

	temp = np.array(temp) ## to be able to use int()
	for i in range(Len):
		out_x[i] = x_s[int(temp[i,0])]
		out_psi[i] = psi[int(temp[i,0])]

	return out_x,out_psi


def compare2TF(args_syst,n0,n_TF):

	U = args_syst['U']

	x_s = extract_coord_from_dict(args_syst['sites_dic'])[0]
	x_s,n0 = reorder(x_s,n0)
	x_s = extract_coord_from_dict(args_syst['sites_dic'])[0]
	x_s,n_TF = reorder(x_s,n_TF)

	pyplot.suptitle('OL Honey Yper %s, %s, %s,\
		 Nx = %.1i, N = %.3e, J = %.2f, \n V = %.3e, U = %.3e' % \
		 (args_syst['Method'],args_syst['Trap'],args_syst['Symm'],\
		args_syst['Nx'],args_syst['N'],args_syst['J'],args_syst['V'],\
		args_syst['U']))

	pyplot.plot(x_s,n_TF,'--',label="$n_{TF}$")
	pyplot.plot(x_s,n0,'-',label="$n_0, U = {}$".format(U),fillstyle='none')
	pyplot.xlabel('x')
	pyplot.grid(axis='both')
	pyplot.legend();

	return 0

def get_eigVect(filename_ferm,filename_bos,my_ky,my_eigVal):

	data_ferm = np.load(filename_ferm+'.npy',allow_pickle=True)
	args_syst = data_ferm[0]
	eigVal_all = data_ferm[1]
	eigVal_TF_all = data_ferm[2]
	kya_list = data_ferm[3]

	J = args_syst['J']
	Nx = args_syst['Nx']
	V = args_syst['V']
	U = args_syst['U']
	N = args_syst['N']
	alpha = args_syst['alpha']
	sites_dic = args_syst['sites_dic']

	data_bos = np.load(filename_bos+'.npy',allow_pickle=True)
	mu = data_bos[2]+3*J
	psi0 = data_bos[3]
	n0 = np.abs(psi0)**2
	n_TF = (mu-trap(args_syst))/U/N
	n_TF = remove_under0(n_TF) # negative values are put to 0

	ky,index_ky = GP.nearest_value(kya_list,my_ky)

	H_c = H_honey_c(n0*N,args_syst,my_ky)
	#H_c = H_honey_c(n_TF*N,args_syst,my_ky)

	eigVal,eigVect = linalg.eigh(H_c)

	E,index_eigVal = GP.nearest_value(eigVal,my_eigVal)
	my_eigVect = eigVect[:,index_eigVal] # eigenvector computed for specified ky and energy level

	#p = 1
	#flag = 0
	# while flag==0:
	# 	for j in range(len(eigVect_TF[0])):
	# 		x_s = OHL.extract_coord_from_dict(args_syst['sites_dic'])[0]
	# 		x_s,psi = OHL.reorder(x_s,eigVect_TF[:,j])

	# 		if np.abs(GP.mean_x(psi,x_s))>p*np.max(x_s) or\
	# 				np.abs(GP.mean_x(psi,x_s))<(1-p)*np.max(x_s):
	# 			ax2.plot(x_s,np.abs(psi),label='number = {}'.format(j))
	# 			#ax1.plot(kya_list,eigVal_TF_all[:,j],'-',label='number = {}'.format(j))
	# 			#print(GP.mean_x(psi,x_s))
	# 			flag = 1
	# 	#print(p)
	# 	#print((1-p)*np.max(x_s))
	# 	p -= 0.001

	return my_eigVect,E,ky	

def LL_eigVect(x_list,N,x0,alpha):
	out = np.zeros(len(x_list))

	N = abs(N) # N could be negative 
	x = sy.symbols('x')
	eigVect = sy.exp(-((x-x0)/alpha)**2/2)*(-1)**N*sy.exp(((x-x0)/alpha)**2)*\
				sy.diff(sy.exp(-((x-x0)/alpha)**2),x,N)
	eigVect = lambdify(x,eigVect,"numpy")

	i = 0
	for x in x_list:
		out[i] = eigVect(x)
		i += 1
	return out/np.sqrt(np.sum(np.abs(out)**2))

def normalize(vector):
	return vector/np.sqrt(np.sum(np.abs(vector)**2))

def relat2max(vector):
	return vector/np.max(np.abs(vector))

def color_position(positions,phi,R): 

	''' Picks a color in function of the mean position of the distribution with respect to
	the center.	 
	'''
	x_c = np.max(positions)/2
	moy_phi = np.sum((positions-x_c)*np.abs(phi)**2/np.sqrt(np.sum(np.abs(phi)**2))) # mean distance from x_c

	## moy_phi goes from -x_c to x_c

	color = moy_phi/x_c 
	## gives a number between 0 and 1 that will correspond to a color (not with RGB)
	## but with a colormap 

	return color

############################################
## Test honey_lattice:
# sites_dic = honey_lattice({'BC' : 'Yper', 'Nx' : 10, 'Ny' : 2})

# for n in range(len(sites_dic)):
# 	[x,y] = sites_dic[n]

# 	pyplot.plot(x,y,'hb')
# pyplot.axis('scaled')
# pyplot.show()
# exit()

## Test (reorder + extract_coord_from_dict):
# args_syst = {
# 			'J' : 1,
# 			'N' : 1e4,
# 			'V' : 1e-4,
# 			'Nx' : 3,
# 			'Ny' : 2,
# 			'U' : 0,
# 			'Method' : 'SC',
# 			'Trap' : 'Harmonic',
# 			'Symm' : 'Isotropic',
# 			'BC' : 'Yper'
# 			}
# args_syst.update({'sites_dic' : honey_lattice(args_syst)})

# sites_dic = args_syst['sites_dic']
# x,y = extract_coord_from_dict(sites_dic)
# print(x,y)

# print(reorder(x))
# exit()


############################################
'''
sites_dic = args_syst['sites_dic']
x_s = extract_coord_from_dict(sites_dic)[0]
x_s,psi0 = reorder(x_s,psi0)

#print(H_KV)
#pyplot.pcolormesh(abs(H_KV))
#pyplot.colorbar()
pyplot.plot(x_s,abs(psi0)**2)
pyplot.show()
'''


