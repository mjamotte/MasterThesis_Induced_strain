import numpy as np
#import sympy as sp
#from sympy import *
#from sympy import Matrix
import matplotlib.pyplot as pyplot
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
import lib_Honeycomb_Lattice as FHoney


##################################### TEST HONEY #####################################
rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

Nx = 0
Ny = 0
t = 1
tau = 1
sites_dic = []

My = 400
Mx = 400
kya_list = np.linspace(-5*np.pi/3/np.sqrt(3),5*np.pi/3/np.sqrt(3),My)
kxa_list = np.linspace(-3*np.pi/3,3*np.pi/3,Mx)

system = ['FullPeriodic',Nx,Ny,t,tau]

eigValues = 1j*np.zeros((len(kya_list),len(kxa_list),2))

i = 0
j = 0
for kxa in kxa_list:
	j = 0
	for kya in kya_list:
		H = FHoney.H_honey(system,kxa,kya,sites_dic)
		eigValues[j,i] = linalg.eigh(H)[0]
		j += 1
	i += 1

if (len(kya_list)>1 and len(kxa_list)>1):
    [Kxa,Kya] = np.meshgrid(kxa_list,kya_list)
    figE = pyplot.figure(figsize=(6,5))
    ax = pyplot.axes(projection='3d')
    
    for l in range(2):
        ax.plot_surface(Kxa, Kya, eigValues[:,:,l].real, rstride=1, cstride=1,
                        edgecolor='none',rasterized=True)
        ax.view_init(elev=10, azim=135)

    ax.set_xlabel('$k_x$',fontsize=25)
    ax.set_ylabel('$k_y$',fontsize=25)
    ax.set_zlabel('$E$',fontsize=25)

pyplot.show()

#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
temp = '%s' %(system[0])
#figE.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows

path_above_Maxime = (sys.path[0]).split('maxime')[0]
figE.savefig(path_above_Maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf")

## Zoom on a Dirac cone

My = 50
Mx = 50
kya_list = np.linspace(4*np.pi/3/np.sqrt(3)-0.3,4*np.pi/3/np.sqrt(3)+0.3,My)
kxa_list = np.linspace(-0.3,0.3,Mx)

system = ['FullPeriodic',Nx,Ny,t,tau]

eigValues = 1j*np.zeros((len(kya_list),len(kxa_list),2))

i = 0
j = 0
for kxa in kxa_list:
	j = 0
	for kya in kya_list:
		H = FHoney.H_honey(system,kxa,kya,sites_dic)
		eigValues[j,i] = linalg.eigh(H)[0]
		j += 1
	i += 1

if (len(kya_list)>1 and len(kxa_list)>1):
    [Kxa,Kya] = np.meshgrid(kxa_list,kya_list)
    figE = pyplot.figure(figsize=(5,5))
    ax = pyplot.axes(projection='3d')
    
    for l in range(2):
        ax.plot_surface(Kxa, Kya, eigValues[:,:,l].real, rstride=1, cstride=1,
                    edgecolor='none',rasterized=True)
        ax.view_init(elev=10, azim=45)

    ax.set_xlabel('$k_x$',fontsize=25)
    ax.set_ylabel('$k_y$',fontsize=25)
    ax.set_zlabel('$E$',fontsize=25)
    pyplot.axis('off')

pyplot.show()

#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
temp = '%sZoom' %(system[0])
#figE.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows

path_above_Maxime = (sys.path[0]).split('maxime')[0]
figE.savefig(path_above_Maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf")