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
from PIL import Image, ImageDraw

import sys
path_above_Codes = (sys.path[0]).split('Codes')[0]
sys.path.append(path_above_Codes+'Codes/')
import lib_GrossPitaevskii as GP
import lib_Manage_files as Mf
import lib_Optical_Honeycomb_Lattice as OHL

rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})
rc('text', usetex=True)

args_syst = {
	'Nx':3,
	'Ny':5,
	'BC':'Open'
}

params = {
	'fontsize':25,
	'linewidth':3
}

drawAB = True
put_delta = False
Brillouin = False
put_t = False

fig_Honey = pyplot.figure(figsize=(args_syst['Nx']+4,args_syst['Ny']+2))

ax = OHL.draw_Honey(args_syst,drawAB,put_delta,put_t,params,Brillouin)

pyplot.axis('off') #scaled

pyplot.show()

#path_above_Maxime = (sys.path[0]).split('Maxime')[0] # Windows
temp = '%s' %(args_syst['BC'])
#figHoney.savefig(path_above_Maxime+'\\Maxime\\ownCloud\\MasterThesis\\figures\\fig_Honey_'+temp+".pdf") # Windows
#path_above_Maxime = (sys.path[0]).split('maxime')[0]
#fig_Honey.savefig(path_above_Maxime+'/maxime/ownCloud/MasterThesis/figures/figHoney'+temp+".pdf")