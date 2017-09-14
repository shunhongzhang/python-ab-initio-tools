#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       plot_kpdos.py                                                #
#  Dependence: none                                                         #
#  Usage:      plot k-resolved-pdos from input data                         #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Dec 15, 2016                                                 #
#                                                                           #
#===========================================================================#

#under development
#import xplot
import os
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import numpy as np
import argparse
import arguments
import parse
import bzkpt


desc_str = '''plot the k-resolved DOS of crystals.'''

parser = argparse.ArgumentParser(prog='plot_kpdos.py', description = desc_str)
arguments.add_io_arguments(parser)
arguments.add_fig_arguments(parser)
arguments.add_plot_arguments(parser)

parser.add_argument('--sigma',type=float,default=0.001,help='smearing width of the ldos')
parser.add_argument('--nedos',type=float,default= 5000,help='number of grid of the ldos')
parser.add_argument('--contrast',type=float,default=0.07,help='magnify the lods by multiplying the factor')
parser.add_argument('--cmap',type=str,default='hot',help='color map for the plot')
parser.add_argument('--nnkpt',type=int,default=1,help='number of interpolated kpts on the path')

args=parser.parse_args()


if args.efermi==0:
   try:
      line=os.popen("grep fermi OUTCAR|tail -1").readline()
      efermi=float(line.split()[2])
      print 'fermi energy is',efermi
   except:
      print "Note: The fermi energy is 0!"
      efermi=args.efermi


prefix='feo'
orb_type={0:'s',1:'p',2:'d'}
efermi=0
filename=prefix+'pdos_atm#2*_wfc#'+str(irob)+'('+orb_type[iorb]+')'
nedos=int(os.popen('wc -l '+filename).readline())-1
f=open(filename,'r')
print f.file()
headline=f.readline().split()
ncolumn=int(len(headline))-2
nspin=(ncolumn-1)/2
data = np.fromfile(f,sep=' ',dtype=float).reshape[nedos,ncolumn]


kpath=np.zeros(nedos)
dos=np.zeros(nedos,nspin)
pdos=np.zeros(nedos,nspin)
kpath=data[:,0]
for ispin in range(nspin):
    dos[:,ispin]=data[:,ispin+1]
    pdos[:,ispin]=data[:,ispin+2,]



fig, ax = plot.subplots(figsize=args.figsize)
cax = ax.imshow(ldos.T, extent=(kpath[0], kpath[-1], ymin, ymax), cmap=plot.get_cmap(args.cmap),aspect='auto')
fig.tight_layout()
cbar = fig.colorbar(cax,ticks=[np.amin(ldos), (np.amin(ldos)+np.amax(ldos))/2, np.amax(ldos)], orientation='horizontal')
cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar

plot.show()


