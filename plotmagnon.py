#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       plotmagnon.py                                                #
#  Dependence: none                                                         #
#  Usage:      calculate and plot magnon dispersion                         #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Dec 15, 2016                                                 #
#                                                                           #
#===========================================================================#



import os
import numpy as np
from math import *
import argparse
import parse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches


desc_str = '''calculate and plot the magnon energy dispersion vs. k vectors.'''

parser = argparse.ArgumentParser(prog='plotmagnon', description = desc_str)
parser.add_argument('--output', type=str, default="Magnon_dispersion",
                    help='Filename for the plot. It accepts '
                    'all image formats supported by matplotlib.')
parser.add_argument('--figsize', type=eval, default=(6, 4),
                    help='Figure size in inches formatted as width,height '
                    '(no spaces allowed). Default is 6,4')
parser.add_argument('--dpi', type=float, default=300,
                    help='Plot resolution in dots per inch. Default is 300.')
parser.add_argument('--marker', type=str, default='o', 
                    help='Marker for the fatband plot. Default is o.')
parser.add_argument('--markersize', type=float, default=20,
                    help='Marker size. Default is 20.')
parser.add_argument('--color', type=str, default="blue",
                    help='Color for the marker. It accepts any color specification '
                    'accepted by matplotlib. Color specified as r,g,b tuple should ' 
                    'not have any spaces.')
parser.add_argument('--elim', type=eval, default=(1,-1),
                    help='Energy range for the plot, specified as emin,emax '
                    '(no spaces allowed). Default is entire band range.')
parser.add_argument('--label_k', type=str, default='$\\barX(\\barY)$ $\Gamma$ $X(Y)$',
		    help='labels for high symmetry k points'
		    'Use the form like \Gamma to represent Greek characters'
		    "example: --label_k 'X \Gamma Y'")
                    
args = parser.parse_args()

# Handle color arguments
def color(string):
    try:
        # Try to convert to tuple
        return eval(string)
    except:
        # If cannot convert to tuple just
        # return the original string
        return string

def distance(k1,k2):
    return sqrt((k1[0]-k2[0])**2+(k1[1]-k2[1])**2)

def Magnon_E(J1,J2,J3,S,a,kx,ky):
    return 4*J1*S*abs(sin(kx*a/2)) - 4*J2*S*(1-cos(kx*a)) - 4*J3*S*(1-cos(ky*a/2))


# Parameters
# Exchange integrals, in unit of meV
J1 = 22.22
J2 = -5 
J3 = -1

# Spin per site
S=3/2

# Lattice constant
a=5.76
b=2*np.pi/a

print "lattice constant:",a,", reciprocal lattice constant:",b

# Number of kpoints per path
nkpt=400

# High symmetry k point path
K=[[-0.5,0],[0,0],[0.5,0]]

k_dist=0
kpath=[]
Mag_E=[]
color_patch=[]
for ipath in range(len(K)-1):
   kx_set = np.linspace (K[ipath][0]*b, K[ipath+1][0]*b, nkpt)
   ky_set = np.linspace (K[ipath][1]*b, K[ipath+1][1]*b, nkpt)
   path_length = distance(K[ipath],K[ipath+1])
   for kx,ky in zip(kx_set,ky_set):
       Mag_E.append(Magnon_E(J1,J2,J3,S,a,kx,ky))
   for ikp in range(nkpt):
       kpath.append(k_dist)
       k_dist=k_dist+path_length/(nkpt-1)

xsym = [kpath[ipath*nkpt] for ipath in range(len(K)-1)]
xsym.append(kpath[-1])
print "high symmetry k points in k-path:",xsym
   
nkp=len(Mag_E)
fig = plot.figure(figsize=args.figsize)
ax = fig.add_subplot(1,1,1)
ax.plot(kpath[0:nkp], Mag_E[0:nkp], color=args.color, ls='-',lw=1.5)
color_patch.append(mpatches.Patch(color=args.color, label='$\\barX-\Gamma-X$'))


e_low=min(Mag_E)
e_high=max(Mag_E)
label_height=e_low-(e_high-e_low)*0.1
#plot.title(args.output)

if args.label_k:
   for i in range(len(xsym)):
       label=args.label_k.split()[i]
       plot.text(xsym[i],label_height,label,ha='center')
   label=args.label_k.split()[-1]
   plot.text(xsym[-1],label_height,label,ha='center')

plot.ylabel("$\hbar\omega_k$ ($meV$)")

K=[[0,-0.5],[0,0],[0,0.5]]

k_dist=0
kpath=[]
Mag_E=[]
for ipath in range(len(K)-1):
   kx_set = np.linspace (K[ipath][0]*b, K[ipath+1][0]*b, nkpt)
   ky_set = np.linspace (K[ipath][1]*b, K[ipath+1][1]*b, nkpt)
   path_length = distance(K[ipath],K[ipath+1])
   for kx,ky in zip(kx_set,ky_set):
       Mag_E.append(Magnon_E(J1,J2,J3,S,a,kx,ky))
   for ikp in range(nkpt):
       kpath.append(k_dist)
       k_dist=k_dist+path_length/(nkpt-1)

xsym = [kpath[ipath*nkpt] for ipath in range(len(K)-1)]
xsym.append(kpath[-1])
print "high symmetry k points in k-path:",xsym

nkp=len(Mag_E)
#fig = plot.figure(figsize=args.figsize)
#ax = fig.add_subplot(1,1,1)
ax.plot(kpath[0:nkp], Mag_E[0:nkp], color='red', ls='-',lw=1.5)
color_patch.append(mpatches.Patch(color='red', label='$\\barY-\Gamma-Y$'))


plot.legend(bbox_to_anchor=(0.2,1.0),handles=[item for item in color_patch],loc=0)


# Set x-axis ticks
plot.xlim(0,kpath[-1])
plot.xticks(xsym, ['']*len(xsym))

# Save the figure
plot.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
