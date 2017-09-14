#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       plot_kpdos.py                                                #
#  Dependence: none                                                         #
#  Usage:      plot k-resolved-pdos from input data                         #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Apr 16, 2017                                                 #
#                                                                           #
#===========================================================================#


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

def mesh_energy(args,eigenval):
    nkpt=eigenval.shape[0]
    if args.elim[1]>args.elim[0]:
       emin,emax=float(args.elim[0]),float(args.elim[1])
    else:
       emin, emax=np.amin(eigenval),np.amax(eigenval)
    deltaE=(emax-emin)/(args.nedos-1)
    print 'nedos= {0:5d},deltaE={1:12.9f}'.format(args.nedos,deltaE)
    energy=np.zeros((args.nedos),float)
    for ie in range(args.nedos):
        energy[ie]=emax-ie*deltaE
        #energy[ie]=emin+ie*deltaE
    return energy

def gaussian(x,miu,sigma2):
    return np.exp(-(x-miu)**2/sigma2/2)/np.sqrt(2*np.pi*sigma2)

def interpolate_ldos(nnkpt,energy,eigenval,weight,sigma,contrast=0):
    sigma2=sigma**2
    nedos=energy.shape[0]
    nkpt=eigenval.shape[0]
    ldos=np.zeros((nkpt,nedos),float)
    nband=eigenval.shape[1]
    print 'interpolating...'
    print 'ldos mesh grid: {0:5d} k-pionts x {1:6d} energy points'.format(nkpt,nedos)
    print 'interpolated from ldos of {0:5d} k-points x {1:5d} bands'.format(nkpt,nband)
    for ikpt in range(nkpt):
        for ie in range(nedos):
            ldos[ikpt*nnkpt,ie]=sum(weight[ikpt,:]*gaussian(energy[ie],eigenval[ikpt,:],sigma2))
        if np.mod(ikpt+1,nkpt/10)==0:
           print '{0:3d} of {1:5d} points finished'.format(ikpt+1,nkpt)
    print 'done'
    if contrast:
       ldos=ldos**contrast
    return ldos


if args.efermi==0:
   try:
      line=os.popen("grep fermi OUTCAR|tail -1").readline()
      efermi=float(line.split()[2])
      print 'fermi energy is',efermi
   except:
      print "Note: The fermi energy is 0!"
      efermi=args.efermi

struct = parse.parse_poscar(args.poscar)
orbitals, kpt, kweights, eigenval, occupancies, weights, phases = parse.parse_procar(args.procar,efermi=efermi)
kpt=np.matrix(kpt)
kpt=kpt*struct._reciprocal_cell()


'''
efermi=0
filename='band.dat'
v = np.fromfile(filename,sep=' ')
ncol=4
v = v.reshape(v.shape[0]/ncol,ncol)
nkpt=len(set(v[:,0]))
nrow=v.shape[0]/nkpt
print " {0:8d} bands, {1:8d} k-points".format(nrow,nkpt)

kpt = v[:,0].reshape(nkpt,nrow)
eigenval = v[:,1].reshape(nkpt,nrow)-efermi
weight = v[:,2].reshape(nkpt,nrow)

print "mesh grid:",nkpt,nedos

#For data from projwfc.x of Quantum ESPRESSO
#kpt = kpt.T[:][::-1]
#eigenval = eigenval.T[:][::-1]
#weight = weight.T[:][::-1]

kpath=np.zeros((nkpt,nedos),float)
for ikpt in range(nkpt):
    kpath[ikpt,:]=ikpt
'''



kpath = bzkpt.get_path(kpt)
energy=mesh_energy(args,eigenval[:,:,0])
if args.elim[1]>args.elim[0]:
   ymin,ymax=float(args.elim[0]),float(args.elim[1])
else:
   ymin,ymax=np.amin(energy),np.amax(energy)
xsym = bzkpt.guess_xsym(kpath)

norb=9
nkpt=kpt.shape[0]
nband=eigenval.shape[1]
res_weight=np.zeros((nkpt,nband),float)
for iat in range(2,8):
    res_weight[:,:]=res_weight[:,:] + weights[:,iat,norb,:,0,0]

'''
for iband in range(nband):
    plt.plot(kpath,eigenval[:,iband])
plt.show()
'''
ldos=interpolate_ldos(args.nnkpt,energy,eigenval[:,:,0],res_weight,args.sigma,contrast=args.contrast)


fig, ax = plot.subplots(figsize=args.figsize)
cax = ax.imshow(ldos.T, extent=(kpath[0], kpath[-1], ymin, ymax), cmap=plot.get_cmap(args.cmap),aspect='auto')
fig.tight_layout()
cbar = fig.colorbar(cax,ticks=[np.amin(ldos), (np.amin(ldos)+np.amax(ldos))/2, np.amax(ldos)], orientation='horizontal')
cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])  # horizontal colorbar

plot.show()


#if args.label_k:
#   xplot.plot_high_sym(xsym,args.label_k,ymin,ymax)
