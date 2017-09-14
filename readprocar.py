#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       readprocar.py                                                #
#  Dependence: arguments.py,parse.py,bzkpt.py,xplot.py                      #
#  Usage:      read the PROCAR file, plot weighted bands                    #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       May 22, 2017                                                 #
#                                                                           #
#===========================================================================#


import os
import numpy as np
import argparse
import arguments
from crystal_structure import *
import bzkpt
import parse
import write
import xplot

def get_bandgap(nspin,energy,nelect):
    if energy.shape[2]==1 and nspin==2:
       e1=energy[:,0::2,:]
       e2=energy[:,1::2,:]
    energy=np.vstack([e1,e2])
    val_idx = nelect/2-1
    con_idx = val_idx+1
    VBM = np.zeros((nspin),float)
    CBM = np.zeros((nspin),float)
    for ispin in range(nspin):
       VBM[ispin] = max(energy[:,val_idx,ispin])
       CBM[ispin] = min(energy[:,con_idx,ispin])
       print "spin channel {0:1d}".format(ispin)
       print "VBM is {0:10.6f} eV, CBM is {1:10.6f} eV".format(VBM[ispin],CBM[ispin])
       if CBM[ispin] - VBM[ispin] >=0:
          print "band gap is {0:10.6f} eV".format(CBM[ispin]-VBM[ispin])
       else:
          print "metallic!"

desc_str = '''Simple program used to quickly plot the orbital weights
of the band structure contained in the specified PROCAR file.
'''

parser = argparse.ArgumentParser(prog='plotprocar', description = desc_str)

arguments.add_io_arguments(parser)
arguments.add_fig_arguments(parser)
arguments.add_plot_arguments(parser)
                    
args = parser.parse_args()

struct = parse.parse_poscar(args.poscar)
orbitals, kpt, kweights, energy, occupancies, weights, phases = parse.parse_procar(args.procar)
kpt=np.matrix(kpt)
kpt=kpt*struct._reciprocal_cell()

nspin=int(os.popen('grep ISPIN OUTCAR').readline().split()[2])
if os.popen('grep LNONCOLLINEAR OUTCAR').readline().split()[2]=='T':
   print 'noncollinear spin polarized calculations'
   nspin=2
nelect=float(os.popen('grep NELECT OUTCAR').readline().split()[2])
nelect=int(nelect)
print 'energy.shape:',energy.shape
print 'weights.shape:',weights.shape
#get_bandgap(nspin,energy,nelect)

if args.efermi==0:
   try:
      line=os.popen("grep fermi OUTCAR|tail -1").readline()
      efermi=float(line.split(":")[1].split()[0])
   except:
      print "Note: The fermi energy is 0!"
      efermi=args.efermi
else:
   efermi=args.efermi
energy = energy - efermi
path=bzkpt.get_path(kpt)
xsym=bzkpt.guess_xsym(path)
xplot.plot_weighted_band(path,energy,weights,struct,orbitals,xsym,args)
write.write_weighted_band(path,energy,weights,struct,orbitals,xsym,args)
