#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       readband.py                                                  #
#  Dependence: arguments.py,parse.py,bzkpt.py,write.py, xplot.py            #
#  Usage:      convert the EIGENVAL file to a file for band structure plot  #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Jul 14, 2017                                                 #
#                                                                           #
#===========================================================================#

import os
import math
import numpy as np
import sys
import argparse
import arguments
import bzkpt
import parse
import write
import xplot


desc_str = '''plot the band structure of crystals.'''

parser = argparse.ArgumentParser(prog='plotband.py', description = desc_str)
arguments.add_io_arguments(parser)
arguments.add_fig_arguments(parser)
arguments.add_plot_arguments(parser)
args = parser.parse_args()
#argument.check_args()


def get_bandgap(nspin,energy,nelect):
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

nspin,nelect,efermi,kpt,kpt_cart,kweight,energy = parse.parse_eigenval(efermi=args.efermi)
get_bandgap(nspin,energy,nelect)


segment_nkpt = int(open("KPOINTS",'r').readlines()[1])
if segment_nkpt:
   get_mode=open("KPOINTS",'r').readlines()[2][0]
   if get_mode=="L" or get_mode=="l":
      print "line mode for band structure calculations"
      kp_mode='line'
      path = bzkpt.get_path(kpt_cart)
      xsym = bzkpt.guess_xsym(path)
      print 'high sym on kpath:',','.join(["{0:10.6f} ".format(ixsym) for ixsym in xsym])
      #xsym = bzkpt.get_xsym(path,segment_nkpt)
      #print 'high sym on kpath:',','.join(["{0:10.6f} ".format(ixsym) for ixsym in xsym])
      write.write_band(path,segment_nkpt,energy,nelect)
      write.write_specified_band(kpt,kweight,path,energy,nelect)
   else:
      print "K points specified explicitly in the KPINTS file"
      kp_mode="specified"
      write.write_specified_band(kpt_cart,kweight,path,energy,nelect)
else:
   print "Monkhort Pack k points grid"
   kp_mode='mp'
   grid=[int(item) for item in open("KPOINTS","r").readlines()[3].split()]
   #write.write_specified_band(kpt_cart,energy,nspin,nelect)
   write.write_mesh_band(kpt_cart,energy,nspin,nelect,grid,isoenergy=-0.7)


if kp_mode=='line':
   struct=parse.parse_poscar()
   ibrav,brav,center=struct._get_ibrav()
   label_k = bzkpt.get_label_k(4,kpt)
   print label_k
   struct._is_slab()
   xplot.plot_band(args,path,energy,xsym)
elif kp_mode == 'specified':
   for ikpt in range(nkpt):
       if kweight[ikpt]!=0:
          ikpt+=1
       else:
          start_kpt=ikpt
          break
   path = bzkpt.get_path(kpt_cart)
   xsym = bzkpt.guess_xsym(path)
   kpt=kpt[start_kpt:,:]
   energy=energy[start_kpt:,:,:]
   xplot.plot_band(args,path,energy,xsym)

