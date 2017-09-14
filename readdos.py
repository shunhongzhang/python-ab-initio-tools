#!/usr/bin/env python

#===========================================================================#
#                                                                           #
#  File:       readdos.py                                                   #
#  Dependence: util.py,parse.py,arguments.py,write.py,xplot.py              #
#  Usage:      read the DOSCAR and write files, plot the DOS                #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Apr 22, 2017                                                 #
#                                                                           #
#===========================================================================#


desc_str='''
In the VASP output DOSCAR, we have density of states (DOS):
1. total DOS, in order of
energy total_dos integrated_dos
2. projected DOS (PDOS) if the LORBIT tag is set as 11 or 12 in the INCAR file for the DOS calculation
energy s py pz px dxy dyz dz2 dzx dx2-dy2
In spin-polarized calculations the PDOS corresponds to each orbital will be splitted into to spin channles
energy s_up py_up pz_up px_up dxy_up dyz_up dz2_up dzx_up dx2-dy2_up s_dn py_dn pz_dn px_dn dxy_dn dyz_dn dz2_dn dzx_dn dx2-dy2_dn
'''

note="All indice start from 0!"

import os
from termcolor import colored
import numpy as np
import parse
import arguments
import write
import xplot
import argparse
import dos

parser = argparse.ArgumentParser(prog='readdos.py', description = desc_str)
arguments.add_io_arguments(parser)
arguments.add_fig_arguments(parser)
arguments.add_plot_arguments(parser)
args = parser.parse_args()

print colored(desc_str,'green')
print colored(note,'red')

try:
   struct = parse.parse_poscar('POSCAR')
except:
   print "POSCAR not found!Please input the following information manually:"
   species=raw_input("Please input the chemical species in the same order as that in POSCAR, separated by space:")
   counts=raw_input("please input the number of each kind of atoms, in the same order as that in POSCAR, separated by space:")


energy,emin,emax,dos_total,dos_integral,dos_orb=parse.parse_doscar(args.doscar)
dos_atom = np.sum(dos_orb,axis=2)

print "writing dat files..."
write.write_dos_total(energy,dos_total,dos_integral)
write.write_pdos_atom(energy,dos_atom,struct._symbols)
#write.write_pdos_species(energy,dos_species,dos_species_orb,struct._symbols,struct._counts)
print "done"

os.system("rm -r dat 2>/dev/null")
os.system("mkdir dat")
os.system("mv *.dat dat")

xplot.plot_total_dos(args,energy,dos_total,dos_integral)
ax = xplot.plot_dos(args,energy,dos_orb,struct)


'''
dos_obj=dos.electronic_dos()
at_list  = [int(projection.split('_')[0]) for projection in args.proj_index.split(',')]
for iat in at_list:
    print struct._symbols[iat],iat
    for iorb in range(0,9):
        print '{0:15s} {1:10.7f} {2:10.7f}'.format(xplot.orb_dic[iorb],dos_obj._integrated_dos(iat,iorb,0)[-1], dos_obj._integrated_dos(iat,iorb,1)[-1])
'''
