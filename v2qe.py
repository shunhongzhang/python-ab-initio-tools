#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       v2qe.py                                                      #
#  Dependence: parse.py,write.py (required)                                 #
#  Usage:      convert the POSCAR file to part of input file for PWscf(QE)  #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Aug 20, 2017                                                 #
#                                                                           #
#===========================================================================#

import numpy as np
from math import *
from termcolor import colored
import os
import parse
import write
import argparse
import crystal_structure

Note='''
This file can be used to generate input file for PWscf (Quantum ESPRESSO) by using VASP-POSCAR as input.
The definition of primitive cell basis vectors follows the Quantum ESPRESSO code, please refer to:
http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html#idm6425376
'''
Usage='''
Usage: Please prepare the POSCAR file in the conventional cell form, use direct (fractional) coordinates to indicate atomic positions.
Then run this script by type the command: python v2qe.py
'''
Alert="NOTE: for body/face/base-centered structures we can only process POSCAR of the CONVENTIONAL CELL!\n"


def get_connect(ibrav,spg,latt):
    connect_dic={
                  2 : np.matrix(([-1./2.,  0.0,1./2.],[   0.0,1./2.,1./2.],[-1./2., 1./2.,  0.0]),float),
                  3 : np.matrix(([1./2., 1./2.,1./2.],[-1./2.,1./2.,1./2.],[-1./2.,-1./2.,1./2.]),float),
                  7 : np.matrix(([1./2.,-1./2.,1./2.],[ 1./2.,1./2.,1./2.],[-1./2.,-1./2.,1./2.]),float),
                  9 : np.matrix(([1./2., 1./2.,  0.0],[-1./2.,1./2.,  0.0],[   0.0,   0.0,  1.0]),float),
                 -9 : np.matrix(([1./2.,-1./2.,  0.0],[ 1./2.,1./2.,1./2.],[   0.0,   0.0,  1.0]),float),
                 10 : np.matrix(([1./2.,   0.0,1./2.],[ 1./2.,1./2.,  0.0],[   0.0, 1./2.,1./2.]),float),
                 11 : np.matrix(([1./2., 1./2.,1./2.],[-1./2.,1./2.,1./2.],[-1./2.,-1./2.,1./2.]),float),
                 13 : np.matrix(([1./2,    0.0,-1./2],[ 0.0,    1.0,  0.0],[ 1./2.,  0.0,  1./2]),float)}
    if abs(ibrav)==5:
       a=latt['a'];c=latt['c']
       R_cosalpha=(2*c**2-3*a**2)/(2*c**2+6*a**2)
       print "cosalpha(celldm4)=",R_cosalpha
       tx=sqrt((1-R_cosalpha)/2)
       ty=sqrt((1-R_cosalpha)/6)
       tz=sqrt((1+2*R_cosalpha)/3)
       A = np.matrix(([[1,0,0],[-1./2.,sqrt(3)/2,0],[0,0,c/a]]),float)   #basis vectors of the hexagonal representation
       if ibrav==5 or ibrav==-5:
          B = 1/sqrt(2-2*R_cosalpha)*np.matrix(([[tx,-ty,tz],[0,2*ty,tz],[-tx,-ty,tz]]),float)
       #elif ibrav==-5:
       #   u = tz + 2*sqrt(2)*ty
       #   v = tz - sqrt(2)*ty
       #   B = 1/sqrt(6-6*R_cosalpha)*np.matrix(([[u,v,v],[v,u,v],[v,v,u]]),float)
       #print np.linalg.det(B*A**-1)
       connect_dic.setdefault(ibrav,B*np.linalg.inv(A))
    else:
       connect_dic.setdefault(ibrav,np.matrix(([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]),float))

    connect=connect_dic[ibrav]
    if np.linalg.det(connect)==0:
      print "The connect np.matrix is:\n",connect
      exit("The connect np.matrix is singular! Check your structure carefully!")
    return connect

def standardize_poscar(ibrav,spg,poscar='POSCAR'):
    struct = parse.parse_poscar(poscar)
    if abs(ibrav)==9 and spg.split()[0]=="A":
       sc=np.matrix(([[0,1,0],[0,0,1],[1,0,0]]),float)
    elif ibrav==13 and latt['beta']!=90:
       sc=np.matrix(([[1,0,0],[0,0,1],[0,-1,0]]),float)
    else:
       sc=np.matrix(([1,0,0],[0,1,0],[0,0,1]),float)
    struct_std = struct.build_supercell(sc)
    write.write_poscar_head(1,struct_std._cell,system='standardized poscar',filename='POSCAR_standardized')
    write.write_poscar_atoms(struct_std,filename='POSCAR_standardized')


if __name__=='__main__':
   parser = argparse.ArgumentParser(prog='v2qe.py', description = Note)
   parser.add_argument('--poscar',type=str,default='POSCAR',help='name of the POSCAR file')
   parser.add_argument('--symmprec', type=float, default=5e-4, help='deviation tolerance for finding crystal symmetry, in angstrom')
   parser.add_argument('--ecutwfc', type=float, default=100, help='plane wave cutoff')
   parser.add_argument('--ecutrho', type=float, default=500, help='charge density cutoff')
   parser.add_argument('--prefix', type=str, default="pw", help='prefix for the pw calculation')
   parser.add_argument('--pseudo_dir',type=str,default="'/home/hydrogen/pseudo/pbe'",help="directory for pseudopotential files")
   parser.add_argument('--upf',    type=str, default=".pbe-mt_fhi.UPF",help="type of pseudopotentail")
   parser.add_argument('--kmesh', type=str, default="5 5 5", help='k point mesh using the Monkhorst Pack scheme')
   parser.add_argument('--kshift', type=str, default="0 0 0", help='k point mesh shift from the Gamma point')
   args=   parser.parse_args()
   print colored(Note,'blue')
   print colored(Usage,'blue')
   print colored(Alert,"red")

   struct = parse.parse_poscar(args.poscar)
   spg, spg_no = struct._find_symmetry()
   ibrav,brav,center = struct._get_ibrav()
   print "space group:",spg,", No.",spg_no,", Lattice type:",brav,center

   print colored(crystal_structure.def_ibrav,'green')
   print colored("ibrav=","red"),ibrav
   print colored("Please be cautious with the space group and ibrav found by this code if you are dealing with the following systems",'red')
   print colored("Low dimensional materials: The periodicity in the vacuum direction(s) are inrealistic, so the 'spacegroup' may be wrong.","green")
   print colored("Magnetic materials: The spin polarization may adds extra properties to the atoms, so the magnetic unit cell may differs from the chemical primitive cell.","green")
   print colored("Base-centered structures: The choice of the base plane is alternative, so please check the generated structure carefully using xcrysden before calculation.","green")
   ldef=raw_input(colored("Do you want to define ibrav manually?(y/n)","red"))
   if ldef=="y":
      ibrav=input(colored("ibrav = ","red"))

   connect=get_connect(ibrav,spg,struct._latt_param())
   print "connect\n",connect
   standardize_poscar(ibrav,spg)
   struct_std = parse.parse_poscar('POSCAR_standardized')
   if ldef=='y':
      struct_std=struct
   celldm=struct_std._find_celldm(ibrav=ibrav)
   print "Lattice Parameters of standardized POSCAR:"
   print "a=",struct_std._latt_param()['a'],"angstrom, b=",struct_std._latt_param()['b'],"angstrom, c=",struct_std._latt_param()['c'],"angstrom"
   print "alpha=",struct_std._latt_param()['alpha'],"beta=",struct_std._latt_param()['beta'],"gamma=",struct_std._latt_param()['gamma']
   struct_pm = struct_std.build_supercell(connect)
   write.write_poscar_head(1,struct_pm._cell,system='primitive cell',filename="POSCAR_Primitive")
   write.write_poscar_atoms(struct_pm,filename="POSCAR_Primitive")
   setup_dic={"prefix" : args.prefix,
              "ecutwfc": args.ecutwfc,
              "ecutrho": args.ecutrho,
           "pseudo_dir": args.pseudo_dir,
              "upf"    : args.upf,
              "kmesh"  : args.kmesh,
              "kshift" : args.kshift}
   print colored("\n\nSample input for PWscf(Quantum ESPRESSO),Start","green")
   write.write_pwi(setup_dic,ibrav,celldm,struct_pm,filename=None)
   write.write_pwi(setup_dic,ibrav,celldm,struct_pm,filename="scf.in")
   print colored("Sample input for PWscf(Quantum ESPRESSO),End\n\n","green")

   #struct._visualize_struct()
