#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       match.py                                                     #
#  Dependence: parse.py, write.py (required)                                #
#  Usage:      find matched superlattices for two given structures          #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       May 01, 2017                                                 #
#                                                                           #
#===========================================================================#


import numpy as np
from termcolor import colored
from math import *
import crystal_structure as cs
import argparse
import parse
import write


desc_str = '''Objective: find match superlattice for two two-dimensional (2-D) structures'''

notes='''
Note 1: Before usage prepare the two 2-D structures in POSCAR1, POSCAR2, respectively, in VASP5-POSCAR format
Note 2: Set the vacuum layer along the c-axis of the unit cell
Note 3: The vertical positiions of atoms will not be changed, please pre-arrange them properly in the two POSCAR files
'''

alert='''
================================WARNING!===================================
   Benchmark for this code has been carried out for some specific cases.
   But the author cannot warrant full correctness
   Be cautious when you use this code for real research or pulication! 
   Check your output structures carefully.
   Bug report or improvement suggestions are welcome via email.
===========================================================================
'''

parser = argparse.ArgumentParser(prog='match_latt.py', description = desc_str)
parser.add_argument('--poscar_1', type=str, default='POSCAR1', help='POSCAR file for the  first structure')
parser.add_argument('--poscar_2', type=str, default='POSCAR2', help='POSCAR file for the second structure')
parser.add_argument('--maxarea', type=float, default=200, help='maximum area for the superlattice, in unit of angstrom square')
parser.add_argument('--tolerance', type=float, default=0.02, help='tolerance for lattice mismatch')
parser.add_argument('--min_sc_angle',type=float,default=30,help='minimum angle for supercells, in degree, in rangg of (0,180)')
parser.add_argument('--max_sc_angle',type=float,default=150,help='maximum angle for supercells, in degree, in range of (min_sc_angle, 180)')

def gcd(a, b): return gcd(b, a % b) if b else a
def area_cell(cell): return abs(np.linalg.det(np.matrix([cell[0],cell[1],(0,0,1)],float)))

#find the commensurate supercell area for two 2-D structures
def find_match_area(cell_1,cell_2,maxarea,tolerance):
    area_1=area_cell(cell_1)
    area_2=area_cell(cell_2)
    n1=area_1/area_2
    scale_list=[]
    for i in range(1,int(maxarea/area_1)):
        for j in range(1,int(maxarea/area_2)):
           n2=float(i)/float(j)
           if abs(n1*n2-1)<tolerance and gcd(i,j)==1:
              scale_list.append((i,j))
    return scale_list

#find possible shpaes for a supercell with spefified scale relative to the unit cell
def find_supercell(scale):
    trans=[]
    for i in range(1,scale+1):
        for j in range(scale+1):
            for m in range(-scale-1,scale+1):
               for n in range(1,scale+1):
                   if i*n-j*m == scale:
                      trans.append([i,j,m,n])
    return trans

#Provide basic infomation of a supercell built from a unit cell by a rotation matrix [[i,j],[m,n]]
def supercell_info(trans,cell):
    i,j,m,n=trans
    sc_cell=[cell[0]*i+cell[1]*j,cell[0]*m+cell[1]*n,cell[2]]
    sc_a=np.linalg.norm(sc_cell[0])
    sc_b=np.linalg.norm(sc_cell[1])
    sc_angle=acos(np.dot(sc_cell[0],sc_cell[1])/sc_a/sc_b)*180/pi               #angle in degree
    return [sc_a,sc_b,sc_angle]

def err(x,y):
    if x==y:
       return 0
    elif y!=0:
       return abs(x/y-1)
    else:
       return abs(y/x-1)

#Judge whether two supercells are matched within specified tolerance
def match_supercell(cell_1,trans_1,cell_2,trans_2,tolerance,min_sc_angle,max_sc_angle):
    sa1,sb1,sangle1 = supercell_info(trans_1,cell_1)
    sa2,sb2,sangle2 = supercell_info(trans_2,cell_2)
    if err(sangle1, sangle2)<tolerance and sangle1>=min_sc_angle and sangle2>=min_sc_angle and sangle1<=max_sc_angle and sangle2<=max_sc_angle:
       if (err(sa1,sa2)<tolerance and err(sb1,sb2)<tolerance):
          #print 'angle1=',sangle1,'angle2=',sangle2,', angle matched!'
          #print 'a1=',sa1,', b1=',sb1,', a2=',sa2,', b2=',sb2,', lattice length matched!'
          return 'T'
       elif (err(sa1,sb2)<tolerance and err(sa2,sb1)<tolerance):
          #print 'angle1=',info1[2],'angle2=',info2[2],', angle matched!'
          #print 'a1=',info1[0],', b1=',info1[1],', a2=',info2[0],', b2=',info2[1],', lattice length matched! But a rotation is needed!'
          return 'R'
    else:
       return 'F'

def find_match_sc(scale_set,cell_1,cell_2):
    match_sc_set=[]
    for scale_1,scale_2 in scale_set:
        sc_set_1=find_supercell(scale_1)
        sc_set_2=find_supercell(scale_2)
        for sc_1 in sc_set_1:
            for sc_2 in sc_set_2:
               if match_supercell(cell_1,sc_1,cell_2,sc_2,args.tolerance,args.min_sc_angle,args.max_sc_angle) == 'T':
                  print supercell_info(sc_1,cell_1)
                  print supercell_info(sc_2,cell_2)
                  match_sc_1 = np.matrix(([[sc_1[0],sc_1[1],0],[sc_1[2],sc_1[3],0],[0,0,1]]),float)
                  match_sc_2 = np.matrix(([[sc_2[0],sc_2[1],0],[sc_2[2],sc_2[3],0],[0,0,1]]),float)
                  match_sc_set.append((scale_1,match_sc_1,scale_2,match_sc_2))
                  print "Struct 1: ", scale_1, sc_1, ", Struct 2: ",scale_2,sc_2
                  print match_sc_2
    return match_sc_set


if __name__=='__main__':
   print ''
   print colored(desc_str,'blue')
   print colored(notes,'green')
   print colored(alert,'red')

   args = parser.parse_args()
   struct_1 = parse.parse_poscar(args.poscar_1)
   struct_2 = parse.parse_poscar(args.poscar_2)

   print "Finding matched area for superlattices...",
   scale_set = find_match_area(struct_1._cell,struct_2._cell,args.maxarea,args.tolerance)
   print "done"

   print "Finding matched superlattices..."
   match_sc_set=find_match_sc(scale_set,struct_1._cell,struct_2._cell)
   print 'done\n {0} matched supercells found in total.'.format(len(match_sc_set))

   print "Building superlattices...",
   superlatt_list=[]
   for scale_1,sc_1,scale_2,sc_2 in match_sc_set:
       sc_struct_1 = struct_1.build_supercell(sc_1)
       sc_struct_2 = struct_2.build_supercell(sc_2)
       sc_cell    = sc_struct_2._cell
       sc_species = sc_struct_1._species + sc_struct_2._species
       sc_symbols = sc_struct_1._symbols + sc_struct_2._symbols
       sc_counts  = sc_struct_1._counts  + sc_struct_2._counts
       sc_pos     = np.vstack([sc_struct_1._pos,sc_struct_2._pos])
       superlatt = {'scale_1':scale_1,'sc_1':sc_1,'scale_2':scale_2,'sc_2':sc_2,
                    'struct':cs.cryst_struct(sc_cell,sc_species,sc_symbols,sc_counts,sc_pos)}
       superlatt_list.append(superlatt)
   print "done"

   print "Writing structures...",
   for sc_index,superlatt in enumerate(superlatt_list):
       flpos='POSCAR_'+str(superlatt['scale_1'])+'_'+str(superlatt['scale_2'])+'_'+str(superlatt['struct']._natom)+'atoms_sc'+str(sc_index)
       write.write_poscar_head(1,superlatt['struct']._cell,flpos,flpos)
       write.write_poscar_atoms(superlatt['struct'],flpos)
   print "done"
