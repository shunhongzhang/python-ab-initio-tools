#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       v2openmx.py                                                  #
#  Dependence: parse.py,write.py                                            #
#  Usage:      convert the POSCAR file to the dat file for openmx           #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Apr 21, 2017                                                 #
#                                                                           #
#===========================================================================#

import os
import parse
import argparse
import write

if __name__=='__main__':
   desc_str='convert the POSCAR to the dat file for OPENMX'
   parser = argparse.ArgumentParser(prog='readdos.py', description = desc_str)
   parser.add_argument('--prefix',type=str,default='system',help='the prefix for the input file')
   args = parser.parse_args()

   datapath=os.popen('locate DFT_DATA13').readlines()[0].rstrip('\n')
   print 'find DATA.PATH: ',datapath,'\nplease confirm that the PAO and VPS files are in this directory!\n'
   struct = parse.parse_poscar('POSCAR')
   write.write_openmx_dat(args.prefix,struct,datapath)

