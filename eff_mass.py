#!/usr/bin/python
#===========================================================================#
#                                                                           #
#  File:       eff_mass.py                                                  #
#  Dependence: none                                                         #
#  Usage:      calculate the effective mass of carriers from bands          #
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       May 01, 2017                                                 #
#                                                                           #
#===========================================================================#
import os
import numpy as np
from termcolor import colored

def get_rec_cell(filename):
    try:
        gl = open(filename)
    except:
        exit("fail to open ",filename)
    gl.readline()
    scale = float(gl.readline())
    cell=scale*np.fromfile(gl,sep=' ',count=9).reshape(3,3)
    rec_cell = np.zeros((3,3),float)
    real_cell_volume = abs(np.linalg.det(cell))
    for i in range(3):
        rec_cell[i] = 2*np.pi*np.cross(cell[np.mod(i+1,3)],cell[np.mod(np.mod(i+1,3)+1,3)])/real_cell_volume
    return rec_cell

def parse_eigenval(filename='EIGENVAL',efermi=0):
    rec_cell = get_rec_cell('POSCAR')
    if not efermi:
       efermi = float(os.popen("grep E-fermi OUTCAR 2>/dev/null").readline().split()[2])
    nspin = int(os.popen("grep SPIN OUTCAR 2>/dev/null").readline().split()[2])
    print "Fermi energy = ",efermi,"eV\nISPIN = ",nspin
    print colored("Warning: The Fermi energy is read from the band structure calculation and may be incorrect!",'red')
    f = open(filename,'r')
    [f.readline() for i in range(4)]
    system = f.readline()
    line = f.readline().split()
    nelect,nkpt,nband=[int(item) for item in line]
    val_idx,con_idx = nelect/2-1,nelect/2+1
    kpt = np.zeros((nkpt,3),float)
    kweight = np.zeros(nkpt,float)
    energy = np.zeros((nkpt,nband,nspin),float)
    for ikpt in range(nkpt):
       f.readline()
       line = f.readline().split()
       kpt[ikpt] = ([float(line[i]) for i in range(3)])
       kweight[ikpt] = (float(line[3]))
       for iband in range(nband):
          line = f.readline().split()
          energy[ikpt,iband,0] = float(line[1]) - efermi
          if nspin == 2:
              energy[ikpt,iband,1] = float(line[2]) - efermi
    f.close()
    #tranfer the K points coordinates into cartesian format in the reciprocal space
    kpt = np.matrix(kpt,float)*rec_cell
    print "\nNumber of Electrons = ",nelect,"\nNumber of K points =",nkpt,"\nNumber of bands = ",nband,"\n"
    return nspin,nelect,efermi,kpt,kweight,energy

def get_path(kpt):
    dk = [np.linalg.norm(kpt[ikpt] - kpt[ikpt-1]) for ikpt in range(1,len(kpt))]
    return np.concatenate(([0],[np.cumsum(dk)[i] for i in range(len(np.cumsum(dk)))]),axis=0)

def get_eff_mass(kpath,energy):
    eVtoHartree=27.21138602
    AngstromtoBohr=0.52917721067
    if len(kpath)!=len(energy):
       exit( 'The numbher of k-points and band energys should be the same!')
    coeff= np.polyfit(kpath, energy, 2)
    eff_mass = 1/(coeff[0]/eVtoHartree/AngstromtoBohr**2)
    data=open('band_data','w')
    for i in range(len(kpath)):
        print>>data,'{0:10.7f} {1:10.7f}'.format(kpath[i],energy[i])
    data.close()
    return eff_mass  #in unit of the rest mass of electron

#======================
#main body of the code
#======================
nspin,nelect,efermi,kpt,kweight,energy = parse_eigenval()
path=get_path(kpt)
val_idx=nelect/2-1
print colored('path length:','red'),path[-1],colored('angstrom**-1','red')
print colored('carrier effective mass at k-point (certesian coord., in unit of Angstrom**-1): ','red'),kpt[0]
print colored('transport direction','red'),(kpt[-1]-kpt[0])/np.linalg.norm(kpt[-1]-kpt[0])
print 'valance band is band #',val_idx,' (band counting starts from 0)'
for ispin in range(nspin):
   print colored('\nspin channel {0:2d}:'.format(ispin),'red')
   for iband in range(val_idx-1,val_idx+3):
       print "band index:",iband,
       if iband<=val_idx:
          print ', occupied  ',
       else:
          print ',    empty  ',
       emass = get_eff_mass(np.array(path),np.array(energy[:,iband,ispin]))
       if emass<0:
          print '    hole eff_mass: {0:10.6f} m_e0'.format(emass)
       else:
          print 'electron eff_mass: {0:10.6f} m_e0'.format(emass)
