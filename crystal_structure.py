#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       crystal_structure.py                                         #
#  Dependence:                                                              #
#  Usage:      define a class for crystal structures                        #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Aug 27, 2017                                                 #
#                                                                           #
#===========================================================================#

import os
import numpy as np
from math import *
from termcolor import colored


def parse_cif(fil):
    f=open(fil)
    title=f.readline()
    a,b,c = tuple([float(item.split()[1]) for item in os.popen('grep _cell_length '+fil).readlines()])
    alpha,beta,gamma = tuple([float(item.split()[1]) for item in os.popen('grep _cell_angle '+fil).readlines()])
    alpha_r, beta_r,gamma_r = alpha*np.pi/180,beta*np.pi/180,gamma*np.pi/180
    v1=[a,0,0]
    v2=[b*cos(gamma_r),b*sin(gamma_r),0]
    v3=[c*cos(beta_r),c*cos(beta_r)*cos(gamma_r)/sin(gamma_r),
        c*sqrt( 1 + 2*cos(alpha_r)*cos(beta_r)*cos(gamma_r)- cos(alpha_r)**2-cos(beta_r)**2-cos(gamma_r)**2)/sin(gamma_r)]
    cell = np.array([v1,v2,v3],float)
    print cell
    atoms =  np.array([item.split() for item in os.popen('grep Uiso '+fil).readlines()])
    sites=atoms[:,0]
    symbols=atoms[:,1]
    pos=np.array(atoms[:,2:5],float)

def_ibrav='''
ibrav =  0: Free lattice (The card CELL_PARAMETERS required in this case)
ibrav =  1: Simple cubic lattice (sc),        v1=a(1,0,0),       v2=a(0,1,0),                       v3=a(0,0,1)
ibrav =  2: Face-centered cubic (fcc),        v1=(a/2)(-1,0,1),  v2=(a/2)(0,1,1),                   v3=(a/2)(-1,1,0)
ibrav =  3: Body-centered cubic (bcc),        v1=(a/2)(1,1,1),   v2=(a/2)(-1,1,1),                  v3=(a/2)(-1,-1,1)
ibrav =  4: Hexagonal and Trigonal P,         v1=a(1,0,0),       v2=a(-1/2,sqrt(3)/2,0),            v3=a(0,0,c/a)
ibrav =  5: Trigonal R, 3 fold axis c,        v1=a(-tx,ty,tz),   v2=a(0,2ty,tz),                    v3=a(-tx,-ty,tz)
ibrav = -5: Trigonal R, 3 fold axis <111>,    v1=(u,v,v),        v2=(v,u,v),                        v3=(v,v,u)        (under test)
ibrav =  6: Tetragonal P (st),                v1=a(1,0,0),       v2=a(0,1,0),                       v3=a(0,0,c/a)
ibrav =  7: Body-centered tetragonal (bct),   v1=a/2*(1,-1,c/a), v2=a/2*(1,1,c/a),                  v3=a/2*(-1,-1,c/a)
ibrav =  8: Orthorhombic P,                   v1=a(1,0,0),       v2=a(0,b/a,0),                     v3=a(0,0,c/a)
ibrav =  9: Base-centered orthorhombic,       v1=(a/2,b/2,  0),  v2=(-a/2,b/2,  0),                 v3=(   0,   0,  c)  (Note: C as the base face)
ibrav = -9: Base-centered orthorhombic,       v1=(a/2,-b/2, 0)   v2=( a/2,b/2,  0),                 v3=(   0,   0,  c)  (Note: C as the base face)
ibrav = 10: Face-centered orghogonal (fco),   v1=(a/2, 0, c/2),  v2=( a/2,b/2,  0),                 v3=(   0, b/2,c/2)
ibrav = 11: Body-centered orthorhombic (bco), v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),                 v3=(-a/2,-b/2,c/2)
ibrav = 12: Monoclinic P,                     v1=(a,0,0),        v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
ibrav =-12: Monoclinic P,                     v1 = (a,0,0),      v2 = (0,b,0),                      v3 = (c*cos(beta),0,c*sin(beta))
ibrav = 13: Based-centered monoclinic,        v1=(a/2, 0,-c/2),  v2=(b*cos(gamma),b*sin(gamma),0),  v3=(a/2,0,c/2)      (Note: B as the base face)
ibrav = 14: Triclinic,                        v1=(a,0,0),        v2=(b*cos(gamma),b*sin(gamma),0),  v3=(c*cos(beta),c*cos(beta)cos(gamma)/sin(gamma),c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
                     - cos(alpha)^2-cos(beta)^2-cos(gamma)^2)/sin(gamma))
'''

class cryst_struct(object):
      def __init__(self,cell,species,symbols,counts,pos):
          self._cell=cell
          self._species=species
          self._symbols=symbols
          self._counts=counts
          self._pos=pos                    # atomic positions in fracitonal coordinates
          self._natom=sum(counts)
          self._center=np.average(pos,axis=0)
          #self._pos_cart = np.dot(self._pos-origin,cell)
          self._pos_cart = np.dot(self._pos,cell)
      
      def _latt_param(self):
          a,b,c = np.linalg.norm(self._cell,axis=1)
          alpha = acos(np.dot(self._cell[1],self._cell[2].T)/b/c)/pi*180
          beta  = acos(np.dot(self._cell[2],self._cell[0].T)/c/a)/pi*180
          gamma = acos(np.dot(self._cell[0],self._cell[1].T)/a/b)/pi*180
          latt  = {'a':a,'b':b,'c':c,'alpha':alpha,'beta':beta,'gamma':gamma}
          return latt

      def _real_cell_volume(self): return abs(np.linalg.det(self._cell))

      def _reciprocal_cell(self):           # Calcualte the reciprocal lattice vectors
         rec_cell = np.zeros((3,3),float)
         for i in range(3):
             rec_cell[i] = 2*np.pi*np.cross(self._cell[np.mod(i+1,3)],self._cell[np.mod(np.mod(i+1,3)+1,3)])/self._real_cell_volume()
         return rec_cell


      def _bond_length(self,i,j):
          dist=[]
          for ii in range(-1,2):
              for jj in range(-1,2):
                  for kk in range(-1,2):
                      vec_R=ii*self._cell[0]+jj*self._cell[1]+kk*self._cell[2]
                      dist.append(np.linalg.norm(self._pos_cart[i]-self._pos_cart[j]-vec_R))
          return min(dist)
      
      def _is_slab(self,vacuum=10):
          max_diff = np.max(self._pos_cart,axis=0)-np.min(self._pos_cart,axis=0)
          slab=True
          if self._latt_param()['a'] - max_diff[0]>=vacuum:
             print 'slab structure, vacuum along crystal axis a'
          elif self._latt_param()['b'] - max_diff[1]>=vacuum:
             print 'slab structure, vacuum along crystal axis b'
          elif self._latt_param()['c'] - max_diff[2]>=vacuum:
             print 'slab structure, vacuum along crystal axis c'
          else:
             slab=False
          return slab

      def _find_symmetry(self,symprec=1e-4):
          if self._is_slab:
             print 'slab structure, be cautious when using the 3D space group for symmetry analysis!'
          find_sym_cmd = os.popen('which phonopy').read().rstrip('\n')+' --symmetry --tolerance={0}'.format(symprec) 
          if find_sym_cmd:
             get_sym=os.popen(find_sym_cmd).readlines()
             spg=get_sym[1].split()[1].strip("'")
             #spg_no = int(get_sym[1].split()[2].lstrip('(').rstrip(')')) # for old version
             spg_no = int(get_sym[2].split(':')[1])
          elif os.popen('which findsym').read() and os.popen('which pos2find').read():
             get_sym=os.popen('pos2find {}|grep "Space Group"'.format(np.log(symprec))).readline().split()
             spg_no = int(get_sym[2])
             spg=get_sym[4]
          else:
             exit('we need phonopy or ISOTROPY to help find symmetry currently, please install either of them')
          return spg,spg_no

      def _get_ibrav(self):
          spg,spg_no=self._find_symmetry()
          center_dic={"F":("Face-centered"),
                      "I":("Body-centered"),
                      "A":("Base-centered"),
                      "B":("Base-centered"),
                      "C":("Base-centered"),
                      "P":("P"            ),
                      "R":("Rhombohedral"  )}
          brav_dic={}
          [brav_dic.setdefault(num,"Triclinic"   ) for num in range(  1,  3)]
          [brav_dic.setdefault(num,"Monoclinic"  ) for num in range(  3, 15)]
          [brav_dic.setdefault(num,"Orthorhombic") for num in range( 15, 75)]
          [brav_dic.setdefault(num,"Tetragonal"  ) for num in range( 75,143)]
          [brav_dic.setdefault(num,"Trigonal"    ) for num in range(143,168)]
          [brav_dic.setdefault(num,"Hexagonal"   ) for num in range(168,195)]
          [brav_dic.setdefault(num,"Cubic"       ) for num in range(195,230)]

          ibrav_dic={("Free","P"):0,
               ("Cubic","P"):1, ("Cubic","Face-centered"):2, ("Cubic","Body-centered"):3,
               ("Hexagonal","P"):4,
               ("Trigonal","P"):4,("Trigonal","Rhombohedral"):5,
               ("Tetragonal","P"):6,("Tetragonal","Body-centered"):7,
               ("Orthorhombic","P"):8,("Orthorhombic","Base-centered"):9,("Orthorhombic","Face-centered"):10,("Orthorhombic","Body-centered"):11,
               ("Monoclinic","P"):12,("Monoclinic","Base-centered"):13,
               ("Triclinic","P"):14}
          brav   = brav_dic[spg_no]
          center = center_dic[spg[0]]
          ibrav=ibrav_dic[brav,center]
          return ibrav,brav,center


      def _find_celldm(self,ibrav=None):
          bohr_to_angs=0.52917720859
          if not ibrav:
             ibrav,brav,center=self._get_ibrav()
          print 'ibrav=',ibrav
          latt=self._latt_param()
          celldm=np.zeros(7,float)
          celldm[1]=latt['a']/bohr_to_angs     #in unit of Bohr
          if ibrav>=8:
             celldm[2]=latt['b']/latt['a']
          if ibrav==4 or ibrav>=6:
             celldm[3]=latt['c']/latt['a']
          if abs(ibrav)==5:
             a=latt['a']
             c=latt['c']
             R_cosalpha=(2*c**2-3*a**2)/(2*c**2+6*a**2)
             celldm[4]=R_cosalpha
             celldm[1]=a/sqrt(2-2*R_cosalpha)/bohr_to_angs
          if ibrav==12 or ibrav==13:
             celldm[4]=abs(cos(latt['gamma']/180*pi))
          if ibrav==14:
             celldm[4]=abs(cos(latt['alpha']/180*pi))
             celldm[5]=abs(cos(latt['beta']/180*pi))
             celldm[6]=abs(cos(latt['gamma']/180*pi))
          return celldm


      def _get_label_k(self):
          if self._is_slab==False:
             ibrav,brav,center=self._get_ibrav()
             label_dic = {'G':[0.0, 0.0, 0.0],
                          'K':[1./3.,1./3.,0.0]}
             if ibrav==4:
                label_dic.setdefault('M',[0.5,  0.0,  0.0])
                label_dic.setdefault('K',[1./3.,1./3.,0.0])
                label_dic.setdefault('A',[0.0,  0.0,  0.5])
             elif ibrav==6 or ibrav==8:
                label_dic.setdefault('X',[0.5, 0.0, 0.0])
                label_dic.setdefault('M',[0.5, 0.5, 0.0])
                label_dic.setdefault('Z',[0.0, 0.0, 0.5])
                if ibrav==8:
                   label_dic.setdefault('Y',[0.0, 0.5, 0.0])
             else:
                  print ('sorry! we do not have the label lib for this symmetry now!')


      def _visualize_struct(self):
          try: import matplotlib.pyplot as plt
          except: exit('cannot import matplotlib!')
          try: from mpl_toolkits.mplot3d import Axes3D
          except: exit('cannot import mpl_toolkits!')
          fig = plt.figure()
          ax = fig.add_subplot(111, projection='3d')
          markersize=40
          colors=['r','b','g']
          orig=np.array([0,0,0])
          pp=(orig,)+tuple(self._cell)
          for j in range(4):
              lines = np.array([[pp[j],pp[i]+pp[j]] for i in range(1,4)])
              print lines
              [ax.plot3D(*(tuple(points.T)),color='g',ls='--') for points in lines]
 
          for iat in range(self._natom):
              #if magmom[iat]>0:
              ax.scatter(self._pos_cart[iat,0],self._pos_cart[iat,1],self._pos_cart[iat,2],s=50,edgecolor='blue',facecolor='r')
              #ax.scatter(tuple(self._pos[iat,:]),s=markersize,edgecolor='blue',facecolor='r')
              #elif magmom[iat]<0:
              #   ax.scatter(struct._pos[iat,0],struct._pos[iat,1],struc.t_pos[iat,2],s=markersize,edgecolor='red',facecolor='r')
          plt.show()

      def _shift_pos(self,idir,shift): # shift atoms in crystal along idir (crystal axis) with a distance $shift (in Angstreom)
          for iat in range(self._natom):
              latt=self._latt_param()
              if idir==0: length=shift/latt['a']
              elif idir==1: length=shift/latt['b']
              elif idir==2: length=shift/latt['c']
              else: exit('invalid idir for shift! valide values: 0/1/2')
              self._pos[iat,idir] += shift/length
              while self._pos[iat,idir]>1: self._pos[iat,idir] += -1
             
      def build_supercell(self,sc,eps=-1e-4):
          def image_atom(rx,ry,rz,pos):
              return np.matrix(([[pos[0]+rx,pos[1]+ry,pos[2]+rz]]),float)
          def is_pos_at_home(pos,eps=1e-4):
              return (pos[0]>=eps and pos[0]-1.0<eps and pos[1]>=eps and pos[1]-1.0<eps and pos[2]>=eps and pos[2]-1.0<eps)

          #np.linalg.det(sc)
          sc_scale=np.linalg.det(sc)
          if sc_scale==0:
             print colored('Warning: The supercell size is zero!','red')
          elif sc_scale<0:
             print colored("The supercell basis vectors are not right-handed!",'red')

          sc_scale=abs(sc_scale)
          sc_counts=[int(round(item*sc_scale)) for item in self._counts]

          sc_species=self._species
          sc_symbols=[]
          sc_cell=sc*self._cell
          sc_pos=np.zeros((0,3),float)
          nimage_atom=np.zeros(self._natom,float)

          rb1 = int(round(abs(sc[0,0])+abs(sc[1,0])+abs(sc[2,0])))+1
          rb2 = int(round(abs(sc[0,1])+abs(sc[1,1])+abs(sc[2,1])))+1
          rb3 = int(round(abs(sc[0,2])+abs(sc[1,2])+abs(sc[2,2])))+1
 
          for iat in range(self._natom):
              for r1 in range(-rb1,rb1+1):
                for r2 in range(-rb2,rb2+1):
                    for r3 in range(-rb3,rb3+1):
                       if nimage_atom[iat]==sc_scale:
                           break
                       else:
                           pos_home = self._pos[iat,:]        # atom in the home unit cell
                           pos_image = image_atom(-r1,-r2,-r3,pos_home)       
                           # periodic image, using -r1 instead of r1 for possible improvement of efficiency
                           pos_image = (np.linalg.inv(sc.T)*pos_image.T).T
                           pos_image = [pos_image[0,i] for i in range(3)]
                           if is_pos_at_home(pos_image):
                              sc_pos=np.vstack([sc_pos,pos_image])
                              nimage_atom[iat]+=1
                              sc_symbols.append(self._symbols[iat])
                           #else:
                               #print '{0:4d} {1:5d} {2:5d} {3:5d}'.format(iat,r1,r2,r3),' '.join(['{0:8.5f}'.format(item) for item in pos_image])
          print sc
          if abs(sc_scale*self._natom -  float(len(sc_pos))) > abs(eps):
             print colored('=============================================WARNING!===========================================','red')
             print colored("The number of atoms in the supercell is incorrect, print it out and check it carefully!","red")
             print colored('supercell size = ',"red"),sc_scale
             print colored(', No. of atoms in unit cell = ',"red"),self._natom,
             print colored(', No. of atoms in supercell=',"red"),len(sc_pos)
             print colored('images for each atom:','red'),nimage_atom
          sc_struct=cryst_struct(sc_cell,sc_species,sc_symbols,sc_counts,sc_pos)
          return sc_struct   
          #return sc_cell,sc_pos,species,sc_symbols,sc_counts

      def redefine_lattice(self,sc1,sc2,sc3):
          sc_redef=np.matrix([sc1,sc2,sc3],dtype=float)
          struct_redef = self.build_supercell(sc_redef)
          latt=struct_redef._latt_param()
          struct_redef._cell[0] = np.array([latt['a'],0,0])
          struct_redef._cell[1] = np.array([latt['b']*cos(latt['gamma']/180*np.pi),latt['b']*sin(latt['gamma']/180*np.pi),0])
          cx=latt['c']*cos(latt['beta']/180*np.pi)
          cy=latt['c']*cos(latt['alpha']/180*np.pi)*sin(latt['gamma']/180*np.pi)
          cz=sqrt(latt['c']**2-cx**2-cy**2)
          struct_redef._cell[2] = np.array([cx,cy,cz])
          print struct_redef._cell
          return struct_redef

      def build_slab(self,h,k,l,thickness,vacuum,atom_shift):
          self._shift_pos(2,atom_shift*self._latt_param()['c'])
          print self._pos
          print 'build slab model, thickness=',thickness,'UCs, vacuum distance=',vacuum,'Angstroem'
          def gcd(a,b): return gcd(b, a % b) if b else a
          def lcm(a,b): return max(a,b)*int(a*b==0) + a*b/gcd(a,b)*int(a*b!=0)
          if abs(h)+abs(k)+abs(l)==0:
             exit('Error! The miller index should not be (000)!')
          if abs(h)+abs(k)==0:
             u = np.array(([1,0,0]),float)
             v = np.array(([0,1,0]),float)
          elif abs(k)+abs(l)==0:
             u = np.array(([0,1,0]),float)
             v = np.array(([0,0,1]),float)
          elif abs(h)+abs(l)==0:
             u = np.array(([0,0,1]),float)
             v = np.array(([1,0,0]),float)
          else:
             #intercept on the crystal axis
             hkl_lcm=lcm(h,lcm(k,l))
             v = np.zeros(3,float)
             try:
               hh = hkl_lcm/h
             except:
               hh = 0
               v=np.array(([1,0,0]),float)
             try:
               kk = hkl_lcm/k
             except:
               kk = 0
               v=np.array(([0,1,0]),float)
             try:
               ll = hkl_lcm/l
             except:
               ll=0
               v=np.array(([0,0,1]),float)

             print 'intercept on the crystal axis',hh,kk,ll
             u=np.array(([hh,-kk,0]),float)/gcd(hh,kk)
             if np.linalg.norm(v)==0:
                v=np.array(([hh,0,-ll]),float)/gcd(hh,ll)
          w=np.array(([h,k,l]),float)/gcd(gcd(h,k),l)
          print 'u= {0}\nv= {1}\nw= {2}'.format(u,v,w)
          sc = np.matrix(([u,v,w*thickness]),float)
          sc_size=np.linalg.det(sc)
          if sc_size==0:
             exit('error! supercell sizse is 0')
          elif sc_size<0:
             sc = np.matrix(([v,u,w*thickness]),float)
          struct_slab = self.build_supercell(sc)
          # rescale the lattice constant and the fractional coordinates in the vacuum direction
          cc = np.linalg.norm(struct_slab._cell[2])
          struct_slab._pos[:,2]=struct_slab._pos[:,2]*cc/(cc+vacuum)+vacuum/2/(cc+vacuum)
          struct_slab._cell[2]=struct_slab._cell[2]/cc*(cc+vacuum)
          return struct_slab


      def build_tube(self,n1,m1):
          def find_normal_vector(u,v,n1,m1,eps=1e-4,max_int=100):
              # Asume u and v are to noncollinear vectors, now we have A=n1*u+m1*v.
              # We want to find B=n2*u+m2*v which is (almost) perpendicular to A. 
              # Here n1,n2,m1,m2 are all integers.
              A = n1*u + m1*v
              LA = np.linalg.norm(A)
              for nn in range(max_int+1):
                  for mm in range(max_int+1):
                      vecs=[(nn,mm),(-nn,mm),(nn,-mm),(-nn,-mm)]
                      for (n2,m2) in vecs:
                          connect=np.matrix(([A,n2*u+m2*v,[0,0,1]]),float)
                          if abs(np.linalg.det(connect))>0 and abs(np.dot(A,n2*u+m2*v)/LA/np.linalg.norm(n2*u+m2*v))<=eps:
                             return n2,m2
                          else:
                             pass
              if n2==max_int or m2==max_int:
                 exit("cannot find the normal vector!")

          n2,m2 = find_normal_vector(self._cell[0],self._cell[1],n1,m1)
          sc = np.matrix(([[n1,m1,0],[n2,m2,0],[0,0,1]]),float)
          struct_rect = self.build_supercell(sc)

          alat,blat=30,30
          clat=np.linalg.norm(struct_rect._cell[1])
          L = np.linalg.norm(struct_rect._cell[0])  # The circumference of the tube
          R0 = L/2/np.pi                  # The radius of the tube
          ave_z=sum(struct_rect._pos[i,2]*struct_rect._cell[2,2] for i in range(struct_rect._natom))/struct_rect._natom
          Origin=[0.5,0.5,0.0]
          print "the tube radius is {0:10.7f} Angstrom".format(R0)

          tube_pos=np.zeros((struct_rect._natom,3),float)
          for iat in range(struct_rect._natom):
              theta = struct_rect._pos[iat,0]*2*np.pi
              R = R0+(struct_rect._pos[iat,2]*struct_rect._cell[2,2]-ave_z)
              tube_pos[iat,2] = struct_rect._pos[iat,1]+Origin[2]
              tube_pos[iat,0] = R/alat*np.cos(theta)+Origin[0]
              tube_pos[iat,1] = R/blat*np.sin(theta)+Origin[1]

          tube_cell = np.matrix([[alat,0,0],[0,blat,0],[0,0,clat]])
          tube_symbols=struct_rect._symbols
          tube_counts=struct_rect._counts
          struct_tube = cryst_struct(tube_cell,struct_rect._species,tube_symbols,tube_counts,tube_pos)
          return struct_tube

      def _write_poscar_head(self,scale=1,system='system',filename = 'POSCAR'):
          fpos=open(filename,'w')
          print >>fpos,system
          print >>fpos,scale
          for i in range(3):
              for j in range(3):
                  print >>fpos,'{0:15.12f}'.format(struct._cell[i,j]),
              print >>fpos,''
          fpos.close()

      def _write_poscar_atoms(self,filename='POSCAR'):
          fpos=open(filename,'a')
          for sym in self._species:
              print >>fpos,sym,
          print>>fpos,''
          for i in range(len(self._counts)):
              print >>fpos,self._counts[i],
          print >>fpos,'\nDirect'
          print >>fpos,'\n'.join('{0:15.12f} {1:15.12f} {2:15.12f}'.format(pos[0],pos[1],pos[2]) for pos  in self._pos)
          print >>fpos,''
          fpos.close()


if __name__=='__main__':
   import parse
   import termcolor
   import write
   import argparse

   desc_str='''
   input exmaple: 
   crystal_structure.py --task=redefine --sc1=1,-1,0 --sc2=1,1,0 sc3=0,0,1
   crystal_structure.py --task=slab --hkl=121
   crystal_structure.py --task=tube --chiral_num=2,4
   crystal_structure.py --task=cif
   crystal_structure.py --task=bond --atom_index=0,1
   '''
   print colored(desc_str,'green')

   parser = argparse.ArgumentParser(prog='struct_modify.py', description = desc_str)
   parser.add_argument('--poscar',type=str,default='POSCAR',help='the POSCAR file')
   parser.add_argument('--task',type=str,default=None,help='task to modify the structure')
   parser.add_argument('--sc1',type=eval,default=(1,0,0), help='sc1')
   parser.add_argument('--sc2',type=eval,default=(0,1,0), help='sc2')
   parser.add_argument('--sc3',type=eval,default=(0,0,1), help='sc3')
   parser.add_argument('--hkl',type=str,default='111',help='the Miller index to build slab model')
   parser.add_argument('--thickness',type=int,default=2,help='the thickness of the slab')
   parser.add_argument('--vacuum',type=float,default=20,help='the vacuum distance of slab model in Angstrom')
   parser.add_argument('--chiral_num',type=eval,default=(3,3),help='number to define the chiral vector of nanotube')
   parser.add_argument('--atom_index',type=eval,default=None,help='indices for two ions to calculate the distance')
   parser.add_argument('--atom_shift',type=float,default=0,help='shift atoms for building slab in layered structures')
   args = parser.parse_args()

   struct = parse.parse_poscar(args.poscar)
   flpos='POSCAR_'+args.task

   if args.task=='redefine':
      redef_struct = struct.redefine_lattice(args.sc1,args.sc2,args.sc3)
      write.write_poscar_head(1,redef_struct._cell,flpos,flpos)
      write.write_poscar_atoms(redef_struct,flpos)

   elif args.task=='slab':
      hkl=args.hkl
      if len(hkl)>3:
         exit('we do not allow negative index!')
      h=int(hkl[0])
      k=int(hkl[1])
      l=int(hkl[2])
      print 'hkl=',h,k,l
      struct_slab = struct.build_slab(h,k,l,args.thickness,args.vacuum,args.atom_shift)
      write.write_poscar_head(1,struct_slab._cell,system='slab',filename=flpos)
      write.write_poscar_atoms(struct_slab,filename=flpos)

   elif args.task=='tube':
      n,m=args.chiral_num
      print 'chiral vector:',n,m
      struct_tube = struct.build_tube(n,m)
      fltube = "POSCAR_tube_"+str(n)+"_"+str(m)
      write.write_poscar_head(1,struct_tube._cell,system=str(n)+'_'+str(m)+'_nanotube',filename=flpos)
      write.write_poscar_atoms(struct_tube,filename=flpos)

   elif args.task==None:
      struct=parse.parse_poscar()
      struct._visualize_struct()
   elif args.task=='cif':
      parse_cif('stru.cif')
   elif args.task=='bond':
      i,j=args.atom_index
      print 'min distance between {0}{1} and {2}{3} is: {4:8.5f} Angstreom'.format(struct._symbols[i],i,struct._symbols[j],j,struct._bond_length(i,j))
   else:
      exit('we currently do not accept tasks other than redefine/slab/tube!')


