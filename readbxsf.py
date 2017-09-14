#!/usr/bin/python

import numpy as np
import os
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D


#Energy resolution of the Fermi surface, eigenvalues within efermi+/-eps are included
eps=1e-3
#kpt resolution
keps=1e-3

def is_in_BZ(kpt,rec_cell):
    for i in range(-1,2):
      for j in range(-1,2):
        for k in range(-1,2):
            vec_g=np.array([i,j,k])*rec_cell  #reciprocal lattice vector
            kpt_cart=kpt*rec_cell
            vec_k=kpt_cart-vec_g/2
            vec_g=[vec_g[0,0],vec_g[0,1],vec_g[0,2]]
            vec_k=[vec_k[0,0],vec_k[0,1],vec_k[0,2]]
            if np.dot(vec_g,vec_k)>0:
               return False
    return True

def is_BZ_border(kpt,rec_cell):
    for i in range(-1,2):
      for j in range(-1,2):
        for k in range(-1,2):
            vec_g=np.array([i,j,k])*rec_cell  #reciprocal lattice vector
            vec_k=kpt-vec_g/2
            vec_g=[vec_g[0,0],vec_g[0,1],vec_g[0,2]]
            vec_k=[vec_k[0,0],vec_k[0,1],vec_k[0,2]]
            if abs(np.dot(vec_g,vec_k))<keps:
               return True
    return False


alat=float(os.popen("grep 'celldm(1)' scf.in").readline().split('=')[1])

prefix=os.popen('grep prefix fs.in').readline().split('=')[1].rstrip().lstrip('"').rstrip('"')
f=open(prefix+'_fs.bxsf','r')
[f.readline() for i in range(7)]
efermi=float(f.readline().split()[2])
[f.readline() for i in range(4)]
nbands=int(f.readline())
grid=[int(item) for item in f.readline().split()]
ngrid=grid[0]*grid[1]*grid[2]
origin=[float(item) for item in f.readline().split()]
rec_cell=np.fromfile(f,count=9,sep=' ',dtype=float).reshape(3,3)
rec_cell=np.matrix(rec_cell)
print 'reciprocal lattice vectors:'
print "in unit of 2pi/alat, alat = {0:8.5f} angstrom\n".format(alat),rec_cell
print "volume of reciprocal unit cell:",np.linalg.det(rec_cell),'\n'

eigenval=np.zeros((nbands,grid[0],grid[1],grid[2]),float)
for iband in range(nbands):
    band_index=int(f.readline().split()[1])
    data = np.fromfile(f,sep=' ',count=ngrid,dtype=float).reshape(grid)
    eigenval[iband]=data
f.close()

border=1
lborder=-border
rborder=border
gridx=grid[0]*2-1
gridy=grid[1]*2-1
gridz=grid[2]
xx=np.linspace(lborder,rborder,gridx)
yy=np.linspace(lborder,rborder,gridy)
zz=np.linspace(lborder,rborder,gridz)
kx,ky,kz=np.meshgrid(xx,yy,zz)

print is_in_BZ([0.7,0,0],rec_cell)
print "extracting data..."
kz0=0.0
kk=[]
fw=open('log','w')
for iband in range(nbands):
  print "band #{0:<2d}".format(iband)
  kk.append([])
  for i in range(gridx):
    for j in range(gridy):
      for k in range(gridz):
          if abs(eigenval[iband,abs(grid[0]-i-1),abs(grid[1]-j-1),k]-efermi)<eps:
             kpt=np.array([kx[i,j,k],ky[i,j,k],kz[i,j,k]])
             kpt_cart=(kpt*rec_cell)
             if is_in_BZ(kpt,rec_cell):
             #if is_BZ_border(kpt_cart,rec_cell):
                #if abs(kpt_cart[0,2]-kz0)<keps:
                   print>>fw,kpt
                   kk[iband].append(kpt_cart)
                   kk[iband].append(-kpt_cart) #use time reversal symmetry
fw.close()
print "done"
print "plotting Fermi surface.."
fig = plot.figure()
ax=fig.add_subplot(111)
#ax=fig.add_subpplot(111,projection='3d')
color_list=['violet','red','green','blue']
for iband in range(nbands):
    for ikk in range(len(kk[iband])):
      ax.scatter(kk[iband][ikk][0,0],kk[iband][ikk][0,1],color=color_list[iband])
      #ax.scatter(kk[iband][ikk][0],kk[iband][ikk][1],kk[iband][ikk][2],color=color_list[iband])
'''
for iax in range(3):
    plot.arrow(0,0, rec_cell[iax,0], rec_cell[iax,1], color='orange', head_width=0.05, head_length=0.02)
    #ax.annotate('b1', xy=(1, -1), xytext=(3, 1.5), arrowprops=dict(facecolor='black', shrink=0.05))
plot.xlim(-1.1,1.1)
plot.ylim(-1.1,1.1)
'''
plot.show()
print "done"
