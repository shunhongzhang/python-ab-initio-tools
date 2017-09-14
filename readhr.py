#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       readhr.py                                                    #
#  Dependence: none                                                         #
#  Usage:      read the hamiltonian matrix elements from the _hr.dat file   #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Apr 28, 2017                                                 #
#                                                                           #
#===========================================================================#


import numpy as np
import os

filename=os.popen('ls *_hr.dat').readline().rstrip()
f=open(filename,'r')
f.readline()
nwann=int(f.readline())
nrpt=int(f.readline())

degeneracy=np.fromfile(f,sep=' ',count=nrpt,dtype=int)
data = np.fromfile(f,sep=' ',count=nwann**2*nrpt*7,dtype=float).reshape(nrpt,nwann,nwann,7)
R=np.zeros((nrpt,nwann,nwann,3),int)
site_i=np.zeros((nrpt,nwann,nwann),int)
site_j=np.zeros((nrpt,nwann,nwann),int)
hopping_t=np.zeros((nrpt,nwann,nwann),complex)
for irpt in range(nrpt):
    for i in range(nwann):
        for j in range(nwann):
            R[irpt,i,j]=np.array(data[irpt,i,j,0:3],int)
            site_i[irpt,i,j]=int(data[irpt,i,j,3])
            site_j[irpt,i,j]=int(data[irpt,i,j,4])
            hopping_t[irpt,i,j]=complex(data[irpt,i,j,5],data[irpt,i,j,6])

#print 'mymodel.set_hop(t{0:<3d}, {1:2d}, {2:2d}, [ {3:2d}, {4:2d}])'.format(2,site_i[0,1,1],site_j[0,1,1],R[0,1,1,0],R[0,1,1,1]),hopping_t[0,1,1]

N=nrpt*nwann**2
# avoid double-counting by only take the first half of the data
R=R[0:(nrpt+1)/2,:,:,:].reshape((nrpt+1)/2*nwann**2,3)
site_i=site_i[0:(nrpt+1)/2,:,:].flatten()
site_j=site_j[0:(nrpt+1)/2,:,:].flatten()
t=hopping_t[0:(nrpt+1)/2,:,:].flatten()
abs_t=[abs(item) for item in t]

Rx=R[:,0]
Ry=R[:,1]
Rz=R[:,2]
print Rx.shape,site_i.shape,t.shape
get_sorted_hop = sorted(zip(abs_t,site_i,site_j,Rx,Ry,Rz,t))[::-1]
sorted_hop = [(site_i,site_j,[Rx,Ry,Rz],t) for (abs_t,site_i,site_j,Rx,Ry,Rz,t) in get_sorted_hop]

onsite=np.zeros(nwann,complex)
hop=open('hop_for_pythtb','w')
for i in range((nrpt+1)/2*nwann**2):
    if sorted_hop[i][0]==sorted_hop[i][1] and (abs(sorted_hop[i][2][0])+abs(sorted_hop[i][2][1])+abs(sorted_hop[i][2][2])==0): 
       print sorted_hop[i][0],sorted_hop[i][1],sorted_hop[i][2]
       onsite[sorted_hop[i][0]-1]=sorted_hop[i][3]
    else:
       if abs(sorted_hop[i][2][0])+abs(sorted_hop[i][2][1])+abs(sorted_hop[i][2][2])==0 and sorted_hop[i][0]>sorted_hop[i][1]:
          pass
       else:
          print >>hop,'my_model.set_hop({0:8.5f}+{1:8.5f}j, {2:2d}, {3:2d}, [ {4:2d}, {5:2d}, {6:2d}])'.format(sorted_hop[i][3].real,sorted_hop[i][3].imag,sorted_hop[i][0]-1,sorted_hop[i][1]-1,sorted_hop[i][2][0],sorted_hop[i][2][1],sorted_hop[i][2][2])
hop.close()
print 'my_model.set_onsite(['
for i in range(nwann):
    print '{0:8.5f}+{1:8.5f}j,'.format(onsite[i].real,onsite[i].imag)
print '])'

get_real_vec=os.popen("grep -A 3 'Lattice Vectors' *wout").readlines()[1:]
real_vec=[[float(get_real_vec[i].split()[j]) for j in range(1,4)]for i in range(3)]
real_vec=np.matrix(real_vec,float)
trans=np.linalg.inv(real_vec.T)


get_WF_centres=os.popen('grep X *centres.xyz').readlines()
WF_centres=np.zeros((nwann,3),float)

for i in range(nwann):
    WF_centres[i] = np.array([float(get_WF_centres[i].split()[j+1]) for j in range(3)])
    WF_centres[i] = WF_centres[i]*trans
print 'WF_centers in fracational coordinates:\n',WF_centres
