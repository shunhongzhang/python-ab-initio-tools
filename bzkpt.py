#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       bzkpt.py                                                     #
#  Dependence: none                                                         #
#  Usage:      proceed with the k points in the BZ                          #
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       May 03, 2017                                                 #
#                                                                           #
#===========================================================================#



from math import *
import numpy as np


def get_ksym_dic(ibrav):
    if ibrav==4:
       return {
       'G':[0.0, 0.0, 0.0],
       'M':[0.5, 0.0, 0.0],
       'K':[1./3.,1./3.,0.0],
       'A':[0.0, 0.0, 0.5] }
    elif ibrav==6:
       return {
       'G':[0.0, 0.0, 0.0],
       'X':[0.5, 0.0, 0.0],
       'M':[0.5, 0.5, 0.0],
       'Y':[0.0, 0.5, 0.0],
       'Z':[0.0, 0.0, 0.5] }
    else:
       print ('sorry! we do not have the label lib for this symmetry now!')
       return None

      
def get_label_k(ibrav,kpt,eps=1e-4):
    ksym_dic = get_ksym_dic(ibrav)
    ksym_list=[kpt[0]]+[x for n, x in enumerate(kpt[:-1]) if np.linalg.norm(x-kpt[n+1])<eps]+[kpt[-1]]
    label_k_list=[]
    for ksym in ksym_list:
        for key in ksym_dic.keys():
            if np.linalg.norm(np.array(ksym)-np.array(ksym_dic[key]))<eps:
               label_k_list.append(key)
    return label_k_list
    #exit('The k point is not in the high symmetry k points set!')

def get_path(kpt):
    dk = [np.linalg.norm(kpt[ikpt] - kpt[ikpt-1]) for ikpt in range(1,len(kpt))]
    return np.concatenate(([0],np.cumsum(dk)),axis=0)

def get_xsym(path,segment_nkpt): return  path[0::segment_nkpt].tolist()+[path[-1]]

def guess_xsym(path):  return [path[0]]+[x for n, x in enumerate(path) if x in path[:n]]+[path[-1]]

def gen_kp(weight=None):
    f=open("KPOINTS")
    f.readline()
    segment_nkpt=int(f.readline())
    line = f.readline()
    if line[0]!='l' and line[0]!='L':
       exit('please use line mode to generate the kpath')
    ctype = f.readline()
    if ctype[0]=='D' or ctype[0]=='d':
       print 'direct coordinates for k points'
    elif ctype[0]=='C' or ctype[0]=='c':
       print 'cartesian coordinates for k points'
    data = np.fromfile(f,sep=' ',dtype=float)
    f.close()
    data = data.reshape(data.shape[0]/3,3)
    ksym = np.vstack([data[::2],data[-1]])
    print "high symm kpts:\n",ksym
    npath=len(ksym)-1
    weight = [1.0/npath/segment_nkpt for i in range(npath*segment_nkpt) if not weight]
    kpt=np.zeros((npath,segment_nkpt,3),float)
    for ipath in range(npath):
        kpt[ipath,] = np.array([np.linspace(ksym[ipath,i], ksym[ipath+1,i], segment_nkpt) for i in range(3)]).T
    fw=open('kpath','w')
    for ipath in range(npath):
        print>>fw,'\n'.join(["{0:12.8f} {1:12.8f} {2:12.8f} ".format(kpt[ipath,ikpt,0],kpt[ipath,ikpt,1],kpt[ipath,ikpt,2]) for ikpt in range(segment_nkpt)])
    fw.close()
    return kpt


def first_BZ(sym):
    if sym==1:
       print "rectangle lattice"
    elif sym==2:
       print "hexagonal lattice,please choose type:\n"
       print "type 1: b1 = (1/2,-sqrt(3)/2), b2 = (1/2,sqrt(3)/2)"
       print "type 2: b1 = (sqrt(3)/2,-1/2), b2 = (sqrt(3)/2,1/2)"
       type=input("type=")
    x_min=-1.0
    y_min=-1.0
    x_max=1.0
    y_max=1.0
    nx=10
    ny=10
    xgrid=np.linspace(x_min,x_max,nx)
    ygrid=np.linspace(y_min,y_max,nx)
    xx,yy=np.meshgrid(xgrid,ygrid)
    kp_d=open("k-points-fractional","w")
    kp_c=open("k-points-cartesian","w")  #Note: here we assume the length of the reciprocal lattice vecotrs are 1

    r3=sqrt(3) 
    for x,y in zip(xx.flatten(),yy.flatten()):
        if sym==1:
           print>>kp_d, "{0:8.5f} {1:8.5f} {2:8.5f} {3:8.5f}".format(x,y,0,1)
        elif sym==2:
           if type==1 and x+y+1>=0 and x+y-1<=0 and x-2*y-1<=0 and x-2*y+1>=0 and 2*x-y+1>=0 and 2*x-y-1<=0:
              print >>kp_d,"{0:8.5f} {1:8.5f} {2:8.5f} {3:8.5f}".format(x,y,0,1)
              print >>kp_c,"{0:8.5f} {1:8.5f}".format(x/2+y/2,-x*r3/2+y*r3/2)
           elif type==2 and x-y+1>=0 and x-y-1<=0 and 2*x+y+1>=0 and 2*x+y-1<=0 and x+2*y+1>=0 and x+2*y-1<=0:
              print >>kp_d,"{0:8.5f} {1:8.5f} {2:8.5f} {3:8.5f}".format(x,y,0,1)
              print >>kp_c,"{0:8.5f} {1:8.5f}".format(x/2+y/2,-x*r3/2+y*r3/2)
    kp_d.close()
    kp_c.close()

def plot_kpath(kpt):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    kpt=np.array(kpt)
    ax.plot(kpt[:,0],kpt[:,1],kpt[:,2],color='green',linewidth=1.5)
    plt.show()


if __name__ == "__main__":
   #cell, latt, rec_vec, pos, species, symbols, counts=parse.parse_poscar('POSCAR')
   first_BZ(2)
   #get_label_k(6,[0,0,0])
   #gen_kp()
