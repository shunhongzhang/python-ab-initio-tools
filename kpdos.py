#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       kpdos.py                                                     #
#  Dependence: none                                                         #
#  Usage:      three classes: kpdos, electronic_dos, phdos                  #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Sep 08, 2017                                                 #
#                                                                           #
#===========================================================================#

import os
import numpy as np
import parse

orb_dic={0:'s',1:'p_y',2:'p_z',3:'p_x',
         4:'d_{xy}',5:'d_{yz}',6:'d_{z^2}',7:'d_{xz}',8:'d_{x^2-y^2}'}

spin_dic={0:'up',1:'dn'}

# parse structure from OUTCAR
def parse_outcar_struct(fil='OUTCAR'):
    counts=[int(item) for item in os.popen('grep "ions per type" {0}'.format(fil)).read().split()[4:]]
    species=[item.split('=')[1].split(':')[0] for item in os.popen('grep VRH {0}'.format(fil)).readlines()]
    symbols=[item for ispec in range(len(species)) for item in counts[ispec]*[species[ispec]]]
    get_rec_cell=os.popen('grep -A 3 "reciprocal lattice vectors" OUTCAR|tail -3').readlines()
    cell=np.matrix(([[float(coord) for coord in item.split()[:3]] for item in get_rec_cell]),float)
    rec_cell=np.matrix(([[float(coord) for coord in item.split()[3:]] for item in get_rec_cell]),float)
    get_pos = os.popen('grep -A {0} "position of ions in fractional" {1}'.format(sum(counts),fil)).readlines()[1:]
    pos=np.array([[float(item) for item in line.split()] for line in get_pos])
    try:
       import crystal_structure as cs
       return cs.cryst_struct(cell,species,symbols,counts,pos)
    except: return (cell,rec_cell,species,symbols,counts,pos)



# this function is under test
def refine_xsym_xlabels(xsym,xlabels,xlabel_thr):
    klength=xsym[-1]-xsym[0]
    new_xsym=[]
    new_labels=[]
    i=len(xsym)-1
    while i>0:
        dk=xsym[i]-xsym[i-1]
        if dk<klength/args.xlabel_thr:
           new_xsym.append(xsym[i-1]/2+xsym[i]/2)
           new_labels.append(xlabels[i-1]+xlabels[i])
           i=i-1
        else:
           new_xsym.append(xsym[i])
           new_labels.append(xlabels[i])
        i=i-1
    if i==0:
       new_xsym.append(xsym[0])
       new_labels.append(xlabels[0])
    return [item for item in reversed(new_xsym)],[item for item in reversed(new_labels)]

def parse_outcar_band(filename='OUTCAR'):
    lsoc  = os.popen('grep LSORBIT {0}'.format(filename)).read().split()[2]
    lncl  = os.popen('grep LNONCOLLINEAR {0}'.format(filename)).read().split()[2]
    lhf   = os.popen('grep LHFCALC {0}'.format(filename)).read().split()[2]
    nspin = int(os.popen('grep ISPIN {0}'.format(filename)).read().split()[2])
    nkpt  = int(os.popen('grep NKPTS {0}'.format(filename)).read().split()[3])
    nband = int(os.popen('grep NBAND {0}|head -1'.format(filename)).read().split()[-1])
    nelect= float(os.popen('grep NELECT {0}'.format(filename)).read().split()[2])
    efermi= float(os.popen('grep E-fermi {0}'.format(filename)).read().split()[2])
    startline=int(os.popen('grep -n E-fermi {0}'.format(filename)).read().split(':')[0])+2
    f=open(filename)
    [f.readline() for i in range(startline)]
    energy=[];occ=[]
    for ispin in range(nspin):
        if nspin==2: [f.readline() for i in range(2)]
        energy.append([])
        occ.append([])
        for ikpt in range(nkpt):
           [f.readline() for i in range(2)]
           data=np.fromfile(f,dtype=float,count=nband*3,sep=' ').reshape(nband,3)
           energy[ispin].append(data[:,1].tolist())
           occ[ispin].append(data[:,2].tolist())
    energy=np.array(energy)
    get_kpts=os.popen('grep -A {0} "in reciprocal lattice and weight" {1}'.format(nkpt,filename)).readlines()[1:]
    kpt = np.array([[float(kcoord.rstrip('\),')) for kcoord in kk.split()[3:6]] for kk in  get_kpts])
    alat=float(os.popen('grep SCALE {0}'.format(filename)).read().split()[2])
    get_kpts_cart=os.popen('grep -A {0} "2pi/SCALE" {1}'.format(nkpt,filename)).readlines()[1:]
    kpt_cart = np.array([[float(kcoord) for kcoord in kk.split()[:-1]] for kk in  get_kpts_cart])*2*np.pi/alat
    kweight = np.array([float(item.split()[-1]) for item in get_kpts])

    startk=0
    if lhf=='T':
       startk=kweight.tolist().index(0)
       kpt=kpt[startk:]
       kpt_cart=kpt_cart[startk:]
       kweight=kweight[startk:]
       energy=energy[:,startk:,:]
    return lsoc,lncl,lhf,nelect,efermi,kpt,kpt_cart,kweight,startk,energy
        


def plot_band(args,path,eigenval,xsym,label_k=None,output='band_structure'):
    import matplotlib.pyplot as plot
    import matplotlib.patches as mpatches


    print "Plotting band structure ..."
    nspin=eigenval.shape[0]
    nband=eigenval.shape[2]
    fig = plot.figure(figsize=args.figsize)
    ax = fig.add_subplot(1,1,1)
    paths=np.array([path for iband in range(nband)])
    for ispin in range(nspin):
        if args.style == 'line':
           [ax.plot(path,eigenval[ispin,:,iband],color=args.color.split()[ispin],lw=args.linewidth,label='spin '+spin_dic[ispin]) for iband in range(nband)]
        elif args.style == 'dot':
           ax.scatter(paths, eigenval[ispin,:,:].T, s=args.markersize, marker=args.marker, facecolor='none',edgecolor=args.color.split()[ispin],label='spin '+spin_dic[ispin])
    color_patch=[]
    if nspin>1:
       color_patch.append(mpatches.Patch(color='blue', label='spin up'))
       color_patch.append(mpatches.Patch(color= 'red', label='spin dn'))
       if args.legend_switch=='on':
          leg = ax.legend(handles=color_patch,loc='upper center',fontsize=10,numpoints=1)

    if args.elim[1]>args.elim[0]:
       plot.ylim(args.elim[0],args.elim[1])
    ymin,ymax = plot.ylim()
    if args.title:  ax.set_title(output)
    if not label_k: label_k=['$'+label+'$' for label in args.label_k.split()]
    [ax.plot([xi, xi], [ymin, ymax], color='gray', zorder=-1) for xi in xsym]
    ax.plot([path[0],path[-1]], [0, 0], color='green', ls='dashed', zorder=-1)

    if label_k:
       ax.set_xticks(xsym)
       ax.set_xticklabels(label_k)
    plot.xlim(0,path[-1])
    if args.lylabel==True: ax.set_ylabel(r'$E-E_f$ ($eV$)')
    plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
    print "Done"


class kpdos(object):
      def __init__(self,filename='PROCAR'):
          if not os.path.isfile(filename):
             exit('cannot find the '+filename+' file! put it in the current directory!')
          else:
             #fix possible format errors in the procar file
             fix='sed -i "s/\([0-9]\)\-\([0-9]\)/\\1 -\\2/g" '+filename
             os.system(fix)
          if not os.path.isfile('OUTCAR'):
             exit('cannot find the OUTCAR file! put it in the current directory!')
          if 'phase' in open(filename).readline().split(): self._lphase=True
          else: self._lphase=False
          allbands=parse_outcar_band()
          self._soc,self._lncl,self._lhf,self._nelect,self._efermi,self._kpt,self._kpt_cart,self._kweight,self._startk,self._energy=allbands
          numbers=os.popen('grep \# PROCAR 2>/dev/null|head -1').read().split()
          self._nspin = int(os.popen('grep ISPIN OUTCAR 2>/dev/null').read().split()[2])
          #self._ncl = (os.popen('grep LNONCOLLINEAR OUTCAR 2>/dev/null').read().split()[2]=='T')
          self._filename = filename
          self._nkpt = int(numbers[3])
          self._nband = int(numbers[7])
          self._nat = int(numbers[11])
          self._norb = 9
          self._ndim=1+self._nspin*(self._lphase)*2+(self._lncl=='T')*3
          #get_energy = os.popen('grep energy '+filename).readlines()
          #self._energy = np.array([float(item.split()[4]) for item in get_energy]).reshape(self._nspin,self._nkpt,self._nband)
          #self._efermi = float(os.popen('grep fermi OUTCAR|tail -1').read().split()[2])
          print 'fermi energy = {0:10.7f} eV, read from OUTCAR, please make sure it is correct\n'.format(self._efermi)
          if self._lncl=='T':
             self._ndim=4
             if self._lphase: self._ndim=6
             print 'noncollinear spin polarized calculation,ISPIN=',
          elif self._nspin==2:
             print 'collinear spin polarized calculation,ISPIN=',
          elif self._nspin==1:
             print 'non spin polarized calculation,ISPIN=',
          print self._nspin
          print '\n------------------band data summary-----------------\n'
          if self._nspin==2: print 'spin polarized ',
          print 'band structure calculated by {0}'.format('VASP')
          if self._lncl=='T': print 'spin-orbit coupling is included'
          if self._lhf=='T':  print 'Hartree-Fock type functional used?\nplease make sure the result is what you want.'
          print 'nspin = {0}'.format(self._nspin)
          print 'number of kpts  = {0}'.format(self._nkpt)
          print 'number of bands = {0}'.format(self._nband)
          print 'number of electrons = {0}'.format(self._nelect)
          print 'femri energy = {0} eV, read from OUTCAR \nplease make sure it is correct'.format(self._efermi)
          print '\n------------------band data summary-----------------\n'

      def _get_kpts(self):
          try:
             import parse
             struct = parse.parse_poscar('POSCAR')
          except:
             struct = parse_outcar_struct()
          get_kpt=os.popen('grep k-point PROCAR').readlines()[1:self._nkpt+1]
          kpt=np.matrix([[float(item.split()[3]),float(item.split()[4]),float(item.split()[5])] for item in get_kpt])
          rec_cell=struct._reciprocal_cell()
          kpt = np.matrix(kpt,float)*rec_cell
          return kpt

      def _get_spin_weight(self,args):
          if not self._lncl: exit('error! non collinear calculation is required to get the sz component')
          f=open(self._filename)
          lines=f.readlines()
          data=[]
          # Note: the index of iatom start from 0, but iat starts from 1 (in PROCAR)
          nn=3
          fmt='{0:'+str(nn)+'d}'
          # Note: the index of iatom start from 0, but iat starts from 1 (in PROCAR)
          if args.proj_spinor_index==-1: iat='tot'
          else: iat=str(fmt.format(args.proj_spinor_index+1))
          for line in lines:
              if line.startswith(iat):
                 data.append(float(line.rstrip('\n').split()[-1]))
          data=np.array(data)
          if self._lhf=='T':
             startk=self._kweight.tolist().index(0)*self._ndim
             data=data[startk:]
          sx_weight = data[1::self._ndim].reshape(self._nspin,self._nkpt,self._nband)
          sy_weight = data[2::self._ndim].reshape(self._nspin,self._nkpt,self._nband)
          sz_weight = data[3::self._ndim].reshape(self._nspin,self._nkpt,self._nband)
          return sx_weight,sy_weight,sz_weight

      def _plot_spinor_band(self,args):
          import matplotlib
          matplotlib.use('Agg')
          import matplotlib.pyplot as plot
          import matplotlib.patches as mpatches

          if args.proj_spinor==False: return 0
          dk = [np.linalg.norm(self._kpt_cart[ikpt] - self._kpt_cart[ikpt-1]) for ikpt in range(1,self._nkpt)]
          path = np.concatenate(([0],[np.cumsum(dk)[i] for i in range(len(np.cumsum(dk)))]),axis=0)
          xsym=[path[0]]+[x for n, x in enumerate(path) if x in path[:n]]+[path[-1]]
          xticks,xlabels=refine_xsym_xlabels(xsym,args.label_k.split(),args.xlabel_thr)

          sx_weight,sy_weight,sz_weight=self._get_spin_weight(args)
          print np.max(sx_weight),np.max(sy_weight),np.max(sz_weight)
          fig = plot.figure(figsize=args.figsize)
          paths=np.array([path for iband in range(self._nband)])
          if args.spinor_dir==0 :cc=(sx_weight+1)*10;prefix='$S_x$'
          if args.spinor_dir==1 :cc=(sy_weight+1)*10;prefix='$S_y$'
          if args.spinor_dir==2 :cc=(sz_weight+1)*10;prefix='$S_z$'
          if args.nplot==(1,2):
             ax1=fig.add_subplot(131)
             cax1 = ax1.scatter(paths,self._energy[0].T,s=args.markersize,facecolor=np.where(cc>10, cc, 0),edgecolor='None',cmap='cool')
             #cbar = fig.colorbar(cax,ticks=[np.min(cc),(np.min(cc)+np.max(cc))/2,np.max(cc)])
             #cbar.ax.set_yticklabels([prefix+'='+(item) for item in ['-0.5','0','0.5']])
             ax2=fig.add_subplot(132)
             cax2 = ax2.scatter(paths,self._energy[0].T,s=args.markersize,facecolor=cc,edgecolor='None',cmap='cool')
             #cbar = fig.colorbar(cax,ticks=[np.min(cc),(np.min(cc)+np.max(cc))/2,np.max(cc)])
             #cbar.ax.set_yticklabels([prefix+'='+(item) for item in ['-0.5','0','0.5']])
             if args.label_k:
                [ax.set_xticks(xsym) for ax in (ax1,ax2)]
                [ax.set_xticklabels(['$'+label+'$' for label in args.label_k.split()]) for ax in (ax1,ax2)]
                [ax.set_ylim(args.elim[0],args.elim[1]) for ax in (ax1,ax2) if args.elim[1]>args.elim[0]]
                ymin,ymax=plot.ylim()
                [[ax.plot([xi, xi], [ymin, ymax], color='gray', zorder=-1) for xi in xsym] for ax in (ax1,ax2)]
                [ax.plot([path[0],path[-1]], [0, 0], color='green', ls='dashed', zorder=-1) for ax in (ax1,ax2)]
                [ax.set_xlim(path[0],path[-1]) for ax in (ax1,ax2)]

          else:
             ax=fig.add_subplot(111)
             cax = ax.scatter(paths,self._energy[0].T,s=args.markersize,facecolor=cc,edgecolor='None',cmap='cool')
             cbar = fig.colorbar(cax,ticks=[np.min(cc),(np.min(cc)+np.max(cc))/2,np.max(cc)])
             cbar.ax.set_yticklabels([prefix+'='+(item) for item in ['-0.5','0','0.5']])

             if args.elim[1]>args.elim[0]:  plot.ylim(args.elim[0],args.elim[1])
             ymin,ymax = plot.ylim()
             if args.title:  ax.set_title(output)
             label_k=['$'+label+'$' for label in args.label_k.split()]
             [ax.plot([xi, xi], [ymin, ymax], color='gray', zorder=-1) for xi in xsym]
             ax.plot([path[0],path[-1]], [0, 0], color='green', ls='dashed', zorder=-1)

             if label_k:
                ax.set_xticks(xsym)
                ax.set_xticklabels(label_k)
             plot.xlim(0,path[-1])
             if args.lylabel==True: ax.set_ylabel(r'$E-E_f$ ($eV$)')
          plot.savefig('sz_band.png', dpi=args.dpi, bbox_inches='tight')
 


      def _get_kpdos(self,iatom):
          f=open(self._filename)
          lines=f.readlines()
          data=[]
          nn=3
          fmt='{0:'+str(nn)+'d}'
          # Note: the index of iatom start from 0, but iat starts from 1 (in PROCAR)
          iat=str(fmt.format(iatom+1))
          for line in lines:
              if line.startswith(iat):
                 data.append(float(line.rstrip('\n').split()[-1]))
          data=np.array(data)
          if self._lhf=='T':
             startk=self._kweight.tolist().index(0)*self._ndim
          data = data[::self._ndim]
          kpdos = data.reshape(self._nspin,self._nkpt,self._nband)
          if self._lhf=='T':
             startk=self._kweight.tolist().index(0)*self._ndim
             kpdos = kpdos[:,startk:,:]
          return kpdos

      def _get_kpdos_orb(self,iatom):
          f=open(self._filename)
          lines=f.readlines()
          data=[]
          nn=3
          fmt='{0:'+str(nn)+'d}'
          # Note: the index of iatom starts from 0, but iat starts from 1 (in PROCAR)
          iat=str(fmt.format(iatom+1))
          for line in lines:
              if line.startswith(iat):
                 data.append([float(item) for item in line.rstrip('\n').split()[1:10]])
          data=np.array(data)
          kpdos_orb = data[::self._ndim].reshape(self._nspin,self._nkpt,self._nband,self._norb)
          if self._lhf=='T':
             startk=self._kweight.tolist().index(0)*self._ndim
             kpdos = kpdos[:,sta:,startk:,:,:]
          return kpdos_orb


      def _get_bandgap_old(self):
          nelect=float(os.popen('grep NELECT OUTCAR').read().split()[2])
          if nelect-int(nelect)!=0:
             print "Warning: NELECT is non-ingeter, the band gap reported here may be incorrect!"
          else:
             nelect=int(nelect)
             print 'number of electron: {0}'.format(nelect)
          val_idx = nelect/2-1
          con_idx = val_idx+1
          print 'assign the {0} band and {1} band as valance and conduction band respectively'.format(val_idx+1,con_idx+1)
          VBM = np.zeros((self._nspin),float)
          CBM = np.zeros((self._nspin),float)
          for ispin in range(self._nspin):
              VBM[ispin] = max(self._energy[ispin,:,val_idx]-self._efermi)
              CBM[ispin] = min(self._energy[ispin,:,con_idx]-self._efermi)
              print "spin channel {0:1d}".format(ispin)
              print "VBM is {0:10.6f} eV, CBM is {1:10.6f} eV".format(VBM[ispin],CBM[ispin])
              if CBM[ispin] - VBM[ispin] >=0:
                 print "band gap is {0:10.6f} eV\n".format(CBM[ispin]-VBM[ispin])
              else:  print "metallic!\n"


      def _get_bandgap(self):
          print '\n--------------band gap estimated from bands along high symmetry kpath------------\n'
          if self._nspin==2: print '{0:4s}'.format('spin'),
          print ' '.join(['{0:12s}'.format(item.rjust(12)) for item in ['VBM(eV)','VBM_index','CBM(eV)','CBN_index','band gap']])
          for ispin in range(self._nspin):
              energy = self._energy-self._efermi
              VBM=max([item for item in energy[ispin].flatten().tolist() if item <0])
              CBM=min([item for item in energy[ispin].flatten().tolist() if item >0])
              iv=np.mod(energy[ispin].flatten().tolist().index(VBM),self._nband)
              ic=np.mod(energy[ispin].flatten().tolist().index(CBM),self._nband)
              if self._nspin==2: print '{0:4s}'.format(spin_dic[ispin]),
              if ic<=iv:
                 print '{0:12.6f} {1:11d} {2:12.6f} {3:11d} {4:14.6f}'.format(VBM,iv,CBM,ic,VBM-CBM),
                 print ' metallic ? \n\n  we found (VBM index >= CBM index, please check band structure'
              else:
                 print '{0:12.6f} {1:11d} {2:12.6f} {3:11d} {4:14.6f}'.format(VBM,iv,CBM,ic,CBM-VBM)
          print '\n--------------band gap estimated from bands along high symmetry kpath-------------\n'


      def _plot_kpdos(self,args):
          import matplotlib
          matplotlib.use('Agg')
          import matplotlib.pyplot as plot
          import matplotlib.patches as mpatches
          import matplotlib.gridspec as gridspec

          try:
             import parse 
             struct = parse.parse_poscar('POSCAR')
          except:
             struct = parse_outcar_struct()
          dk = [np.linalg.norm(self._kpt_cart[ikpt] - self._kpt_cart[ikpt-1]) for ikpt in range(1,len(self._kpt))]
          path = np.concatenate(([0],[np.cumsum(dk)[i] for i in range(len(np.cumsum(dk)))]),axis=0)
          xsym=[path[0]]+[x for n, x in enumerate(path) if x in path[:n]]+[path[-1]]
          xticks,xlabels=refine_xsym_xlabels(xsym,args.label_k.split(),args.xlabel_thr)
          print '\n---------high symmetry k points-----------\n'
          print '{0:13s} {1:12s}'.format('xticks'.center(13),'xlabels'.center(20))
          print '\n'.join(['{0:12.6f} {1:20s}'.format(tick,label.strip('$').center(20)) for tick,label in zip(xticks,xlabels)])
          print '\n---------high symmetry k points-----------\n'
          print 'plot kpdos, proj_type=',args.proj_type,', proj_index=',args.proj_index
          if args.shift_fermi==True:
             self._energy += -self._efermi

          projections = args.proj_index.split(',')
          for proj in projections:
              if '_d' in proj or '_p' in proj:
                  idx = projections.index(proj)
                  head = proj.split('_')[0]
                  projections.remove(proj)
                  if '_p' in proj:
                     for i in range(1,4): projections.insert(idx+i-1,head+'_'+str(i))
                  if '_d' in proj:
                     for i in range(4,9): projections.insert(idx+i-4,head+'_'+str(i))

          cc = np.concatenate((np.array([0]),np.cumsum(struct._counts)),axis=0)
          if args.proj_type == 'species':
             #print [int(ispec) for ispec in args.proj_index.split(',')]
             #print [[iat for iat in range(cc[int(ispec)],cc[int(ispec)+1])] for ispec in args.proj_index.split(',')]
             weights_list = [sum([self._get_kpdos(iat) for iat in range(cc[int(ispec)],cc[int(ispec)+1])]) for ispec in args.proj_index.split(',')]
             legend_list  = [struct._symbols[cc[int(ispec)]] for ispec in args.proj_index.split(',')]
          elif args.proj_type == 'atom':
             weights_list = [self._get_kpdos(int(iat)) for iat in args.proj_index.split(',')]
             legend_list  = [struct._symbols[int(iat)]+str(filter(lambda x: x>=0, [int(iat)-ic for ic in cc])[-1]+1) for iat in args.proj_index.split(',')]
          elif args.proj_type == 'atom_subshell':
             def start_orb(isubshell):
                 if isubshell == 's': return 0
                 if isubshell == 'p': return 1
                 if isubshell == 'd': return 4
             def last_orb(isubshell):
                 if isubshell == 's': return 0
                 if isubshell == 'p': return 3
                 if isubshell == 'd': return 8
             at_list = [int(projection.split("_")[0]) for projection in projections]
             subshell_list  = [int(projection.split("_")[1]) for projection in projections]
             weights_list = [sum(self._get_kpdos_orb(int(iat))[start_orb(isubshell):last_orb(isubshell)+1],axis=3) for iat,isubshell in zip(at_list,subshell_list)]
             legend_list  = [struct._symbols[int(iat)]+str(filter(lambda x: x>=0, [int(iat)-ic for ic in cc])[-1]+1)+isubshell for iat,isubshell in zip(at_list,subshell_list)]
          elif args.proj_type=='orbital':
             at_list  = [int(projection.split('_')[0]) for projection in projections]
             orb_list = [int(projection.split('_')[1]) for projection in projections]
             weights_list = [self._get_kpdos_orb(iat)[:,:,:,iorb] for iat,iorb in zip(at_list,orb_list)]
             legend_list = [struct._symbols[iat]+str(iat)+'_$'+orb_dic[iorb]+'$' for iat,iorb in zip(at_list,orb_list)]
          elif args.proj_type == 'species_orbital':
             cc = np.concatenate((np.array([0]),np.cumsum(struct._counts)),axis=0)
             spec_list = [int(projection.split("_")[0]) for projection in projections]
             orb_list  = [int(projection.split("_")[1]) for projection in projections]
             weights_list = [sum([self._get_kpdos_orb(iat)[:,:,:,iorb] for iat in range(cc[int(ispec)],cc[int(ispec)+1])]) for ispec,iorb in zip(spec_list,orb_list)]
             legend_list = [struct._species[ispec]+'_$'+orb_dic[iorb]+'$' for ispec,iorb in zip(spec_list,orb_list)]
          elif args.proj_type == 'None':
               plot_band(args,path,self._energy,xsym,label_k=None,output='band_structure')
               exit('band structure without weights plotted')


          def plot_color_band(ax,ispin,weights_list,legend_list,args):
              color_patch=[]
              fw=open('kpdos.dat','w')
              for icolor,(weights,legend) in enumerate(zip(weights_list,legend_list)):
                  print '{0:20s}, spin {1}'.format(legend,spin_dic[ispin]),', nspin,nkpt,nband=',weights.shape
                  if args.pow != 1: weights = np.power(weights, args.pow)
                  for bi, wi in zip(self._energy[ispin,:,:].T, weights[ispin,:,:].T):
                      print >> fw, '\n'.join('{0:10.7f} {1:10.7f} {2:10.7f}'.format(path[ikpt],bi[ikpt],wi[ikpt]) for ikpt in range(self._nkpt,self._startk))
                      scolor=args.color.split()[icolor]
                      if args.spin_color:
                         scolor=args.color.split()[icolor*2+ispin]
                      ax.scatter(path, bi, s=args.markersize*wi, marker=args.marker,
                      facecolor='none', edgecolor=scolor, lw=args.linewidth,alpha=0.7)
                  color_patch.append(mpatches.Patch(color=scolor, label=legend))

              fw.close()
              ax.set_xticks(xticks)
              if args.yticks:
                 ax.set_yticks(args.yticks)
              [ax.plot([xi, xi], [args.elim[0],args.elim[1]], color='gray', zorder=-1) for xi in xsym]
              ax.set_xticklabels(['$'+ilabel+'$' for ilabel in xlabels],fontsize=args.label_fontsize,visible=(args.nplot[0]==1 or (not args.merge_spin or not(self._nspin-ispin-1))))
              ax.plot([np.min(path),np.max(path)], [(args.shift_fermi==False)*self._efermi, (args.shift_fermi==False)*self._efermi], ls='--',color='grey', zorder=-1)
              ax.set_xlim(np.min(path),np.max(path))
              ax.set_ylim(args.elim)
              ax.tick_params(top='off',right='off',direction='out')
              if args.yticks:
                 ax.set_yticks(args.yticks)
              [label.set_fontsize(args.label_fontsize) for label in ax.get_yticklabels()]
              ax.get_yaxis().set_tick_params(direction='out')
              ax.get_xaxis().set_tick_params(length=0.0)
              #ax.set_ylabel(args.ylabel,fontsize=14,visible=not args.merge_spin)
              if args.subtitles:
                 ax.set_title(args.subtitles.split(',')[ispin])
              return color_patch


          if args.legend_content:
             legend_list=[item for item in args.legend_content.split(',')]

          if args.merge_spin==False or self._nspin==1:
             print 'plot each spin band in a single file'
             for ispin in range(self._nspin):
                 fig = plot.figure(figsize=args.figsize)
                 ax=fig.add_subplot(111)
                 color_patch = plot_color_band(ax,ispin,weights_list,legend_list,args)
                 if args.legend_switch=='on':
                    plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],prop={'size':args.legend_fontsize})
                 ax.set_ylabel(args.ylabel,fontsize=args.label_fontsize)
                 output=args.proj_type+'_projected_bands'
                 if self._nspin>1: output = output + '_{0}'.format(spin_dic[ispin])
                 plot.tight_layout(pad=0.7,w_pad=1.5,h_pad=1.2)
                 plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
          elif args.merge_spin==True and self._nspin>1:
             print 'merge spin bands in one plot'
             fig = plot.figure(figsize=args.figsize)
             nrow,ncol=args.nplot
             gs=gridspec.GridSpec(nrow,ncol,height_ratios=[item for item in args.hratio],width_ratios=[item for item in args.wratio])
             ax=[]
             for ispin in range(self._nspin):
                 if nrow==1 and ispin>0:
                    ax.append(fig.add_subplot(gs[ispin],sharey=ax[0]))
                 else:
                    ax.append(fig.add_subplot(gs[ispin]))
                 color_patch=plot_color_band(ax[ispin],ispin,weights_list,legend_list,args)
                 ax[ispin].set_ylim(args.elim[0],args.elim[1])
             ax[0].set_ylabel(args.ylabel)
             for label in ax[1].get_yticklabels():
                 label.set_visible(args.nplot[0]!=1)

             #if nrow>1:
             #   plot.text(args.ylabel_pos[0],args.ylabel_pos[1],args.ylabel,fontsize=14,ha='center',va='center',rotation='vertical')
             if args.legend_switch=='on':
                plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],prop={'size':args.legend_fontsize})

             output=args.proj_type+'_projected_bands'
             plot.tight_layout(pad=0.7,w_pad=1.5,h_pad=1.2)
             plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
          print "Done"




if __name__=='__main__':
  import argparse
  import arguments
  desc_str='kpdos'
  parser = argparse.ArgumentParser(prog='plotprocar', description = desc_str)

  arguments.add_io_arguments(parser)
  arguments.add_fig_arguments(parser)
  arguments.add_plot_arguments(parser)
  args = parser.parse_args()
  parse_outcar_struct()
  #struct=parse.parse_poscar()
  arguments.check_args(args)

  kpdos=kpdos()
  kpdos._get_bandgap()
  kpdos._plot_kpdos(args)
  kpdos._plot_spinor_band(args)
