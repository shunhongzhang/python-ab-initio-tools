#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       bands.py                                                     #
#  Dependence: parse.py                                                     #
#  Usage:      convert the EIGENVAL file to a file for band structure plot  #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Sep 13, 2017                                                 #
#                                                                           #
#===========================================================================#

import os
import numpy as np


spin_dic={0:'up',1:'dn'}

def refine_xsym_xlabels(xsym,xlabels,xlabel_thr):
    klength=xsym[-1]-xsym[0]
    new_xsym=[]
    new_labels=[]
    i=len(xsym)-1
    while i>0:
        dk=xsym[i]-xsym[i-1]
        if dk<klength/xlabel_thr:
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


def parse_qe_band(filename='band.out'):
    lsoc  = 'T' if os.popen('grep spin-orbit {0}'.format(filename)).read() else 'F'
    lhf   = 'T' if 'HSE' in os.popen('grep HSE {0}|head -1'.format(filename)).read().split() else 'F'
    nspin = len(os.popen('grep SPIN {0}'.format(filename)).readlines())
    if nspin ==0: nspin=1
    nkpt  = int(os.popen('grep "k points" {0}'.format(filename)).read().split()[4])
    nband = int(os.popen('grep states {0} '.format(filename)).read().split()[4])
    nelect= float(os.popen('grep "number of electrons" {0}'.format(filename)).read().split()[-1])
    # efermi= float(os.popen('grep Fermi scf.out').read().split()[2])
    startline=int(os.popen('grep -n "End of band structure calculation" {0}'.format(filename)).read().split(':')[0])+1
    err_info='''set verbosity='high' to print the bands.'''
    if os.popen('grep "{0}" {1}'.format(err_info,filename)).read(): exit(err_info)
    f=open(filename)
    [f.readline() for i in range(startline)]
    energy=[]
    for ispin in range(nspin):
        if nspin==2: [f.readline() for i in range(3)]
        energy.append([])
        for ikpt in range(nkpt):
           f.readline()
           data=np.fromfile(f,dtype=float,count=nband,sep=' ')
           energy[ispin].append(data)
    energy=np.array(energy)
    alat=float(os.popen('grep alat {0}|grep lattice'.format(filename)).read().split()[4])
    get_kpts=os.popen('grep -A {0} "2pi/alat" {1}'.format(nkpt,filename)).readlines()[1:]
    kpt_cart = np.array([[float(kcoord.rstrip('\),')) for kcoord in kk.split()[4:7]] for kk in  get_kpts])*2*np.pi/alat
    kweight = np.array([float(item.split()[-1]) for item in get_kpts])
    return lsoc,lhf,nelect,kpt_cart,kweight,energy


def parse_outcar_band(filename='OUTCAR'):
    lsoc  = os.popen('grep LSORBIT {0}'.format(filename)).read().split()[2]
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

    if lhf=='T':
       startk=kweight.tolist().index(0) 
       kpt=kpt[startk:]
       kpt_cart=kpt_cart[startk:]
       kweight=kweight[startk:]
       energy=energy[:,startk:,:]
    return lsoc,lhf,nelect,efermi,kpt,kpt_cart,kweight,energy


    

class ebands(object):
      def __init__(self,fil='EIGENVAL',source='VASP'):
          self._source=source
          if source=='VASP':
             try:
                import parse
                self._lsoc,nspin,self._nelect,self._efermi,self._kpt,self._kweight,self._energy=parse.parse_eigenval(filename=fil)
                if self._lsoc=='T': self._energy=energy
                else: self._energy=np.rollaxis(energy,2)
                try: 
                   struct = parse.parse_poscar('POSCAR')
                   rec_cell=struct._reciprocal_cell()
                   print 'reciprocal lattice calculated from POSCAR'
                except:
                   get_rec_cell=os.popen('grep -A 3 "reciprocal lattice vectors" OUTCAR|tail -3').readlines()
                   rec_cell=np.matrix(([[float(coord) for coord in item.split()[3:]] for item in get_rec_cell]),float)
                   print 'reciprocal lattice extracted from OUTCAR'
                   self._kpt_cart=kpt*rec_cell
             except:
                self._lsoc,self._lhf,self._nelect,self._efermi,self._kpt,self._kpt_cart,self._kweight,self._energy=parse_outcar_band()
                #except: exit('we need parse to parse the EIGENVAL file')
          elif source=='QE':
             self._lsoc,self._lhf,self._nelect,self._kpt_cart,self._kweight,self._energy = parse_qe_band(fil)
             kweight=np.array([1/len(self._kpt_cart) for i in range(len(self._kpt_cart))])
             ffs=['scf.out','relax.out','rx.out']
             for ff in ffs:
                try:
                   self._efermi=float(os.popen('grep Fermi {0}|tail -1'.format(ff)).read().rstrip('\n').split()[4])
                   print 'Fermi energy = {0} eV, read from {1}, make sure it is right'.format(self._efermi,ff)
                   break
                except:
                   pass
          elif source=='Abinit':
             kpt,kweight,energy=qe_parse.parse_abinit_self._energy(fil) 
          else:
             exit('can not get data')
             
          self._nspin,self._nkpt,self._nband=self._energy.shape
          dk = [np.linalg.norm(self._kpt_cart[ikpt] - self._kpt_cart[ikpt-1]) for ikpt in range(1,self._nkpt)]
          self._path = np.concatenate(([0],[np.cumsum(dk)[i] for i in range(len(np.cumsum(dk)))]),axis=0)
          if source=='VASP':
             self._xsym = [self._path[0]]+[x for n, x in enumerate(self._path) if x in self._path[:n]]+[self._path[-1]]
          elif source=='QE':
             get_sym_kpts=[item for item in os.popen('grep -A 100 crystal_b band.in').readlines()][1:]
             nxsym=int(get_sym_kpts[0])
             segment_nkpt=[int(get_sym_kpts[i].split()[3]) for i in range(1,nxsym+1)][:-1]
             self._xsym = [self._path[0]]+[self._path[isym] for isym in np.cumsum(segment_nkpt)]
          print '\n------------------band data summary-----------------\n'
          if self._nspin==2: print 'spin polarized ',
          print 'band structure calculated by {0}'.format(self._source)
          if self._lsoc=='T': print 'spin-orbit coupling is included'
          if self._lhf=='T':  print 'Hartree-Fock type functional used?\nplease make sure the result is what you want.'
          print 'nspin = {0}'.format(self._nspin)
          print 'number of kpts  = {0}'.format(self._nkpt)
          print 'number of bands = {0}'.format(self._nband)
          print 'number of electrons = {0}'.format(self._nelect)
          print 'femri energy = {0} eV, read from OUTCAR \nplease make sure it is correct'.format(self._efermi)
          print '\n------------------band data summary-----------------\n'



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
              print '{0:12.6f} {1:11d} {2:12.6f} {3:11d} {4:14.6f}'.format(VBM,iv,CBM,ic,CBM-VBM),
              if iv-ic>=0: print ' metallic ? (VBM index >= CBM index)'
              else: print ''
          print '\n--------------band gap estimated from bands along high symmetry kpath-------------\n'


      def _write_eig(self,args):
          fw=open(args.wan_seedname+'.eig','w')
          for ispin in range(self._nspin):
              for ikpt in range(self._nkpt):
                  print>>fw, '\n'.join(['{0:12d}{1:12d}{2:22.12f}'.format(iband+1,ikpt+1,self._energy[ispin,ikpt,iband]) for iband in range(self._nband)])
          fw.close()

      def _write_dat(self):
          fw=open('band.dat','w')
          for iband in range(self._nband):
              for ikpt in range(self._nkpt):
                  print >>fw, '{0:12.7f}'.format(self._path[ikpt]),
                  print >>fw,' '.join(['{0:12.7f}'.format(en) for en in self._energy[:,ikpt,iband]])
              print >>fw,''
          fw.close()

      def _plot_band(self,args):
          import matplotlib.pyplot as plot
          import matplotlib.patches as mpatches
          from collections import Counter

          if args.plus_wan:
             data=np.fromfile(open(args.wan_seedname+'_band.dat'),sep=' ',dtype=float)
             data = data.reshape(len(data)/2,2).T
             self._wan_path = data[0]
             self._wan_energy = data[1] - args.wan_efermi
             self._wan_nband = len(Counter(self._wan_path).values())
             self._wan_nkpt = self._wan_energy.shape[0]/self._wan_nband
             if args.wan_bandpoint!=100:
                step=100/args.wan_bandpoint
                self._wan_path=self._wan_path[0::step]
                self._wan_energy=self._wan_energy[0::step]
          
          if args.shift_fermi: energy = self._energy - self._efermi
          else: energy = self._energy
          fig = plot.figure(figsize=args.figsize)
          ax = fig.add_subplot(1,1,1)
          paths=np.array([self._path for iband in range(self._nband)])
          for ispin in range(self._nspin):
              if args.style == 'line':
                 [ax.plot(self._path,energy[ispin,:,iband],color=args.color.split()[ispin],lw=args.linewidth,label='spin '+spin_dic[ispin]) for iband in range(self._nband)]
                 if args.plus_wan==True:
                    ax.scatter(self._wan_path,self._wan_energy,s=5,color=args.color.split()[ispin+2],lw=args.linewidth,label='wan spin '+spin_dic[ispin])
              if args.style == 'dot':
                 ax.scatter(paths, energy[ispin,:,:].T, s=args.markersize, marker=args.marker, facecolor='none',edgecolor=args.color.split()[ispin],label='spin '+spin_dic[ispin])
              color_patch=[]
              if self._nspin>1:
                 color_patch.append(mpatches.Patch(color=args.color.split()[0], label='spin up'))
                 color_patch.append(mpatches.Patch(color=args.color.split()[1], label='spin dn'))
              if args.legend_switch=='on': leg = ax.legend(handles=color_patch,loc='upper center',fontsize=10,numpoints=1)

          if args.elim[1]>args.elim[0]: plot.ylim(args.elim[0],args.elim[1])
          ymin,ymax = plot.ylim()
          if args.title:  ax.set_title(output)
          [ax.plot([xi, xi], [ymin, ymax], color='gray', zorder=-1) for xi in self._xsym]
          ax.plot([self._path[0],self._path[-1]], [0, 0], color='green', ls='dashed', zorder=-1)


          xlabels=['$'+label+'$' for label in args.label_k.split()]
          if len(xlabels)<len(self._xsym): exit('Error! No. of xlabels (={0}) < No. of xsym (={1})'.format(len(xlabels),len(self._xsym)))
          xticks,xlabels = refine_xsym_xlabels(self._xsym,xlabels,args.xlabel_thr)
          ax.set_xticks(xticks)
          ax.set_xticklabels(xlabels)
          [tick.set_visible(False) for tick in ax.get_xticklines()]

          print '\n--------------high symmetry k points--------------\n'
          print '{0:13s} {1:12s}'.format('xticks'.center(13),'xlabels'.center(20))
          print '\n'.join(['{0:12.6f} {1:20s}'.format(tick,label.strip('$').center(20)) for tick,label in zip(xticks,xlabels)])
          print '\n--------------high symmetry k points--------------\n'
          plot.xlim(0,self._path[-1])
          if args.lylabel==True: ax.set_ylabel(r'$E-E_f$ ($eV$)')
          plot.savefig('band_structure.png', dpi=args.dpi, bbox_inches='tight')



if __name__ == "__main__":
   import argparse
   import arguments

   desc_str = '''band structure of crystals.'''
   parser = argparse.ArgumentParser(prog='plotband.py', description = desc_str)
   parser.add_argument('--source',type=str,default='VASP',help='source file type: VASP/QE/Abinit')
   arguments.add_io_arguments(parser)
   arguments.add_fig_arguments(parser)
   arguments.add_plot_arguments(parser)
   arguments.add_wan_arguments(parser)
   args = parser.parse_args()
   #arguments.check_args(args)

   if args.source=='VASP': filename='EIGENVAL'
   #if args.source=='QE': filename=os.popen('grep filband *.in').read().rstrip('\n').split('=')[1].strip("'")
   elif args.source=='QE': filename='band.out'
   elif args.source=='Abinit': filename=os.popen('ls *EIG').read().split()[-1]
   eband=ebands(fil=filename,source=args.source)
   eband._plot_band(args)
   #parse.parse_abinit_struct()
   eband._write_eig(args)
   eband._get_bandgap()
   eband._write_dat()
