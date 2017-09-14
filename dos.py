#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       dos.py                                                       #
#  Dependence: none                                                         #
#  Usage:      electronic density of states                                 #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Sep 06, 2017                                                 #
#                                                                           #
#===========================================================================#


orb_dic={0:'s',
         1:'p_y',2:'p_z',3:'p_x',
         4:'d_{xy}',5:'d_{yz}',6:'d_{z^2}',7:'d_{xz}',8:'d_{x^2-y^2}'}

spin_dic={0:'up',1:'dn'}

def parse_projection(projection):
    print 'orbital index:'
    print orb_dic
    for proj in projection:
        for orb in orb_dic.values():
            if orb in proj:
               head=proj.split('_')[0]
               idx=orb_dic.values().index(orb)
               projection.remove(proj)
               projection.insert(idx,head+'_'+str(idx))
    for proj in projection:
        if ('_d' in proj and not '_d_'in proj) or ('_p' in proj and not '_p_' in proj):
           idx = projection.index(proj)
           head = proj.split('_')[0]
           projection.remove(proj)
           if '_p' in proj:
              for i in range(1,4): projection.insert(idx+i-1,head+'_'+str(i))
           if '_d' in proj:
              for i in range(4,9): projection.insert(idx+i-4,head+'_'+str(i))
    return projection

def parse_outcar_struct(fil='OUTCAR'):
    import os
    import numpy as np
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


def parse_doscar(filename='DOSCAR'):
    import os
    import numpy as np
    "Parsing {0} ...".format(filename)
    try:
       nspin = int(os.popen("grep ISPIN OUTCAR 2>/dev/null").readline().split()[2])
    except:
       nspin=input("nspin=")
    if nspin==2: print "Spin polarized DOS read"
    elif nspin==1: print "Spin unpolarized DOS read"
    lncl  = os.popen('grep LNONCOLLINEAR {0}'.format('OUTCAR')).read().split()[2]
    ndim=(lncl=='T')*3+1


    f=open(filename)
    natom=int(f.readline().split()[0])
    [f.readline() for i in range(3)]
    system=f.readline()
    line = f.readline().split()
    emax,emin,nedos,efermi=[float(item) for item in line[:-1]]
    nedos=int(nedos)
    norb=9       #Note: for elements that contains f orbitals this number should be changed accordingly (typically it should be 16)

    energy=np.zeros(nedos,float)
    dos_total=np.zeros((nedos,nspin),float)
    dos_integral=np.zeros((nedos,nspin),float)
    dos_orb=np.zeros((natom,nedos,norb,nspin),float)

    print "NEDOS=",nedos
    print "E-fermi=",efermi,"eV (shifted to 0 in the plot)"
    for ie in range(nedos):
       line=f.readline().split()
       energy[ie]=float(line[0])
       dos_total[ie,0]=float(line[1])
       try:
          dos_total[ie,1]=float(line[2])
          dos_integral[ie,0]=float(line[3])
          dos_integral[ie,1]=float(line[4])
       except:
          dos_integral[ie,0]=float(line[2])

    try:
       print "Parsing PDOS..."
       for iat in range(natom):
           line=f.readline()
           data = np.fromfile(f,sep=' ',count=nedos*(ndim*norb*nspin+1)).reshape(nedos,ndim*norb*nspin+1)
           for ispin in range(nspin):
               dos_orb[iat,:,:,ispin] = data[:,ispin+1+(ndim-1)::nspin*ndim]
    except:
       print "PDOS not read!"
    f.close()

    energy=energy-efermi
    return lncl,energy,emin,emax,dos_total,dos_integral,dos_orb


class edos(object):
      def __init__(self,source='VASP'):
         try:
            import parse
         except:
            exit('we need the parse.py file!')
         filename = 'DOSCAR'
         if source == 'QE':
            filename = prefix+'.dos'
            exit('under development')
         elif source == 'OPENMX':
            filename = 'xxx'
            exit('under development')
         self._lncl,energy,emin,emax,dos_total,dos_integral,dos_orb = parse_doscar(filename)
         self._energy =  energy
         self._emin   =  emin
         self._emax   =  emax
         self._tdos   =  dos_total
         self._nspin  =self._tdos.shape[1]
         self._intdos =  dos_integral
         self._pdos   =  dos_orb
         self._nedos=energy.shape[0]

      def _species_pdos(self):
          return 1

      def _atom_pdos(self):
          return np.sum(self._pdos,axis=2)

      def _species_orb_pdos(self):
          import numpy as np
          try:
             import parse
             struct = parse.parse_poscar('POSCAR')
             species=struct._species
             counts=struct._counts
          except:
             print "POSCAR not found!Please input the following information manually:"
             species=raw_input("Please input the chemical species in the same order as that in POSCAR, separated by space:")
             counts=raw_input("please input the number of each kind of atoms, in the same order as that in POSCAR, separated by space:")
          cc = np.concatenate((np.array([0]),np.cumsum(counts)),axis=0)
          dos_species_orb = np.array([np.sum(self._pdos[cc[ispec]:cc[ispec+1],:,:,:],axis=0) for ispec in range(len(species))])
          return dos_species_orb

      def _integrated_dos(self,iat,iorb,ispin):
          intdos_calc=[]  # test int dos
          for ie in range(1,self._nedos-1):
              if self._energy[ie] <= 0:
                 intdos_calc.append(np.trapz(self._pdos[iat,1:ie,iorb,ispin],x=self._energy[1:ie]))
              else:
                 break
          intdos_calc=np.array(intdos_calc)
          return intdos_calc

      def _write_tdos(self):
          fw=open('Total_DOS.dat','w')
          print >>fw, '\n'.join(['{0:12.6f}  '.format(e)+' '.join(['{0:12.6f}'.format(t) for t in tdos]) 
                + ' '.join(['{0:12.6f}'.format(it) for it in intdos]) for e,tdos,intdos in zip(self._energy,self._tdos,self._intdos)])
          fw.close()

      def _plot_total_dos(self,args):
          import numpy as np
          import matplotlib.pyplot as plot
          import matplotlib.patches as mpatches
          nspin=self._tdos.shape[1]
          print 'plotting PDOS'
          fig=plot.figure(figsize=args.figsize)
          ax = fig.add_subplot(1,1,1)
          color_patch=[]
          [ax.plot(self._energy[1:-1], -np.sign(ispin-0.5)*self._tdos[1:-1,ispin], color='black', ls='-',lw=args.linewidth) for ispin in range(self._nspin)]
          color_patch.append(mpatches.Patch(color='black', label="Total"))
          if args.dos_intgl==True:
             [ax.plot(energy[1:-1], -np.sign(ispin-0.5)*self._intdos[1:-1,0], color='violet', ls='-',lw=args.linewidth) for ispin in range(self._nspin)]
          axes=plot.gca()
          ymin, ymax = axes.get_ylim()
          plot.ylim(ymin,ymax)
          plot.plot([0,0], [ymin,ymax], color='gray', ls='dashed', zorder=-1)
          plot.plot([plot.xlim()[0],plot.xlim()[1]], [0,0], color='gray', lw=0.5, zorder=-1)
          if args.elim[0] < args.elim[1]:
             plot.xlim(args.elim[0], args.elim[1])
          else:
             plot.xlim(min(energy),max(energy))
          if nspin==1:
             plot.ylim(0)
          ax.set_xlabel('$E-E_f$ ($eV$)')
          ax.set_ylabel('$Total DOS$ ($states/eV/u.c.$)')
          plot.savefig('total_dos', dpi=args.dpi, bbox_inches='tight')
          fig=plot.figure(figsize=args.figsize)
          ax = fig.add_subplot(111)
          ax.plot(self._energy[1:-1],np.sum(self._tdos[1:-1,:],axis=1))
          ax.set_xlabel('$E-E_f$ ($eV$)')
          ax.set_ylabel('$Total DOS$ ($eV^{-1}$)')
          ax.set_xlim(args.elim[0],args.elim[1])
          plot.savefig('total_dos_sum_spin', dpi=args.dpi, bbox_inches='tight')
          return ax

      def _plot_pdos(self,args):
          import numpy as np
          import matplotlib.pyplot as plot
          import matplotlib.patches as mpatches
          try:
             import parse
             struct=parse.parse_poscar()
          except:
             struct=parse_outcar_struct()

          nedos=self._energy.shape[0]
          dos_orb=self._pdos
          nspin=dos_orb.shape[3]
          dos_list=[]
          intdos_list=[]
          legend_list=[]

          projections = parse_projection(args.proj_index.split(','))
          cc = np.concatenate((np.array([0]),np.cumsum(struct._counts)),axis=0)

          def start_orb(isubshell):
              if isubshell == 's': return 0
              if isubshell == 'p': return 1
              if isubshell == 'd': return 4
          def last_orb(isubshell):
              if isubshell == 's': return 0
              if isubshell == 'p': return 3
              if isubshell == 'd': return 8


          if args.proj_type == 'species':
             dos_atom = np.sum(dos_orb,axis=2)
             dos_species = np.array([np.sum(dos_atom[cc[ispec]:cc[ispec+1],:,:],axis=0) for ispec in range(len(struct._species))],float)
             dos_list = [dos_species[int(ispec),:,:] for ispec in args.proj_index.split(',')]
             legend_list = [struct._species[int(ispec)] for ispec in args.proj_index.split(',')]
          elif args.proj_type == 'atom':
             dos_atom = np.sum(dos_orb,axis=2)
             dos_list = [dos_atom[int(iat),:,:] for iat in args.proj_index.split(',')]
             legend_list = [struct._symbols[int(iat)]+str(filter(lambda x: x>=0, [int(iat)-ic for ic in cc])[-1]+1) for iat in args.proj_index.split(',')]
          elif args.proj_type == 'atom_subshell':
             at_list = [int(projection.split("_")[0]) for projection in args.proj_index.split(',')]
             subshell_list  = [projection.split("_")[1] for projection in args.proj_index.split(',')]
             dos_list = [np.sum(self._pdos[int(iat),:,start_orb(isubshell):last_orb(isubshell)+1,:],axis=1) for iat,isubshell in zip(at_list,subshell_list)]
             legend_list  = [struct._symbols[int(iat)]+str(filter(lambda x: x>=0, [int(iat)-ic for ic in cc])[-1]+1)+'_'+isubshell for iat,isubshell in zip(at_list,subshell_list)]
          elif args.proj_type == 'orbital':      #under test
             at_list  = [int(projection.split('_')[0]) for projection in projections]
             orb_list = [int(projection.split('_')[1]) for projection in projections]
             dos_list = [dos_orb[iat,:,iorb,:] for iat,iorb in zip(at_list,orb_list)]
             legend_list = [struct._symbols[int(iat)]+str(filter(lambda x: x>=0, [int(iat)-ic for ic in cc])[-1]+1)+'_$'+orb_dic[iorb]+'$' for iat,iorb in zip(at_list,orb_list)]
          elif args.proj_type == 'species_orbital':
             spec_list = [int(projection.split("_")[0]) for projection in projections]
             orb_list  = [int(projection.split("_")[1]) for projection in projections]
             dos_list = [self._species_orb_pdos()[ispec,:,iorb,:] for ispec,iorb in zip(spec_list,orb_list)]
             legend_list = [struct._species[ispec]+'_$'+orb_dic[iorb]+'$' for ispec,iorb in zip(spec_list,orb_list)]
          elif args.proj_type=='species_subshell':
             spec_list  = [int(projection.split('_')[0]) for projection in args.proj_index.split(',')]
             subshell_list = [projection.split('_')[1] for projection in args.proj_index.split(',')]
             dos_list = [np.sum(self._species_orb_pdos()[ispec,:,start_orb(isubshell):last_orb(isubshell),:],axis=1) for ispec,isubshell in zip(spec_list,subshell_list)]
             legend_list = [struct._species[ispec]+'_'+isubshell for ispec,isubshell in zip(spec_list,subshell_list)]
          if args.legend_content: legend_list=args.legend_content.split(',')
          if args.legend_switch=='off': legend_list=[]

          for ilist in range(len(dos_list)):
              intdos_list.append([])
              if args.int_pdos:
                 intdos_list[ilist] = np.array([[np.trapz(dos_list[ilist][1:ie,ispin],x=self._energy[1:ie]) for ispin in range(self._nspin)] for ie in range(0,nedos)])

          fw=open('pdos.dat','w')
          for ie in range(self._nedos):
              print >>fw,'{0:10.6f}'.format(self._energy[ie]),
              for item in dos_list:
                  print >>fw, ' '.join(['{0:10.6f}'.format(item[ie,ispin]) for ispin in range(self._nspin)]),
              print >>fw,''
          fw.close()

          fig=plot.figure(figsize=args.figsize)
          ax = fig.add_subplot(1,1,1)
          color_patch=[]
          output=args.proj_type+'_projected_DOS'
          print 'plotting '+output+'...'
          for icolor,(dos,intdos,legend) in enumerate(zip(dos_list,intdos_list,legend_list)):
              print 'plotting {0:20s}'.format(legend)
              if self._lncl=='T':
                 [ax.plot(self._energy[1:-1], dos[1:-1,ispin], color=args.color.split()[icolor], ls='-',lw=args.linewidth,alpha=0.7) for ispin in range(self._nspin)]
              else:
                 [ax.plot(self._energy[1:-1], -np.sign(ispin-0.5)*dos[1:-1,ispin], color=args.color.split()[icolor], ls='-',lw=args.linewidth) for ispin in range(self._nspin)]
                 if args.int_pdos:
                    [ax.plot(self._energy[1:-1], -np.sign(ispin-0.5)*intdos[1:-1,ispin], color=args.color.split()[icolor], ls='--',lw=args.linewidth) for ispin in range(self._nspin)]
              color_patch.append(mpatches.Patch(color=args.color.split()[icolor], label=legend))
          if args.yticks: ax.set_yticks(args.yticks)
          if args.legend_switch=='on':  plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],loc=0,fontsize=args.legend_fontsize)

          if args.yannotates:
             delta=0.3
             for yannotate in args.yannotates:
                 ax.annotate(str(yannotate),xy=(0,yannotate),xytext=(0+delta,yannotate), arrowprops=dict(arrowstyle='-'),ha='center',va='center',fontsize=9)

          axes=plot.gca()
          ymin, ymax = axes.get_ylim()
          if args.ylim:  ymin,ymax=args.ylim
          plot.ylim(ymin,ymax)
          plot.plot([0,0], [ymin,ymax], color='gray', ls='dashed', zorder=-1)
          if args.elim[0] < args.elim[1]:
             plot.xlim(args.elim[0], args.elim[1])
          plot.plot([plot.xlim()[0],plot.xlim()[1]], [0,0], color='gray', lw=0.5, zorder=-1)
          if self._nspin==1 and self._lncl=='F':  plot.ylim(0,ymax)
          ax.set_xlabel('$E-E_f$ ($eV$)',fontsize=args.label_fontsize)
          ax.set_ylabel('$DOS$ ($states/eV/u.c.$)',fontsize=args.label_fontsize)
          plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
          return ax




class ph_dos(object):
      def __init__(self,source='phonopy'):
          filename=prefix+'.DOS'
          if source == 'QE':
             filename = 'xxx.dos'
          elif source == 'CASTEP':
             filename = 'xxx'
             exit('under development')


if __name__=='__main__':
   import argparse
   import arguments
   desc_str='dos'
   parser = argparse.ArgumentParser(prog='readdos.py', description = desc_str)
   arguments.add_io_arguments(parser)
   arguments.add_fig_arguments(parser)
   arguments.add_plot_arguments(parser)
   args = parser.parse_args()
   arguments.check_args(args)

   edos=edos()
   edos._plot_total_dos(args)
   edos._write_tdos()
   edos._plot_pdos(args)
   #dos=phdos()
