#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       grid.py                                                      #
#  Dependence: arguments.py, parse.py, crystal_structure.py                 #
#  Usage:      define a class for real space grid in a crystal              #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Aug 20, 2017                                                 #
#                                                                           #
#===========================================================================#

import parse
import numpy as np
import argparse
import arguments


class grid(object):
      def __init__(self,fil='CHG'):
          self._type='charge'
          if fil=='ELFCAR': self._type='ELF'
          elif fil=='LOCPOT': self._type='potential'
          print 'grid file type: {0}'.format(self._type)
          self._struct,NGXF,NGYF,NGZF,self._data=parse.parse_chg(filename=fil)
          self._grid_den = np.array([NGXF,NGYF,NGZF])
          self._ngrid = np.prod(self._grid_den)
          if self._type=='charge': self._data=self._data/self._struct._real_cell_volume()
          self._grid_data = np.array(self._data).reshape(tuple(reversed(self._grid_den))).T

      def _avrg(self,idir=2):
          dim = [0,1,2]
          try:    dim.remove(idir)
          except: exit('Valid idir values: 0/1/2')
          latt=self._struct._latt_param()
          if idir==0: xlen=latt['a']
          elif idir==1: xlen=latt['b']
          elif idir==2: xlen=latt['c']
          xstep=xlen/self._grid_den[idir]
          xticks=[xstep*item for item in range(self._grid_den[idir])]
          return xticks,np.average(self._grid_data,axis=tuple(dim))

      def _write_avrg(self,idir=2):
          xticks,avrg_data=self._avrg(idir=idir)
          fw=open('average_'+fil+'_'+str(idir)+'.dat')
          print '\n'.join(['{0} {1}'.format(xx,yy) for (xx,yy) in zip(xticks,avrg_data)])
          fw.close()

      def _plot_avrg(self,args,idir=2):
          try: import matplotlib.pyplot as plt
          except: exit('cannot import matplotlib')
          xticks,avrg_data = self._avrg(idir)
          fig = plt.figure(figsize=args.figsize)
          ax = fig.add_subplot(111)
          ax.plot(xticks,avrg_data,ls='-',lw=args.linewidth,c='r')
          fig.tight_layout()
          fig.savefig(self._type+'_avrg_'+str(idir)+'.png',dpi=args.dpi)

      def _slice(self,idir,position):
          # currently can only plot slices spanned by two lattice vectors normal to each other
          try: import matplotlib.pyplot as plt
          except: exit('cannot import matplotlib')
          dim = [0,1,2]
          try:    dim.remove(idir)
          except: exit('Valid idir values: 0/1/2')
          for i in range(self._grid_den[idir]):
              if float(i)/self._grid_den[idir]>position:
                 cc=i
                 break
          if idir==0: slice_data = self._grid_data[cc-1:cc+1,:,:]
          elif idir==1: slice_data = self._grid_data[:,cc-1:cc+1,:]
          elif idir==2: slice_data = self._grid_data[:,:,cc-1:cc+1]
          w2 = (cc-position)/self._grid_den[idir]
          w1 = 1-w2
          slice_data = np.average(slice_data, axis=idir,weights=[w1,w2])
	  fig = plt.figure(figsize=args.figsize)
          ax = fig.add_subplot(111)
          ax.imshow(slice_data, extent=(0,10,0,10), cmap='hot', aspect='auto')
          fig.tight_layout()
          fig.savefig(self._type+'_slice_'+str(idir)+'.png',dpi=args.dpi)

      def _gen_grd(self):
          filename=self._type+'.grd'
          latt=self._struct._latt_param()
          fw=open(filename,"w")
          print>>fw,'grd for visualization in Materials Studio\n(1p,e12.5)'
          print>>fw,latt['a'],latt['b'],latt['c'],latt['alpha'],latt['beta'],latt['gamma']
          print>>fw,' '.join(['{0:5d}'.format(item-1) for item in self._grid_den])
          print>>fw,'1 0 {0:5d} 0 {1:5d} 0 {2:5d}'.format(*tuple(self._grid_den))
          print>>fw,'\n'.join(['{0:10.5f}'.format(item) for item in self._data])


if __name__=='__main__':

  desc_str = 'real space grid'
  parser = argparse.ArgumentParser(prog='readdos.py', description = desc_str)
  arguments.add_io_arguments(parser)
  arguments.add_fig_arguments(parser)
  arguments.add_plot_arguments(parser)
  parser.add_argument('--filgrid',type=str,default='CHG',help='the grid file name')
  parser.add_argument('--idir',type=int,default=2,help='direction normal to which the grid average is calculated, 0/1/2 for a/b/c crystal axes respectively')
  args = parser.parse_args()

  grd = grid(fil=args.filgrid)
  grd._gen_grd()
  grd._plot_avrg(args,idir=args.idir)
  grd._slice(args.idir,0.5)
