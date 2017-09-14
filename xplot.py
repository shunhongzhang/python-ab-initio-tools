#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       xplot.py                                                     #
#  Dependence: none                                                         #
#  Usage:      plot figures from input data                                 #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Jul 12, 2017                                                 #
#                                                                           #
#===========================================================================#

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches

orb_dic={0:'s',1:'p_y',2:'p_z',3:'p_x',
         4:'d_{xy}',5:'d_{yz}',6:'d_{z^2}',7:'d_{xz}',8:'d_{x^2-y^2}'}

spin_dic={0:'up',1:'dn'}

# Handle color arguments
def color(string):
    try:
        # Try to convert to tuple
        return eval(string)
    except:
        # If cannot convert to tuple just
        # return the original string
        return string

def plot_total_dos(args,energy,dos_total,dos_integral):
    nedos=energy.shape[0]
    nspin=dos_total.shape[1]
    fig=plot.figure(figsize=args.figsize)
    ax = fig.add_subplot(1,1,1)
    color_patch=[]
    #intdos_calc=[]  # test int dos
    #for i in range(1,nedos-1):
    #    intdos_calc.append(np.trapz(dos_total[1:i,0],x=energy[1:i]))
    #intdos_calc=np.array(intdos_calc)
    ax.plot(energy[1:-1], dos_total[1:-1,0], color='black', ls='-',lw=args.linewidth)
    if nspin==2:
       dos_total[:,1] = - dos_total[:,1]
       ax.plot(energy[1:-1], dos_total[1:-1,1], color='black', ls='-',lw=args.linewidth)
    color_patch.append(mpatches.Patch(color='black', label="Total"))
    if args.dos_intgl==True:
       ax.plot(energy[1:-1], dos_integral[1:-1,0], color='violet', ls='-',lw=args.linewidth)
       if nspin==2:
          dos_integral[:,1] = - dos_integral[:,1]
          ax.plot(energy[1:-1], dos_integral[1:-1,1], color='violet')
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
    print "Done"
    return ax

def plot_dos(args,energy,dos_orb,struct):
    nedos=energy.shape[0]
    nspin=dos_orb.shape[3]
    dos_list=[]
    legend_list=[]
    if args.proj_type == 'species':
       dos_atom = np.sum(dos_orb,axis=2)
       cc = np.concatenate((np.array([0]),np.cumsum(struct._counts)),axis=0)
       dos_species = np.array([np.sum(dos_atom[cc[ispec]:cc[ispec+1],:,:],axis=0) for ispec in range(len(struct._species))],float)
       dos_list = [dos_species[int(ispec),:,:] for ispec in args.proj_index.split(',')]
       legend_list = [struct._species[int(ispec)] for ispec in args.proj_index.split(',')]
    elif args.proj_type == 'atom':
       dos_atom = np.sum(dos_orb,axis=2)
       dos_list = [dos_atom[int(iat),:,:] for iat in args.proj_index.split(',')]
       legend_list = [struct._symbols[int(iat)]+iat for iat in args.proj_index.split(',')]
    elif args.proj_type == 'orbital':      #under test
       at_list  = [int(projection.split('_')[0]) for projection in args.proj_index.split(',')]
       orb_list = [int(projection.split('_')[1]) for projection in args.proj_index.split(',')]
       dos_list = [dos_orb[iat,:,iorb,:] for iat,iorb in zip(at_list,orb_list)]
       legend_list = [struct._symbols[iat]+str(iat)+'_$'+orb_dic[iorb]+'$' for iat,iorb in zip(at_list,orb_list)]
    elif args.proj_type == 'species_orbital':
       cc = np.concatenate((np.array([0]),np.cumsum(struct._counts)),axis=0)
       dos_species_orb = np.array([np.sum(dos_orb[cc[ispec]:cc[ispec+1],:,:,:],axis=0) for ispec in range(len(struct._species))])
       spec_list = [int(projection.split("_")[0]) for projection in args.proj_index.split(',')]
       orb_list  = [int(projection.split("_")[1]) for projection in args.proj_index.split(',')]
       dos_list = [dos_species_orb[ispec,:,iorb,:] for ispec,iorb in zip(spec_list,orb_list)]
       legend_list = [struct._species[ispec]+'_$'+orb_dic[iorb]+'$' for ispec,iorb in zip(spec_list,orb_list)]
    elif argfs.proj_type=='subshell':
       exit('uncoded')
       at_list  = [int(projection.split('_')[0]) for projection in args.proj_index.split(',')]
       orb_list = [int(projection.split('_')[1]) for projection in args.proj_index.split(',')]
       dos_list = [dos_species_orb[ispec,:,iorb,:] for ispec,iorb in zip(spec_list,orb_list)]
       legend_list = [struct._species[ispec]+'_$'+orb_dic[iorb]+'$' for ispec,iorb in zip(spec_list,orb_list)]


    if args.legend_content:
       legend_list=args.legend_content.split(',')

    if args.legend_switch=='off':
       legend_list=[]

    fig=plot.figure(figsize=args.figsize)
    ax = fig.add_subplot(1,1,1)
    color_patch=[]
    output=args.proj_type+'_projected_DOS'
    print 'plotting '+output+'...'
    for ispin in range(nspin):
        icolor=0
        for icolor,(dos,legend) in enumerate(zip(dos_list,legend_list)):
            sign=-2*ispin+1
            ax.plot(energy[1:-1], sign*dos[1:-1,ispin], color=args.color.split()[icolor], ls='-',lw=args.linewidth)
            if ispin==0:
               color_patch.append(mpatches.Patch(color=args.color.split()[icolor], label=legend))
    if args.legend_switch=='on': 
       plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],loc=0,fontsize=args.legend_fontsize)
    axes=plot.gca()
    ymin, ymax = axes.get_ylim()
    if args.ylim:
       ymin,ymax=args.ylim
    plot.ylim(ymin,ymax)
    plot.plot([0,0], [ymin,ymax], color='gray', ls='dashed', zorder=-1)
    # Fix x and y-axis boundaries
    if args.elim[0] < args.elim[1]:
       plot.xlim(args.elim[0], args.elim[1])
    plot.plot([plot.xlim()[0],plot.xlim()[1]], [0,0], color='gray', lw=0.5, zorder=-1)
    if nspin==1:
       plot.ylim(0,ymax)
    ax.set_xlabel('$E-E_f$ ($eV$)',fontsize=args.label_fontsize)
    ax.set_ylabel('$DOS$ ($states/eV/u.c.$)',fontsize=args.label_fontsize)
    plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
    print "Done"
    return ax


# Plot band dispersion along high symmetry k paths
def plot_band(args,path,eigenval,xsym,label_k=None,output='band_structure'):
    print "Plotting band structure ..."
    nband,nspin=eigenval.shape[1:]
    fig = plot.figure(figsize=args.figsize)
    ax = fig.add_subplot(1,1,1)
    for ispin in range(nspin):
        for iband in range(nband):
            if args.style == 'line':
               ax.plot(path[:], eigenval[:,iband,ispin], color=args.color.split()[ispin], ls='-',lw=args.linewidth)
            elif args.style == 'dot':
             ax.scatter(path[:], eigenval[:,iband,ispin], s=args.markersize, marker=args.marker, facecolor='none',edgecolor=args.color.split()[ispin])
    color_patch=[]
    if nspin>1:
       color_patch.append(mpatches.Patch(color='blue', label='spin up'))
       color_patch.append(mpatches.Patch(color= 'red', label='spin dn'))
       if args.legend_switch=='on':
          plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],loc=0,prop={'size':8})

    if args.elim[1]>args.elim[0]:
       plot.ylim(args.elim[0],args.elim[1])
    ymin,ymax = plot.ylim()
    if args.title:
       plot.title(output)
    if not label_k:
       label_k=['$'+label+'$' for label in args.label_k.split()]
    [ax.plot([xi, xi], [ymin, ymax], color='gray', zorder=-1) for xi in xsym]
    if label_k:
       ax.set_xticks(xsym)
       ax.set_xticklabels(label_k)
       print "high symmetric k points:",label_k
    plot.plot([path[0],path[-1]], [0, 0], color='green', ls='dashed', zorder=-1)
    plot.xlim(0,path[-1])
    if args.lylabel==True:
       if args.ylabel:
          ax.set_ylabel(args.ylabel)
       else:
          ax.set_ylabel(r"Energy (eV)")

    plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
    print "Done"


def plot_weighted_band(path,energy,weights,struct,orbitals,xsym,args,title=None):
    nkpt,natom,norb,nband,ndim,nspin=weights.shape
    norb=norb-1
    nspec=len(struct._species)
    if args.pow != 1:
       weights = np.power(weights, args.pow)

    if args.proj_type == 'species':
       cc = np.concatenate((np.array([0]),np.cumsum(struct._counts)),axis=0)
       weights_list = [np.sum(weights[:,cc[int(ispec)]:cc[int(ispec)+1],norb,:,0,:],axis=1) for ispec in args.proj_index.split(',')]
       legend_list  = [struct._symbols[cc[int(ispec)]] for ispec in args.proj_index.split(',')]
    elif args.proj_type == 'atom':
       weights_list = [weights[:,int(iat),norb,:,0,:] for iat in args.proj_index.split(',')]
       legend_list  = [struct._symbols[int(iat)]+str(iat) for iat in args.proj_index.split(',')]
    elif args.proj_type=='orbital':
       iat_list  = [int(projection.split('_')[0]) for projection in args.proj_index.split(',')]
       iorb_list = [int(projection.split('_')[1]) for projection in args.proj_index.split(',')]
       weights_list = [weights[:,iat,iorb,:,0,:] for iat,iorb in zip(iat_list,iorb_list)]
       legend_list = ['$'+struct._symbols[iat]+str(iat)+'$_$'+orb_dic[iorb]+'$' for iat,iorb in zip(iat_list,iorb_list)]

    '''
    for ispin in range(nspin):
        fig=plot.figure(figsize=args.figsize)
        ax = fig.add_subplot(1,1,1)
        color_patch=[]
        output=args.proj_type+'_projected_bands'
        if ndim==1:
           output+='_'+spin_dic[ispin]
        print 'plotting {0} ...'.format(output)
        for icolor,(weights,legend) in enumerate(zip(weights_list,legend_list)):
            for bi, wi in zip(energy[:,:,ispin].T, weights[:,:,ispin].T):
                ax.scatter(path, bi, s=args.markersize*wi, marker=args.marker,
                facecolor='none', edgecolor=args.color.split()[icolor], lw=args.linewidth)
            color_patch.append(mpatches.Patch(color=args.color.split()[icolor], label=legend))

        if args.legend_switch=='on':
           plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],loc=0,prop={'size':8})
        if args.elim[1]>args.elim[0]:
            plot.ylim(args.elim[0],args.elim[1])
        if args.label_k:
           ax.set_xticklabels(args.label_k)
        ymin,ymax = plot.ylim()
        plot.ylim(ymin,ymax)
        if title:
           plot.title(output)
        if args.label_k:
           plot_high_sym(ax,args,xsym,args.label_k.split(),ymin,ymax)
        ax.get_yaxis().set_tick_params(direction='out')
        plot.ylabel(args.ylabel)
        plot.plot([path[0],path[-1]], [0, 0], color='gray', ls='dashed', zorder=-1)
        plot.xlim(0,path[-1])
        plot.xticks(xsym, ['']*len(xsym))
        plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
        print "Done"

    '''
    nrows=2
    ncols=1
    if nspin==1:
       nrows=1
       ncols=1
    fig = plot.figure(figsize=args.figsize)
    '''
    ax = fig.add_subplot(111)
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
    ax.set_ylabel(args.ylabel)
    '''

    ymin,ymax=args.elim
    if args.legend_content:
       legend_list=[item for item in args.legend_content.split(',')]

    for ispin in range(nspin):
        ax = fig.add_subplot(2,1,ispin+1)
        color_patch=[]
        for icolor,(weights,legend) in enumerate(zip(weights_list,legend_list)):
            for bi, wi in zip(energy[:,:,ispin].T, weights[:,:,ispin].T):
                ax.scatter(path, bi, s=args.markersize*wi, marker=args.marker,
                facecolor='none', edgecolor=args.color.split()[icolor], lw=args.linewidth)
            color_patch.append(mpatches.Patch(color=args.color.split()[icolor], label=legend))

        ax.set_xticks(xsym)
        if args.yticks:
           ax.set_yticks(args.yticks)
        [ax.plot([xi, xi], [ymin, ymax], color='gray', zorder=-1) for xi in xsym]
        ax.set_xticklabels(['$'+ilabel+'$' for ilabel in args.label_k.split()],visible=ispin)
        ax.plot([np.min(path),np.max(path)], [0, 0], ls='--',color='grey', zorder=-1)
        ax.set_xlim(np.min(path),np.max(path))
        ax.set_ylim(args.elim)
        ax.get_yaxis().set_tick_params(direction='out')
        '''
        if ispin==0 and args.lylabel==True:
           if args.ylabel:
              plot.ylabel(args.ylabel)
           else:
              plot.ylabel(r"$E-E_f$ ($eV$)")
         '''

    plot.text(-0.8,0.2,args.ylabel,fontsize=14,ha='center',va='center',rotation='vertical')
    if args.legend_switch=='on':
       plot.legend(bbox_to_anchor=(args.legend_pos),handles=[item for item in color_patch],prop={'size':args.legend_fontsize})

    output=args.proj_type+'_projected_bands'
    plot.tight_layout(pad=0.7,w_pad=1.5,h_pad=1.2)
    plot.savefig(output, dpi=args.dpi, bbox_inches='tight')
    print "Done"
