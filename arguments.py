#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       arguments.py                                                 #
#  Dependence: none                                                         #
#  Usage:      process with arguments                                       #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Jul 31, 2017                                                 #
#                                                                           #
#===========================================================================#


import argparse

def add_io_arguments(parser):
    parser.add_argument('--poscar', type=str, default='POSCAR', help='POSCAR file')
    parser.add_argument('--procar', type=str, default='PROCAR', help='POSCAR file')
    parser.add_argument('--doscar', type=str, default='DOSCAR', help='DOSCAR file')
    parser.add_argument('--output', type=str, default='output_plot',
                        help='Filename for the plot. It accepts all image formats supported by matplotlib.')
    return parser


def add_fig_arguments(parser):
   parser.add_argument('--figsize', type=eval, default=(6, 4),
                       help='Figure size in inches formatted as width,height '
                       '(no spaces allowed). Default is 6,4')
   parser.add_argument('--dpi', type=float, default=600,
                       help='Plot resolution in dots per inch. Default is 300.')
   parser.add_argument('--subplot',type=bool,default=False,
                       help='whether or not to plot in subplots')
   parser.add_argument('--marker', type=str, default='o',
                       help='Marker for the fatband plot. Default is o.')
   parser.add_argument('--markersize', type=float, default=20,
                       help='Marker size. Default is 20.')
   parser.add_argument('--linewidth', type=float, default=0.5,
                       help='lien width of the DOS plot')
   parser.add_argument('--nplot',type=eval,default=(2,1),help='nrow and ncol of subplots') 
   parser.add_argument('--hratio',type=eval,default=(1,1),help='height ratios of subplots')
   parser.add_argument('--wratio',type=eval,default=(1,1),help='width ratios of subplots')

   return parser


def add_plot_arguments(parser):
    parser.add_argument('--title', type=str, default=None, help="title of plot")
    parser.add_argument('--subtitles',type=str,default=None,help='title for each subplot')
    parser.add_argument('--style', type=str, default='line', help='style to plot the data points')
    parser.add_argument('--color', type=str, default="red blue green orange violet cyan magenta black yellow",
                        help='Color for the marker. It accepts any color specification '
                        'accepted by matplotlib. Color specified as r,g,b tuple should ' 
                        'not have any spaces.')
    parser.add_argument('--elim', type=eval, default=(1,-1),
                        help='Energy range for the plot, specified as emin,emax '
                        '(no spaces allowed). Default is entire band range.')
    parser.add_argument('--efermi', type=float, default=0,
                        help='Fermi energy. Default is 0.')
    parser.add_argument('--shift_fermi',type=eval,default=True,help='shift the fermi energy to 0 or not?')
    parser.add_argument('--proj_type',type=str,default='species',
                        help='prjection type, the value can be: species, atoms, orbitals')
    parser.add_argument('--proj_index',type=str,default="0,1",
                        help='index of species/atoms/orbitals for projection')
    parser.add_argument('--proj_spec_orb_index',type=str,default=None,
                        help='for species_orbtials projection')
    parser.add_argument('--label_k', type=str, default=None,
                        help='labels for high symmetry k points'
                        'Use the form like \Gamma to represent Greek characters'
                        "example: --label_k '\Gamma X M \Gamma Y M'")
    parser.add_argument('--legend_switch',type=str,default='on',help='switch to control legends')
    parser.add_argument('--legend_content',type=str,default=None,help='customized legend from user')
    parser.add_argument('--legend_pos',type=eval,default=(0.9,0.9),help='left up position of legend box')
    parser.add_argument('--legend_fontsize',type=int,default=10,help='fontsize for legend texts')
    parser.add_argument('--pow', type=float, default=1,
                        help='Raise orbital weights to the specified integer power. '
                        'Powers larger than 1 help to filter out the ghost energy '
                        'in the unfolded band structures.')
    parser.add_argument('--dos_intgl',type=bool,default=False,
                        help='whether or not to plot the integrated DOS')
    parser.add_argument('--y_nticks',type=int,default=0,help='# of y ticks')
    parser.add_argument('--yticks',type=eval,default=None,help='costomerized yticks')
    parser.add_argument('--ylim',type=eval,default=None,
                        help='limit of the DOS in y axis')
    parser.add_argument('--lylabel',type=eval,default=True,help='write ylabel or not')
    parser.add_argument('--ylabel',type=str,default='$E-E_f$ ($eV$)',help='ylabel')
    parser.add_argument('--ylabel_pos',type=eval,default=(-0.5,0.5),help='ylabel position')
    parser.add_argument('--label_fontsize',type=int,default=10,help='font size for labels')
    parser.add_argument('--merge_spin',type=eval,default=True,help='merge plots of different spins or not?')
    parser.add_argument('--spin_color',type=eval,default=False,help='use different colors for spins')
    parser.add_argument('--proj_spinor',type=eval,default=False,help='projected on the sz weights, only used for noncollinear spin polarized calculations')
    parser.add_argument('--proj_spinor_index',type=int,default=-1,help='spinors of single atom or total,-1 for total')
    parser.add_argument('--spinor_dir',type=int,default=2,help='primary direction of spinor bands')
    parser.add_argument('--xlabel_thr',type=int,default=10,help='the threhold for refine the xlabels, read the code or consult the author for details')
    parser.add_argument('--int_pdos',type=eval,default=False,help='integrate pdos or not')
    parser.add_argument('--yannotates',type=eval,default=None,help='annotation for yaxis, read the code or consult the author for details')

    return parser

def add_wan_arguments(parser):
    parser.add_argument('--plus_wan',type=eval,default=False,help='plot the bands from wanniersation for comparison')
    parser.add_argument('--wan_seedname',type=str,default='wannier90',help='seedname for wannierisation')
    parser.add_argument('--wan_efermi',type=float,default=0,help='fermi energy for wannierisation')
    parser.add_argument('--wan_bandpoint',type=int,default=100,help='num of band points for wann band plot')


def plot_color_gradients(cmap_category, cmap_list):
    cmaps = [('Perceptually Uniform Sequential',
                                ['viridis', 'inferno', 'plasma', 'magma']),
             ('Sequential',     ['Blues', 'BuGn', 'BuPu',
                                 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                                 'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                                 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']),
             ('Sequential (2)', ['afmhot', 'autumn', 'bone', 'cool',
                                 'copper', 'gist_heat', 'gray', 'hot',
                                 'pink', 'spring', 'summer', 'winter']),
             ('Diverging',      ['BrBG', 'bwr', 'coolwarm', 'PiYG', 'PRGn', 'PuOr',
                                 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn', 'Spectral',
                                 'seismic']),
             ('Qualitative',    ['Accent', 'Dark2', 'Paired', 'Pastel1',
                                 'Pastel2', 'Set1', 'Set2', 'Set3']),
             ('Miscellaneous',  ['gist_earth', 'terrain', 'ocean', 'gist_stern',
                                 'brg', 'CMRmap', 'cubehelix',
                                 'gnuplot', 'gnuplot2', 'gist_ncar',
                                 'nipy_spectral', 'jet', 'rainbow',
                                 'gist_rainbow', 'hsv', 'flag', 'prism'])]


    nrows = max(len(cmap_list) for cmap_category, cmap_list in cmaps)
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    fig, axes = plt.subplots(nrows=nrows)
    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
    axes[0].set_title(cmap_category + ' colormaps', fontsize=14)

    for ax, name in zip(axes, cmap_list):
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axes:
        ax.set_axis_off()

def show_colormap(cmps):
    for cmap_category, cmap_list in cmaps:
       plot_color_gradients(cmap_category, cmap_list)
       plt.show()



def check_args(args):
    import os
    try: 
        fil='OUTCAR'
        counts=[int(item) for item in os.popen('grep "ions per type" {0}'.format(fil)).read().split()[4:]]
        species=[item.split('=')[1].split(':')[0] for item in os.popen('grep VRH {0}'.format(fil)).readlines()]
    except:
        species=raw_input('input name of each species:')
        counts=input('input number of ions for each species:')
    proj_type_list=['species','atom','orbital','species_orbital','atom_subshell','species_subshell','None']
    if args.proj_type not in proj_type_list:
       exit('proj_type {0} currently not supported! \
             \nValid proj_type are {1}'.format(args.proj_type,proj_type_list))
    if args.proj_type=='species':
       nspec=len(counts)
       for ispec in args.proj_index.split(','):
           if int(ispec) >= nspec:
              exit('error! the crystal has only {0} species, valid proj_indice are 0 to {1}!'.format(nspec,nspec-1))
    if args.elim != (1,-1) and args.elim[0]>=args.elim[1]:
       exit('error! elim[0] must be smaller than elim[1]!')
