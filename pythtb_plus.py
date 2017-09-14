#!/usr/bin/python

#------------------------------------------------------------------------------
# PythTB_plus python tight binding module.
# Version 0.0.1, Sep 14, 2017
# by Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>, Tsinghua University
#
# This file is an extension module of PythTB (v1.7.0), 
# a free python code by by Sinisa Coh and David Vanderbilt
#------------------------------------------------------------------------------


#---------------------------------------------------------------------------------
# This module is a modification of pythtb under the terms of the GNU General
# Public License as published by the Free Software Foundation.
#
# This module is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# A copy of the GNU General Public License should be available
# alongside this source in a file named gpl-3.0.txt.  If not,
# see <http://www.gnu.org/licenses/>.
#
# PythTB is availabe at http://www.physics.rutgers.edu/pythtb/
#---------------------------------------------------------------------------------



try: 
   from pythtb import *
except: 
   ImportError 
   exit('you need to install pythtb first!')
import sys
import copy
from math import *

def display(my_model):
        r"""
        Prints on the screen some information about this tight-binding
        model. This function doesn't take any parameters.
        """
        print '---------------------------------------'
        print 'report of tight-binding model'
        print '---------------------------------------'
        print 'k-space dimension           =',my_model._dim_k
        print 'r-space dimension           =',my_model._dim_r
        print 'number of spin components   =',my_model._nspin
        print 'periodic directions         =',my_model._per
        print 'number of orbitals          =',my_model._norb
        print 'number of electronic states =',my_model._nsta
        print 'lattice vectors:'
        for i,o in enumerate(my_model._lat):
            print " #",_nice_int(i,2)," ===>  [",
            for j,v in enumerate(o):
                print _nice_float(v,7,4),
                if j!=len(o)-1:
                    print ",",
            print "]"
        print 'positions of orbitals:'
        for i,o in enumerate(my_model._orb):
            print " #",_nice_int(i,2)," ===>  [",
            for j,v in enumerate(o):
                print _nice_float(v,7,4),
                if j!=len(o)-1:
                    print ",",
            print "]"
        print 'site energies:'
        for i,site in enumerate(my_model._site_energies):
            print " #",_nice_int(i,2)," ===>  ",
            if my_model._nspin==1:
                print _nice_float(site,7,4)
            elif my_model._nspin==2:
                print str(site).replace("\n"," ")
        print 'hoppings:'
        for i,hopping in enumerate(my_model._hoppings):
            print "<",_nice_int(hopping[1],2),"| H |",_nice_int(hopping[2],2),
            if len(hopping)==4:
                print "+ [",
                for j,v in enumerate(hopping[3]):
                    print _nice_int(v,2),
                    if j!=len(hopping[3])-1:
                        print ",",
                    else:
                        print "]",
            print ">     ===> ",
            if my_model._nspin==1:
                print _nice_complex(hopping[0],7,4)
            elif my_model._nspin==2:
                #print str(hopping[0]).replace("\n"," ")
                print ' '.join([_nice_complex(item,8,4) for item in hopping[0].flatten()])
        print

def visualize(my_model,dir_first,dir_second=None,eig_dr=None,draw_hoppings=True,ph_color="black"):
        r"""

        Rudimentary function for visualizing tight-binding model geometry,
        hopping between tight-binding orbitals, and electron eigenstates.

        If eigenvector is not drawn, then orbitals in home cell are drawn
        as red circles, and those in neighboring cells are drawn with
        different shade of red. Hopping term directions are drawn with
        green lines connecting two orbitals. Origin of unit cell is
        indicated with blue dot, while real space unit vectors are drawn
        with blue lines.

        If eigenvector is drawn, then electron eigenstate on each orbital
        is drawn with a circle whose size is proportional to wavefunction
        amplitude while its color depends on the phase. There are various
        coloring schemes for the phase factor; see more details under
        *ph_color* parameter. If eigenvector is drawn and coloring scheme
        is "red-blue" or "wheel", all other elements of the picture are
        drawn in gray or black.

        :param dir_first: First index of Cartesian coordinates used for
          plotting.

        :param dir_second: Second index of Cartesian coordinates used for
          plotting. For example if dir_first=0 and dir_second=2, and
          Cartesian coordinates of some orbital is [2.0,4.0,6.0] then it
          will be drawn at coordinate [2.0,6.0]. If dimensionality of real
          space (*dim_r*) is zero or one then dir_second should not be
          specified.

        :param eig_dr: Optional parameter specifying eigenstate to
          plot. If specified, this should be one-dimensional array of
          complex numbers specifying wavefunction at each orbital in
          the tight-binding basis. If not specified, eigenstate is not
          drawn.

        :param draw_hoppings: Optional parameter specifying whether to
          draw all allowed hopping terms in the tight-binding
          model. Default value is True.

        :param ph_color: Optional parameter determining the way
          eigenvector phase factors are translated into color. Default
          value is "black". Convention of the wavefunction phase is as
          in convention 1 in section 3.1 of :download:`notes on
          tight-binding formalism  <misc/pythtb-formalism.pdf>`.  In
          other words, these wavefunction phases are in correspondence
          with cell-periodic functions :math:`u_{n {\bf k}} ({\bf r})`
          not :math:`\Psi_{n {\bf k}} ({\bf r})`.

          * "black" -- phase of eigenvectors are ignored and wavefunction
            is always colored in black.

          * "red-blue" -- zero phase is drawn red, while phases or pi or
            -pi are drawn blue. Phases in between are interpolated between
            red and blue. Some phase information is lost in this coloring
            becase phase of +phi and -phi have same color.

          * "wheel" -- each phase is given unique color. In steps of pi/3
            starting from 0, colors are assigned (in increasing hue) as:
            red, yellow, green, cyan, blue, magenta, red.

        :returns:
          * **fig** -- Figure object from matplotlib.pyplot module
            that can be used to save the figure in PDF, EPS or similar
            format, for example using fig.savefig("name.pdf") command.
          * **ax** -- Axes object from matplotlib.pyplot module that can be
            used to tweak the plot, for example by adding a plot title
            ax.set_title("Title goes here").

        Example usage::

          # Draws x-y projection of tight-binding model
          # tweaks figure and saves it as a PDF.
          (fig, ax) = tb.visualize(0, 1)
          ax.set_title("Title goes here")
          fig.savefig("model.pdf")

        See also these examples: :ref:`edge-example`,
        :ref:`visualize-example`.

        """

        # check the format of eig_dr
        if not (eig_dr is None):
            if eig_dr.shape!=(my_model._norb,):
                raise Exception("\n\nWrong format of eig_dr! Must be array of size norb.")
        
        # check that ph_color is correct
        if ph_color not in ["black","red-blue","wheel"]:
            raise Exception("\n\nWrong value of ph_color parameter!")

        # check if dir_second had to be specified
        if dir_second==None and my_model._dim_r>1:
            raise Exception("\n\nNeed to specify index of second coordinate for projection!")

        # start a new figure
        import pylab as plt
        fig=plt.figure(figsize=[plt.rcParams["figure.figsize"][0],
                                plt.rcParams["figure.figsize"][0]])
        ax=fig.add_subplot(111, aspect='equal')

        def proj(v):
            "Project vector onto drawing plane"
            coord_x=v[dir_first]
            if dir_second==None:
                coord_y=0.0
            else:
                coord_y=v[dir_second]
            return [coord_x,coord_y]

        def to_cart(red):
            "Convert reduced to Cartesian coordinates"
            return np.dot(red,my_model._lat)

        # define colors to be used in plotting everything
        # except eigenvectors
        if (eig_dr is None) or ph_color=="black":
            c_cell="b"
            c_orb="r"
            c_nei=[0.85,0.65,0.65]
            c_hop="g"
        else:
            c_cell=[0.4,0.4,0.4]
            c_orb=[0.0,0.0,0.0]
            c_nei=[0.6,0.6,0.6]
            c_hop=[0.0,0.0,0.0]
        # determine color scheme for eigenvectors
        def color_to_phase(ph):
            if ph_color=="black":
                return "k"
            if ph_color=="red-blue":
                ph=np.abs(ph/np.pi)
                return [1.0-ph,0.0,ph]
            if ph_color=="wheel":
                if ph<0.0:
                    ph=ph+2.0*np.pi
                ph=6.0*ph/(2.0*np.pi)
                x_ph=1.0-np.abs(ph%2.0-1.0)
                if ph>=0.0 and ph<1.0: ret_col=[1.0 ,x_ph,0.0 ]
                if ph>=1.0 and ph<2.0: ret_col=[x_ph,1.0 ,0.0 ]
                if ph>=2.0 and ph<3.0: ret_col=[0.0 ,1.0 ,x_ph]
                if ph>=3.0 and ph<4.0: ret_col=[0.0 ,x_ph,1.0 ]
                if ph>=4.0 and ph<5.0: ret_col=[x_ph,0.0 ,1.0 ]
                if ph>=5.0 and ph<=6.0: ret_col=[1.0 ,0.0 ,x_ph]
                return ret_col

        # draw origin
        ax.plot([0.0],[0.0],"o",c=c_cell,mec="w",mew=0.0,zorder=7,ms=4.5)

        # first draw unit cell vectors which are considered to be periodic
        for i in my_model._per:
            # pick a unit cell vector and project it down to the drawing plane
            vec=proj(my_model._lat[i])
            ax.plot([0.0,vec[0]],[0.0,vec[1]],"--",c=c_cell,lw=1.5,zorder=7)


        # now draw all orbitals
        for i in range(my_model._norb):
            # find position of orbital in cartesian coordinates
            pos=to_cart(my_model._orb[i])
            pos=proj(pos)
            ax.plot([pos[0]],[pos[1]],"o",c=c_orb,mec="w",mew=0.0,zorder=100,ms=4.0)

        # draw hopping terms
        if draw_hoppings==True:
            for h in my_model._hoppings:
                # draw both i->j+R and i-R->j hop
                for s in range(2):
                    # get "from" and "to" coordinates
                    pos_i=np.copy(my_model._orb[h[1]])
                    pos_j=np.copy(my_model._orb[h[2]])
                    # add also lattice vector if not 0-dim
                    if my_model._dim_k!=0:
                        if s==0:
                            pos_j[my_model._per]=pos_j[my_model._per]+h[3][my_model._per]
                        if s==1:
                            pos_i[my_model._per]=pos_i[my_model._per]-h[3][my_model._per]
                    # project down vector to the plane
                    pos_i=np.array(proj(to_cart(pos_i)))
                    pos_j=np.array(proj(to_cart(pos_j)))
                    # add also one point in the middle to bend the curve
                    prcnt=0.05 # bend always by this ammount
                    pos_mid=(pos_i+pos_j)*0.5
                    dif=pos_j-pos_i # difference vector
                    orth=np.array([dif[1],-1.0*dif[0]]) # orthogonal to difference vector
                    orth=orth/np.sqrt(np.dot(orth,orth)) # normalize
                    pos_mid=pos_mid+orth*prcnt*np.sqrt(np.dot(dif,dif)) # shift mid point in orthogonal direction
                    # draw hopping
                    all_pnts=np.array([pos_i,pos_mid,pos_j]).T
                    ax.plot(all_pnts[0],all_pnts[1],"-",c=c_hop,lw=0.75,zorder=8)
                    # draw "from" and "to" sites
                    ax.plot([pos_i[0]],[pos_i[1]],"o",c=c_nei,zorder=9,mew=0.0,ms=4.0,mec="w")
                    ax.plot([pos_j[0]],[pos_j[1]],"o",c=c_nei,zorder=9,mew=0.0,ms=4.0,mec="w")

        # now draw the eigenstate
        if not (eig_dr is None):
            for i in range(my_model._norb):
                # find position of orbital in cartesian coordinates
                pos=to_cart(my_model._orb[i])
                pos=proj(pos)
                # find norm of eigenfunction at this point
                nrm=(eig_dr[i]*eig_dr[i].conjugate()).real
                # rescale and get size of circle
                nrm_rad=2.0*nrm*float(my_model._norb)
                # get color based on the phase of the eigenstate
                phase=np.angle(eig_dr[i])
                c_ph=color_to_phase(phase)
                ax.plot([pos[0]],[pos[1]],"o",c=c_ph,mec="w",mew=0.0,ms=nrm_rad,zorder=11,alpha=0.8)

        # center the image
        #  first get the current limit, which is probably tight
        xl=ax.set_xlim()
        yl=ax.set_ylim()
        # now get the center of current limit
        centx=(xl[1]+xl[0])*0.5
        centy=(yl[1]+yl[0])*0.5
        # now get the maximal size (lengthwise or heightwise)
        mx=max([xl[1]-xl[0],yl[1]-yl[0]])
        # set new limits
        extr=0.05 # add some boundary as well
        ax.set_xlim(centx-mx*(0.5+extr),centx+mx*(0.5+extr))
        ax.set_ylim(centy-mx*(0.5+extr),centy+mx*(0.5+extr))

        # return a figure and axes to the user
        return (fig,ax)


def dist_hop(w90_model, dim=3):
        """

        This is one of the diagnostic tools that can be used to help
        in determining *min_hopping_norm* and *max_distance* parameter in
        :func:`pythtb.w90.model` function call.

        This function returns all hopping terms (from orbital *i* to
        *j+R*) as well as the distances between the *i* and *j+R*
        orbitals.  For well localized Wannier functions hopping term
        should decay exponentially with distance.

        :param: w90_model, an object subjected to the w90 class
        :param: dim, the real dimensionality of the system

        :returns:
           * **dist** --  Distances between Wannier function centers (*i* and *j+R*) in Angstroms.

           * **ham** --  Corresponding hopping terms in eV.

        Example usage::

          # get distances and hopping terms
          (dist,ham)=silicon.dist_hop()

          # plot logarithm of the hopping term as a function of distance
          import pylab as plt
          fig, ax = plt.subplots()
          ax.scatter(dist,np.log(np.abs(ham)))
          fig.savefig("localization.pdf")

        """

        if dim==2: w90_model.xyz_cen[:,2]=0  # important change to reduce the dimensionality

        ret_ham=[]
        ret_dist=[]
        site_idx=[]
        site_pairs=[]
        for R in w90_model.ham_r:
            # treat diagonal terms differently
            if R[0]==0 and R[1]==0 and R[2]==0:
                avoid_diagonal=True
            else:
                avoid_diagonal=False

            # get R vector
            vecR=_red_to_cart((w90_model.lat[0],w90_model.lat[1],w90_model.lat[2]),[R])[0]
            for i in range(w90_model.num_wan):
                vec_i=w90_model.xyz_cen[i]
                for j in range(w90_model.num_wan):
                    vec_j=w90_model.xyz_cen[j]
                    # diagonal terms
                    if not (avoid_diagonal==True and i==j):

                        # divide the matrix element from w90 with the degeneracy
                        ret_ham.append(w90_model.ham_r[R]["h"][i,j]/float(w90_model.ham_r[R]["deg"]))

                        # get distance between orbitals
                        ret_dist.append(np.sqrt(np.dot(-vec_i+vec_j+vecR,-vec_i+vec_j+vecR)))

                        #get the site pairs involved in the hopping
                        site_pairs.append((vec_i,vec_j+R))
                        site_idx.append((i,j,R))


        return (np.array(ret_dist),np.array(ret_ham),np.array(site_pairs),site_idx)

def _cart_to_red((a1,a2,a3),cart):
    "Convert cartesian vectors cart to reduced coordinates of a1,a2,a3 vectors"
    # matrix with lattice vectors
    cnv=np.array([a1,a2,a3])
    # transpose a matrix
    cnv=cnv.T
    # invert a matrix
    cnv=np.linalg.inv(cnv)
    # reduced coordinates
    red=np.zeros_like(cart,dtype=float)
    for i in range(0,len(cart)):
        red[i]=np.dot(cnv,cart[i])
    return red

def _red_to_cart((a1,a2,a3),red):
    "Convert reduced to cartesian vectors."
    # cartesian coordinates
    cart=np.zeros_like(red,dtype=float)
    for i in range(0,len(cart)):
        cart[i,:]=a1*red[i][0]+a2*red[i][1]+a3*red[i][2]
    return cart

def calc_projection(my_model,k_vec):
    dim=[item for item in k_vec.shape[:-1]]
    ndim=len(dim)
    nkpt=np.prod(np.array(dim))
    nband,norb,nspin=my_model._nsta,my_model._norb,my_model._nspin
    if dim[0]==1:  energy,wfc_single_k = my_model.solve_one(k_vec[0],eig_vectors=True)
    else:
       energy=[]
       my_array_1=wf_array(my_model,[nkpt])
       for ikpt,kpt in enumerate(k_vec):
           evalue,evec = my_model.solve_one([kk for kk in tuple(kpt)],eig_vectors=True)
           energy.append(evalue)
           my_array_1[ikpt]=evec
    energy=np.array(energy)
    print 'calculating spin components, dimension of k points: ',dim
    print 'number of bands:',nband,', number of orbitals:',norb
    ldos=np.zeros((nkpt,nband,norb*nspin),float)

    fw=open('PROCAR','w')
    print >>fw, 'PROCAR lm decomposed'
    print >>fw, '# of k-points: {0:<5d},    # of bands: {1:<5d},    # of orbitals: {2:<5d}'.format(nkpt,nband,norb)
    for ikpt,kpt in enumerate(k_vec):
        print >>fw, '\nk-point {0}'.format(ikpt+1)+' '.join(['{0:12.7f}'.format(component) for component in kpt])
        for iband in range(nband):
            print >>fw,'\nband {0} # energy = {1:12.7f}'.format(iband+1,energy[ikpt,iband])
            # wavefunction in format of (a1_up,a1_dn,...;an_up,an_dn,...)
            if dim[0]==1: wfc=wfc_single_k[iband].flatten()
            else: wfc = my_array_1[ikpt][iband].flatten()
            for ii in range(norb*nspin): 
                ldos[ikpt,iband,ii] = (sum(wfc[ii].conjugate()*wfc)).real
            print >>fw, '\n'.join(['{0:3d} {1:12.7f}'.format(ii+1,ldos[ikpt,iband,ii]) for ii in range(norb*nspin)])
    return ldos


def plot_proj_band(my_model,path,label):
    import matplotlib.pyplot as plt
    (k_vec,k_dist,k_node)=my_model.k_path(path,101,report=False)
    evals=my_model.solve_all(k_vec)
    ldos=calc_projection(my_model,k_vec)
    color_list=['red','blue','green','violet','orange','yellow']
    fig, ax = plt.subplots()
    for iband in range(len(evals)):
        for ii in range(my_model._norb*my_model._nspin):
            ax.scatter(k_dist,evals[iband],s=100*abs(ldos[:,iband,ii]),edgecolor=color_list[ii],facecolor='None')
    ax.set_xlim([0,k_node[-1]])
    ax.set_xticks(k_node)
    ax.set_xticklabels(label)
    for n in range(len(k_node)): ax.axvline(x=k_node[n],linewidth=0.5, color='k')
    ax.set_ylabel('$Energy$ ($eV$)')
    #[ax.plot(k_dist,evals[i]) for i in range(len(evals))]
    fig.tight_layout()
    fig.savefig("proj_band.png",dpi=600)


# calculate the spin components at each k-point in a spin-orbti tight-binding model
# the projection on each orbital is also calculated
# suitbale for both k_list for band and k_mesh for grid plot
# return: spin components -- sx, sy, sz
#         norm -- the norm of each wavefunction
def calc_spin_band(my_model,k_vec):
    # calculate the spin component S_z of each eigenstate
    # added by Shunhong Zhang, Aug 9, 2017
    dim=[item for item in k_vec.shape[:-1]]
    ndim=len(dim)
    nkpt=np.prod(np.array(dim))
    nband,norb,nspin=my_model._nsta,my_model._norb,my_model._nspin
    if dim[0]==1: 
       energy,wfc_single_k = my_model.solve_one(k_vec[0],eig_vectors=True)
    else:
       if ndim>1: k_vec = k_vec.reshape(nkpt,k_vec.shape[-1])
       energy=[]
       my_array_1=wf_array(my_model,[nkpt])
       for ikpt,kpt in enumerate(k_vec):
           evalue,evec = my_model.solve_one([kk for kk in tuple(kpt)],eig_vectors=True)
           energy.append(evalue)
           my_array_1[ikpt]=evec
    energy=np.array(energy)
    print 'calculating spin components, dimension of k points: ',dim
    print 'number of bands:',nband,', number of orbitals:',norb

    pauli_sigma_x=np.zeros((norb*nspin,norb*nspin),complex)
    pauli_sigma_y=np.zeros((norb*nspin,norb*nspin),complex)
    pauli_sigma_z=np.zeros((norb*nspin,norb*nspin),complex)
    pauli_sigma_i=np.diag(np.ones(norb*nspin,complex))
    #for iorb in range(norb):
    pauli_sigma_x[0::2,1::2]=1
    pauli_sigma_x[1::2,0::2]=1
    pauli_sigma_y[0::2,1::2]=-1.j
    pauli_sigma_y[1::2,0::2]=1.j
    pauli_sigma_z[0::2,0::2]=1
    pauli_sigma_z[1::2,1::2]=-1
    pauli_sigma_i[0::2,0::2]=1
    pauli_sigma_i[1::2,1::2]=1

    sx=np.zeros((nkpt,nband),float)
    sy=np.zeros((nkpt,nband),float)
    sz=np.zeros((nkpt,nband),float)
    norm=np.zeros((nkpt,nband),float)
    ss=np.zeros((nkpt,nband),float)

    fw=open('spin_texture','w')
    for ikpt,kpt in enumerate(k_vec):
        for iband in range(nband):
            if dim[0]==1: wfc=np.matrix(wfc_single_k[iband].flatten())
            else: wfc = np.matrix(my_array_1[ikpt][iband].flatten())
            norm[ikpt,iband] = (wfc.conjugate()*pauli_sigma_i*wfc.T).real
            sx[ikpt,iband]=(wfc.conjugate()*pauli_sigma_x/2*wfc.T).real/norm[ikpt,iband]
            sy[ikpt,iband]=(wfc.conjugate()*pauli_sigma_y/2*wfc.T).real/norm[ikpt,iband]
            sz[ikpt,iband]=(wfc.conjugate()*pauli_sigma_z/2*wfc.T).real/norm[ikpt,iband]
            #ss[ikpt,iband]=np.linalg.norm([sx[ikpt,iband],sy[ikpt,iband],sz[ikpt,iband]])
            print >>fw,' '.join([' '.join(['{0:9.6f}'.format(item) for item in (sx[ikpt,iband],sy[ikpt,iband],sz[ikpt,iband])]) for iband in range(nband)])
    if ndim>1: 
       energy=energy.reshape(tuple(dim+[nband]))
       sx = sx.reshape(tuple(dim+[nband]))
       sy = sy.reshape(tuple(dim+[nband]))
       sz = sz.reshape(tuple(dim+[nband]))
       norm = norm.reshape(tuple(dim+[nband]))
    fw.close()
    return energy,sx,sy,sz,norm

def plot_spin_texture(my_model,center=[0,0,0],bound=[0.6,0.6,0],nkx=15,nky=15,nkz=1,ktype='frac',startband=None,lastband=None):
    kx=np.linspace(center[0]-bound[0],center[0]+bound[0],num=nkx)
    ky=np.linspace(center[1]-bound[1],center[1]+bound[1],num=nky)
    kz=np.linspace(center[2]-bound[2],center[2]+bound[2],num=nkz)
    center=[center[i] for i in my_model._per]
    bound=[bound[i] for i in my_model._per]
    lat_per=np.copy(my_model._lat)
    #lat_per=lat_per[my_model._per]

    if (lat_per.shape[0]==lat_per.shape[1]):
       # lat_per is invertible
       lat_per_inv=2*np.pi*np.linalg.inv(lat_per).T
       print 'reciprocal-space lattice vectors\n', lat_per_inv

    if ktype=='cart': print 'cartersian k-point coordinates used'
    elif ktype=='frac': print 'fractional k-point coordinates used'
    else: exit('ktype need to be specified as frac (Default) or cart!')
    try: import matplotlib.pyplot as plt
    except:
         print ('skipping plotting because failure in loading python-matplotlib')
         return 1
    # calculate the spin texture, added by Shunhong Zhang, Sep 2, 2017
    print ('Calculating spin texture...')
    if ktype=='frac': 
       k_vec=np.array(np.meshgrid(kx,ky,kz))
       tmp=np.rollaxis(k_vec,3)
       tmp=np.rollaxis(tmp,3)
       k_vec=np.rollaxis(tmp,3)
       k_vec=k_vec[:,:,:,my_model._per]
       k_vec_cart=np.array([[[np.dot(k_vec[i,j,k],lat_per_inv) for k in range(nkz)] for j in range(nky)] for i in range(nkx)])
    if ktype=='cart': 
       k_vec_cart = np.array(np.meshgrid(kx,ky,kz))
       tmp=np.rollaxis(k_vec_cart,3)
       tmp=np.rollaxis(tmp,3)
       k_vec_cart=np.rollaxis(tmp,3)
       k_vec_cart=k_vec_cart[:,:,:,my_model._per]
       k_vec = np.array([[[np.dot(k_vec_cart[i,j,k],np.linalg.inv(lat_per_inv)) for k in range(nkz)] for j in range(nky)] for i in range(nkx)])
    lim=max(np.max(k_vec_cart),abs(np.min(k_vec_cart)),np.max(lat_per_inv),abs(np.min(lat_per_inv)))*1.2
    energy_c,sx_c,sy_c,sz_c,norm_c=calc_spin_band(my_model,np.array([center]))
    print 'plot range in reciprocal lattice coordinates:'
    print center[0]-bound[0],center[0]+bound[0],center[1]-bound[1],center[1]+bound[1]
    print 'plot center: ('+' '.join(['{0:8.5f}'.format(item) for item in center]),')'
    print 'plot center(cart):',np.dot(lat_per_inv,center)
    print 'the spin component is:\n'+' '.join(['{0:11s}'.format(item) for item in ['band','energy','sx','sy','sz']])
    for iband in range(energy_c.shape[-1]):
        print '{0:3d}'.format(iband).center(5)+' '.join(['{0:12.6f}'.format(item).center(12) for item in [energy_c[iband],sx_c[0,iband],sy_c[0,iband],sz_c[0,iband]]])
    energy,sx,sy,sz,norm=calc_spin_band(my_model,k_vec)
    fw=open('kpts','w')
    for i in range(nkx):
        for j in range(nky):
            print >>fw,''.join(['{0:9.5f}'.format(item) for item in np.concatenate((k_vec[i,j,0],k_vec_cart[i,j,0]),axis=0)])
    fw.close()
    if not startband: startband=0
    if not lastband:  lastband=energy.shape[-1]
    for iband in range(startband,lastband):
        print 'calculate and plot spin texture of band # {0}'.format(iband)
        fig = plt.figure(figsize=(8,6))
        plt.title('spin texture of band {0}'.format(iband))
        colormap='cool'
        Q = plt.quiver(k_vec_cart[:,:,0,0], k_vec_cart[:,:,0,1], sx[:,:,0,iband], sy[:,:,0,iband], sz[:,:,0,iband], 
            cmap=colormap, units='inches', pivot='mid', width=0.01, scale=3,headlength=10,headwidth=5)
        cb = plt.colorbar(Q)
        plt.xlim(-lim,lim)
        plt.ylim(-lim,lim)
        plt.quiver(center[0],center[1],lat_per_inv[0,0],lat_per_inv[0,1],edgecolor='None',facecolor='g',units='x',scale=1,pivot='tail',linestyle='dashed',width=0.005*lim,headwidth=4)
        plt.quiver(center[0],center[1],lat_per_inv[1,0],lat_per_inv[1,1],edgecolor='None',facecolor='g',units='x',scale=1,pivot='tail',linestyle='dashed',width=0.005*lim,headwidth=4)
        fig.tight_layout()
        fig.savefig('band_'+str(iband)+"_spin_texture.png",bbox_inches='tight',dpi=600)
    print 'Done\n'



def _nice_float(x,just,rnd):
    return str(round(x,rnd)).rjust(just)
def _nice_int(x,just):
    return str(x).rjust(just)
def _nice_complex(x,just,rnd):
    ret=""
    ret+=_nice_float(complex(x).real,just,rnd)
    if complex(x).imag<0.0:
        ret+=" - "
    else:
        ret+=" + "
    ret+=_nice_float(abs(complex(x).imag),just,rnd)
    ret+=" i"
    return ret



if __name__=='__main__':
   print 'this is an extension python code for PythTB'
   print 'a code for tight-binding modeling'
   print 'please refer to '
   print 'http://www.physics.rutgers.edu/pythtb/'
   print 'for more information about PythTB'
