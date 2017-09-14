#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       parse.py                                                     #
#  Dependence: crystal_structure.py                                         #
#  Usage:      convert the EIGENVAL file to a file for band structure plot  #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Aug 02, 2017                                                 #
#  Note:       coded partially based on the vasp_unfold code                #
#              by Milan Tomic <tomic@th.physik.uni-frankfurt.de>            #
#                                                                           #
#===========================================================================#


import numpy as np
from crystal_structure import *
from math import *
from utils import *
import os


#orb_dic={0:'s',1:'py',2:'pz',3:'px',
#         4:'dxy',5:'dyz',6:'dz2',7:'dxz',8:'dx2-y2'}


def parse_poscar(filename='POSCAR'):
    '''Parses POSCAR file. Returns 3x3 float array
    containing unit cell vectors as rows, Nx3
    array containing fractional positions of atoms,
    N being number of atoms and a list of chemical
    symbols, and a list of number of atoms for each chemical species.
    '''
    gl=open(filename) 
    gl.readline()
    scale = float(gl.readline())
    cell=scale*np.fromfile(gl,sep=' ',count=9).reshape(3,3)
    line=gl.readline().split()
    try:
        counts = np.array(line, int)                    #VASP 4
        syms = raw_input("input chemical species:")
        addline="sed -i '5a\\"+syms+"' "+filename
        print addline
        os.system(addline)
        syms=syms.split()
    except:
        syms =  line                                    #VASP 5
        counts = np.array(gl.readline().split(), int)

    # Rearrange into the long list of chemical symbols
    symbols = []
    for s, c in zip(syms, counts):
        symbols += [s]*c
    
    # Cartesian or fractional coordinates?
    ctype = gl.readline()[0].lower()
    
    if ctype == 'c':
        mult = np.linalg.inv(cell)
    elif ctype == 'd':
        mult = np.eye(3)
    else:
        exit('"{0}" is unknown POSCAR option'.format(plines[7].strip()))
    
    # Allocate storage for positions and read positions from file
    spos = np.zeros((len(symbols), 3))
    for i in xrange(len(symbols)):
        spos[i] = np.array(gl.readline().split()[:3], float)
    # If necessary, this will convert from Cartesian to fractional coordinates
    spos = np.dot(spos, mult)

    sym_seen=set()
    species=[]
    for sym in symbols:
        if sym not in sym_seen:
           sym_seen.add(sym)
           species.append(sym)

    struct=cryst_struct(cell,species,symbols,counts,spos)
    return struct
    
def parse_doscar(filename='DOSCAR'):
    "Parsing {0} ...".format(filename)
    try:
       nspin = int(os.popen("grep ISPIN OUTCAR 2>/dev/null").readline().split()[2])
    except:
       nspin=input("nspin=")
    if nspin==2: print "Spin polarized DOS read"
    elif nspin==1: print "Spin unpolarized DOS read"

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
           data = np.fromfile(f,sep=' ',count=nedos*(norb*nspin+1)).reshape((nedos,norb*nspin+1))
           for ispin in range(nspin):
               dos_orb[iat,:,:,ispin] = data[:,ispin+1::nspin]
    except:
       print "PDOS not read!"
    print 'done'
    f.close()
    
    energy=energy-efermi
    return energy,emin,emax,dos_total,dos_integral,dos_orb

def parse_eigenval(filename='EIGENVAL',efermi=0):
    struct = parse_poscar('POSCAR')
    
    #read the specified path in the reciprocal space from the KPOINTS file
    if not efermi:
       efermi = float(os.popen("grep E-fermi OUTCAR 2>/dev/null").readline().split()[2])
       print 'we read fermi energy from the OUTCAR file of the band structure calculation'
       print 'plase do a dense-k-mesh nscf calculation if you want to get a more accurate fermi energy'
    nspin = int(os.popen("grep SPIN OUTCAR 2>/dev/null").readline().split()[2])
    print "Fermi energy = ",efermi,"eV\nISPIN = ",nspin
    lsorbit=os.popen('grep LSORBIT OUTCAR 2>/dev/null').readline().split()[2]

    f = open(filename,'r')
    head = np.fromfile(f,sep=' ',count=5)
    [f.readline() for i in range(3)]
    system = f.readline()[4]
    nelect, nkpt, nband =np.fromfile(f, sep=' ',count=3, dtype=int)
    val_idx = nelect/2-1
    con_idx = val_idx+1

    kpt = np.zeros((nkpt,3),float)
    kweight = np.zeros(nkpt,float)
    if lsorbit=='F':
       energy = np.zeros((nkpt,nband,nspin),float)
       for ikpt in range(nkpt):
           kpt[ikpt] = np.fromfile(f, sep=' ', count=3, dtype=float)
           kweight[ikpt] = np.fromfile(f, sep=' ', count=1, dtype=float)
           for iband in range(nband):
               line = np.array([item for item in f.readline().split()],float)
               energy[ikpt,iband,:] = line[1:nspin+1]
    elif lsorbit=='T':
       energy = np.zeros((nspin,nkpt,nband),float)
       print 'band structure with spin-orbit coupling'
       for ikpt in range(nkpt):
           line=f.readline().split()
           kpt[ikpt] = np.array([float(item) for item in line[:-1]])
           kweight[ikpt] = float(line[-1])
           for iband in range(nband):
               line=f.readline().split()
               for ispin in range(nspin): energy[ispin,ikpt,iband]=float(line[1])
           f.readline()
    f.close()
    kpt = np.matrix(kpt,float)
    #kpt_cart = kpt*struct._reciprocal_cell()
    print "Number of Electrons = ",nelect,"\nNumber of K points =",nkpt,"\nNumber of bands = ",nband
    return lsorbit,nspin,nelect,efermi,kpt,kweight,energy


def parse_procar(filename,efermi=0):
    '''This function parses a PROCAR file. It returns a tuple
    consisting of following elements:
    
    orbitals    - (norbs) string array of orbital labels (s, px, py etc...)
    kpoints     - (nkpt,3) float array of k-point coordinates
    kweights    - (nkpt) float array of k-point weights
    energy      - (nkpt,nband,nspin) float array of band energies
    occupancies - (nkpt,nband,nspin) float array of band occupancies
    weights     - (nkpt,natom*norbs,nband,ndim,nspin) float array
                  of orbital weights
    phases      - (nkpt,natom*norbs,nband,nspin) complex array of
                  phases of orbital weights if LORBIT=12, otherwise None
                  
    Where:
    
    norbs   - number of orbitals (It can be 9 or 16 with f orbitals)
    nkpt    - number of k-points
    nband   - number of energy
    nspin   - number of spins (1 for non spin-polarized, 2 otherwise)
    natom   - number of atoms
    ndim    - orbital weight dimensionality (1 for collinear, 4 otherwise)
    '''

    print "Parsing {0} ...".format(filename)
    #fix possible format errors in the procar file
    fix='sed -i "s/\([0-9]\)\-\([0-9]\)/\\1 -\\2/g" '+filename
    os.system(fix)


    try:
        gl = Getlines(filename)
    except:
        post_error('Unable to open "{0}" for reading.'.format(filename))


    header_1 = gl.readline()
    header_2 = gl.readline().split()
    
    nkpt = int(header_2[3])
    nband = int(header_2[7])
    natom = int(header_2[-1])
    
    print "nkpt=",nkpt,"nband=",nband,"natom=",natom
    
    # Remember the position in file
    start = gl.tell()
    
    # Skip two lines containing first k-point and band
    gl.readline()
    gl.readline()
    
    # Determine the number of orbitals
    orbitals = gl.readline().split()[1:-1]
    
    norbs = len(orbitals)
    print 'number of orbitals:',norbs
    
    # Allocate maximal storage. In case calculation was
    # non spin-polarized or collinear we can just trim
    # the excess components at the end
    kpoints = np.zeros((nkpt, 3), float)
    kweights = np.zeros(nkpt, float)
    energy = np.zeros((nkpt, nband, 2), float)
    occupancies = np.zeros((nkpt, nband, 2), float)
    weights = np.zeros((nkpt, natom, (norbs+1), nband, 4, 2), float)
    
    dim = 0
    
    # Determine if the calculation was non-collinear
    # by counting how many lines in the first band
    # block begin with tot. That number will be equal
    # to the number of sub-blocks for orbital weights
    # (1 in case of collinear and 4 otherwise)
    while True:
        line = gl.readline()
        
        if line.startswith('tot'):
            dim += 1
        elif line.startswith('band'):
            break
    
    # This function will read block of absolute weights
    # for i-th k-point, j-th band and s-th spin component
    def get_absweights(ikpt, iband, ispin):
        # Skip line with orbital names        
        gl.readline()
        
        for idim in xrange(dim):
            # Fetch entiret orbital weight block
            data = np.fromfile(gl, sep=" ", count=natom*(norbs+2))
            # Cast it into tabular shape
            data = data.reshape((natom, norbs+2))
            
            # Discard first column and store weights
            #weights[ikpt,:,:,iband,idim,s] = data[:,1:].flatten()
            weights[ikpt,:,:,iband,idim,ispin] = data[:,1:]

            
            # Skip line with the totals
            gl.readline()
            
    # Check whether phase information is included
    if '+ phase' in header_1:
        # Allocate storage for phases
        phases = np.zeros((nkpt, natom, norbs, nband, 2), complex)

        # Declare nested function that handles 
        # parsing of complex weights
        def get_weights(ikpt, iband, ispin):
            # Read abs values of weights
            get_absweights(ikpt, iband, ispin)
            
            # Skip line with orbital names
            gl.readline()
            
            # Fetch entire phase block
            data = np.fromfile(gl, sep=" ", count=2*natom*(norbs+1))
            # Cast it into tabular shape
            data = data.reshape((2*natom, norbs+1))

            # Discard first column and store real and imaginary
            # parts respectively
            phases[ikpt,:,:,iband,ispin] =  data[::2,1:]     #The real part
            phases[ikpt,:,:,iband,ispin] += 1j*data[1::2,1:]   #The complex
    else:
        # Phases are None in this case
        phases = None
        
        # In this case we just have absolutes of weights
        # No need for a new function
        get_weights = get_absweights
    
    # Go back to the beginning of the first k-point
    gl.seek(start)
    
    for ikpt in xrange(nkpt):
        # Parse k-point coordinates
        k_line = gl.readline().split()
        
        kpoints[ikpt] = [float(k_line[c]) for c in [3, 4, 5]]
        kweights[ikpt] = float(k_line[-1])

        
        for iband in xrange(nband):
            # Parse band energy
            band_line = gl.readline().split()

            energy[ikpt, iband, 0] = float(band_line[4])
            occupancies[ikpt, iband, 0] = float(band_line[-1])
            
            # Parse orbital weights
            get_weights(ikpt, iband, 0)
    
    # Seek now for the second spin component
    res = gl.readline(False)
    
    # If there is no second spin component, finish by
    # returning just the first component
    if res is None:
        # Trim the second spin component
        energy = energy[:,:,:1]
        occupancies = occupancies[:,:,:1]
        weights = weights[:,:,:,:,:dim,:1]
        
        if phases is not None:
            phases = phases[:,:,:,:1]
        
        if efermi:
           energy=energy-efermi
        print "Done"
        return [orbitals, kpoints, kweights, energy, occupancies, \
            weights, phases]
    
    # Otherwise, read energy and weights for the second component
    for ikpt in xrange(nkpt):
        # Skip k-point coordinates
        gl.readline()
        
        for iband in xrange(nband):
            # Parse band energy
            band_line = gl.readline().split()

            energy[ikpt, iband, 1] = float(band_line[4])
            occupancies[ikpt, iband, 1] = float(band_line[-1])
            
            # Parse orbital weights
            get_weights(ikpt, iband, 1)
    if efermi:
       energy=energy-efermi 
    print "Done"
    return [orbitals, kpoints, kweights, energy, occupancies, \
        weights[:,:,:,:,:dim,:], phases]


def parse_chg(filename='CHG'):
    print "Parsing {0} ...".format(filename)
    gl=open(filename,'r')
    gl.readline()
    scale = float(gl.readline())
    cell=scale*np.fromfile(gl,sep=' ',count=9).reshape(3,3)
    line=gl.readline().split()
    try:
        counts = np.array(line, int)                    #VASP 4
        syms = raw_input("input chemical species:")
        addline="sed -i '5a\\"+syms+"' "+filename
        print addline
        os.system(addline)
        syms=syms.split()
    except:
        syms =  line                                    #VASP 5
        counts = np.array(gl.readline().split(), int)
    symbols = []
    for s, c in zip(syms, counts):
        symbols += [s]*c
    ctype = gl.readline()[0].lower()
    if ctype == 'c':
        mult = np.linalg.inv(cell)
    elif ctype == 'd':
        mult = np.eye(3)
    else:
        exit('"{0}" is unknown POSCAR option'.format(plines[7].strip()))
    spos = np.zeros((len(symbols), 3))
    for i in xrange(len(symbols)):
        spos[i] = np.array(gl.readline().split()[:3], float)
    spos = np.dot(spos, mult)
    sym_seen=set()
    species=[]
    for sym in symbols:
        if sym not in sym_seen:
           sym_seen.add(sym)
           species.append(sym)
    struct=cryst_struct(cell,species,symbols,counts,spos)

    gl.readline() 
    get_grid=gl.readline().split()
    NGXF=int(get_grid[0])
    NGYF=int(get_grid[1])
    NGZF=int(get_grid[2])
    ngrid=NGXF*NGYF*NGZF
    print "real space grid:", NGXF, NGYF, NGZF
    print "total grid points:",ngrid
    data = np.fromfile(gl, sep=" ", count = ngrid)
    chg = [float(item) for item in data]
    print "Done"
    return struct, NGXF, NGYF, NGZF, chg


#----------------------#
# parse abinit files   #
#----------------------#

def parse_abinit_struct(filename=None):
    if filename:
       fil=filename
    else:
       fil=os.popen('ls *.out').read()
    get_latt=os.popen('grep -A 3 "space primitive vectors" '+fil).readlines()[-3:]
    cell=np.array([[float(item) for item in line.split()[1:4]] for line in get_latt])
    nspec=int(os.popen('grep ntypat '+fil).readlines()[-1].split()[1])
    count=[int(item) for item in os.popen('grep typat '+fil).readlines()[-1].split()[1:]]
    znucl=[int(float(item)) for item in os.popen('grep znucl '+fil).readlines()[-1].split()[1:]]
    symbols=[open('/home/hydrogen/app/tools/elements').readlines()[item-1].split()[3] for item in znucl]
    print symbols
