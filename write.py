#!/usr/bin/python

#===========================================================================#
#                                                                           #
#  File:       write.py                                                     #
#  Dependence: mass.py,element (required)                                   #
#  Usage:      write files for VASP or Quantum ESPRESSO                     #      
#  Author:     Shunhong Zhang <zhangshunhong@tsinghua.edu.cn>               #
#  Date:       Jul 14, 2017                                                 #
#                                                                           #
#===========================================================================#


import os
import numpy as np
from math import *
from termcolor import colored


#==========================================#
# input and post-processing file for VASP  #
#==========================================#

def write_poscar_head(scale,cell,system='system',filename = 'POSCAR'):
    fpos=open(filename,'w')
    print >>fpos,system
    print >>fpos,scale
    cell=np.array(cell)
    print >>fpos,'\n'.join('{0:15.12f} {1:15.12f} {2:15.12f}'.format(ce[0],ce[1],ce[2]) for ce in cell)
    fpos.close()

def write_poscar_atoms(struct,filename='POSCAR'):
    fpos=open(filename,'a')
    for sym in struct._species:
        print >>fpos,sym,
    print>>fpos,''
    for i in range(len(struct._counts)):
        print >>fpos,struct._counts[i],
    print >>fpos,'\nDirect'
    print >>fpos,'\n'.join('{0:15.12f} {1:15.12f} {2:15.12f}'.format(pos[0],pos[1],pos[2]) for pos  in struct._pos)
    fpos.close()

def write_dos_total(energy,dos_total,dos_integral):
    nedos=len(energy)
    nspin=dos_total.shape[1]
    fw_total=open("Total_DOS.dat","w")
    for ie in range(nedos):
        print >>fw_total,"{0:9.5f} ".format(energy[ie]),
        for ispin in range(nspin):
            print >>fw_total,'{0:12.6f} {1:12.6f}'.format(dos_total[ie,ispin],dos_integral[ie,ispin]),
        print>>fw_total,''
    fw_total.close()

def write_pdos_atom(energy,dos_atom,symbols):
    natom,nedos,nspin=dos_atom.shape
    for iatom in range(natom):
        fw_atom=open(symbols[iatom]+str(iatom)+'_pdos.dat','w')
        for ie in range(nedos):
            if nspin==1:
               print>>fw_atom,"{0:9.5f} {1:10.6f}".format(energy[ie],dos_atom[iatom,ie,0])
            elif nspin==2:
               print>>fw_atom,"{0:9.5f} {1:10.6f} {2:10.6f}".format(energy[ie],dos_atom[iatom,ie,0],dos_atom[iatom,ie,1])
        fw_atom.close()
    

def write_pdos_species(energy,dos_species,dos_species_orb,symbols,counts):
    nspecies=dos_species.shape[0]
    nedos=len(energy)
    norb=dos_species_orb.shape[2]
    nspin=dos_species.shape[2]
    start_count=0
    for i in range(nspecies):
        fw_spec=open(symbols[start_count]+'_pdos.dat','w')
        for ie in range(nedos):
            print>>fw_spec,"{0:9.5f}".format(energy[ie]),
            for ispin in range(nspin):
                for iorb in range(norb):
                    print>>fw_spec," {0:10.6f} ".format(dos_species_orb[i,ie,iorb,0]),
                    print>>fw_spec,"{0:10.6f}".format(dos_species[i,ie,0]),
        start_count=np.cumsum(counts)[i]
        fw_spec.close()

#write the file for drawing the band structure
def write_band(path,segment_nkpt,energy,nelect):
    nkpt,nband,nspin =energy.shape
    fw = open('Band.dat','w')
    for iband in range(nband):
        for ikpt in range(nkpt):
            if ikpt%segment_nkpt != 0 or ikpt == 0:
               print>>fw,' {0:8.5f}  {1:8.5f}'.format(path[ikpt],energy[ikpt,iband,0]),
               if nspin == 2:
                  print>>fw,' {0:8.5f}'.format(energy[ikpt,iband,1]),
               print >>fw,''
        print>>fw,''
    fw.close()

def write_specified_band(kpt,kweight,path,energy,nelect):
    nkpt,nband,nspin =energy.shape
    val_idx=nelect/2-1    #Note: band index starts from 0
    start_band=val_idx-1
    last_band=val_idx+1
    start_kpt=0
    for ikpt in range(nkpt):
        if kweight[ikpt]!=0:
           ikpt+=1
        else:
           start_kpt=ikpt
           break

    fw = open('specified_band.dat','w')
    for ikpt in range(start_kpt,nkpt):
        #print >>fw, '{0:>4d}  {1:8.5f}  {2:8.5f}  {3:8.5f}'.format(ikpt,kpt[ikpt,0],kpt[ikpt,1],kpt[ikpt,2]),
        print >>fw, '{0:8.5f}'.format(path[ikpt]),
        for ispin in range(nspin):
            for iband in range(start_band,last_band+1):
                print >>fw, '{0:10.5f} '.format(energy[ikpt,iband,ispin]),
        print>>fw,""
    fw.close()


def write_mesh_band(kpt,energy,nspin,nelect,grid,isoenergy=0):
    nkpt =energy.shape[0]
    nband=energy.shape[1]
    nspin=energy.shape[2]
    val_idx=nelect/2-1    #Note: band index starts from 0
    start_band=val_idx-1
    last_band=val_idx+1

    deltaE=0.002
    fs = open('Fermi_surface','w')
    for ispin in range(nspin):
        print >>fs,"spin channel ",ispin,"(0 and 1 represent two spin channels)"
        for ikpt in range(nkpt):
            for iband in range(nband):
                if   energy[ikpt,iband,ispin]>=0 and energy[ikpt,iband,ispin]<=deltaE:
                     print >>fs,"{0:10.5f} {1:10.5f} {2:2d}".format(kpt[ikpt,0],kpt[ikpt,1],1)
                elif energy[ikpt,iband,ispin]<=0 and energy[ikpt,iband,ispin]>=-deltaE:
                     print >>fs,"{0:10.5f} {1:10.5f} {2:2d}".format(kpt[ikpt,0],kpt[ikpt,1],-1)
    fs.close()
    grid_x=grid[0]
    if np.mod(nkpt,grid_x)!=0:
       print "the k-point grid is irregular, please switch off symmetry in the calculation (ISYM=0)"
       exit() 

    fm=open("matrix_band","w")
    for idim in range(2):
       for ikpt in range(nkpt):
            print >>fm,"{0:10.6f}".format(kpt[ikpt,idim]),
            if np.mod(ikpt+1,grid_x)==0:
               print >>fm,''
    print >>fm,"energy"
    for ispin in range(nspin):
       for iband in range(start_band,last_band+1):
           for ikpt in range(nkpt):
               print >>fm,"{0:10.6f}".format(energy[ikpt,iband,ispin]),
               if np.mod(ikpt+1,grid_x)==0:
                  print >>fm,'' 
    fm.close()

def write_weighted_band(path,energy,weights,struct,orbitals,xsym,args):
    nkpt,natom,norb,nband,ndim,nspin=weights.shape
    if args.proj_type=='orbitals':
       for ispin in range(nspin):
           for projection in args.proj_index.split(","):
               composition=projection.split('_')
               iat=int(projection.split('_')[0])
               iorb=int(projection.split('_')[1])
               fw=open('proj_band_'+struct._symbols[iat]+str(iat)+'_'+str(iorb)+'.dat','w')
               for iband in range(nband):
                   for ikpt in range(nkpt):
                       print >>fw,'{0:8.5f} {1:10.6f} {2:10.6f}'.format(path[ikpt],energy[ikpt,iband,ispin],weights[ikpt,iat,iorb,iband,0,ispin])
                   print>>fw,''
               fw.close()

#============================================================#
# input file for quantum ESPRESSO (only for SCF calculation) #
#============================================================#

def write_pwi(setup_dic,ibrav,celldm,struct,filename="scf.in"):
    if filename:
       filename=open(filename,"w")
    mass_path=os.popen('which mass.py').readline().rstrip('\n')+' '
    atomic_mass=[(float(os.popen(mass_path+sym).readline().split()[1])) for sym in struct._species]
    print>>filename,"&CONTROL"
    print>>filename,"calculation = ","'scf'"
    print>>filename,"restart_mode = ","'from_scratch'"
    print>>filename,"outdir = './tmp/'"
    print>>filename,"pseudo_dir = ",setup_dic['pseudo_dir']
    print>>filename,"prefix = ",setup_dic['prefix']
    print>>filename,"/"
    print>>filename,"&SYSTEM"
    print>>filename,"ibrav=",ibrav
    #celldm=find_celldm(ibrav,latt)
    for i in range(1,7):
       if celldm[i]!=0:
          print>>filename,"celldm("+str(i)+") =",celldm[i]
    print>>filename,"nat =",struct._natom
    print>>filename,"ntyp =",len(struct._species)
    print>>filename,"ecutwfc = ",setup_dic['ecutwfc']
    print>>filename,"ecutrho = ",setup_dic['ecutrho']
    print>>filename,"occupations = smearing"
    print>>filename,"smearing ='gaussian'"
    print>>filename,"degauss = 0.001"
    print>>filename,"/"
    print>>filename,"&ELECTRONS"
    print>>filename,"conv_thr=1.d-8"
    print>>filename,"/"
    print>>filename,"ATOMIC_SPECIES"
    for item,mass in zip(struct._species,atomic_mass):
        print>>filename, "{0:2s} {1:10.5f} {2:>2s}".format(item,mass,item+setup_dic['upf'])
    print>>filename,"ATOMIC_POSITIONS crystal"
    for sym,atom in zip(struct._symbols,struct._pos):
        print>>filename, "{0:2s} {1:15.14f} {2:15.14f} {3:15.14f}".format(sym,atom[0],atom[1],atom[2])
    print>>filename,"K_POINTS automatic"
    print>>filename, setup_dic['kmesh'],setup_dic['kshift']
    try:
        filename.close()
    except:
        pass 

#==================================================#
# input file for Materials Studio (visualization ) #
#==================================================#

def write_grid(latt,NGXF,NGYF,NGZF,data,filename='density.grd'):
    fw=open(filename,"w")
    print>>fw,"DMol3 total electron density"
    print>>fw,'(1p,e12.5)'
    print>>fw,latt['a'],latt['b'],latt['c'],latt['alpha'],latt['beta'],latt['gamma']
    print>>fw,NGXF-1,NGYF-1,NGZF-1
    print>>fw,1,0,NGXF,0,NGYF,0,NGZF
    for i in range(len(data)):
        print>>fw,'{0:10.5f}'.format(data[i])

#==================================================#
# input file for OPENMX (only for SCF calculation) #
#==================================================#

def write_openmx_dat(prefix,struct,datapath):
    import pao_set
    print "the generated input file should be checked carefully before real calculation!"
    print "make sure that it is consistent with what you want!"
    fw=open(prefix+'.dat','w')
    print>>fw,'#\n# .dat file for OPENMX, generated by v2openmx\n#'
    print>>fw,'System.CurrrentDirectory         ./    # default=./'
    print>>fw,'System.Name                   '+ prefix
    print>>fw,'level.of.stdout                  1     # default=1 (1-3)'
    print>>fw,'level.of.fileout                 1     # default=1 (0-2)\n'
    print>>fw,'\n#\n# Definition of Atomic Species\n#\n'
    print>>fw,'Species.Number',len(struct._species) 
    print>>fw,'<Definition.of.Atomic.Species'
    for ispec in struct._species:
        get_vps='locate '+datapath+'/VPS/'+ispec+'_'
        if ispec=='E':
           vps='E.vps'
        else:
           vps=os.popen(get_vps+'|grep PBE13S').readline().split('/')[-1].rstrip('\n')
           if not vps:
              vps=os.popen(get_vps+'|grep PBE13').readline().split('/')[-1].rstrip('\n')
        print ispec, vps
        print >>fw,'{0:8s} {1:30s} {2:20s}'.format(ispec,pao_set.pao_set[ispec],vps.rstrip('.vps'))
    print>>fw,'Definition.of.Atomic.Species>\n'
    print>>fw,'Atoms.SpeciesAndCoordinates.Unit  FRAC    # Ang|AU|FRAC '
    print>>fw,'Atoms.Number ','{0:4d}'.format(struct._natom)
    print>>fw,'<Atoms.SpeciesAndCoordinates'
    for iat in range(struct._natom):
        if struct._symbols[iat]=='E':
           charge_up,charge_dn=0,0
        else:
           get_ve='grep valence.electron '+datapath+'/VPS/'+struct._symbols[iat]
           try:
              ve=float(os.popen(get_ve+'_PBE13.vps 2>/dev/null').readline().split()[1])
           except:
              ve=float(os.popen(get_ve+'_PBE13S.vps').readline().split()[1])
           charge_up=ve/2.0
           charge_dn=ve/2.0
        print >>fw, '{0:4d} {1:4<s}'.format(iat+1, struct._symbols[iat]),("%20.16f %20.16f %20.16f" % (tuple(struct._pos[iat]))),
        print >>fw, '{0:4.2f} {1:4.2f}'.format(charge_up, charge_dn)
    print>>fw,'Atoms.SpeciesAndCoordinates>\n'
    print>>fw,'Atoms.UnitVectors.Unit            Ang # Ang|AU'
    print>>fw,'<Atoms.UnitVectors'
    for i in range(3):
        print >>fw, ("%20.16f %20.16f %20.16f" % (tuple(struct._cell[i])))
    print>>fw,'Atoms.UnitVectors>\n'
    print>>fw,'\n#\n# SCF or Electronic System\n#\n'
    print>>fw,'scf.XcType                 GGA-PBE     # LDA|LSDA-CA|LSDA-PW|GGA-PBE'
    print>>fw,'scf.SpinPolarization        off        # On|Off|NC'
    print>>fw,'scf.SpinOrbit.Coupling      off        # On|Off, default=off'       
    print>>fw,'scf.ElectronicTemperature  300.0       # default=300 (K)'
    print>>fw,'scf.energycutoff           200.0       # default=150 (Ry)'
    print>>fw,'scf.maxIter                  40        # default=40'
    print>>fw,'scf.EigenvalueSolver       band        # DC|GDC|Cluster|Band'
    print>>fw,'scf.Kgrid                  1 1 1       # means n1 x n2 x n3'
    print>>fw,'scf.Mixing.Type           rmm-diisk    # Simple|Rmm-Diis|Gr-Pulay|Kerker|Rmm-Diisk'
    print>>fw,'scf.Init.Mixing.Weight     0.20        # default=0.30 '
    print>>fw,'scf.Min.Mixing.Weight      0.001       # default=0.001 '
    print>>fw,'scf.Max.Mixing.Weight      0.500       # default=0.40 '
    print>>fw,'scf.Mixing.History          7          # default=5'
    print>>fw,'scf.Mixing.StartPulay       7          # default=6'
    print>>fw,'scf.Mixing.EveryPulay       1          # default=6'
    print>>fw,'scf.criterion             1.0e-12      # default=1.0e-6 (Hartree) '
    print>>fw,'scf.lapack.dste            dstevx      # dstevx|dstedc|dstegr,default=dstevx'
    print>>fw,'scf.restart               off'
    print>>fw,'orderN.HoppingRanges        6.5        # default=5.0 (Ang) '
    print>>fw,'orderN.NumHoppings           2         # default=2'
    print>>fw,'orderN.KrylovH.order        400'
    print>>fw,'orderN.Expand.Core           on'
    print>>fw,'orderN.Recalc.Buffer         on'
    print>>fw,'orderN.Exact.Inverse.S       on'
    print>>fw,'\nDATA.PATH ',datapath
    fw.close()

