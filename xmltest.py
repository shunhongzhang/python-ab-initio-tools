#!/usr/bin/python

import sys
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import io
import numpy as np

# vasprun.xml handling
#
class VasprunWrapper(object):
    """VasprunWrapper class
    This is used to avoid VASP 5.2.8 vasprun.xml defect at PRECFOCK,
    xml parser stops with error.
    """
    def __init__(self, fileptr):
        self._fileptr = fileptr

    def read(self, size=None):
        element = self._fileptr.next()
        if element.find("PRECFOCK") == -1:
            return element
        else:
            return "<i type=\"string\" name=\"PRECFOCK\"></i>"


class Vasprun(object):
    def __init__(self, fileptr, use_expat=False):
        self._fileptr = fileptr
        self._use_expat = use_expat

    def read_forces(self):
        if self._use_expat:
            return self._parse_expat_vasprun_xml()
        else:
            vasprun_etree = self._parse_etree_vasprun_xml(tag='varray')
            return self._get_forces(vasprun_etree)

    def read_force_constants(self):
        vasprun = self._parse_etree_vasprun_xml()
        return self._get_force_constants(vasprun)


    def _get_forces(self, vasprun_etree):
        """
        vasprun_etree = etree.iterparse(fileptr, tag='varray')
        """
        forces = []
        for event, element in vasprun_etree:
            if element.attrib['name'] == 'forces':
                for v in element:
                    forces.append([float(x) for x in v.text.split()])
        return np.array(forces)

    def _get_force_constants(self, vasprun_etree):
        fc_tmp = None
        num_atom = 0
        for event, element in vasprun_etree:
            if num_atom == 0:
                atomtypes = self._get_atomtypes(element)
                if atomtypes:
                    num_atoms, elements, elem_masses = atomtypes[:3]
                    num_atom = np.sum(num_atoms)
                    masses = []
                    for n, m in zip(num_atoms, elem_masses):
                        masses += [m] * n

            # Get Hessian matrix (normalized by masses)
            if element.tag == 'varray':
                if element.attrib['name'] == 'hessian':
                    fc_tmp = []
                    for v in element.findall('./v'):
                        fc_tmp.append([float(x) for x in v.text.strip().split()])

        if fc_tmp is None:
            return False
        else:
            fc_tmp = np.array(fc_tmp)
            if fc_tmp.shape != (num_atom * 3, num_atom * 3):
                return False
            # num_atom = fc_tmp.shape[0] / 3
            force_constants = np.zeros((num_atom, num_atom, 3, 3), dtype='double')

            for i in range(num_atom):
                for j in range(num_atom):
                    force_constants[i, j] = fc_tmp[i*3:(i+1)*3, j*3:(j+1)*3]

            # Inverse normalization by atomic weights
            for i in range(num_atom):
                for j in range(num_atom):
                    force_constants[i, j] *= -np.sqrt(masses[i] * masses[j])

            return force_constants, elements

    def _get_atomtypes(self, element):
        atom_types = []
        masses = []
        valences = []
        num_atoms = []

        if element.tag == 'array':
            if 'name' in element.attrib:
                if element.attrib['name'] == 'atomtypes':
                    for rc in element.findall('./set/rc'):
                        atom_info = [x.text for x in rc.findall('./c')]
                        num_atoms.append(int(atom_info[0]))
                        atom_types.append(atom_info[1].strip())
                        masses.append(float(atom_info[2]))
                        valences.append(float(atom_info[3]))
                    return num_atoms, atom_types, masses, valences

        return None

    def _parse_etree_vasprun_xml(self, tag=None):
        if self._is_version528():
            return self._parse_by_etree(VasprunWrapper(self._fileptr), tag=tag)
        else:
            return self._parse_by_etree(self._fileptr, tag=tag)

    def _parse_by_etree(self, fileptr, tag=None):
        try:
            import xml.etree.ElementTree as etree
            for event, elem in etree.iterparse(fileptr):
                if tag is None or elem.tag == tag:
                    yield event, elem
        except ImportError:
            print("Python 2.5 or later is needed.")
            print("For creating FORCE_SETS file with Python 2.4, you can use "
                  "phonopy 1.8.5.1 with python-lxml .")

    def _parse_by_expat(self, fileptr):
        vasprun = VasprunxmlExpat(fileptr)
        vasprun.parse()
        return vasprun.get_forces()[-1]

    def _is_version528(self):
        for line in self._fileptr:
            if '\"version\"' in str(line):
                self._fileptr.seek(0)
                if '5.2.8' in str(line):
                    sys.stdout.write(
                        "\n"
                        "**********************************************\n"
                        "* A special routine was used for VASP 5.2.8. *\n"
                        "**********************************************\n")
                    return True
                else:
                    return False


if __name__=='__main__':
   import xml.etree.ElementTree as ET
   vpr=Vasprun('vasprun.xml')
   vpr._is_version528()
   print vpr.read_forces()
   #print vpr.read_force_constants()
   #print vpr.read_atomtypes()
   vpr_etree = vpr._parse_etree_vasprun_xml()
   for event, element in vpr_etree:
       if element.tag=='incar':
          print '\n'.join(['{0:12s} {1:16s}'.format(a.attrib['name'],a.text.rstrip().rjust(16)) for a in element.findall('./i')])
       if element.tag=='eigenvalues':
          for nn in element.iter('kpoint'):
              print nn
          print element.attrib
          print [a.text for a in element.findall('kpoint')]

   vtree = ET.parse('vasprun.xml').getroot()
   for child in vtree:
       if child.tag =='structure' and child.attrib['name']=='initialpos':
          for gchild in child:
              if 'name' in gchild.attrib and gchild.attrib['name']=='positions':
                 print 'get it'
                 print [item.text for item in gchild]
