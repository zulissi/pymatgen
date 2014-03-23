#!/usr/bin/env python
from __future__ import division
import tempfile
from monty.tempfile import ScratchDir
from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
import pybel as pb
import os
from subprocess import Popen, PIPE
import numpy as np


PACKMOL_DEBUG = False

class PackmolRunner(object):
    """
    Create MD simulation box using packmol.
    """

    def __init__(self, mols, param_list):
        """
        Create PackMolRunner
        :param mols: [Molecule] - list of Molecules to pack
        :param param_list: [{}, {}] - array of parameters containing dicts for each Structure
        """
        self.mols = mols
        self.param_list = param_list

    def _format_packmol_str(self, some_obj):
        """
        Internal method to format packmol strings
        :param some_obj: Some object to turn into String
        :return:
        """
        if isinstance(some_obj,list):
            return ' '.join(str(x) for x in some_obj)
        else:
            return some_obj


    def run(self):
        """
        Runs packmol
        :return: a Molecule object
        """

        # open temp files
        scratch = tempfile.gettempdir()
        with ScratchDir(scratch, copy_to_current_on_exit=PACKMOL_DEBUG) as d:
        #with ScratchDir(scratch, copy_to_current_on_exit=True) as d:

            # calculate box size for mixture
            volume = 0.0
            for idx, m in enumerate(self.mols):
                lx,ly,lz = get_auto_boxsize(m,self.param_list[idx]['number'])
                volume = volume + lx*ly*lz
            boxlength = volume**(1/3)

            for idx, m in enumerate(self.mols):
                use_auto_box = True
                for key in self.param_list[idx]:
                    if key.endswith(' box'):
                        use_auto_box = False

                if use_auto_box:
                    #self.param_list[idx]['inside box'] = [0, 0, 0] + get_auto_boxsize(m, self.param_list[idx]['number'])
                    # use box size for mixture to make sure box is large enough
                    self.param_list[idx]['inside box'] = [0,0,0,boxlength,boxlength,boxlength]

            # convert mols to pdb files
            for idx, m in enumerate(self.mols):
                a = BabelMolAdaptor(m)
                pm = pb.Molecule(a.openbabel_mol)
                pm.write("pdb", filename=os.path.join(d, '{}.pdb'.format(idx)), overwrite=True)

            # write packmol input file
            with open(os.path.join(d, 'pack.inp'), 'w') as inp:
                inp.write('tolerance 2.0\n')
                inp.write('filetype pdb\n')
                inp.write('output {}\n'.format(os.path.join(d, "box.pdb")))
                for idx, m in enumerate(self.mols):
                    inp.write('\n')
                    inp.write('structure {}.pdb\n'.format(os.path.join(d, str(idx))))

                    for k, v in self.param_list[idx].iteritems():
                        inp.write('  {} {}\n'.format(k, self._format_packmol_str(v)))
                    inp.write('end structure\n')

            # run packmol
            proc = Popen(['packmol'], stdin=open(os.path.join(d, 'pack.inp'), 'r'),stdout=PIPE)
            (stdout, stderr) = proc.communicate()
#            print stdout
#            print stderr


            a = BabelMolAdaptor.from_file(os.path.join(d, "box.pdb"), "pdb")
            return a.pymatgen_mol

def get_auto_boxsize(molecule, number, tolerance=2):
    """
    Calculate box length
    :param molecule - a Molecule
    :param number - number of molecules to pack (?)
    :param tolerance - tolerance for length
    :return box lengths as 3-element array
    """

    lx, ly, lz = np.max(molecule.cart_coords, 0)-np.min(molecule.cart_coords, 0)
    max_length = max(lx, ly, lz) + tolerance
    return [(max_length**3*number)**(1/3)] * 3


if __name__ == '__main__':
    coords = [[0.000000, 0.000000, 0.000000],
           [0.000000, 0.000000,  1.089000],
           [1.026719, 0.000000, -0.363000],
           [-0.513360, -0.889165, -0.363000],
           [-0.513360, 0.889165, -0.363000]]
    mol = Molecule(["C", "H", "H", "H", "H"], coords)
    #pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5, "inside box":[0.,0.,0.,40.,40.,40.]}])
    pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5}])
    #pmr = PackmolRunner([mol, mol], [{"number":4},{"number":5}])
    s = pmr.run()
#    print s