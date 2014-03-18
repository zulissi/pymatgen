#!/usr/bin/env python
import tempfile
from monty.tempfile import ScratchDir
from pymatgen import Molecule
from pymatgen.io.babelio import BabelMolAdaptor
import pybel as pb
import os
from subprocess import Popen, PIPE
import numpy as np


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


#    def _get_auto_boxsize(self):
#        """
#        TODO: Put docs here!!
#        :param idx:
#        :raise NotImplementedError:
#        """
#
#        volume=0.0
#        for idx,mol in enumerate(self.mols):
#            lx,ly,lz = np.max(mol.cart_coords,0)-np.min(mol.cart_coords,0)
##            print lx,ly,lz
#            lx += 2.0  # assume tolerance = 2.0
#            ly += 2.0
#            lz += 2.0
#            length=max(lx,ly,lz)
#            volume += length**(3.0)*float(self.param_list[idx]['number'])
#        length = volume**(1.0/3.0)
##        print length
#
#        for idx,mol in enumerate(self.mols):
#            self.param_list[idx]['inside box']='0.0 0.0 0.0 {} {} {}'.format(length,length,length)
##        raise NotImplementedError('Auto box size is not implemented yet!')


    def run(self):
        """
        Runs packmol
        :return: a Molecule object
        """

        scratch = tempfile.gettempdir()
        with ScratchDir(scratch,copy_to_current_on_exit=True) as d:
        #with ScratchDir(scratch) as d:
            # convert mols to pdb files
            for idx, mol in enumerate(self.mols):
                a = BabelMolAdaptor(mol)
                pm = pb.Molecule(a.openbabel_mol)
                pm.write("pdb", filename=os.path.join(d, '{}.pdb'.format(idx)), overwrite=True)

            # TODO: also check if user specified outside box, etc.
            # Do not use auto mode if user specified any type of box

            # ignore user input for now, always calculate box size
            #self._get_auto_boxsize()

            with open(os.path.join(d, 'pack.inp'), 'w') as inp:
                # create packmol control file
                inp.write('tolerance 2.0\n')
                inp.write('filetype pdb\n')
                inp.write('output {}\n'.format(os.path.join(d, "box.pdb")))
                for idx, mol in enumerate(self.mols):
                    inp.write('\n')
                    inp.write('structure {}.pdb\n'.format(os.path.join(d, str(idx))))

                    #ignore box size from user input for now
                    boxlength = get_auto_boxsize(mol,self.param_list[idx]['number'])
                    #print boxlength
                    self.param_list[idx]['inside box']='0.0 0.0 0.0 {} {} {}'.format(boxlength,boxlength,boxlength)
                    for k, v in self.param_list[idx].iteritems():
                        inp.write('  {} {}\n'.format(k, self._format_packmol_str(v)))
                    inp.write('end structure\n')

            proc = Popen(['packmol'], stdin=open(os.path.join(d, 'pack.inp'), 'r'),stdout=PIPE)
            (stdout, stderr) = proc.communicate()
#            print stdout
#            print stderr


            a = BabelMolAdaptor.from_file(os.path.join(d, "box.pdb"), "pdb")
            return a.pymatgen_mol

def get_auto_boxsize(mol,number):
    """
    TODO: Put docs here!!
    :param mol:
    :raise NotImplementedError:
    """

    #print number
    volume=0.0
    lx,ly,lz = np.max(mol.cart_coords,0)-np.min(mol.cart_coords,0)
    #print lx,ly,lz
    lx += 2.0  # assume tolerance = 2.0
    ly += 2.0
    lz += 2.0
    length=max(lx,ly,lz)
    volume += length**(3.0)*float(number)
    length = volume**(1.0/3.0)
    #print length
    return length


if __name__ == '__main__':
    coords = [[0.000000, 0.000000, 0.000000],
           [0.000000, 0.000000,  1.089000],
           [1.026719, 0.000000, -0.363000],
           [-0.513360, -0.889165, -0.363000],
           [-0.513360, 0.889165, -0.363000]]
    mol = Molecule(["C", "H", "H", "H", "H"], coords) 
#    pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5, "inside box":[0.,0.,0.,40.,40.,40.]}])
    pmr = PackmolRunner([mol, mol], [{"number":4,"inside box":[0.,0.,0.,40.,40.,40.]}, {"number":5}])
    s = pmr.run()
#    print s


