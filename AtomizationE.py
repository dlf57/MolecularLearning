import sys, os
from os.path import basename
import pybel
# don't need to reimport openbabel
ob = pybel.ob

# syntax:
# molml.py [files]

def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()

# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders
for argument in sys.argv[1:]:
    filename, extension = os.path.splitext(argument)
    name = os.path.basename(argument)
    print name
    #print(os.path.basename(argument))

    # read the molecule from the supplied file
    mol = next(pybel.readfile(extension[1:], argument))
    #mol = next(pybel.readfile("out", "rmsd28_opt.out"))
    print mol.energy # in kcal/mol
    print #should probably take this out later but for now gives spacing