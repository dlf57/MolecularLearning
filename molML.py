#!/usr/bin/env python

import sys, os
import numpy as np
import pybel
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.kernel_ridge import KernelRidge
# don't need to reimport openbabel
ob = pybel.ob

# syntax:
# molml.py

def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()

row_list = []

# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders
for argument in sys.argv[1:]:
    filename, extension = os.path.splitext(argument)
    # Include the filename as to know which file is being read
    name = os.path.basename(argument)

    # read the molecule from the supplied file
    try:
        # Use this for Python 2.7
        # mol = pybel.readfile(extension[1:], argument).next()
        # Use this for Python 3.6
        mol = next(pybel.readfile(extension[1:], argument))
    except:
        pass


    # In this case I do not include Energy as that is our dependent variable
    #print mol.energy # in kcal/mol
    # ideally, we should turn this into an atomization energy\
    energy = [mol.energy]
    # energy = np.array(mol.energy)


    # iterate through all atoms
    #  .. this is commented out because Bag Of Bonds doesn't use atomic charges
    # for atom in mol.atoms:
    #    print "Atom %d, %8.4f" % (atom.type, atom.partialcharge)

    # iterate through all bonds
    bonds = []
    for bond in ob.OBMolBondIter(mol.OBMol):
        begin = atomType(mol, bond.GetBeginAtomIdx())
        end = atomType(mol, bond.GetEndAtomIdx())
        if (end < begin):
            # swap them for lexographic order
            begin, end = end, begin
        # bonds.append("Bond %s-%s, %8.4f" % (begin, end, bond.GetLength()) )
        # bonds.append("%s-%s, %8.4f" % (begin, end, bond.GetLength()) )
        bonds.append("%8.4f" % (bond.GetLength()))
        # print (bonds[-1])
        # np.array(bonds)
    bonds = bonds + ([None] * 98)
    bonds = np.asarray(bonds, dtype=float)
    bonds[np.isnan(bonds)] = -99999


    # iterate through all angles
    angles = []
    for angle in ob.OBMolAngleIter(mol.OBMol):
        a = (angle[0] + 1)
        b = mol.OBMol.GetAtom(angle[1] + 1)
        c = (angle[2] + 1)

        aType = atomType(mol, a)
        cType = atomType(mol, c)
        if (cType < aType):
            # swap them for lexographic order
            aType, cType = cType, aType
        # angles.append( "Angle %s-%s-%s, %8.3f" % (aType, b.GetType(), cType, b.GetAngle(a, c)) )
        # angles.append( "%s-%s-%s, %8.3f" % (aType, b.GetType(), cType, b.GetAngle(a, c)) )
        angles.append("%8.3f" % (b.GetAngle(a, c)))
        #print (angles[-1])


    # iterate through all torsions
    torsions = []
    for torsion in ob.OBMolTorsionIter(mol.OBMol):
        a = (torsion[0] + 1)
        b = (torsion[1] + 1)
        c = (torsion[2] + 1)
        d = (torsion[3] + 1)

        aType = atomType(mol, a)
        bType = atomType(mol, b)
        cType = atomType(mol, c)
        dType = atomType(mol, d)

        # output in lexographic order
        if (aType < dType):
            # torsions.append( "Torsion %s-%s-%s-%s, %8.3f" % (aType, bType, cType, dType, mol.OBMol.GetTorsion(a, b, c, d)) )
            # torsions.append( "%s-%s-%s-%s, %8.3f" % (aType, bType, cType, dType, mol.OBMol.GetTorsion(a, b, c, d)) )
            torsions.append( "%8.3f" % (mol.OBMol.GetTorsion(a, b, c, d)) )
        else:
            # torsions.append( "Torsion %s-%s-%s-%s, %8.3f" % (dType, cType, bType, aType, mol.OBMol.GetTorsion(a, b, c, d)) )
            # torsions.append( "%s-%s-%s-%s, %8.3f" % (dType, cType, bType, aType, mol.OBMol.GetTorsion(a, b, c, d)) )
            torsions.append("%8.3f" % (mol.OBMol.GetTorsion(a, b, c, d)))
            #print (torsions[-1])
    torsions = np.asarray(torsions, dtype=float)


    dict1 = {}
    dict1.update({'Name': name})
    dict1.update({'Bonds': bonds})
    dict1.update({'Angles': angles})
    dict1.update({'Torsions': torsions})
    dict1.update({'Energy': energy})
    row_list.append(dict1)

df = pd.DataFrame(row_list, columns=['Name','Bonds','Angles','Torsions', 'Energy'])

molecules = df['Bonds'] #+ df['Torsions']
Energy = df['Energy']


X =  np.asarray(list(molecules), dtype=np.float)
y = np.asarray(list(Energy), dtype=np.float)


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

clf = KernelRidge()
clf.fit(X_train, y_train)
accuracy = clf.score(X_test,y_test)
print(accuracy)
