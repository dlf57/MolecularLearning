#!/usr/bin/env python

import sys, os
import numpy as np
import pybel
import pandas as pd
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.kernel_ridge import KernelRidge
from sklearn.feature_extraction import DictVectorizer
from sklearn.feature_extraction import DictVectorizer
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
    # print name

    # read the molecule from the supplied file
    mol = pybel.readfile(extension[1:], argument).next()

    # In this case I do not include Energy as that is our dependent variable
    #print mol.energy # in kcal/mol
    # ideally, we should turn this into an atomization energy\
    energy = mol.energy
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
        # print bonds[-1]
        # np.array(bonds)


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
        #print angles[-1]


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
            #print torsions[-1]


    dict1 = {}
    dict1.update({'Name' : name})
    dict1.update({'Bonds' : bonds})
    dict1.update({'Angles' : angles})
    dict1.update({'Torsions' : torsions})
    dict1.update({'Energy' : energy})
    row_list.append(dict1)

df = pd.DataFrame(row_list, columns=['Name','Bonds','Angles','Torsions', 'Energy'])
df1 = pd.DataFrame(row_list, columns=['Bonds','Angles','Torsions'])
# molecule = np.array(df['Bonds'] + df['Angles'] + df['Torsions'])
# df['Molecule'] = df['Bonds'] + df['Angles'] + df['Torsions']
# molecule = df['Molecule']
# molecule = np.array(molecule)
# bond = df['Bonds']
# bond = np.array(bond)
# print df['Bonds'] #molecule
Energy = df['Energy']
# Energy = np.array(Energy)
# print Energy

X = df1.as_matrix()
y = Energy.as_matrix()
# print len(X), len(y)
# print X, y

# df.dropna(inplace=True)
# y = np.array(df['Energy'])


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)


clf = KernelRidge()

clf.fit(X_train, y_train)
accuracy = clf.score(X_test,y_test)

print(accuracy)
