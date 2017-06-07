#!/usr/bin/env python

import sys, os
import sqlite3
import numpy as np
import pybel
# don't need to reimport openbabel
ob = pybel.ob

# syntax:
# molml.py

conn = sqlite3.connect('MolecularData.db')
c = conn.cursor()

def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


# repeat through all the files on the command-line
# we can change this to use the glob module as well
#  e.g., find all the files in a set of folders
for argument in sys.argv[1:]:
    c.execute('CREATE TABLE IF NOT EXISTS MoleculeData (Name TEXT, Bonds REAL, Angles REAL, Torsions REAL, Energy REAL)')

    filename, extension = os.path.splitext(argument)
    # print #should probably take this out later but for now gives spacing
    # Include the filename as to know which file is being read
    name = os.path.basename(argument)
    # print name


    # read the molecule from the supplied file
    mol = pybel.readfile(extension[1:], argument).next()
    #mol = next(pybel.readfile(extension[1:], argument))
    # DLF used in order to view initial outputs (mostly used in Jupyter for quick looks)
    # mol = next(pybel.readfile("out", "rmsd28_opt.out"))



    # In this case I do not include Energy as that is our dependent variable
    #print mol.energy # in kcal/mol
    # ideally, we should turn this into an atomization energy
    energy = mol.energy

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
        bonds.append("Bond %s-%s, %8.4f" % (begin, end, bond.GetLength()) )
        #print bonds[-1]
    # Would representing the data as an array make it easier to manipulate?
    # Bonds = np.asarray([bonds], dtype=object)
    bond = ' , '.join(bonds)
    # print bond

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
        angles.append( "Angle %s-%s-%s, %8.3f" % (aType, b.GetType(), cType, b.GetAngle(a, c)) )
        #print angles[-1]
    #Angles = np.asarray([angles], dtype=object)
    angle = ' , '.join(angles)
    # print angle

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
            torsions.append( "Torsion %s-%s-%s-%s, %8.3f" % (aType, bType, cType, dType, mol.OBMol.GetTorsion(a, b, c, d)) )
        else:
            torsions.append( "Torsion %s-%s-%s-%s, %8.3f" % (dType, cType, bType, aType, mol.OBMol.GetTorsion(a, b, c, d)) )
            #print torsions[-1]
    #Torsions = np.asarray([torsions], dtype=object)
    torsion = ' , '.join(torsions)
    # print torsion


    c.execute("INSERT INTO MoleculeData (Name, Bonds, Angles, Torsions, Energy) VALUES (?, ?, ?, ?)",
              (name, bond, angle, torsion, energy))
    conn.commit()
    #
    # data = c.fetchall()
    # print data

c.close()
conn.close()
