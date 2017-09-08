#!/usr/bin/env python

import glob
import sqlite3
import numpy as np
import pybel
# don't need to reimport openbabel
ob = pybel.ob

# syntax:
# molml.py

conn = sqlite3.connect('MolecularData.db')
cursor = conn.cursor()


def create_table():
    cursor.execute(
        'CREATE TABLE IF NOT EXISTS MoleculeData (Name TEXT, Conformer TEXT, Bonds REAL, Angles REAL, Torsions REAL)')


def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


# Read through all the files in the folder of this directory
for directory in glob.iglob("/Users/dakota/Documents/Research/conformers/*jobs/*"):
    name = "/".join(directory.split('/')[6:10])  # name of the entry

    for files in glob.iglob(directory + "/rmsd*.mol"):
        conf = files.split('/')[-1]  # conformer name/number

        try:
            # Use this for Python 2.7
            # mol = pybel.readfile('format', argument).next()
            # Use this for Python 3.6
            # readfile(format, filename)
            mol = next(pybel.readfile('mol', files))
        except:
            pass

        # In this case I do not include Energy as that is our dependent variable
        # print mol.energy # in kcal/mol
        # ideally, we should turn this into an atomization energy
        # energy = [mol.energy]

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
            # bonds.append("Bond %s-%s, %8.4f" %
            #              (begin, end, bond.GetLength()))
            bonds.append("%s-%s, %8.4f" %
                         (begin[0], end[0], bond.GetLength()))
            # bonds.append("%8.4f" % (bond.GetLength()))
            # print(bonds[-1])
        bond = '; '.join(bonds)

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
            # angles.append("Angle %s-%s-%s, %8.3f" %
            #                (aType, b.GetType(), cType, b.GetAngle(a, c)))
            angles.append("%s-%s-%s, %8.3f" %
                          (aType[0], b.GetType()[0], cType[0], b.GetAngle(a, c)))
            # angles.append("%8.3f" % (b.GetAngle(a, c)))
            # print(angles[-1])
        angle = '; '.join(angles)

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
                # torsions.append( "Torsion %s-%s-%s-%s, %8.3f" %
                #                  (aType, bType, cType, dType,
                #                   mol.OBMol.GetTorsion(a, b, c, d)) )
                torsions.append("%s-%s-%s-%s, %8.3f" %
                                (aType[0], bType[0], cType[0], dType[0],
                                 mol.OBMol.GetTorsion(a, b, c, d)))
                # torsions.append("%8.3f" % (mol.OBMol.GetTorsion(a, b, c, d)))
            else:
                # torsions.append( "Torsion %s-%s-%s-%s, %8.3f" %
                #                  (dType, cType, bType, aType,
                #                   mol.OBMol.GetTorsion(a, b, c, d)))
                torsions.append("%s-%s-%s-%s, %8.3f" %
                                (dType[0], cType[0], bType[0], aType[0],
                                 mol.OBMol.GetTorsion(a, b, c, d)))
                # torsions.append("%8.3f" % (mol.OBMol.GetTorsion(a, b, c, d)))
                # print(torsions[-1])
        torsion = '; '.join(torsions)

        cursor.execute("INSERT INTO MoleculeData (Name, Conformer, Bonds, Angles, Torsions) VALUES (?, ?, ?, ?, ?)",
                       (name, conf, bond, angle, torsion))
        conn.commit()

    # data = c.fetchall()
    # print data

create_table()
cursor.close()
conn.close()
