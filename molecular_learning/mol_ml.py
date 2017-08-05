#!/usr/bin/env python

import glob
import numpy as np
import pybel
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.kernel_ridge import KernelRidge
import matplotlib.pyplot as plt
from matplotlib import style
# don't need to reimport openbabel
ob = pybel.ob

# syntax:
# molml.py


def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


row_list = []

# Read through all the files in the folder of this directory
for d in glob.iglob("*omegacsd_YAGWOS/*"):
    name = "/".join(d.split('/')[0:2]) # name of the entry

    for f in glob.iglob(d + "/rmsd*.out"):
        # all the base files
        conf = f.split('/')[-1] # conformer name/number

        try:
            # Use this for Python 2.7
            # mol = pybel.readfile('format', argument).next()
            # Use this for Python 3.6
            # readfile(format, filename)
            mol = next(pybel.readfile('out', f))
        except:
            pass

        # In this case I do not include Energy as that is our dependent variable
        # print mol.energy # in kcal/mol
        # ideally, we should turn this into an atomization energy
        energy = [mol.energy]

        # iterate through all atoms
        #  .. this is commented out because Bag Of Bonds doesn't use atomic charges
        # for atom in mol.atoms:
        #    print "Atom %d, %8.4f" % (atom.type, atom.partialcharge)

        bonds = []
        for bond in ob.OBMolBondIter(mol.OBMol):
            begin = atomType(mol, bond.GetBeginAtomIdx())
            end = atomType(mol, bond.GetEndAtomIdx())
            if (end < begin):
                # swap them for lexographic order
                begin, end = end, begin
            # bonds.append("Bond %s-%s, %8.4f" %
            #              (begin, end, bond.GetLength()))
            # bonds.append("%s-%s, %8.4f" %
            #              (begin, end, bond.GetLength()))
            bonds.append("%8.4f" % (bond.GetLength()))
            # print(bonds[-1])

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
            # angles.append("%s-%s-%s, %8.3f" %
            #                (aType, b.GetType(), cType, b.GetAngle(a, c)))
            angles.append("%8.3f" % (b.GetAngle(a, c)))
            # print(angles[-1])

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
                # torsions.append( "%s-%s-%s-%s, %8.3f" %
                #                  (aType, bType, cType, dType,
                #                   mol.OBMol.GetTorsion(a, b, c, d)))
                torsions.append("%8.3f" % (mol.OBMol.GetTorsion(a, b, c, d)))
            else:
                # torsions.append( "Torsion %s-%s-%s-%s, %8.3f" %
                #                  (dType, cType, bType, aType,
                #                   mol.OBMol.GetTorsion(a, b, c, d)))
                # torsions.append( "%s-%s-%s-%s, %8.3f" %
                #                  (dType, cType, bType, aType,
                #                   mol.OBMol.GetTorsion(a, b, c, d)))
                torsions.append("%8.3f" % (mol.OBMol.GetTorsion(a, b, c, d)))
                # print(torsions[-1])
        torsions = np.asarray(torsions, dtype=float)

        bond_difference = len(torsions) - len(bonds)
        bonds = bonds + ([None] * bond_difference)
        bonds = np.asarray(bonds, dtype=float)
        bonds[np.isnan(bonds)] = -99999

        angle_difference = len(torsions) - len(angles)
        angles = angles + ([None] * angle_difference)
        angles = np.asarray(angles, dtype=float)
        angles[np.isnan(angles)] = -99999

        # store the bonds angles and torsions as molecule descriptors
        molecule_descriptor = np.concatenate((bonds, angles, torsions))

        # store the descriptors for the molecules in a dictionary
        dict1 = {}
        dict1.update({'Molecule': molecule_descriptor})
        dict1.update({'Name': name})
        dict1.update({'Conformer': conf})
        dict1.update({'Energy': energy})
        row_list.append(dict1)


# make a dataframe of all the molecules, their descriptors and energies
df = pd.DataFrame(row_list, columns=['Name', 'Conformer', 'Molecule', 'Energy'])
molecules = df['Molecule']
Energy = df['Energy']

X = np.asarray(list(molecules), dtype=np.float)
y = np.asarray(list(Energy), dtype=np.float)

# split the molecules up into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
# print (X_train.shape)
# print (y_train.shape)

clf = KernelRidge()
clf.fit(X_train, y_train)
accuracy = clf.score(X_test, y_test)
print(accuracy)

predicted = clf.predict(X_test)

# dataframe to organize the predicted and actual values
df_predicted = pd.DataFrame(predicted, columns=['Predicted'])
df_actual = pd.DataFrame(y_test, columns=['Actual'])
df_compare = pd.concat([df_predicted, df_actual], axis=1)
print(df_compare)

# plotting the actual vs. the predicted to get a visual representation
plt.figure(figsize=(8, 6))
plt.scatter(y_test, predicted)
plt.title('Comparison of Energies (Kcal/mol)', fontsize=16)
plt.xlabel('Actual', fontsize=12)
plt.ylabel('Predicted', fontsize=12)
plt.grid(True)
plt.show()
