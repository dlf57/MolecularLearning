#!/usr/bin/env python

import glob
import numpy as np
import pybel
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from sklearn.linear_model import BayesianRidge
from sklearn.metrics import mean_absolute_error
import matplotlib.pyplot as plt
from matplotlib import style

# don't need to reimport openbabel
ob = pybel.ob

# syntax:
# molml.py

# organization of all paths for the code
csv_file = '/Users/dakota/Documents/Research/data-all_ordered.csv'
csv_file_columns = ['Name', 'Conformer', 'dftE', 'pm7E', 'mmffE']
jobs_directory = "/Users/dakota/Documents/Research/conformers/*jobs/*"
jobs_file = "/rmsd*.mol"
jobs_file_format = "mol"
figure_path = '/Users/dakota/Documents/Research/data/Figures/'
figure_name = 'figure_name.png'
save_figure = figure_path + figure_name

# Import the energy values
df_csv = pd.read_csv(csv_file)
df_csv.columns = csv_file_columns
df_csv = df_csv.replace('nan', np.nan, regex=True)
# If there is not a dftE than it is dropped so that
# . there are no errors of nan when training
df_csv = df_csv.dropna(axis=0, how='any')


def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


def encoding(atom):
    '''Atom encoding
    Each atom encoding is represented by a string of numbers
    consisting of an element encoding and a hybrid encoding.

    Element encoding is done based off of the periodic table
    where the element that will be encoded is encoded to the
    corresponding element number on the periodic table.

    Type encoding is done based off of the type of atom it is.
     _________________________
    | Element (XX) | Type (X) |

    Ex.
     ____________       ________
    | Carbon | 3 | --> | 06 | 3 | --> 063

    '''
    element_encoding = {'B': '05', 'C': '06', 'F': '09',
                        'I': '53', 'N': '07', 'O': '08',
                        'P': '15', 'S': '16', 'H': '01',
                        'Cl': '17', 'Br': '35'}
    type_encoding = {'1': '1', '2': '2', 'a': '2', '3': '3', '+': '4'}
    for atom_type in atom:
        atom_type = atom[0:2]
        if atom_type == (atom_type[0] + '1'):
            atom_encoding = (str(element_encoding[atom[0]]) +
                             str(type_encoding[atom[1]]))
        elif atom_type == (atom_type[0] + '2'):
            atom_encoding = (str(element_encoding[atom[0]]) +
                             str(type_encoding[atom[1]]))
        elif atom_type == (atom_type[0] + 'a'):
            atom_encoding = (str(element_encoding[atom[0]]) +
                             str(type_encoding[atom[1]]))
        elif atom_type == (atom_type[0] + '3'):
            atom_encoding = (str(element_encoding[atom[0]]) +
                             str(type_encoding[atom[1]]))
        elif atom_type == (atom_type[0] + '+'):
            atom_encoding = (str(element_encoding[atom[0]]) +
                             str(type_encoding[atom[1]]))
        elif atom_type == (atom_type[0] + 'l'):
            atom_encoding = str(element_encoding[atom[0:2]]) + '0'
        elif atom_type == (atom_type[0] + 'r'):
            atom_encoding = str(element_encoding[atom[0:2]]) + '0'
        else:
            atom_encoding = str(element_encoding[atom[0]]) + '0'
    return atom_encoding


# Pad lists to ensure same length
bond_len = 200
angle_len = 400
torsion_len = 600
nb_len = 7500
pm7_len = 2
fill_val = -99999

row_list = []
# Read through all the files in the folder of this directory
for directory in glob.iglob(jobs_directory):
    name = "/".join(directory.split('/')[6:10])  # name of the entry

    for files in glob.iglob(directory + jobs_file):
        conformer = files.split('/')[-1]  # conformer name/number
        conf = conformer.split('.')[0]  # this splits off conformer name
        name = str(name)
        conf = str(' ' + conf)

        # Use this for Python 2.7
        # mol = pybel.readfile('format', argument).next()
        # Use this for Python 3.6
        mol = next(pybel.readfile(jobs_file_format, files))

        # iterate through all atoms
        #  .. this is commented out because Bag Of Bonds doesn't use atomic charges
        # for atom in mol.atoms:
        #    print "Atom %d, %8.4f" % (atom.type, atom.partialcharge)

        # Running this try to ensure files and csv line up
        # . If the molecule file does not exist in the csv
        # . then it is skipped for reading in
        try:
            # pull energy out of csv
            energy = (df_csv[(df_csv['Name'] == name)
                             & (df_csv['Conformer'] == conf)].dftE.item())

            # pull pm7E info out of csv
            pm7e = (df_csv[(df_csv['Name'] == name)
                           & (df_csv['Conformer'] == conf)].pm7E.item())
            pm7e = [pm7e]
            pm7e.extend([fill_val])

            # iterate through all bonds
            bonds = []
            for bond in ob.OBMolBondIter(mol.OBMol):
                begin = atomType(mol, bond.GetBeginAtomIdx())
                end = atomType(mol, bond.GetEndAtomIdx())

                if (end < begin):
                    # swap them for lexographic order
                    begin, end = end, begin

                bond_type = encoding(begin) + encoding(end)
                bond_length = ("%8.4f" % (bond.GetLength()))
                bonds.append(float(bond_type))
                bonds.append(float(bond_length))

            # check to see if length does not exceed index
            if (len(bonds) > bond_len):
                raise Exception('{} bonds. Change padding.'.format(len(bonds)))

            bonds.extend([fill_val] * (bond_len - len(bonds)))

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

                angle_type = (encoding(aType) + encoding(b.GetType())
                              + encoding(cType))
                angle_angle = ("%8.3f" % (b.GetAngle(a, c)))
                angles.append(float(angle_type))
                angles.append(float(angle_angle))

            # check to see if length does not exceed index
            if (len(angles) > angle_len):
                raise Exception(
                    '{} angles. Change padding.'.format(len(angles)))

            angles.extend([fill_val] * (angle_len - len(angles)))

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
                    torsion_type = (encoding(aType) + encoding(bType)
                                    + encoding(cType) + encoding(dType))
                    torsion_angle = ("%8.3f" %
                                     ((mol.OBMol.GetTorsion(a, b, c, d))))
                    torsions.append(float(torsion_type))
                    torsions.append(float(torsion_angle))
                else:
                    torsion_type = (encoding(dType) + encoding(cType)
                                    + encoding(bType) + encoding(aType))
                    torsion_angle = ("%8.3f" %
                                     ((mol.OBMol.GetTorsion(a, b, c, d))))
                    torsions.append(float(torsion_type))
                    torsions.append(float(torsion_angle))

            # check to see if length does not exceed index
            if (len(torsions) > torsion_len):
                raise Exception(
                    '{} torsions. Change padding.'.format(len(torsions)))

            torsions.extend([fill_val] * (torsion_len - len(torsions)))

            # iterate through all non-bonded
            nb = []
            for pair in ob.OBMolPairIter(mol.OBMol):
                (first, second) = pair
                begin_nb = atomType(mol, first)
                end_nb = atomType(mol, second)
                dist = mol.OBMol.GetAtom(first).GetDistance(second)
                if (end_nb[0] < begin_nb[0]):
                    # swap them for lexographic order
                    begin_nb, end_nb = end_nb, begin_nb

                nonbond_type = encoding(begin_nb) + encoding(end_nb)
                nonbond_distance = ("%8.4f" % (dist))
                nb.append(float(nonbond_type))
                nb.append(float(nonbond_distance))

            # check to see if length does not exceed index
            if (len(nb) > nb_len):
                raise Exception(
                    '{} nonbonding. Change padding.'.format(len(nb)))

            nb.extend([fill_val] * (nb_len - len(nb)))

            dict1 = {}
            dict1.update({'Bond': bonds})
            dict1.update({'Angle': angles})
            dict1.update({'Torsion': torsions})
            dict1.update({'NB': nb})
            dict1.update({'Name': name})
            dict1.update({'Conformer': conf})
            dict1.update({'Energy': energy})
            dict1.update({'PM7E': pm7e})
            row_list.append(dict1)

        except ValueError:
            pass


# make a dataframe of all the molecules, their descriptors and energies
df = pd.DataFrame(row_list, columns=['Name', 'Conformer',
                                     'Energy',
                                     'PM7E',
                                     'Bond',
                                     'Angle',
                                     'Torsion',
                                     'NB'])


# make an array of all the parts of the molecular descriptors
bonds = np.asarray(list(df['Bond']), dtype=np.float)
angles = np.asarray(list(df['Angle']), dtype=np.float)
torsions = np.asarray(list(df['Torsion']), dtype=np.float)
nbs = np.asarray(list(df['NB']), dtype=np.float)
pm7 = np.asarray(list(df['PM7E']), dtype=np.float)
energy = np.asarray(list(df['Energy']), dtype=np.float)

# make the molecular descriptor that will be used as the input
molecule = np.concatenate((bonds, angles, torsions, nbs, pm7), axis=1)
# print(molecule.shape) # if array size errors print to see shape

# test, train, split for data
X_train, X_test, y_train, y_test = train_test_split(
    molecule, energy, test_size=0.2)
# print(X_train.shape) # if array size errors print to see shape
# print(y_train.shape) # if array size errors print to see shape

# specify regressor and run
clf = BayesianRidge()
clf.fit(X_train, y_train)
accuracy = clf.score(X_test, y_test)
r2 = ("%8.3f" % accuracy)
print('r^2 =', accuracy)

predicted = clf.predict(X_test)

mae = mean_absolute_error(y_test, predicted)
print('Mean Absolute Error =', mae)


# plot the actual versus the predicted
plt.figure(figsize=(8, 6), dpi=300)
plt.scatter(predicted, y_test)
plt.annotate('$r^2$=' + str(r2), xy=(0.97, 0.10),
             xycoords='axes fraction', ha='right')
plt.annotate('$MAE$=' + str("%8.4f" % mae), xy=(0.97, 0.05),
             xycoords='axes fraction', ha='right')
plt.title('Comparison of DFT Energies', fontsize=16)
plt.xlabel('Predicted', fontsize=12)
plt.ylabel('Actual', fontsize=12)
plt.grid(False)
plt.savefig(save_figure)
plt.show()
plt.close
