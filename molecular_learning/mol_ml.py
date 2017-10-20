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
# mol_ml.py

# organization of all paths for the code
csv_file = '/Users/dakota/Documents/Research/data-all_ordered.csv'
csv_file_columns = ['Name', 'Conformer', 'dftE', 'pm7E', 'mmffE']
jobs_directory = "/Users/dakota/Documents/Research/conformers/*jobs/*"
jobs_file = "/rmsd*.mol"
jobs_file_format = "mol"
figure_path = '/Users/dakota/Documents/Research/data/Figures/'
figure_name = 'figure_name.png'


# import the energy values
df_csv = pd.read_csv(csv_file)
df_csv.columns = csv_file_columns
df_csv = df_csv.replace('nan', np.nan, regex=True)
df_csv = df_csv.dropna(axis=0, how='any')
path_csv = df_csv['Name'] + df_csv['Conformer']


def atomType(mol, atomIdx):
    # get the atomic type given an atom index
    return mol.OBMol.GetAtom(atomIdx).GetType()


# Atom encoding
molec_encoding = {'B': 1, 'C': 2, 'F': 3, 'I': 4,
                  'N': 5, 'O': 6, 'P': 7, 'S': 8, 'H': 9}

# Indexing for dataframes to ensure same length
index = list(range(1, 501))

row_list = []
# Read through all the files in the folder of this directory
for directory in glob.iglob(jobs_directory):
    name = "/".join(directory.split('/')[6:10])  # name of the entry

    for files in glob.iglob(directory + jobs_file):
        conformer = files.split('/')[-1]  # conformer name/number
        opt_type = files.split('-')[-1]
        conf = conformer.split('.')[0]
        opt = opt_type.split('.')[0]
        path = name + ' ' + conf
        for value in path_csv.isin([path]):
            if value == True:
                try:
                    # Use this for Python 2.7
                    # mol = pybel.readfile('format', argument).next()
                    # Use this for Python 3.6
                    # readfile(format, filename)
                    mol = next(pybel.readfile(jobs_file_format, files))
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

                bonds = []
                for bond in ob.OBMolBondIter(mol.OBMol):
                    begin = atomType(mol, bond.GetBeginAtomIdx())[0]
                    end = atomType(mol, bond.GetEndAtomIdx())[0]
                    if (molec_encoding[end] < molec_encoding[begin]):
                        # swap them for lexographic order
                        begin, end = end, begin
                    dict2 = {}
                    bond_type = (
                        (str(molec_encoding[begin]) + str(molec_encoding[end])))
                    bond_length = ("%8.4f" % (bond.GetLength()))
                    dict2.update({'Bond_Type': bond_type})
                    dict2.update({'Bond_Length': bond_length})
                    bonds.append(dict2)

                dfb = pd.DataFrame(
                    bonds, columns=['Bond_Type', 'Bond_Length'], dtype=float)
                dfb = dfb.reindex(index, fill_value=-99999)

                # iterate through all angles
                angles = []
                for angle in ob.OBMolAngleIter(mol.OBMol):
                    a = (angle[0] + 1)
                    b = mol.OBMol.GetAtom(angle[1] + 1)
                    c = (angle[2] + 1)

                    aType = atomType(mol, a)
                    cType = atomType(mol, c)
                    if (molec_encoding[cType[0]] < molec_encoding[aType[0]]):
                        # swap them for lexographic order
                        aType, cType = cType, aType

                    dict3 = {}
                    angle_type = (str(molec_encoding[aType[0]]) + str(molec_encoding[b.GetType()[0]])
                                  + str(molec_encoding[cType[0]]))
                    angle_angle = ("%8.3f" % (b.GetAngle(a, c)))
                    dict3.update({'Angle_Type': angle_type})
                    dict3.update({'Angle': angle_angle})
                    angles.append(dict3)

                dfa = pd.DataFrame(
                    angles, columns=['Angle_Type', 'Angle'], dtype=float)
                dfa = dfa.reindex(index, fill_value=-99999)

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

                    dict4 = {}
                    # output in lexographic order
                    if (molec_encoding[aType[0]] < molec_encoding[dType[0]]):
                        torsion_type = (str(molec_encoding[aType[0]]) + str(molec_encoding[bType[0]])
                                        + str(molec_encoding[cType[0]]) + str(molec_encoding[dType[0]]))
                        torsion_angle = ("%8.3f" %
                                         ((mol.OBMol.GetTorsion(a, b, c, d))))
                        dict4.update({'Torsion_Type': torsion_type})
                        dict4.update({'Torsion': torsion_angle})
                        torsions.append(dict4)
                    else:
                        torsion_type = (str(molec_encoding[dType[0]]) + str(molec_encoding[cType[0]])
                                        + str(molec_encoding[bType[0]]) + str(molec_encoding[aType[0]]))
                        torsion_angle = ("%8.3f" %
                                         ((mol.OBMol.GetTorsion(a, b, c, d))))
                        torsion_angle = float(torsion_angle)
                        dict4.update({'Torsion_Type': torsion_type})
                        dict4.update({'Torsion': torsion_angle})
                        torsions.append(dict4)

                dftor = pd.DataFrame(torsions, columns=[
                                     'Torsion_Type', 'Torsion'], dtype=float)
                dftor = dftor.reindex(index, fill_value=-99999)

                # store the descriptors for the molecules in a dictionary
                dict1 = {}
                dict1.update({'Bond_Type': dfb['Bond_Type']})
                dict1.update({'Bond_Length': dfb['Bond_Length']})
                dict1.update({'Angle_Type': dfa['Angle_Type']})
                dict1.update({'Angle': dfa['Angle']})
                dict1.update({'Torsion_Type': dftor['Torsion_Type']})
                dict1.update({'Torsion': dftor['Torsion']})
                dict1.update({'Name': name})
                dict1.update({'Conformer': conf})
                dict1.update({'Optimization': opt})
                # Comment below out if energy is not included in file
                # dict1.update({'Energy': energy})
                row_list.append(dict1)
            else:
                pass


# make a dataframe of all the molecules, their descriptors and energies
df = pd.DataFrame(row_list, columns=[
                  'Name', 'Conformer', 'Optimization', 'Bond_Type', 'Bond_Length', 'Angle_Type', 'Angle',
                  'Torsion_Type', 'Torsion'])

# sort the dataframe to match with the csv
df.sort_values(by=['Conformer', 'Name'], ascending=[True, True], inplace=True)

# pull DFT energy out of csv
Energy = df_csv['dftE']
y = np.asarray(list(Energy), dtype=np.float)

# make an array of all the parts of the molecular descriptors
bondstuff = np.asarray(list(df['Bond_Type']))
bondlen = np.asarray(list(df['Bond_Length']))
anglestuff = np.asarray(list(df['Angle_Type']))
anglelen = np.asarray(list(df['Angle']))
torsionstuff = np.asarray(list(df['Torsion_Type']))
torsionlen = np.asarray(list(df['Torsion']))
# print(bondstuff.shape) # if array size errors print to see shape

# extract pm7 information from the csv and store in array
pm7 = []
for x in range(0, 23911):
    name = df_csv['Name'].iloc[x]
    conformer = df_csv['Conformer'].iloc[x]
    pm7E = df_csv['pm7E'].iloc[x]
    pm7E = float(pm7E)
    pm7_value = np.asarray([pm7E])
    fix_len = ([-9999] * 499)
    fix_len = np.asarray(fix_len)
    pm7_value = np.append(pm7_value, fix_len, axis=0)
    dict7 = {}
    dict7.update({'Name': name})
    dict7.update({'Conformer': conformer})
    dict7.update({'pm7E': pm7_value})
    pm7.append(dict7)

df_pm7E = pd.DataFrame(pm7, columns=['Name', 'Conformer', 'pm7E'])
pm7_info = np.asarray(list(df_pm7E['pm7E']))
# print(pm7_info.shape) # if array size errors print to see shape


# make the molecular descriptor that will be used as the input
molecule = np.concatenate((bondstuff, bondlen, anglestuff,
                           anglelen, torsionstuff, torsionlen, pm7_info), axis=1)
# print(molecule.shape) # if array size errors print to see shape

# test, train, split for data
X_train, X_test, y_train, y_test = train_test_split(molecule, y, test_size=0.2)
# print(X_train.shape) # if array size errors print to see shape
# print(y_train.shape) # if array size errors print to see shape

# specify regressor and run
clf = BayesianRidge()  # change to KernelRidge or SVR if need to
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
plt.savefig(figure_path + figure_name)
plt.show()
plt.close
