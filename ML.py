import pickle
import numpy as np
import modified_molml
#import AtomizationE
import matplotlib.pyplot as plt
import itertools
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.feature_extraction import DictVectorizer

pickle_in = open('Molecule_Dataset.pickle', 'rb')
molecules = []
while 1:
    try:
        molecules.append(pickle.load(pickle_in))
    except EOFError:
        break
#print molecules

#pickle_in = open('Molecule_Data.pickle', 'rb')
#molecule_data = pickle.load(pickle_in)
#print molecule_data

#A = modified_molml.molecules
A = molecules
#B = AtomizationE.mol.energy

print A

# Training input values from given data set
#X = A

# Training target values atomization energy
#y = B

# Test inputs ideally from other conformers
#X_test = Will uncomment if we choose to test against conformers

# Might not do this if we use the confomers as a test
#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

#reg = KernelRidge()
# Training of the model
#reg.fit(X_train, Y_train)

# Check to see accuracy of the training
#accuracy = reg.score(X_train, Y_train)
#print(accuracy)

# Predicting the output of the conformers
#predict = clf.predict(X_test, Y_test)
#print(predict)


