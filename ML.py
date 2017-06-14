import sqlite3
import numpy as np
import matplotlib.pyplot as plt
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn.feature_extraction import DictVectorizer

conn = sqlite3.connect('MolecularData.db')
cursor = conn.cursor()

# Need to read through the database
cursor.execute('SELECT Name, Bonds, Angles, Torsions, Energy FROM MoleculeData')
for row in cursor:
    name = row[0]
    bonds = row[1]
    angles = row[2]
    torsions = row[3]
    energy = row[4]
    
    entries = bonds.split(";")
    for entry in entries:
        type, length = entry.split(",")
        print float(length)

    # vec = DictVectorizer()
    # vec.fit_transform(bond).toarray()
    # print vec

    # X = [molecule]
    # X = preprocessing.scale(X)
    # y = [energy]
    #
    # X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    #
    # reg =KernelRidge()
    # reg.fit(X_train, y_train)
    #
    # accuracy = reg.score(X_train, y_train)
    # print accuracy
    #
    # predict = reg.predict(X_test, y_test)
    # print(predict)

cursor.close()
conn.close()


# # Test inputs ideally from other conformers
# #X_test = Will uncomment if we choose to test against conformers
#
# # Might not do this if we use the confomers as a test
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
#
# reg = KernelRidge()
# # Training of the model
# reg.fit(X_train, y_train)
#
# # Check to see accuracy of the training
# accuracy = reg.score(X_train, y_train)
# print(accuracy)
#
# # Predicting the output of the conformers
# #predict = clf.predict(X_test, Y_test)
# #print(predict)
