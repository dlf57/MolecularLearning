'''
Example script showing off how to build molecular representations 
to be used with scikit-learn for machine learning.

Author: Dakota Folmsbee <dfolmsbee@gmail.com>
'''

import glob
import numpy as np
import pybel
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from molreps.reps import Representation


# where to get molecule files
directory = '/path/to/molecules/'  # molecule directory
filetype = 'filetype'  # molecule filetype

# import csv with properties
props = pd.read_csv('/path/to/csv/properties.csv')
prop = np.asarray(list(props['Property']), dtype=np.float)

# build molecular descriptors
rep = Representation('BTY')
molecule = np.array([rep.rep(next(pybel.readfile(filetype, molfile)))
                     for molfile in sorted(glob.iglob(directory + '*.' + filetype))])

# test, train, split for data
X_train, X_test, y_train, y_test = train_test_split(
    molecule, prop, test_size=0.2)

# specify regressor and run
regr = RandomForestRegressor(n_jobs=-1)
regr.fit(X_train, y_train)
accuracy = regr.score(X_test, y_test)
r2 = ("%8.3f" % accuracy)
print('r^2 =', accuracy)
predicted = regr.predict(X_test)
mae = mean_absolute_error(y_test, predicted)
print('Mean Absolute Error =', mae)
