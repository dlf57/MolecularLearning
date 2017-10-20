# MolecularLearning
This is a repository for my graduate studies research.

Currently this repository consists of two working files. The file mol_sql.py can be used to create a database using sqlite format. The file mol_ml.py reads in the molecular information, sorts, and runs a machine learning algorithm in order to predict the desired properties. mol_ml.py uses [scikit-learn](http://scikit-learn.org/stable/) library for regressors and classifiers used in machine learning. 

### Future direction:
 - Test on a larger dataset (currently ~700 molecules and ~23,000 conformers have been tested)
 - Modify parameters to increase prediction accuracy 
 - Increase level of theory the ML is trained on

### Files in this Repository:
 - mol_ml.py (Reads through files and uses Bayesian Ridge Regression)
 - mol_sql.py (Reads through files and stores in sqlite database)

If you have any questions or comments please reach out to me at *dfolmsbee@gmail.com* or *dlf57@pitt.edu*.
