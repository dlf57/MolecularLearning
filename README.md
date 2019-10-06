# MolecularLearning
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![Build Status](https://travis-ci.org/dlf57/MolecularLearning.svg?branch=master)](https://travis-ci.org/dlf57/MolecularLearning)
[![codecov](https://codecov.io/gh/dlf57/MolecularLearning/branch/master/graph/badge.svg)](https://codecov.io/gh/dlf57/MolecularLearning)

This library was created as a graduate student project.

This is a library used to create molecular representations for machine learning applications. There are two types of representations with an organized set similar to standard bag of features representations and an unorganized set.

Organized Represntations:
 - BTY &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds w/ Atom Typing)
 - BTYNB &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds and Nonbonding w/ Atom Typing)
 - BATY &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds and Angles w/ Atom Typing)
 - BATYNB &nbsp;&nbsp;&nbsp;(Bonds, Angles, and Nonbonding w/ Atom Typing)
 - BATTY &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds, Angles, and Torsions w/ Atom Typing)
 - BATTYNB &nbsp;(Bonds, Angles, Torsions, and Nonbonding w/ Atom Typing)


Unorganized Representations
 - BTy &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds w/ Atom Typing)
 - BATy &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds and Angles w/ Atom Typing)
 - BATTy &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Bonds, Angles, and Torsions w/ Atom Typing)
 - BATTyNB &nbsp;(Bonds, Angles, Torsions, and Nonbonding w/ Atom Typing)

## Install molreps 
The latest version can be installed by cloning the repository and running:  
```pip install -e .```

### Dependencies
 - openbabel 

## Contact Information
 - Dakota Folmsbee: dfolmsbee@gmail.com or dlf57@pitt.edu
