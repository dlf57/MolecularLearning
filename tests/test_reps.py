'''
Test script to test representation building
'''
# Author: Dakota Folmsbee <dfolmsbee@gmail.com>
# License: GPLv2

import pybel
from molreps.reps import Representation
import pytest

mol_file = './tests/methane.sdf'


def test_bty():
    # test bond representation on methane
    mol = next(pybel.readfile('sdf', mol_file))
    rep = Representation('BTY', n_bonds=5)
    bonds = [63010.0, 1.0919, 63010.0, 1.0919, 63010.0,
             1.0919, 63010.0, 1.0919, -99999, -99999]
    assert rep.rep(mol) == bonds


def test_baty():
    # test bond/angle representation on methane
    mol = next(pybel.readfile('sdf', mol_file))
    rep = Representation('BATY', n_bonds=5, n_angles=8)
    angles = [63010.0, 1.0919, 63010.0, 1.0919, 63010.0, 1.0919, 63010.0, 1.0919, -99999, -99999, 63010010.0, 109.471, 63010010.0,
              109.469, 63010010.0, 109.473, 63010010.0, 109.469, 63010010.0, 109.474, 63010010.0, 109.472, -99999, -99999, -99999, -99999]
    assert rep.rep(mol) == angles


def test_batty():
    # test bond/angle/torsion representation on methane
    mol = next(pybel.readfile('sdf', mol_file))
    rep = Representation('BATTY', n_bonds=5, n_angles=8, n_torsions=2)
    torsions = [63010.0, 1.0919, 63010.0, 1.0919, 63010.0, 1.0919, 63010.0, 1.0919, -99999, -99999, 63010010.0, 109.471, 63010010.0, 109.469,
                63010010.0, 109.473, 63010010.0, 109.469, 63010010.0, 109.474, 63010010.0, 109.472, -99999, -99999, -99999, -99999, -99999, -99999, -99999, -99999]
    assert rep.rep(mol) == torsions


def test_battynb():
    # test bond/angle/torsion/nonbonding representation on methane
    mol = next(pybel.readfile('sdf', mol_file))
    rep = Representation('BATTYNB', n_bonds=5,
                         n_angles=8, n_torsions=2, n_nb=1)
    nb = [63010.0, 1.0919, 63010.0, 1.0919, 63010.0, 1.0919, 63010.0, 1.0919, -99999, -99999, 63010010.0, 109.471, 63010010.0, 109.469, 63010010.0, 109.473,
          63010010.0, 109.469, 63010010.0, 109.474, 63010010.0, 109.472, -99999, -99999, -99999, -99999, -99999, -99999, -99999, -99999, -99999, -99999]
    assert rep.rep(mol) == nb
