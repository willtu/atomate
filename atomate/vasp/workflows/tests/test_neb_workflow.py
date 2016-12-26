# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

import json
import os
import shutil
import unittest

import numpy as np

from pymongo import MongoClient

from fireworks import LaunchPad, FWorker
from fireworks.core.rocket_launcher import rapidfire

from atomate.vasp.powerups import use_fake_vasp
from atomate.vasp.workflows.presets.core import wf_elastic_constant

from pymatgen import SETTINGS
from pymatgen.util.testing import PymatgenTest
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


__author__ = "Hanmei Tang, Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'
__version__ = "1.0"
__status__ = "Development"
__date__ = "December 25, 2016"


class TestNudgedElasticBandWorkflow(unittest.TestCase):
    def test_wf(self):
        pass

if __name__ == "__main__":
    unittest.main()