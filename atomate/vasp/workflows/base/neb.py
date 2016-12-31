#!/usr/bin/env python
# coding: utf-8

"""
This module defines the Climbing Image Nudged Elastic Band (CI-NEB) workflow.
1) initial relaxation fireworks (pre-defined)
2) Generate endpoints --> Two endpoints relaxation
3) Use two endpoints --> Generate images --> CI-NEB
4)                       Images --> CI-NEB
"""

import logging
import yaml
import os

from monty.serialization import loadfn
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from custodian.vasp.jobs import VaspJob, VaspNEBJob
from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet, MVLCINEBSet
from fireworks.core.firework import Firework, Workflow
from atomate.utils.utils import get_logger
from atomate.utils.utils import get_wf_from_spec_dict
from atomate.vasp.firetasks.write_inputs import WriteNEBFromImages
from atomate.vasp.firetasks.neb_tasks import WriteEndpointInputTask, \
    RunVASPCustodianTask, \
    TransferNEBTask, SetupNEBTask  # InsertNEBDBTask, GetEndpointsFromIndexTask, WriteNEBInputTask

from atomate.vasp.firetasks.glue_tasks import PassStressStrainData
from atomate.vasp.fireworks.core import OptimizeFW, TransmuterFW
from atomate.vasp.fireworks.core_add import *


__author__ = "Hanmei Tang, Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'

logger = get_logger(__name__, level=logging.DEBUG)

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

spec = {"_category": "",  # <-- config
        "neb_id": 0,  # update runtime
        "perfect_dir": "",  # update runtime
        "endpoints_dirs": ["", ""],  # update runtime
        "neb_dirs": {},  # key = 1, 2 ... are path_id
        "pathway_sites": [],  # set outside by user
        "perfect_cell": {},  # from input
        "endpoints_0": {},  # from runtime or input
        "endpoints_1": {},  # from runtime or input
        "images": [],  # list of Structure (ini 01 02...end)
        "ini_user_incar_settings": {},  # config
        "ep_user_incar_settings": {},  # <-- config
        "neb_user_incar_settings_init": {},  # <-- config
        "neb_user_incar_settings_accu": {},  # <-- config
        "neb_vasp_cmd": [],  # <-- config, updated by user
        "neb_vasp_cmd_gam": [],  # <-- config
        "_queueadapter": {"nnodes": 1}  # <-- user
        }

def get_wf_neb(structures, name="neb workflow", db_file=None):
    """Return neb workflow according to inputs."""
    if len(structures) == 1:
        workflow = get_wf_neb_from_structure(structures[0])
    elif len(structures) == 2:
        workflow = get_wf_neb_from_endpoints(structures)
    else:
        workflow = get_wf_neb_from_images(structures)
    return workflow


def get_wf_neb_from_structure(structure, name="initial", db_file=None):
    """
    Returns a CI-NEB workflow from a given perfect structure.
    Args:
        structure:
        name:
        db_file:

    Returns:

    """
    workflow = None
    return workflow


def get_wf_neb_from_endpoints(structures, name="endpoints", db_file=None):
    """
    Returns a CI-NEB workflow from two given endpoints.

    Args:
        structures:
        name:
        db_file:

    Returns:

    """
    workflow = None
    return workflow


def get_wf_neb_from_images(structures, name="neb", db_file=None):
    """
    Returns a CI-NEB workflow from given images.

    Args:
        structures ([structure_0, structure_1, ...]): The image structures.
        name (str): some appropriate name for the workflow.
        db_file (str): path to file containing the database credentials.

    Returns:
        Workflow
    """
    logger.info("Get workflow from images.")

    fw1 = NEBFW(structures=structures, name=name)
    fw2 = NEBFW(structures=structures, name=name)
    workflow = Workflow([fw1, fw2], {fw1: [fw2]}, name="NEB Workflow")

    return workflow


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    structures = [PymatgenTest.get_structure("Si")]
    wf = get_wf_neb(structures)