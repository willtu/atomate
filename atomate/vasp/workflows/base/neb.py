#!/usr/bin/env python
# coding: utf-8

"""
This module defines the Climbing Image Nudged Elastic Band (CI-NEB) workflow.
1) initial relaxation fireworks (pre-defined)
2) Generate endpoints --> Two endpoints relaxation
3) Use two endpoints --> Generate images --> CI-NEB
4)                       Images --> CI-NEB
"""

import yaml
import os
from datetime import datetime

from monty.serialization import loadfn
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet
from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet, MVLCINEBSet, get_endpoints_from_index
from fireworks.core.firework import Firework, Workflow
from atomate.utils.utils import get_logger
from atomate.utils.utils import get_wf_from_spec_dict

from atomate.vasp.fireworks.core import OptimizeFW, NEBFW, NEBRelaxationFW

__author__ = "Hanmei Tang, Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'

logger = get_logger(__name__)

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

spec_orig = {"neb_id": 0,
             "path_id": 0,
             "vasp_cmd": ">>vasp_cmd<<",
             "gamma_vasp_cmd": ">>gamma_vasp_cmd<<",
             "neb_vasp_cmd": [],
             "neb_gamma_vasp_cmd": ">>gamma_vasp_cmd<<",
             "db_file": ">>db_file<<",  # TODO: May remove this
             "_category": "",
             "_queueadapter": {"nnodes": 1},
             "calc_locs": "",
             "rlx_dir": "",
             "ep0_dir": "",  # TODO: This is enveloped in CalLocs
             "ep1_dir": "",
             "neb_dir": {},  # key = "1", "2" ... are path_id
             "rlx_st": {},
             "ep0_st": {},
             "ep1_st": {},
             "path_sites": [],
             "images": []}

# TODO: Double check this is useful or not
# def _update_spec_from_config(spec, configfile):
#     """
#     Update spec from config file.
#     Args:
#         spec (dict): input spec.
#         configfile (str): neb workflow config yaml file.
#
#     Returns:
#         spec (dict)
#     """
#     with open(configfile, "r") as yaml_file:
#         config = yaml.load(yaml_file)
#         if "common_params" in config:
#             spec.update(config["common_params"])
#         return spec


def _update_spec_from_inputs(spec, path_sites=None,
                             endpoints=None, images=None):
    """
    Update spec according to inputs.

    Args:
        path_sites ([int, int]): Indicating pathway site indexes.
        endpoints ([ep0_dict, ep1_dict]): The two endpoints structures.
        images ([s0_dict, s1_dict, ...]): The image structures,
            including the two endpoints.
    """
    if path_sites is not None:
        spec["path_sites"] = path_sites

    if endpoints is not None:
        spec["ep0_st"] = endpoints[0].as_dict()
        spec["ep1_st"] = endpoints[1].as_dict()

    if images is not None:
        if len(images) <= 2:
            raise ValueError("Too few images!")
        spec["images"] = images
        n_images = len(images) - 2
        spec["_queueadapter"].update({"nnodes": n_images})

        # Set -np tag for vasp command
        # TODO: What if neb_vasp_cmd is empty?
        neb_vasp_cmd = spec["neb_vasp_cmd"]

        index = neb_vasp_cmd.index('-np')
        np = n_images * int(neb_vasp_cmd[index + 1])
        neb_vasp_cmd[index + 1] = str(np)
        spec["neb_vasp_cmd"] = neb_vasp_cmd
        return spec


# def get_wf_neb(structures):
#     # This should be a very long input, don't use configfile
#     # fw1: write inputs
#     # fw2: run neb
#     # fw3: transfer result to destination
#     """
#     Return neb workflow according to inputs.
#
#     Args:
#         structures:
#         configfile:
#         transfer_dir (str): list of transfer directories.
#
#     Returns:
#
#     """
#     # spec = spec_orig
#     # if configfile is not None:
#     #     spec = _update_spec_from_config(spec, configfile)
#
#     workflow = get_wf_neb_from_images()
#
#     return workflow


def get_wf_neb_from_structure(structure, path_sites,
                              is_optimized=True,
                              wfname=None, neb_round=1,
                              spec=None, uis_ini=None,
                              uis_ep=None, uis_neb=None):
    """
    Get a CI-NEB workflow from given endpoints.
    Workflow: (Init relax) -- Endpoints relax -- NEB_1 -- NEB_2 - ... - NEB_n

    Args:
        structure (Structure): The perfect structure.
        path_sites (list[int]): The two vacancy site indices.
        is_optimized (bool): True implies the provided structure is optimized,
            otherwise run initial relaxation before generating endpoints.
        wfname (str): some appropriate name for the workflow.
        neb_round (int): Times of NEB calculations.
        spec (dict): user setting spec settings to overwrite spec_orig.
        uis_ini (dict): Additional user_incar_settings for initial relaxations.
        uis_ep (dict): Additional user_incar_settings for endpoint relaxations.
        uis_neb (dict): Additional user_incar_settings for NEB.

    Returns:
        Workflow
    """
    logger.info("Get workflow from perfect structure.")

    # Read in params.
    formula = structure.composition.reduced_formula
    wfname = wfname or datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    spec = spec or {}

    uis_ini = uis_ini or {}
    uis_ep = uis_ep or {}
    uis_neb = uis_neb or {}

    # Get neb fireworks.
    neb_fws = []
    for n in range(neb_round):
        fw = NEBFW(spec=spec,
                   name=formula,
                   neb_label=str(n),
                   from_images=False,
                   vasp_input_set=MVLCINEBSet,
                   user_incar_settings=uis_neb)
        neb_fws.append(fw)

    # Get relaxation fireworks.
    if is_optimized:
        endpoints = get_endpoints_from_index(structure, path_sites)
        ep_dict = [e.as_dict() for e in endpoints]
        spec = _update_spec_from_inputs(spec, endpoints=ep_dict)
        rlx_fws = [NEBRelaxationFW(spec=spec,
                                   st_label="ep{}".format(i),
                                   name=formula,
                                   user_incar_settings=uis_ep) for i in [0, 1]]

        links = {rlx_fws[0]: [neb_fws[0]],
                 rlx_fws[1]: [neb_fws[0]]}
    else:
        spec = _update_spec_from_inputs(spec, path_sites=path_sites)
        rlx_fws = [NEBRelaxationFW(spec=spec,
                                   st_label="{}".format(label),
                                   name=formula,
                                   user_incar_settings=incar)
                   for label, incar in zip(["rlx", "ep0", "ep1"],
                                           [uis_ini, uis_ep, uis_ep])]

        links = {rlx_fws[0]: [rlx_fws[1], rlx_fws[1]],
                 rlx_fws[1]: [neb_fws[0]],
                 rlx_fws[2]: [neb_fws[0]]}

    # Append Firework links.
    fws = rlx_fws + neb_fws
    for r in range(neb_round):
        links[neb_fws[r]] = [neb_fws[r + 1]]

    workflow = Workflow(fws, links_dict=links,
                        name="neb_{}".format(wfname))
    return workflow


def get_wf_neb_from_endpoints(endpoints=None, is_optimized=True,
                              wfname=None, neb_round=1,
                              spec=None, uis_ep=None, uis_neb=None):
    """
    Get a CI-NEB workflow from given endpoints.
    Workflow: (Endpoints relax -- ) NEB_1 -- NEB_2 - ... - NEB_n

    Args:
        endpoints (list[Structure]): The image structures,
            if None then read from spec.
        is_optimized (bool): True implies the provided endpoint structures
            are optimized, otherwise run endpoints relaxation before NEB.
        wfname (str): some appropriate name for the workflow.
        neb_round (int): Times of NEB calculations.
        spec (dict): user setting spec settings to overwrite spec_orig.
        uis_ep (dict): Additional user_incar_settings for endpoint relaxations.
        uis_neb (dict): Additional user_incar_settings for NEB.

    Returns:
        Workflow
    """
    logger.info("Get workflow from endpoints.")

    # Read in params.
    if endpoints is not None:
        endpoints_dict = [s.as_dict() for s in endpoints]
    else:
        endpoints_dict = [spec["ep0_st"], spec["ep1_st"]]
    formula = endpoints[0].composition.reduced_formula
    wfname = wfname or datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    spec = spec or {}
    spec = _update_spec_from_inputs(spec, endpoints=endpoints_dict)
    uis_ep = uis_ep or {}
    uis_neb = uis_neb or {}

    neb_fws = []
    for n in range(neb_round):
        fw = NEBFW(spec=spec,
                   name=formula,
                   neb_label=str(n),
                   from_images=False,
                   vasp_input_set=MVLCINEBSet,
                   user_incar_settings=uis_neb)
        neb_fws.append(fw)
    workflow = Workflow(neb_fws, name="neb_workflow_{}".format(wfname))

    # Add endpoints relaxation if structures not optimized.
    if not is_optimized:
        ep_fws = [NEBRelaxationFW(spec=spec,
                                  st_label="ep{}".format(i),
                                  name=formula,
                                  user_incar_settings=uis_ep)
                  for i in [0, 1]]

        # Create Firework links.
        fws = ep_fws + neb_fws
        links = {ep_fws[0]: [neb_fws[0]],
                 ep_fws[1]: [neb_fws[0]]}
        for r in range(neb_round):
            links[neb_fws[r]] = [neb_fws[r + 1]]

        workflow = Workflow(fws, links_dict=links,
                            name="neb_{}".format(wfname))
        return workflow


def get_wf_neb_from_images(images=None, wfname=None, neb_round=1,
                           spec=None, uis_neb=None):
    """
    Get a CI-NEB workflow from given images.
    Workflow: NEB_1 -- NEB_2 - ... - NEB_n

    Args:
        images (list[Structure]): The image structures,
            if None then read from spec["images"]
        wfname (str): some appropriate name for the workflow.
        neb_round (int): Times of NEB calculations.
        spec (dict): user setting spec settings to overwrite spec_orig.
        uis_neb (dict): Additional user_incar_settings for NEB.

    Returns:
        Workflow
    """
    logger.info("Get workflow from images.")

    # Read in params.
    spec = spec or {}
    formula = images[0].composition.reduced_formula
    if images is not None:
        images_dict = [s.as_dict() for s in images]
    else:
        images_dict = spec["images"]
    wfname = wfname or datetime.utcnow().strftime('%Y-%m-%d-%H-%M-%S-%f')
    spec = _update_spec_from_inputs(spec, images=images_dict)
    uis_neb = uis_neb or {}

    fws = []
    for n in range(neb_round):
        fw = NEBFW(spec=spec,
                   name=formula,
                   neb_label=str(n),
                   from_images=True,
                   vasp_input_set=MVLCINEBSet,
                   user_incar_settings=uis_neb)
        fws.append(fw)

    workflow = Workflow(fws, name="neb_{}".format(wfname))

    return workflow


if __name__ == "__main__":
    from pymatgen.util.testing import PymatgenTest

    pass
    # structures = [PymatgenTest.get_structure("Si")]
    # wf = get_wf_neb(structures)
