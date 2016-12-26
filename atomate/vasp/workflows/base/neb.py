#!/usr/bin/env python
# coding: utf-8

"""
Implements a NEBWorkflowManager for creating NEB workflow using neb_tasks.py.
It is recommended to run some primary tests to determine a reasonable
convergence criterion before carrying out high-throughput simulation.
"""

import logging
import yaml
import os

from monty.serialization import loadfn
from pymatgen.core import Structure
from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet, MVLCINEBSet
from custodian.vasp.jobs import VaspJob
from fireworks.core.firework import Firework, Workflow

from pymacy.fireworks.custodian_test import VaspNEBJob
from pymacy.fireworks.neb_firetasks import GetEndpointsFromIndexTask, \
    WriteEndpointInputTask, WriteNEBInputTask, RunVASPCustodianTask, \
    TransferNEBTask, InsertNEBDBTask, SetupNEBTask

try:
    from atomate.utils.utils import get_wf_from_spec_dict
except:
    get_wf_from_spec_dict = None

__author__ = "Hanmei Tang, Iek-Heng Chu"
__version__ = "1.0"
__status__ = "Development"
__date__ = "December 20, 2016"

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class PreNEBFW(Firework):
    """
    Initial relaxation workflow.
    Once finished, generating Endpoint workflow.
    TODO: NEB_ID, PATH_ID , MECHANISM, USE_DB_OR_NOT
    """

    def __init__(self, structure, name="Initial NEB setup"):
        ep_0 = GetEndpointsFromIndexTask(parent_structure=structure)
        setup = SetupNEBTask()


class EndPointRelaxationFW(Firework):
    """
    Endpoint workflow.
    Once finished, generating NEB workflow.
    """

    def __init__(self, structure, ep_index,
                 name="Endpoint relaxations",
                 vasp_input_set=None,
                 vasp_cmd="vasp", gamma_vasp_cmd="vasp.gamma",
                 parents=parents,
                 override_default_vasp_params=None,
                 db_file=None, **kwargs):
        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or MVLCINEBEndPointSet(structure,
                                                               **override_default_vasp_params)

        defaults = {"ISIF": 2}

        if user_incar_settings:
            defaults.update(user_incar_settings)

        job_rlx = [VaspJob(vasp_cmd=vasp_cmd,
                           gamma_vasp_cmd=gamma_vasp_cmd,
                           auto_npar=False)]

        vasp_rlx = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_rlx],
                                        handlers=[],
                                        custodian_params={})

        ep_write = WriteEndpointInputTask(ep_index=ep_index,
                                          user_incar_settings=defaults)
        ep_transfer = TransferRelaxationTask(type="endpoints", ep_index=ep_index)

        super(EndPointRelaxationFW, self).__init__([ep_write, vasp_rlx, ep_transfer],
                                                   parents=parents, name="{}-{}".
                                                   format(structure.composition.reduced_formula, name),
                                                   **kwargs)


class NEBFW(Firework):
    """
    NEB Workflow.
    """

    def __init__(self, structure, name="CINEB calculation"):
        job_neb = [VaspNEBJob(spec["neb_vasp_cmd"],
                              gamma_vasp_cmd=spec["neb_vasp_cmd_gam"],
                              auto_npar=False)]

        vasp_neb = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_neb],
                                        is_neb=True,
                                        handlers=[],
                                        custodian_params={})
