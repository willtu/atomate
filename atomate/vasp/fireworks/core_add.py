# coding: utf-8

from __future__ import absolute_import, division, print_function, unicode_literals

"""
Defines standardized Fireworks that can be chained easily to perform various
sequences of VASP calculations.
This is add-on part to core.py
Under test, separated to avoid conflict.
"""

from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MITMDSet
from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet, MVLCINEBSet
from custodian.vasp.jobs import VaspJob, VaspNEBJob
from fireworks.core.firework import Firework, Workflow

from atomate.vasp.firetasks.write_inputs import *
from atomate.vasp.firetasks.neb_tasks import WriteEndpointInputTask, \
    RunVASPCustodianTask, TransferNEBTask, SetupNEBTask
# InsertNEBDBTask, GetEndpointsFromIndexTask, WriteNEBInputTask
from atomate.common.firetasks.glue_tasks import PassCalcLocs
from atomate.vasp.firetasks.glue_tasks import CopyVaspOutputs, PassEpsilonTask, PassNormalmodesTask
from atomate.vasp.firetasks.parse_outputs import VaspToDbTask, BoltztrapToDBTask
from atomate.vasp.firetasks.run_calc import RunVaspCustodian, RunBoltztrap


class InitialFW(Firework):
    # TODO: This is in fireworks/core/OptimizeFW
    pass


class EndPointRelaxationFW(Firework):
    """
    Endpoint workflow.
    """

    def __init__(self, structure, endpoint_label="0",
                 name="Endpoint relaxations",
                 vasp_input_set=None, vasp_cmd="vasp",
                 user_incar_settings=None,
                 override_default_vasp_params=None,
                 db_file=None, **kwargs):
        """

        Args:
            structure (Structure): The endpoint structure.
            endpoint_label (str): "0" or "1", denoting endpoint.
            name (str): some appropriate name for the workflow.
            vasp_input_set:
            vasp_cmd:
            override_default_vasp_params:
            db_file (str): path to file containing the database credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """

        override_default_vasp_params = override_default_vasp_params or {}
        vasp_input_set = vasp_input_set or \
                         MVLCINEBEndPointSet(structure,
                                             **override_default_vasp_params)

        defaults = {"ISIF": 2}

        if user_incar_settings:
            defaults.update(user_incar_settings)

        ep_write = WriteEndpointInputTask(ep_index=endpoint_label,
                                          user_incar_settings=defaults)

        job = [VaspJob(vasp_cmd=vasp_cmd,
                       gamma_vasp_cmd=">>gamma_vasp_cmd<<",
                       auto_npar=False)]

        run_rlx_task = RunVASPCustodianTask(jobs=[j.as_dict() for j in job],
                                            handlers=[],
                                            custodian_params={})

        ep_transfer = TransferNEBTask(type="endpoints",
                                      ep_index=endpoint_label)
        tasks = [ep_write, run_rlx_task, ep_transfer]

        super(EndPointRelaxationFW, self).__init__(tasks,
                                                   name="ep-{}".format(name),
                                                   **kwargs)


class NEBFW(Firework):
    """
    NEB Workflow.
    """

    def __init__(self, structures, name="neb",
                 vasp_input_set=None, vasp_cmd="vasp",
                 user_incar_settings=None,
                 dest_dir="dest_dir_neb",
                 db_file=None, **kwargs):
        """

        Args:
            structures ([structure_0, structure_1, ...]): The image structures.
            name (str): some appropriate name for the workflow.

            vasp_input_set:
            vasp_cmd:
            override_default_vasp_params:
            db_file (str): path to file containing the database credentials.
            \*\*kwargs: Other kwargs that are passed to Firework.__init__.
        """
        logger.info("NEBFW")

        user_incar_settings = user_incar_settings or {}
        vi = vasp_input_set or MVLCINEBSet(structures)

        write_neb = WriteNEBFromImages(vasp_input_set=vi,
                                       image_files=structures,
                                       **user_incar_settings)
        job = [VaspNEBJob(vasp_cmd,
                          final=False,
                          auto_npar=False,
                          half_kpts=True,
                          auto_gamma=True,
                          gamma_vasp_cmd=">>gamma_vasp_cmd<<")]

        run_neb_task = RunVASPCustodianTask(jobs=[j.as_dict() for j in job],
                                            handlers=[],
                                            custodian_params={})

        transfer_neb = TransferNEBTask(dest_dir=dest_dir)
        tasks = [write_neb, run_neb_task, transfer_neb]

        super(NEBFW, self).__init__(tasks, name="neb-{}".format(name), **kwargs)