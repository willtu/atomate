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

from pymatgen.core import Structure
from custodian.vasp.jobs import VaspJob, VaspNEBJob
from fireworks.core.firework import Firework, Workflow
from monty.serialization import loadfn

from pymacy.fireworks.neb_tasks import GetEndpointsFromIndexTask, \
    WritePerfectStructureInputTask, WriteEndpointInputTask, WriteNEBInputTask,\
    RunVASPCustodianTask, \
    TransferRelaxationTask, TransferNEBTask, InsertNEBDBTask

from atomate.utils.utils import get_logger
from atomate.vasp.firetasks.neb_tasks import GetEndpointsFromIndexTask, \
    WriteEndpointInputTask, WriteNEBInputTask, RunVASPCustodianTask, \
    TransferNEBTask, SetupNEBTask  # InsertNEBDBTask

try:
    from atomate.utils.utils import get_wf_from_spec_dict
except:
    get_wf_from_spec_dict = None

__author__ = "Hanmei Tang and Iek-Heng Chu"
__email__ = 'hat003@eng.ucsd.edu, ihchu@eng.ucsd.edu'

logger = get_logger(__name__, level=logging.DEBUG)

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))



class NEBWorkflowManager(object):

    # This is a complete spec template.
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

    def __init__(self, configfile):
        """
        Initialization of NEBWorkflowManager.

        Args:
            configfile (str): neb workflow config yaml file.
        """

        logger.info("NEBWorkflowManager initialization.")

        with open(configfile, "r") as yaml_file:
            logger.info("Loading config yaml file...")
            config = yaml.load(yaml_file)

            try:
                images_init = config["neb_user_incar_settings_init"]["IMAGES"]
                images_accu = config["neb_user_incar_settings_accu"]["IMAGES"]
                if images_accu != images_init:
                    raise ValueError("Conflict setting in NEB INCAR")
            except:
                pass

            spec = config
            self.spec.update(spec)

            logger.info("Initialization done.")

    def _get_spec(self, pathway_sites=None, endpoints=None, images=None):
        """
        Read in user's settings and update spec.

        Returns:
            spec (dict): standard spec.

        Args:
            pathway_sites ([int, int]): Indicating pathway site indexes.
            endpoints ([ep0_dict, ep1_dict]): The two endpoints structures.
            images ([structure_0, structure_1, ...]): The image structures.
        """
        spec = self.spec

        if pathway_sites is not None:
            spec["pathway_sites"] = pathway_sites

        if endpoints is not None:
            spec["endpoints_0"] = endpoints[0].as_dict()
            spec["endpoints_1"] = endpoints[1].as_dict()

        if images is not None:
            if len(images) <= 2:
                raise ValueError("Too few images!")
            spec["images"] = images
            n_images = len(images) - 2
            spec["_queueadapter"].update({"nnodes": n_images})

            neb_vasp_cmd = spec["neb_vasp_cmd"]
            index = neb_vasp_cmd.index('-np')
            np = n_images * int(neb_vasp_cmd[index + 1])
            neb_vasp_cmd[index + 1] = str(np)
            spec["neb_vasp_cmd"] = neb_vasp_cmd

        self.spec = spec

        return spec

    def get_wf_from_perfect_structure(self, structure,
                                      pathway_sites,
                                      neb_name=None,
                                      is_relaxed=False,
                                      insert_db=False):
        """
        Get workflow from a perfect structure and two pathway site index.
        This designs for regular vacancy diffusion cases.

        Args:
            structure (Structure): perfect structure.
            pathway_sites ([int, int]): pathway site indexes.
            neb_name (str): The name of NEB run. Set up this to override
                            the default name, which is the reduced formula
                            of perfect structure.

            is_relaxed (bool): A tag indicating the perfect structure
                               is relaxed or not.
            insert_db (bool): Tag enabling insertion to neb collection.
        """

        logger.info("Getting workflow from a perfect structure...")

        spec = self._get_spec(pathway_sites=pathway_sites)
        job_name = structure.composition.reduced_formula
        if neb_name:
            job_name = neb_name

        job_rlx = [VaspJob(vasp_cmd=["mpirun", "-np", "24", "vasp"],
                           gamma_vasp_cmd=["vasp.gam"],
                           auto_npar=False)]
        job_neb = [VaspNEBJob(spec["neb_vasp_cmd"],
                              gamma_vasp_cmd=spec["neb_vasp_cmd_gam"],
                              auto_npar=False)]
        vasp_rlx = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_rlx],
                                        is_neb=False,
                                        handlers=[],
                                        custodian_params={})

        vasp_neb = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_neb],
                                        is_neb=True,
                                        handlers=[],
                                        custodian_params={})

        ini_0 = WritePerfectStructureInputTask(structure=structure)
        ini_1 = TransferRelaxationTask(jobname=job_name,
                                       cal_label="perfect")
        ep_0 = GetEndpointsFromIndexTask(structure=structure)
        ep_10 = WriteEndpointInputTask(endpoint_index=0)
        ep_11 = WriteEndpointInputTask(endpoint_index=1)
        ep_20 = TransferRelaxationTask(jobname=job_name, cal_label="0")
        ep_21 = TransferRelaxationTask(jobname=job_name, cal_label="1")
        neb_00 = WriteNEBInputTask(mode="start", step="initial")
        neb_01 = WriteNEBInputTask(mode="continue", step="accurate")
        neb_1 = TransferNEBTask(jobname=job_name)
        db_ins = InsertNEBDBTask()

        fw1 = Firework([ini_0, vasp_rlx, ini_1],
                       name="{} Geometry".format(job_name),
                       spec=spec, fw_id=1)
        fw2 = Firework([ep_0, ep_10, vasp_rlx, ep_20],
                       name="{} Ep_0".format(job_name),
                       spec=spec, fw_id=2)
        fw3 = Firework([ep_0, ep_11, vasp_rlx, ep_21],
                       name="{} Ep_1".format(job_name),
                       spec=spec, fw_id=3)
        fw4 = Firework([neb_00, vasp_neb, neb_1],
                       name="Initial CI-NEB",
                       spec=spec, fw_id=4)
        fw5 = Firework([neb_01, vasp_neb, neb_1],
                       name="Accurate CI-NEB",
                       spec=spec, fw_id=5)
        if insert_db:
            fw5 = Firework([neb_01, vasp_neb, neb_1, db_ins],
                           name="Accurate CI-NEB",
                           spec=spec, fw_id=5)

        workflow = Workflow([fw2, fw3, fw4, fw5],
                            {fw2: [fw4], fw3: [fw4], fw4: [fw5]},
                            name="1-Ep-NEB Workflow")

        if not is_relaxed:
            workflow = Workflow([fw1, fw2, fw3, fw4, fw5],
                                {fw1: [fw2, fw3],
                                 fw2: [fw4],
                                 fw3: [fw4],
                                 fw4: [fw5]},
                                name="1-Rlx-Ep-NEB Workflow")

        return workflow

    def get_wf_from_endpoints(self, neb_name, endpoint_0, endpoint_1,
                              is_relaxed=False):
        """
        Get workflow from two given endpoints.
        This designs for cases:
            1) having more than one mobile species;
            2) interstitial diffusion mechanism.

        Args:
            neb_name (str): The name of NEB run, which is
                            the reduced formula of perfect structure by default.
            endpoint_0 (Structure): Endpoint 0
            endpoint_1 (Structure): Endpoint 1
            is_relaxed (bool): A tag indicating the endpoint structures
                               are relaxed or not.
        """
        logger.info("Get workflow from two endpoints...")
        spec = self._get_spec(endpoints=[endpoint_0, endpoint_1])

        job_rlx = [VaspJob(vasp_cmd=["mpirun", "-np", "24", "vasp"],
                           gamma_vasp_cmd=["vasp.gam"],
                           auto_npar=False)]
        job_neb = [VaspNEBJob(spec["neb_vasp_cmd"],
                              gamma_vasp_cmd=spec["neb_vasp_cmd_gam"],
                              auto_npar=False)]
        vasp_rlx = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_rlx],
                                        is_neb=False,
                                        handlers=[],
                                        custodian_params={})

        vasp_neb = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_neb],
                                        is_neb=True,
                                        handlers=[],
                                        custodian_params={})

        ep_10 = WriteEndpointInputTask(endpoint_index=0)
        ep_11 = WriteEndpointInputTask(endpoint_index=1)
        ep_20 = TransferRelaxationTask(jobname=neb_name, cal_label="0")
        ep_21 = TransferRelaxationTask(jobname=neb_name, cal_label="1")
        neb_00 = WriteNEBInputTask(mode="start", step="initial")
        neb_01 = WriteNEBInputTask(mode="continue", step="accurate")
        neb_1 = TransferNEBTask(jobname=neb_name)

        fw1 = Firework([ep_10, vasp_rlx, ep_20],
                       name="{} Ep_0".format(neb_name),
                       spec=spec, fw_id=1)
        fw2 = Firework([ep_11, vasp_rlx, ep_21],
                       name="{} Ep_1".format(neb_name),
                       spec=spec, fw_id=2)
        fw3 = Firework([neb_00, vasp_neb, neb_1],
                       name="Initial CI-NEB",
                       spec=spec, fw_id=3)
        fw4 = Firework([neb_01, vasp_neb, neb_1],
                       name="Accurate CI-NEB",
                       spec=spec, fw_id=4)

        workflow = Workflow([fw3, fw4], {fw3: [fw4]}, name="NEB Workflow")

        if not is_relaxed:
            workflow = Workflow([fw1, fw2, fw3, fw4],
                                {fw1: [fw3], fw2: [fw3], fw3: [fw4]},
                                name="Ep-NEB Workflow")

        return workflow

    def get_wf_from_images(self, neb_name, structures):
        """
        Get workflow from given images.

        Args:
            neb_name (str): The name of NEB run, which is
                            the reduced formula of perfect structure by default.
            structures ([structure_0, structure_1, ...]): The image structures.
        """

        logger.info("Get workflow from images...")

        spec = self._get_spec(images=structures)
        job_neb = [VaspNEBJob(spec["neb_vasp_cmd"],
                              gamma_vasp_cmd=spec["neb_vasp_cmd_gam"],
                              auto_npar=False)]

        vasp_neb = RunVASPCustodianTask(jobs=[j.as_dict() for j in job_neb],
                                        is_neb=True,
                                        handlers=[],
                                        custodian_params={})

        neb_0 = WriteNEBInputTask(mode="continue", step="initial")
        neb_1 = WriteNEBInputTask(mode="continue", step="accurate")
        neb_2 = TransferNEBTask(jobname=neb_name)

        fw1 = Firework([neb_0, vasp_neb, neb_2],
                       name="Initial CI-NEB",
                       spec=spec, fw_id=1)
        fw2 = Firework([neb_1, vasp_neb, neb_2],
                       name="Accurate CI-NEB",
                       spec=spec, fw_id=2)

        workflow = Workflow([fw1, fw2], {fw1: [fw2]}, name="NEB Workflow")
        return workflow


# def get_wf(structure, wf_filename, params=None, common_params=None, vis=None):
#     """
#     (***Copied from MatMethods for testing purposes.)
#     A generic function to load generic VASP library workflow, while
#     overriding some of the parameters via the function arguments.
#
#     Args:
#         structure: (Structure) structure to run
#         wf_filename: filename in library subdir, e.g. "neb.yaml"
#         params: (list of dicts) set params for each Firework; format is list
#             that is same length as # of fws in the workflow
#         common_params: (dict) set common params
#         vis: (VaspInputSet) A VaspInputSet to use for the first FW
#     Returns:
#         A Workflow
#     """
#     d = loadfn(os.path.join(module_dir, "library", wf_filename))
#
#     if get_wf_from_spec_dict is None:
#         raise NotImplemented("MatMethods is a required module!")
#
#     if params:
#         if len(params) != len(d["fireworks"]):
#             raise ValueError("The length of the params array must match the"
#                              "length of the Fireworks array!")
#         for idx, v in enumerate(params):
#             if "params" not in d["fireworks"][idx]:
#                 d["fireworks"][idx]["params"] = {}
#             d["fireworks"][idx]["params"].update(v)
#
#     if common_params:
#         if 'common_params' not in d:
#             d["common_params"] = {}
#         d["common_params"].update(common_params)
#
#     if vis:
#         if "params" not in d["fireworks"][0]:
#             d["fireworks"][0]["params"] = {}
#         d["fireworks"][0]["params"]["vasp_input_set"] = vis.as_dict()
#
#     return get_wf_from_spec_dict(structure, d)


def wf_nudged_elastic_band(structure, site_indices, name="CI-NEB Workflow",
                       vasp_input_set="MITRelaxSet", vasp_cmd="vasp", neb_id=None,
                       db_config=None, params=None, debug=False):

    if len(site_indices) != 2 or len(set(site_indices)) != 2:
        raise ValueError("Two distinct site indices are required to "
                         "define a hop.")

    common_params = {"vasp_cmd": vasp_cmd, "site_indices": site_indices,
                     "debug": debug}

    if db_config is not None:
        common_params.update({"db_config": db_config})

    if neb_id is not None:
        common_params.update({"neb_id": neb_id})

    wf = get_wf(structure, "neb.yaml", vis=vasp_input_set,
                common_params=common_params, params=params)

    wf.name = "{}:{}".format(structure.composition.reduced_formula, name)

    return wf


if __name__ == "__main__":
    pass
