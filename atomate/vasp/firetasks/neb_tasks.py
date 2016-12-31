#!/usr/bin/env python
# coding: utf-8

"""
Collection of tasks specific for automating NEB simulations.
Keep separated to avoid conflict.
"""

import os
import logging
import glob
import shutil

from pymatgen.core import Structure
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from custodian import Custodian
from monty.json import MontyDecoder
from pymatgen_diffusion.neb.io import MVLCINEBEndPointSet, MVLCINEBSet

# from pymacy.neb.parse import NEBDBManager  # TODO: fix this
from atomate.utils.utils import get_logger

__author__ = "Hanmei Tang, Iek-Heng Chu"
__version__ = "1.0"

CWD = os.getcwd()
HOSTNAME = os.environ['HOSTNAME'].split('-')[0]
ABC_MIN_VAL = 8.0  # Angstrom
PATH_LEN_MAX_VAL = 7.0  # Angstrom
IMAGE_DIST = 0.7  # Angstrom
# Distance tolerance used to match the atomic
# indices between start and end structures.
SORT_TOL = 0  # Sort off

logger = get_logger(__name__, level=logging.DEBUG)
logger.info("CI-NEB workflow started at %s." % HOSTNAME)
logger.info("Working directory: %s" % CWD)


# @explicit_serialize
# class GetEndpointsFromIndexTask(FireTaskBase):
#     """
#     This class reads in one perfect structure and gets hopping site indices
#     from fw_spec["site_indices"]. The two endpoint structures are generated
#     and stored as dict in fw_spec["endpoints_0"] and fw_spec["endpoints_1"].
#
#     Some basic checks are executed:
#         1. site indexes are exactly 2 (# TODO: Extend to 3 or more later)
#         2. perfect cell all(abc length) > ABC_MIN_VAL (otherwise warning only)
#         3. endpoints distance < PATH_LEN_MAX_VAL (otherwise warning only)
#
#     Required parameters:
#         parent_structure (Structure): A parent structure in supercell.
#
#     fw_spec:
#         site_indices ([int, int]): a two-element list indicating site indices
#         from which a hop is defined.
#     """
#
#     required_params = ["parent_structure"]
#
#     def run_task(self, fw_spec):
#
#         logger.info("GetEndpointsFromIndexTask")
#
#         structure = self["parent_structure"]
#
#         # Site indexes check
#         site_indices = fw_spec["site_indices"]
#         if len(site_indices) != 2 or len(set(site_indices)) != 2:
#             raise ValueError("Two distinct indices are required"
#                              "to define a hop.")
#
#         # Check hopping sites species
#         if structure[site_indices[0]].specie != structure[site_indices[1]].specie:
#             raise ValueError("The site indices must be "
#                              "associated with identical species!")
#
#         # Check cell size.
#         for i, latt in enumerate(structure.lattice.abc):
#             if latt < ABC_MIN_VAL:
#                 logger.warn("Lattice parameter {} along {} is less "
#                             "than {}!".format(latt, i + 1, ABC_MIN_VAL))
#
#         # Endpoints distance check.
#         ep_fcoords = [structure[i].frac_coords for i in site_indices]
#         dist = structure.lattice. \
#             get_distance_and_image(ep_fcoords[0], ep_fcoords[1])[0]
#
#         if dist > PATH_LEN_MAX_VAL:
#             logger.warning("Endpoints distance is too large,"
#                            "exceeding upper limit:{0:.2f} Angstrom!".
#                            format(PATH_LEN_MAX_VAL))
#
#         s = structure.copy()
#         sites = s.sites
#
#         # Move hopping atoms to the beginning of species index.
#         init_site = sites[site_indices[0]]
#         final_site = sites[site_indices[1]]
#         sites.remove(init_site)
#         sites.remove(final_site)
#
#         init_sites = copy.deepcopy(sites)
#         final_sites = copy.deepcopy(sites)
#
#         init_sites.insert(0, final_site)
#         final_sites.insert(0, init_site)
#
#         s_0 = Structure.from_sites(init_sites)
#         s_1 = Structure.from_sites(final_sites)
#
#         endpoints = {0: s_0.as_dict(), 1: s_1.as_dict()}
#
#         return FWAction(update_spec={"endpoints": endpoints})


@explicit_serialize
class WriteEndpointInputTask(FiretaskBase):
    """
    This class calls a defect cell structure from fw_spec using endpoint_index.
    Endpoint input sets will be generated using MITRelaxSet with relatively
    strict convergence criteria.
    INCAR settings can be overwrote by fw_spec["ep_user_incar_settings"]
    and final INCAR will update fw_spec["ep_user_incar_settings"] value.

    Note:
        NELECT is set by user (remind the missing NELECT by throwing a warning).

    Required parameters:
        endpoint_index (int): 0 or 1, indicating the endpoint to be used.

    fw_spec:
        ep_user_incar_settings (dict): Additional INCAR settings.
    """

    required_params = ["ep_index", "user_incar_settings"]

    def run_task(self, fw_spec):
        logger.info("WriteEndpointInputTask")
        index = self["ep_index"]

        endpoints = fw_spec["endpoints"]
        structure = Structure.from_dict(endpoints[index])
        defaults = {"EDIFF": 5e-05, "EDIFFG": -0.05, "ISMEAR": 0,
                    "LCHARG": "False"}

        if self["user_incar_settings"] != {}:
            defaults.update(self["user_incar_settings"])

        m = MVLCINEBEndPointSet(structure,
                                user_incar_settings=defaults)
        output_dir = os.path.abspath(".")
        m.write_input(output_dir=output_dir)

        update_spec = {"endpoint_incar": m.incar.as_dict()}
        return FWAction(update_spec=update_spec)


# @explicit_serialize
# class GetImagesFromEndpointsTask(FireTaskBase):
#     """
#     This class generates images from endpoints.
#     """
#     required_params = ["endpoints"]

# @explicit_serialize
# class WriteNEBInputTask(FireTaskBase):
#     """
#     This class generates CI-NEB input sets.
#     INCAR settings can be overwrote by fw_spec["neb_user_incar_settings_init"]
#     or fw_spec["neb_user_incar_settings_accu"] and final INCAR will be updated
#     to fw_spec accordingly.
#
#     Some basic checks for "start" mode are executed:
#         1. endpoint structures must have the same sequence of species.
#
#     MODE options:
#     "start" mode:
#         Generate CI-NEB input sets using endpoint structures.
#         This will call fw_spec["endpoints_0"] and fw_spec["endpoints_1"]
#         and do linear interpolate. The number of images:
#         1) set by fw_spec - user_incar_settings["IMAGES"];
#         2) if not set by user, remind by throwing a warning and
#            calculate using IMAGE_DIST.
#     "continue" mode:
#         Generate CI-NEB input sets using image structures.
#         This simply uses existing images from fw_spec["images"].
#
#     STEP options:
#     "initial" step:
#         Use INCAR setting from fw_spec["neb_user_incar_settings_init"]
#     "accurate" step:
#         Use INCAR setting from fw_spec["neb_user_incar_settings_accu"]
#
#     Required parameters:
#         from_scratch (boolean): Whether to run a CI-NEB from scratch or not.
#
#     Optional parameters:
#         sort_tol (float): Distance tolerance (in Angstrom) used to match the atomic
#                     indices between start and end structures. If it is set 0, then
#                     no sorting will be performed.
#
#     fw_spec:
#         neb_user_incar_settings_init (dict): Additional INCAR settings.
#         neb_user_incar_settings_accu (dict): Additional INCAR settings.
#     """
#
#     required_params = ["from_scratch"]
#     optional_params = ["sort_tol"]
#
#     def run_task(self, fw_spec):
#
#         logger.info("WriteNEBInputTask")
#
#         from_scratch = self["from_scratch"]
#         nimages = fw_spec["_queueadapter"]["nnodes"]
#         user_incar_settings = fw_spec.get("user_incar_settings", {})
#         sort_tol = self.get("sort_tol", 0)
#
#         if from_scratch:
#             endpoint_0 = Structure.from_dict(fw_spec["relaxed_endpoints"][0])
#             endpoint_1 = Structure.from_dict(fw_spec["relaxed_endpoints"][1])
#
#             # Linear interpolation.
#             try:
#                 images = endpoint_0.interpolate(endpoint_1,
#                                                 nimages=nimages + 1,
#                                                 autosort_tol=sort_tol)
#             except Exception as e:
#
#                 if "Unable to reliably match structures " in str(e):
#                     logger.warn("Auto sorting is turned off because it is "
#                                 "unable to match the end-point structures!")
#                     images = endpoint_0.interpolate(endpoint_1,
#                                                     nimages=nimages + 1,
#                                                     autosort_tol=0)
#                 else:
#                     raise e
#
#         else:
#             images = [Structure.from_dict(d) for d in fw_spec["images"]]
#
#         # Set INCAR.
#         defaults = {"IMAGES": nimages}
#         defaults.update(user_incar_settings)
#         neb = MVLCINEBSet(images, user_incar_settings=defaults)
#         neb.write_input(output_dir=".")
#         update_spec = {"neb_incar": neb.incar.as_dict()}
#
#         return FWAction(update_spec=update_spec)


@explicit_serialize  # TODO: Remove this class
class RunVASPCustodianTask(FiretaskBase):
    """
    This class runs vasp job using custodian.

    Required:
        jobs ([VaspJob.as_dict()] or [VaspNEBJob.as_dict()]): List of VaspJobs.
        handlers ([ErrorHandlers.as_dict()]): List of Handlers.
        # is_neb (bool): CI-NEB run or not.

    Optional:
        validators : Custodian vasp validators.
        custodian_params: Parameters to be supplied to custodian init.
    """

    required_params = ["jobs", "handlers"]
    optional_params = ["validators", "custodian_params"]

    def run_task(self, fw_spec):

        logger.info("RunVASPCustodianTask")

        dec = MontyDecoder()
        jobs = dec.process_decoded(self["jobs"])
        fw_env = fw_spec.get("_fw_env", {})

        # Override VASP and gamma VASP commands using fw_env (not NEB)
        if fw_env.get("vasp_cmd"):
            for j in jobs:
                j.vasp_cmd = os.path.expandvars(fw_env["vasp_cmd"])
                j.gamma_vasp_cmd = j.gamma_vasp_cmd
                logger.info("Vasp command is {}".format(j.vasp_cmd))
        if fw_env.get("gamma_vasp_cmd"):
            for j in jobs:
                j.gamma_vasp_cmd = os.path.expandvars(fw_env["gamma_vasp_cmd"])
                logger.info("Vasp gamma command is {}".format(j.gamma_vasp_cmd))

        # Override custodian scratch dir.
        cust_params = self.get("custodian_params", {})
        if fw_env.get("scratch_root"):
            cust_params["scratch_dir"] = os.path.expandvars(
                fw_env["scratch_root"])

        handlers = dec.process_decoded(self["handlers"])
        validators = dec.process_decoded(self.get("validators", []))

        c = Custodian(handlers, jobs, validators=validators, **cust_params)
        output = c.run()

        logger.info("Custodian job has been executed!")
        return FWAction(stored_data=output)


@explicit_serialize
class SetupNEBTask(FiretaskBase):
    """
    This class generates an id for the given workflow, which will be
    stored in fw_spec.
    """

    def run_task(self, fw_spec):

        neb_id = fw_spec.get("neb_id", None)
        dest_root = fw_spec["_fw_env"]["run_dest_root"]
        dest_dir = os.path.join(dest_root,
                                os.environ["USER"], "neb")

        if neb_id is None:
            dirs = glob.glob(os.path.join(dest_dir, "*"))

            if len(dirs) > 0:
                neb_id = max([int(d.split("/")[-1])
                              for d in dirs]) + 1
            else:
                neb_id = 1

            update_spec = {"neb_id": neb_id, "path_id": 1}

        else:
            dirs = glob.glob(os.path.join(dest_dir,
                                          str(neb_id), "*"))

            if len(dirs) == 0:
                raise ValueError("There is no path subfolders "
                                 "inside {}!".format(neb_id))

            path_id = max([int(d.split("/")[-1]) for d in dirs]) + 1
            update_spec = {"neb_id": neb_id, "path_id": path_id}

        return FWAction(update_spec=update_spec)

# TODO: This is a backup TransferNEBTask, I will make it separated
# @explicit_serialize
# class TransferNEBTask(FireTaskBase):
#     """
#     This class transfers CI-NEB outputs from current directory to
#     destination directory. "type" is used to determine the type of
#     calculation and hence the final path. The corresponding structure
#     will be updated in fw_spec before files transferring.
#
#     Required:
#         - type (str): Type of calculation outputs being transferred.
#     """
#     required_params = ["type", "ep_index"]
#
#     def run_task(self, fw_spec):
#
#         neb_id = fw_spec["neb_id"]
#         path_id = fw_spec["path_id"]
#         type = self["type"]
#
#         if type not in ["parent", "endpoints", "neb"]:
#             raise ValueError("The provided type is not supported!")
#
#         logger.info("Transferring {} outputs to "
#                     "destination...".format(type))
#
#         dest_root = fw_spec["_fw_env"]["run_dest_root"]
#         src = os.path.abspath(".")
#         dest_dir = os.path.join(dest_root, os.environ["USER"],
#                                 "neb", neb_id, path_id, type)
#
#         if type == "endpoints":
#             ep_index = self["ep_index"]
#             ep_dir = "{:02}".format(ep_index)
#             dest_dir = os.path.join(dest_dir, ep_dir)
#
#         shutil.copytree(src, dest_dir)
#
#         # Update fw_spec based on the type of calculations.
#         if type == "neb":
#             # Get all result images.
#             subs = glob.glob("[0-9][0-9]")
#             nimages = len(subs)
#             images = []
#
#             concar_list = ["{:02}/CONTCAR".format(i)
#                            for i in range(nimages)[1: -1]]
#
#             for contcar in concar_list:
#                 images.append(Structure.from_file(contcar))
#
#             # End point structures
#             images.insert(0, Structure.from_file("00/POSCAR"))
#             images.append(Structure.from_file("{:02}/POSCAR".
#                                               format(nimages - 1)))
#             images = [s.as_dict() for s in images]
#             update_spec = {"images": images, "neb_dirs": dest_dir}
#
#         elif type == "endpoints":
#
#             # Update relaxed structure
#             file = glob.glob("CONTCAR*")[0]
#             ep = Structure.from_file(file, False)
#
#             # Check if the other relaxed ep structure is ready.
#             endpoints = fw_spec.get("relaxed_endpoints", {})
#             endpoint_dirs = fw_spec.get("endpoint_dirs", {})
#             _queueadapter = fw_spec["_queueadapter"]
#             vasp_cmd = fw_spec["vasp_cmd"]
#
#             if len(endpoints.keys()) == 1:
#                 k = endpoints.keys()[0]
#                 if k == ep_index:
#                     raise ValueError("Identical endpoint index {} already "
#                                      "exists!".format(ep_index))
#
#                 ep_2 = Structure.from_dict(endpoints[k])
#                 max_dist = max(get_ep_distances(ep, ep_2))
#                 nimages = int(max_dist / IMAGE_DIST) or 4
#
#                 # Update number of nodes used for actual CI-NEB calculation.
#                 _queueadapter.update({"nnodes": nimages})
#                 index = vasp_cmd.index('-np')
#                 np = nimages * int(vasp_cmd[index + 1])
#                 vasp_cmd[index + 1] = str(np)
#
#             endpoints[ep_index] = ep.as_dict()
#             endpoint_dirs[ep_index] = dest_dir
#
#             # Set update_spec
#             update_spec = ({"relaxed_endpoints": endpoints,
#                             "endpoint_dirs": endpoint_dirs,
#                             "_queueadapter": _queueadapter,
#                             "vasp_cmd": vasp_cmd})
#
#         elif type == "parent":
#             f = glob.glob("CONTCAR*")[0]
#             s = Structure.from_file(f, False)
#             update_spec = {"relaxed_parent": s.as_dict(),
#                            "parent_dir": dest_dir}
#
#         # Clear current directory.
#         for d in os.listdir(src):
#             try:
#                 os.remove(os.path.join(src, d))
#                 logger.info("Removing {}...".format(src))
#             except:
#                 pass
#
#         return FWAction(update_spec=update_spec)


@explicit_serialize
class TransferNEBTask(FiretaskBase):
    """
    This class transfers CI-NEB outputs from current directory to
    destination directory. "type" is used to determine the type of
    calculation and hence the final path. The corresponding structure
    will be updated in fw_spec before files transferring.

    Required:
        #- type (str): Type of calculation outputs being transferred.
        dest_dir (str): Destination path
    """
    required_params = ["dest_dir"]

    def run_task(self, fw_spec):
        dest_root = fw_spec["_fw_env"]["run_dest_root"]
        src = os.path.abspath(".")
        dest_dir = os.path.join(dest_root, os.environ["USER"],
                                "neb", self["dest_dir"])
        shutil.copytree(src, dest_dir)

        # Clear current directory.
        for d in os.listdir(src):
            try:
                os.remove(os.path.join(src, d))
                logger.info("Removing {}...".format(src))
            except:
                pass


# @explicit_serialize
# class InsertNEBDBTask(FireTaskBase):
#     """
#         This class inserts a complete CI-NEB calculations to neb database.
#         A perfect cell structure, endpoints_dirs/vasprun.xml* and neb_dirs
#         are required.
#         Currently uses debug_db for test purpose.
#     """
#
#     def run_task(self, fw_spec):
#
#         logger.info("InsertNEBDBTask")
#
#         # Create neb_id
#         for poscar in glob.glob(os.path.join(fw_spec["perfect_dir"], "POSCAR*")):
#             parent = Structure.from_file(poscar)
#             break
#
#         n = NEBDBManager(debug=True)  # TODO: under test
#         neb_id = n.insert_neb(parent)
#
#         # Insert endpoints
#         endpoints_dirs = fw_spec["endpoints_dirs"]
#         filenames = []
#         outcars = []
#         for ep in endpoints_dirs:
#             for vasprun in glob.glob(os.path.join(ep, "vasprun.xml*")):
#                 filenames.append(vasprun)
#                 break
#             for outcar in glob.glob(os.path.join(ep, "OUTCAR")):
#                 outcars.append(outcar)
#                 break
#
#         neb_dirs = fw_spec["neb_dirs"]
#         for path_id in range(1, len(neb_dirs) + 1):
#             # Copy OUTCARs from endpoint folders to neb folders
#             n_images = int(fw_spec["_queueadapter"]["nnodes"]) + 1
#             ini_end = ["00", "{:02}".format(n_images)]
#             for sour, dest in zip(outcars, ini_end):
#                 destination = os.path.join(neb_dirs[path_id], dest, "OUTCAR")
#                 shutil.copyfile(sour, destination)
#
#             # Insert neb path. Sort using path_id from small to large.
#             n.insert_neb_path(directory=neb_dirs[path_id],
#                               neb_id=neb_id, path_id=path_id)
#             n.insert_end_points(filenames=filenames,
#                                 neb_id=neb_id, path_id=path_id)
#
#         update_spec = {"neb_id": neb_id}
#
#         return FWAction(update_spec=update_spec)

# def get_ep_distances(ep_0, ep_1):
#     """
#     Calculate a list of site distances between two endpoints. The atomic
#     sequence in the two endpoints is assumed the same. Confirm data
#     validity in NEB workflow.
#     Args:
#         ep_0 (Structure): the first endpoint structure.
#         ep_1 (Structure): the second endpoint structure.
#     Returns:
#         dist (list): a list of atomic distances between two structures.
#     """
#     if ep_0.formula != ep_1.formula:
#         raise ValueError("Formula mismatch!")
#
#     for i, site in enumerate(ep_0):
#         if site.specie.symbol != ep_1[i].specie.symbol:
#             raise ValueError("Site sequence mismatch!")
#
#     fcoords_cp = [(site0.frac_coords, site1.frac_coords)
#                   for site0, site1 in zip(ep_0, ep_1)]
#
#     return [ep_0.lattice.get_distance_and_image(fc[0], fc[1])[0]
#             for fc in fcoords_cp]

@explicit_serialize
class WriteNEBFromImages(FiretaskBase):
    """
    Generate CI-NEB input sets using image structures.
    A list of POSCAR*/CONTCAR* files are provided, which will be used
    to generate KPOINTS, INCAR, POTCAR.

    Notes:
        'nnodes' (number of nodes) will be passed to next FireTask.

    Required parameters:
        vasp_input_set (AbstractVaspInputSet): full VaspInputSet object
        image_files (list of path str): a list of paths to image structure files
                                        Eg. ["00/POSCAR", "01/POSCAR"]

    Optional parameters:
        user_incar_settings (dict): additional INCAR settings.
    """
    required_params = ["vasp_input_set", "image_files"]
    optional_params = ["user_incar_settings"]

    def run_task(self, fw_spec):
        vis_orig = self["vasp_input_set"]
        image_files = self["image_files"]
        user_incar_settings = self.get("user_incar_settings", {})

        # Check images consistence.
        images = [Structure.from_file(i) for i in image_files]
        atomic_numbers = images[0].atomic_numbers
        for i in images:
            if i.atomic_numbers != atomic_numbers:
                raise ValueError("Images are inconsistent!")

        # Set INCAR.
        nimages = len(images)
        defaults = {"IMAGES": nimages}
        defaults.update(user_incar_settings)

        # Check vasp_input_set
        if hasattr(vis_orig, 'write_input'):
            vis = vis_orig(structures=images, user_incar_settings=defaults)
        else:
            raise ValueError("Unknown vasp_input_set!")
        vis.write_input(output_dir=".")
        update_spec = {"_queueadapter": {"nnodes": nimages}}

        return FWAction(update_spec=update_spec)


@explicit_serialize
class WriteNEBFromEndpoints(FiretaskBase):
    """
    Generate CI-NEB input sets using endpoint structures.
    The number of images:
        1) search in "user_incar_settings";
        2) otherwise, calculate using "image_dist".

    Required parameters:
        vasp_input_set (AbstractVaspInputSet): full VaspInputSet object
        endpoints (list of path str): Eg. ["E0/POSCAR", "E1/POSCAR"]

    Optional parameters:
        user_incar_settings (dict): additional INCAR settings.
        sort_tol (float): Distance tolerance (in Angstrom) used to match the atomic
                    indices between start and end structures. If it is set 0, then
                    no sorting will be performed.
        image_dist (float): distance in Angstrom, used in calculating number of images.
                            Default 0.7 Angstrom.

    fw_spec:
    """

    required_params = ["vasp_input_set", "endpoint_files"]
    optional_params = ["user_incar_settings", "sort_tol", "image_dist"]

    def run_task(self, fw_spec):
        self._set_params()

        # Get number of images.
        user_incar_settings = self.user_incar_settings
        if "IMAGES" in user_incar_settings:
            nimages = user_incar_settings["IMAGES"]
        else:
            nimages = self._get_nimages()

        # Dump images to temporary POSCARs.
        image_files = []
        images = self._get_images_by_linear_interp(nimages=nimages)
        for i in images:
            i.to(fmt='poscar', filename="POSCAR_tmp_{}".format(i))
            image_files.append(os.path.abspath("POSCAR_tmp_{}".format(i)))

        # TODO: Not sure if this works
        WriteNEBFromImages(vasp_input_set=self["vasp_input_set"],
                           image_files=image_files,
                           user_incar_settings=user_incar_settings)

        for i in image_files:
            os.remove(i)

    def _set_params(self):
        self.ep_0 = Structure.from_file(self["endpoint_files"][0])
        self.ep_1 = Structure.from_file(self["endpoint_files"][1])
        self.user_incar_settings = self.get("user_incar_settings", {})
        self.image_dist = self.get("image_dist", 0.7)
        self.sort_tol = self.get("sort_tol", 0)

    def _get_nimages(self):
        """
        Calculate the number of images using "image_dist", which can be
        overwritten in a optional_params list.

        Returns:
            nimages (int): number of images.
        """
        ep_0 = self.ep_0
        ep_1 = self.ep_1
        image_dist = self.image_dist

        # Check endpoints consistence.
        if ep_0.atomic_numbers != ep_1.atomic_numbers:
            raise ValueError("Endpoints are inconsistent!")

        max_dist = 0
        for s_0, s_1 in zip(ep_0, ep_1):
            site_dist = ep_0.lattice.get_distance_and_image(s_0.frac_coords,
                                                            s_1.frac_coords)[0]
            if site_dist > max_dist:
                max_dist = site_dist

        # Number of images must more than one.
        nimages = int(max_dist / image_dist) or 1

        return nimages

    def _get_images_by_linear_interp(self, nimages):
        ep_0 = self.ep_0
        ep_1 = self.ep_1
        sort_tol = self.sort_tol

        try:
            images = ep_0.interpolate(ep_1,
                                      nimages=nimages + 1,
                                      autosort_tol=sort_tol)
        except Exception as e:
            if "Unable to reliably match structures " in str(e):
                logger.warn("Auto sorting is turned off because it is "
                            "unable to match the end-point structures!")
                images = ep_0.interpolate(ep_1,
                                          nimages=nimages + 1,
                                          autosort_tol=0)
            else:
                raise e

        return images
