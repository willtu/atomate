# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks for writing vasp input sets for various types of vasp calculations
"""

import os
from six.moves import range
from importlib import import_module

import numpy as np

from fireworks import FiretaskBase, explicit_serialize
from fireworks.utilities.dict_mods import apply_mod

from pymatgen.core.structure import Structure
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.io.vasp import Incar, Poscar
from pymatgen.io.vasp.sets import MPStaticSet, MPNonSCFSet, MPSOCSet, MPHSEBSSet

from atomate.utils.utils import env_chk, load_class, get_logger

__author__ = 'Anubhav Jain, Shyue Ping Ong, Kiran Mathew'
__email__ = 'ajain@lbl.gov'

logger = get_logger(__name__)


@explicit_serialize
class WriteVaspFromIOSet(FiretaskBase):
    """
    Create VASP input files using implementations of pymatgen's
    AbstractVaspInputSet. An input set can be provided as an object or as a
    String/parameter combo.

    Required params:
        structure (Structure): structure
        vasp_input_set (AbstractVaspInputSet or str): Either a VaspInputSet
            object or a string name for the VASP input set (e.g.,
            "MPRelaxSet").

    Optional params:
        vasp_input_params (dict): When using a string name for VASP input set,
            use this as a dict to specify kwargs for instantiating the input
            set parameters. For example, if you want to change the
            user_incar_settings, you should provide:
            {"user_incar_settings": ...}.
            This setting is ignored if you provide the full object
            representation of a VaspInputSet rather than a String.
    """

    required_params = ["structure", "vasp_input_set"]
    optional_params = ["vasp_input_params"]

    def run_task(self, fw_spec):
        # if a full VaspInputSet object was provided
        if hasattr(self['vasp_input_set'], 'write_input'):
            vis = self['vasp_input_set']

        # if VaspInputSet String + parameters was provided
        else:
            vis_cls = load_class("pymatgen.io.vasp.sets", self["vasp_input_set"])
            vis = vis_cls(self["structure"], **self.get("vasp_input_params", {}))
        vis.write_input(".")


@explicit_serialize
class WriteVaspFromPMGObjects(FiretaskBase):
    """
    Write VASP files using pymatgen objects.

    Required params:
        (none)

    Optional params:
        incar (Incar): pymatgen Incar object
        poscar (Poscar): pymatgen Poscar object
        kpoints (Kpoints): pymatgen Kpoints object
        potcar (Potcar): pymatgen Potcar object
    """

    optional_params = ["incar", "poscar", "kpoints", "potcar"]

    def run_task(self, fw_spec):
        if "incar" in self:
            self["incar"].write_file("INCAR")
        if "poscar" in self:
            self["poscar"].write_file("POSCAR")
        if "kpoints" in self:
            self["kpoints"].write_file("KPOINTS")
        if "potcar" in self:
            self["potcar"].write_file("POTCAR")


@explicit_serialize
class ModifyIncar(FiretaskBase):
    """
    Modify an INCAR file.

    Required params:
        (none)

    Optional params:
        incar_update (dict): overwrite Incar dict key. Supports env_chk.
        incar_multiply ([{<str>:<float>}]) - multiply Incar key by a constant
            factor. Supports env_chk.
        incar_dictmod ([{}]): use DictMod language to change Incar.
            Supports env_chk.
        input_filename (str): Input filename (if not "INCAR")
        output_filename (str): Output filename (if not "INCAR")
    """

    optional_params = ["incar_update", "incar_multiply", "incar_dictmod", "input_filename",
                       "output_filename"]

    def run_task(self, fw_spec):

        # load INCAR
        incar_name = self.get("input_filename", "INCAR")
        incar = Incar.from_file(incar_name)

        # process FireWork env values via env_chk
        incar_update = env_chk(self.get('incar_update'), fw_spec)
        incar_multiply = env_chk(self.get('incar_multiply'), fw_spec)
        incar_dictmod = env_chk(self.get('incar_dictmod'), fw_spec)

        if incar_update:
            incar.update(incar_update)

        if incar_multiply:
            for k in incar_multiply:
                incar[k] = incar[k] * incar_multiply[k]

        if incar_dictmod:
            apply_mod(incar_dictmod, incar)

        # write INCAR
        incar.write_file(self.get("output_filename", "INCAR"))


@explicit_serialize
class WriteVaspStaticFromPrev(FiretaskBase):
    """
    Writes input files for a static run. Assumes that output files from a
    relaxation job can be accessed. Also allows lepsilon calcs.

    Required params:
        (none)

    Optional params:
        (documentation for all optional params can be found in
        MPStaticSet)
    """

    required_params = ["prev_calc_dir"]
    optional_params = ["reciprocal_density", "small_gap_multiply", "standardize", "sym_prec",
                       "international_monoclinic", "lepsilon", "other_params"]

    def run_task(self, fw_spec):
        lepsilon = self.get("lepsilon")

        default_reciprocal_density = 100 if not lepsilon else 200
        other_params = self.get("other_params", {})

        # for lepsilon runs, set EDIFF to 1E-5 unless user says otherwise
        user_incar_settings = self.get("other_params", {}).get("user_incar_settings", {})

        if lepsilon and "EDIFF" not in user_incar_settings and "EDIFF_PER_ATOM" not in user_incar_settings:
            if "user_incar_settings" not in other_params:
                other_params["user_incar_settings"] = {}
            other_params["user_incar_settings"]["EDIFF"] = 1E-5

        vis = MPStaticSet.from_prev_calc(prev_calc_dir=self["prev_calc_dir"],
                                         reciprocal_density=self.get("reciprocal_density", default_reciprocal_density),
                                         small_gap_multiply=self.get("small_gap_multiply", None),
                                         standardize=self.get("standardize", False),
                                         sym_prec=self.get("sym_prec", 0.1),
                                         international_monoclinic=self.get("international_monoclinic", True),
                                         lepsilon=lepsilon, **other_params)
        vis.write_input(".")


@explicit_serialize
class WriteVaspHSEBSFromPrev(FiretaskBase):
    """
    Writes input files for HSE Gap run. Assumes that output files from a
    an NSCF job (for getting VBM/CBM) can be accessed.

    Required params:
        prev_calc_dir

    Optional params:
        (documentation for all optional params can be found in
        MPHSEBSSet)
    """

    required_params = ["prev_calc_dir"]
    optional_params = ["mode", "reciprocal_density", "kpoints_line_density"]

    def run_task(self, fw_spec):
        vis = MPHSEBSSet.from_prev_calc(self["prev_calc_dir"], mode=self.get("mode", "Uniform"),
                                        reciprocal_density=self.get("reciprocal_density", 50),
                                        kpoints_line_density=self.get("kpoints_line_density", 10),
                                        copy_chgcar=False)
        vis.write_input(".")


@explicit_serialize
class WriteVaspNSCFFromPrev(FiretaskBase):
    """
    Writes input files for an NSCF static run. Assumes that output files from an
    scf job can be accessed. There are many options, e.g. uniform mode,
    line mode, adding the optical properties, etc.

    Required params:
        (none)

    Optional params:
        (documentation for all optional params can be found in
        NonSCFVaspInputSet)
    """

    required_params = ["prev_calc_dir"]
    optional_params = ["copy_chgcar", "nbands_factor", "reciprocal_density", "kpoints_line_density",
                       "small_gap_multiply", "standardize", "sym_prec", "international_monoclinic",
                       "mode", "nedos", "optics", "other_params"]

    def run_task(self, fw_spec):
        vis = MPNonSCFSet.from_prev_calc(
            prev_calc_dir=self["prev_calc_dir"],
            copy_chgcar=self.get("copy_chgcar", False),
            nbands_factor=self.get("nbands_factor", 1.2),
            reciprocal_density=self.get("reciprocal_density", 100),
            kpoints_line_density=self.get("kpoints_line_density", 20),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get("international_monoclinic", True),
            mode=self.get("mode", "uniform"),
            nedos=self.get("nedos", 601),
            optics=self.get("optics", False),
            **self.get("other_params", {}))
        vis.write_input(".")


@explicit_serialize
class WriteVaspSOCFromPrev(FiretaskBase):
    """
    Writes input files for a spinorbit coupling calculation.

    Required params:
        prev_calc_dir: path to previous calculation
        magmom (list): magnetic moment values for each site in the structure.
        saxis (list): magnetic field direction

    Optional params:
        (none)
    """
    required_params = ["prev_calc_dir", "magmom", "saxis"]

    optional_params = ["copy_chgcar", "nbands_factor", "reciprocal_density", "small_gap_multiply",
                       "standardize", "sym_prec", "international_monoclinic", "other_params"]

    def run_task(self, fw_spec):
        vis = MPSOCSet.from_prev_calc(
            prev_calc_dir=self["prev_calc_dir"],
            magmom=self["magmom"],
            saxis=self["saxis"],
            copy_chgcar=self.get("copy_chgcar", False),
            nbands_factor=self.get("nbands_factor", 1.2),
            reciprocal_density=self.get("reciprocal_density", 100),
            small_gap_multiply=self.get("small_gap_multiply", None),
            standardize=self.get("standardize", False),
            sym_prec=self.get("sym_prec", 0.1),
            international_monoclinic=self.get("international_monoclinic", True),
            **self.get("other_params", {}))
        vis.write_input(".")


@explicit_serialize
class WriteTransmutedStructureIOSet(FiretaskBase):
    """
    Apply the provided transformations to the input structure and write the
    input set for that structure. Reads structure from POSCAR if no structure provided

    Required params:
        structure (Structure): input structure
        transformations (list): list of names of transformation classes as defined in
            the modules in pymatgen.transformations
        vasp_input_set (VaspInputSet): VASP input set.

    Optional params:
        transformation_params (list): list of dicts where each dict specifies the input parameters
            to instantiate the transformation class in the transformations list.
        override_default_vasp_params (dict): additional user input settings.
        prev_calc_dir: path to previous calculation if using structure from another calculation.
    """

    required_params = ["structure", "transformations", "vasp_input_set"]
    optional_params = ["prev_calc_dir", "transformation_params", "override_default_vasp_params"]

    def run_task(self, fw_spec):

        transformations = []
        transformation_params = self.get("transformation_params",
                                         [{} for i in range(len(self["transformations"]))])
        for t in self["transformations"]:
            for m in ["advanced_transformations", "defect_transformations",
                      "site_transformations", "standard_transformations"]:
                mod = import_module("pymatgen.transformations.{}".format(m))
                try:
                    t_cls = getattr(mod, t)
                except AttributeError:
                    continue
                t_obj = t_cls(**transformation_params.pop(0))
                transformations.append(t_obj)

        structure = self['structure'] if 'prev_calc_dir' not in self else \
                Poscar.from_file(os.path.join(self['prev_calc_dir'], 'POSCAR')).structure
        ts = TransformedStructure(structure)
        transmuter = StandardTransmuter([ts], transformations)
        final_structure = transmuter.transformed_structures[-1].final_structure.copy()

        vis_orig = self["vasp_input_set"]
        vis_dict = vis_orig.as_dict()
        vis_dict["structure"] = final_structure.as_dict()
        vis_dict.update(self.get("override_default_vasp_params", {}) or {})
        vis = vis_orig.__class__.from_dict(vis_dict)
        vis.write_input(".")


@explicit_serialize
class WriteNormalmodeDisplacedPoscar(FiretaskBase):
    """
    Displace the structure from the previous calculation along the provided normal mode by the
    given amount and write the corresponding Poscar file.

    Required params:
        mode (int): normal mode index
        displacement (float): displacement along the normal mode in Angstroms
    """

    required_params = ["mode", "displacement"]

    def run_task(self, fw_spec):
        mode = self["mode"]
        disp = self["displacement"]
        structure = Structure.from_file("POSCAR")
        nm_eigenvecs = np.array(fw_spec["normalmodes"]["eigenvecs"])
        nm_norms = np.array(fw_spec["normalmodes"]["norms"])

        # displace the sites along the given normal mode
        nm_displacement = nm_eigenvecs[mode, :, :] * disp / nm_norms[mode, :, np.newaxis]
        for i, vec in enumerate(nm_displacement):
            structure.translate_sites(i, vec, frac_coords=False)

        # write the modified structure to poscar
        structure.to(fmt="poscar", filename="POSCAR")


@explicit_serialize
class WriteNEBFromImages(FiretaskBase):
    """
    Generate CI-NEB input sets using given image structures.

    Required parameters:
        vasp_input_set (AbstractVaspInputSet): full VaspInputSet object

    Optional parameters:
        user_incar_settings (dict): additional INCAR settings.
        images (list[str]): list of files of image structures,
            including the two endpoints.

    fw_spec:
        images (list[dict]): list of image structure dict,
            including the two endpoints.
    """
    required_params = ["vasp_input_set"]
    optional_params = ["images", "user_incar_settings"]

    def run_task(self, fw_spec):
        logger.info("WriteNEBFromImages")

        vis_orig = self["vasp_input_set"]
        user_incar_settings = self.get("user_incar_settings", {})
        if "images" in self:
            images = [Structure.from_file(i) for i in self["images"]]
        else:
            images = [Structure.from_dict(i) for i in fw_spec["images"]]

        # Check images consistence.
        atomic_numbers = images[0].atomic_numbers
        for i in images:
            if i.atomic_numbers != atomic_numbers:
                raise ValueError("Images are inconsistent!")

        # Set INCAR.
        nimages = fw_spec["_queueadapter"]["nnodes"]
        defaults = {"IMAGES": nimages}
        defaults.update(user_incar_settings)

        # Check vasp_input_set.
        if hasattr(vis_orig, 'write_input'):
            vis = vis_orig(structures=images, user_incar_settings=defaults)
        else:
            raise ValueError("Unknown vasp_input_set!")

        vis.write_input(output_dir=".")


@explicit_serialize
class WriteNEBFromEndpoints(FiretaskBase):
    """
    Generate CI-NEB input sets using endpoint structures.
    The number of images:
        1) search in "user_incar_settings";
        2) otherwise, calculate using "image_dist".

    Required parameters:
        vasp_input_set (AbstractVaspInputSet): full VaspInputSet object

    Optional parameters:
        endpoints (list of path str): Eg. ["E0/POSCAR", "E1/POSCAR"]
        user_incar_settings (dict): additional INCAR settings.
        sort_tol (float): Distance tolerance (in Angstrom) used to match the atomic
                    indices between start and end structures. If it is set 0, then
                    no sorting will be performed.
        image_dist (float): distance in Angstrom, used in calculating number of images.
                            Default 0.7 Angstrom.

    fw_spec:
    """

    required_params = ["vasp_input_set"]
    optional_params = ["endpoint_files", "user_incar_settings",
                       "sort_tol", "image_dist"]

    def run_task(self, fw_spec):
        logger.info("WriteNEBSetFromEndpoints")

        self._set_params(fw_spec)

        # Get number of images.
        user_incar_settings = self.user_incar_settings
        if "IMAGES" in user_incar_settings:
            nimages = user_incar_settings["IMAGES"]
        else:
            nimages = self._get_nimages()

        image_files = []
        images = self._get_images_by_linear_interp(nimages=nimages)
        images = [i.as_dict() for i in images]
        fw_spec.update({"images": images})

        write = WriteNEBFromImages(vasp_input_set=self["vasp_input_set"],
                                   user_incar_settings=user_incar_settings)

        write.run_task(fw_spec=fw_spec)  # TODO: Double check

    def _set_params(self, fw_spec):
        if "endpoint_files" in self:
            self.ep_0 = Structure.from_file(self["endpoint_files"][0])
            self.ep_1 = Structure.from_file(self["endpoint_files"][1])
        else:
            self.ep_0 = Structure.from_dict(fw_spec["ep0_st"])
            self.ep_1 = Structure.from_dict(fw_spec["ep1_st"])

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
