# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module defines tasks that support running vasp in various ways.
A Vasp emulator RunNEBVaspFake, which is isolated to avoid conflict.
"""

import shutil
import shlex
import subprocess
import os
import six
import glob

from pymatgen.io.vasp import Incar, Kpoints, Poscar, Potcar
from pymatgen.io.vasp.sets import get_vasprun_outcar
from pymatgen.electronic_structure.boltztrap import BoltztrapRunner

from custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler, AliasingErrorHandler, MeshSymmetryErrorHandler, \
    UnconvergedErrorHandler, MaxForceErrorHandler, PotimErrorHandler, FrozenJobErrorHandler, \
    NonConvergingErrorHandler, PositiveEnergyErrorHandler, WalltimeHandler
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VasprunXMLValidator, VaspFilesValidator

from fireworks import explicit_serialize, FiretaskBase

from atomate.utils.utils import env_chk, get_logger

__author__ = 'Anubhav Jain <ajain@lbl.gov>'
__credits__ = 'Shyue Ping Ong <ong.sp>'

logger = get_logger(__name__)


# @explicit_serialize
# class RunNEBVaspDirect(FiretaskBase):
#     """
#     Run VASP directly for CI-NEB (no custodian).
#
#     Required params:
#         vasp_cmd (str): the name of the full executable for running VASP.
#         Supports env_chk.
#     """
#     # TODO: Combine this method with RunVaspDirect in run_calc.py.
#     required_params = ["vasp_cmd"]
#
#     def run_task(self, fw_spec):
#         vasp_cmd = env_chk(self["vasp_cmd"], fw_spec)
#         logger.info("Running VASP using exe: {}".format(vasp_cmd))
#         return_code = subprocess.call(vasp_cmd, shell=True)
#         logger.info("VASP finished running with returncode: {}".format(return_code))


@explicit_serialize
class RunNEBVaspFake(FiretaskBase):
    """
     Vasp Emulator for CI-NEB, which has a different file arrangement.

     Required params:
         ref_dir (string): Path to reference vasp run directory with input files
            in the folder named 'inputs' and output files in the folder named 'outputs'.

     Optional params:
         params_to_check (list): optional list of incar parameters to check.
     """
    required_params = ["ref_dir"]
    optional_params = ["params_to_check"]

    def run_task(self, fw_spec):
        self._get_params()
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _get_params(self):
        """Define some convenient variables."""
        self.VASP_NEB_OUTPUT_FILES = ['INCAR', 'KPOINTS', 'POTCAR', 'vasprun.xml']
        self.VASP_NEB_OUTPUT_SUB_FILES = ['CHG', 'CHGCAR', 'CONTCAR', 'DOSCAR',
                                          'EIGENVAL', 'IBZKPT', 'PCDAT', 'POSCAR',
                                          'PROCAR', 'OSZICAR', 'OUTCAR', 'REPORT',
                                          'WAVECAR', 'XDATCAR']
        self.user_dir = os.getcwd()
        self.ref_dir_input = os.path.join(self["ref_dir"], "inputs")
        self.ref_dir_output = os.path.join(self["ref_dir"], "outputs")

        user_sdir = glob.glob(os.path.join(os.getcwd(), "[0-9][0-9]"))
        ref_sdir_input = glob.glob(os.path.join(self["ref_dir"],
                                                "inputs",
                                                "[0-9][0-9]"))
        ref_sdir_output = glob.glob(os.path.join(self["ref_dir"],
                                                 "outputs",
                                                 "[0-9][0-9]"))
        user_sdir.sort()
        ref_sdir_input.sort()
        ref_sdir_output.sort()

        # Check sub-folders consistence.
        if len(user_sdir) != len(ref_sdir_input):
            raise ValueError("Sub-folder numbers are inconsistent!"
                             "Paths are:\n{}\n{}".format(self.user_dir,
                                                         self.ref_dir_input))
        self.user_sdir = user_sdir
        self.ref_sdir_input = ref_sdir_input
        self.ref_sdir_output = ref_sdir_output

    def _verify_inputs(self):
        """Validation of input files under user NEB directory."""
        user_incar = Incar.from_file(os.path.join(self.user_dir, "INCAR"))
        ref_incar = Incar.from_file(os.path.join(self.ref_dir_input, "INCAR"))

        # Check INCAR
        params_to_check = self.get("params_to_check", [])
        defaults = {"ICHAIN": 0, "LCLIMB": True}
        for p in params_to_check:
            if user_incar.get(p, defaults.get(p)) != \
                    ref_incar.get(p, defaults.get(p)):
                raise ValueError("INCAR value of {} is "
                                 "inconsistent!".format(p))

        # Check KPOINTS
        user_kpoints = Kpoints.from_file(os.path.join(self.user_dir,
                                                      "KPOINTS"))
        ref_kpoints = Kpoints.from_file(os.path.join(self.ref_dir_input,
                                                     "KPOINTS"))
        if user_kpoints.style != ref_kpoints.style or \
           user_kpoints.num_kpts != ref_kpoints.num_kpts:
            raise ValueError("KPOINT files are inconsistent! "
                             "Paths are:\n{}\n{}".format(self.user_dir,
                                                         self.ref_dir_input))

        # Check POTCAR
        user_potcar = Potcar.from_file(os.path.join(self.user_dir,
                                                    "POTCAR"))
        ref_potcar = Potcar.from_file(os.path.join(self.ref_dir_input,
                                                   "POTCAR"))
        if user_potcar.symbols != ref_potcar.symbols:
            raise ValueError("POTCAR files are inconsistent! "
                             "Paths are:\n{}\n{}".format(self.user_dir,
                                                         self.ref_dir_input))

        # Check POSCARs
        for u, r in zip(self.user_sdir, self.ref_sdir_input):
            user_poscar = Poscar.from_file(os.path.join(u, "POSCAR"))
            ref_poscar = Poscar.from_file(os.path.join(r, "POSCAR"))
            if user_poscar.natoms != ref_poscar.natoms or \
               user_poscar.site_symbols != ref_poscar.site_symbols:
                raise ValueError("POSCAR files are inconsistent! "
                                 "Paths are:\n{}\n{}".format(u, r))

        logger.info("RunNEBVaspFake: verified inputs successfully.")

    def _clear_inputs(self):
        """Remove all input files from user NEB directory."""
        # Clear neb directory
        for x in self.VASP_NEB_OUTPUT_FILES:
            p = os.path.join(os.getcwd(), x)
            if os.path.exists(p):
                os.remove(p)

        # Clear neb sub-directory
        for d in self.user_sdir:
            for x in self.VASP_NEB_OUTPUT_SUB_FILES:
                p = os.path.join(d, x)
                if os.path.exists(p):
                    os.remove(p)

    def _generate_outputs(self):
        """Copy valid outputs from a reference folder to user NEB directory."""
        # Copy NEB files.
        for file_name in os.listdir(self.ref_dir_output):
            full_file_name = os.path.join(self.ref_dir_output, file_name)
            if os.path.isfile(full_file_name):
                shutil.copy(full_file_name, os.getcwd())

        # Copy NEB sub-files.
        for u_dir, r_dir in zip(self.user_sdir, self.ref_sdir_output):
            for file_name in os.listdir(r_dir):
                full_file_name = os.path.join(r_dir, file_name)
                if os.path.isfile(full_file_name):
                    shutil.copy(full_file_name, u_dir)

        logger.info("RunNEBVaspFake: ran fake VASP, generated outputs.")
