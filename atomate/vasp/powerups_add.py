# coding: utf-8

from __future__ import division, print_function, unicode_literals, absolute_import

"""
This module modifies use_fake_vasp method in powerups.py once finished.
Keeping a copy to avoid conflict. Under testing.
"""

from fireworks import Workflow, FileWriteTask
from fireworks.core.firework import Tracker
from fireworks.utilities.fw_utilities import get_slug

from atomate.utils.utils import get_meta_from_structure, get_fws_and_tasks, update_wf
from atomate.vasp.firetasks.glue_tasks import CheckStability, CheckBandgap
from atomate.vasp.firetasks.run_calc import RunVaspCustodian, RunVaspDirect, RunVaspFake
from atomate.vasp.firetasks.write_inputs import ModifyIncar
from atomate.vasp.config import ADD_NAMEFILE, SCRATCH_DIR, ADD_MODIFY_INCAR
from atomate.vasp.firetasks.run_calc_neb import RunNEBVaspFake  # TODO: Combine source

from pymatgen.core import Structure

__author__ = 'Hanmei Tang'
__email__ = 'hat003@eng.ucsd.edu'


def use_neb_fake_vasp(original_wf, ref_dirs, params_to_check=None, is_neb=False):
    """
    Replaces all tasks with "RunVasp" (e.g. RunVaspDirect) to be
    RunVaspFake. Or replace all "RunNEBVasp" to be "RunNEBVaspFake".
    Thus, we do not actually run VASP but copy
    pre-determined inputs and outputs.

    Args:
        original_wf (Workflow)
        ref_dirs (dict): key=firework name, value=path to the reference vasp calculation directory
        params_to_check (list): optional list of incar parameters to check.
        is_neb (bool): using "RunNEBVaspFake"
    """
    if not params_to_check:  # TODO: Check the default INCAR setting.
        params_to_check = ["ISPIN", "ENCUT", "ISMEAR", "SIGMA", "IBRION", "LORBIT", "NBANDS", "LMAXMIX"]
    wf_dict = original_wf.to_dict()
    for idx_fw, fw in enumerate(original_wf.fws):
        for job_type in ref_dirs.keys():
            if job_type in fw.name:
                for idx_t, t in enumerate(fw.tasks):
                    if "RunVasp" in str(t):
                        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = \
                            RunVaspFake(ref_dir=ref_dirs[job_type], params_to_check=params_to_check).to_dict()
                    if "RunVasp" in str(t) and is_neb:  # TODO: Mark of additional tag.
                        wf_dict["fws"][idx_fw]["spec"]["_tasks"][idx_t] = \
                            RunNEBVaspFake(ref_dir=ref_dirs[job_type], params_to_check=params_to_check).to_dict()
    return Workflow.from_dict(wf_dict)

