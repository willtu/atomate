# coding: utf-8

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import unittest
from collections import defaultdict

from fireworks import FiretaskBase, Firework, Workflow, explicit_serialize, FWAction

from atomate.utils.utils import env_chk, get_logger, get_mongolike
from atomate.utils.utils import append_fw_wf, remove_fws, recursive_get_result

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

logger = get_logger(__name__)


@explicit_serialize
class Task1(FiretaskBase):
    def run_task(self, fw_spec):
        print("task1", fw_spec)
        return FWAction(stored_data={"color": "red"})


@explicit_serialize
class Task2(FiretaskBase):
    def run_task(self, fw_spec):
        print("task2", fw_spec)
        return FWAction(stored_data={"color": "yellow"})


class UtilsTests(unittest.TestCase):

    @classmethod
    def setUp(cls):
        cls.fw1 = Firework(Task1())
        cls.fw2 = Firework([Task2(), Task2()], parents=cls.fw1)
        cls.fw3 = Firework(Task1(), parents=cls.fw1)

    def test_env_chk(self):
        fw_spec_valid = {"_fw_env": {"hello": "there"}}
        fw_spec_invalid = {}

        self.assertEqual(env_chk("hello", fw_spec_valid), "hello")
        self.assertEqual(env_chk([1, 2, 3], fw_spec_valid), [1, 2, 3])
        self.assertEqual(env_chk(defaultdict(int), fw_spec_valid),
                         defaultdict(int))
        self.assertEqual(env_chk(">>hello<<", fw_spec_valid), "there")
        self.assertRaises(KeyError, env_chk, ">>hello1<<", fw_spec_valid)
        self.assertEqual(env_chk(">>hello1<<", fw_spec_valid, False), None)

        self.assertRaises(KeyError, env_chk, ">>hello<<", fw_spec_invalid)
        self.assertEqual(env_chk(">>hello<<", fw_spec_invalid, False), None)
        self.assertEqual(env_chk(">>hello<<", fw_spec_invalid,
                                 False, "fallback"), "fallback")

        self.assertEqual(env_chk(None, fw_spec_valid, False), None)
        self.assertEqual(env_chk(None, fw_spec_valid, False, "fallback"),
                         "fallback")

    def test_get_mongolike(self):
        d = {"a": [{"b": 1}, {"c": {"d": 2}}], "e": {"f": {"g": 3}}, "g": 4, "h":[5, 6]}

        self.assertEqual(get_mongolike(d, "g"), 4)
        self.assertEqual(get_mongolike(d, "e.f.g"), 3)
        self.assertEqual(get_mongolike(d, "a.0.b"), 1)
        self.assertEqual(get_mongolike(d, "a.1.c.d"), 2)
        self.assertEqual(get_mongolike(d, "h.-1"), 6)

    # TODO: @matk86 - remove this once removing append_fw_wf() method -computron
    def test_append_fw(self):
        fw_new = Firework(Task1())
        fws = [self.fw1, self.fw2, self.fw3]
        wflow = Workflow(fws)
        leaf_ids = wflow.leaf_fw_ids
        append_fw_wf(wflow, fw_new)
        new_lead_ids = wflow.leaf_fw_ids
        for i in leaf_ids:
            self.assertEqual(wflow.links[i], [new_lead_ids[0]])

    # TODO: @matk86 - move this test to FireWorks after moving remove_fws() there
    def test_remove_leaf_fws(self):
        fw4 = Firework(Task1(), parents=[self.fw2, self.fw3])
        fws = [self.fw1, self.fw2, self.fw3, fw4]
        wflow = Workflow(fws)
        leaf_ids = wflow.leaf_fw_ids
        parents = []
        for i in leaf_ids:
            parents.extend(wflow.links.parent_links[i])
        new_wf = remove_fws(wflow, wflow.leaf_fw_ids)
        new_leaf_ids = new_wf.leaf_fw_ids
        self.assertEqual(new_leaf_ids, parents)

    # TODO: @matk86 - move this test to FireWorks after moving remove_fws() there
    def test_remove_root_fws(self):
        fw4 = Firework(Task1(), parents=[self.fw2, self.fw3])
        fws = [self.fw1, self.fw2, self.fw3, fw4]
        wflow = Workflow(fws)
        root_ids = wflow.root_fw_ids
        children = []
        for i in root_ids:
            children.extend(wflow.links[i])
        new_wf = remove_fws(wflow, wflow.root_fw_ids)
        new_root_ids = new_wf.root_fw_ids
        self.assertEqual(sorted(new_root_ids), sorted(children))

    def test_recursive_get_result(self):
        # Basic functionality with dictionary key
        result_dict = {"output":[{}, {"data":[0, 1, 2, 3]}]}
        out = recursive_get_result({"my_data":">>output.-1.data.2"}, result_dict)
        self.assertEqual(out["my_data"], 2)
        # Basic functionality with attribute
        task1 = Task1()
        out_attr = recursive_get_result({"fw_name":"a>>_fw_name"}, task1)
        self.assertEqual(out_attr["fw_name"],"{{atomate.utils.tests.test_utils.Task1}}")
        # Testing as_dict functionality
        out_as_dict = recursive_get_result({"fw_name":">>_fw_name"}, task1)
        self.assertEqual(out_as_dict["fw_name"],"{{atomate.utils.tests.test_utils.Task1}}")
