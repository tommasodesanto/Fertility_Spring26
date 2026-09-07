"""Pure orchestration checks for isolated joint-nested policy cases."""
from __future__ import annotations

import os
from pathlib import Path
from queue import Empty
from tempfile import TemporaryDirectory
from unittest import TestCase, main

import run_e5f_joint_nested_finalize as finalizer


class FakeQueue:
    def __init__(self): self.items=[]
    def put(self,item): self.items.append(item)
    def get(self,timeout=None):
        if not self.items: raise Empty
        return self.items.pop(0)
    def close(self): pass
    def join_thread(self): pass


class FakeProcess:
    def __init__(self,target,args,keep_alive=False):
        self.target=target;self.args=args;self.keep_alive=keep_alive
        self.alive=False;self.exitcode=None;self.terminated=False;self.joined=False;self.pid=None
    def start(self):
        self.pid=id(self);self.alive=True;self.target(*self.args)
        if not self.keep_alive:
            self.alive=False;self.exitcode=0
    def is_alive(self): return self.alive
    def terminate(self): self.terminated=True;self.alive=False;self.exitcode=-15
    def kill(self): self.terminate()
    def join(self,timeout=None):
        self.joined=True
        if self.alive and not self.keep_alive:
            self.alive=False;self.exitcode=0


class FakeContext:
    def __init__(self,keep_alive=False): self.keep_alive=keep_alive;self.processes=[]
    def Queue(self): return FakeQueue()
    def Process(self,*,target,args):
        process=FakeProcess(target,args,self.keep_alive);self.processes.append(process)
        return process


def configured_process_runner(name,folder,*_):
    if finalizer.policy.calendar.apply_fertility is not finalizer.policy.transition.apply_sequential_fertility:
        raise RuntimeError('child fertility hook not configured')
    if finalizer.policy.calendar.advance_calendar_distribution is not finalizer.policy.transition.advance_sequential_calendar_distribution:
        raise RuntimeError('child transition hook not configured')
    folder.mkdir(parents=True,exist_ok=True)
    return dict(status='complete',pid=os.getpid(),case=name),None


class FinalizerParallelTests(TestCase):
    def runner(self,calls):
        def run(name,folder,*_):
            folder.mkdir(parents=True,exist_ok=True);calls.append((name,folder))
            return dict(status='complete',dates=2,gates=dict(case=name),source_summary_sha256='selected'),None
        return run

    def test_serial_and_four_process_paths_keep_case_order_and_separate_folders(self):
        with TemporaryDirectory() as tmp:
            root=Path(tmp);calls=[]
            serial=finalizer.run_policy_cases(root/'serial',None,None,None,Path('selected'),1,1,case_runner=self.runner(calls))
            context=FakeContext()
            parallel=finalizer.run_policy_cases(root/'parallel',None,None,None,Path('selected'),1,4,
                process_context=context,case_runner=self.runner(calls))
        self.assertEqual(list(serial[0]),list(finalizer.CASES))
        self.assertEqual(list(parallel[0]),list(finalizer.CASES))
        self.assertEqual(serial,parallel)
        self.assertEqual({folder.name for _,folder in calls},set(finalizer.CASES))
        self.assertEqual(len({folder for _,folder in calls if folder.parent.name == 'parallel'}),4)
        self.assertTrue(all(process.joined for process in context.processes))

    def test_real_spawn_initializes_hooks_in_four_distinct_processes(self):
        with TemporaryDirectory() as tmp:
            receipts,failures=finalizer.run_policy_cases(Path(tmp),None,None,None,Path('selected'),1,4,case_runner=configured_process_runner)
        self.assertFalse(failures)
        self.assertEqual(list(receipts),list(finalizer.CASES))
        pids={row['pid'] for row in receipts.values()}
        self.assertEqual(len(pids),4);self.assertNotIn(os.getpid(),pids)

    def test_failed_process_start_reaps_only_started_children(self):
        class StartFailureContext(FakeContext):
            def Process(self,**kwargs):
                process=super().Process(**kwargs)
                if len(self.processes)==2:
                    def broken_start():raise RuntimeError('cannot pickle worker input')
                    process.start=broken_start
                return process
        with TemporaryDirectory() as tmp:
            context=StartFailureContext(keep_alive=True)
            with self.assertRaisesRegex(RuntimeError,'cannot pickle worker input'):
                finalizer.run_policy_cases(Path(tmp),None,None,None,Path('selected'),1,4,process_context=context,case_runner=self.runner([]))
        self.assertTrue(context.processes[0].terminated and context.processes[0].joined)
        self.assertIsNone(context.processes[1].pid)
        self.assertFalse(context.processes[1].joined)

    def test_contract_allows_only_documented_worker_counts(self):
        self.assertEqual(finalizer.policy_worker_count({}),1)
        self.assertEqual(finalizer.policy_worker_count(dict(policy_workers=1)),1)
        self.assertEqual(finalizer.policy_worker_count(dict(policy_workers=4)),4)
        for value in (0,2,5,True,'4'):
            with self.assertRaisesRegex(RuntimeError,'exactly 1 or 4'):
                finalizer.policy_worker_count(dict(policy_workers=value))

    def test_unexpected_child_failure_cancels_and_reaps_all_children(self):
        def failing(name,*_):
            if name == 'baseline': raise RuntimeError('unexpected policy failure')
            return dict(status='complete'),None
        with TemporaryDirectory() as tmp:
            context=FakeContext(keep_alive=True)
            with self.assertRaisesRegex(RuntimeError,'baseline failed unexpectedly'):
                finalizer.run_policy_cases(Path(tmp),None,None,None,Path('selected'),1,4,
                    process_context=context,case_runner=failing)
        self.assertEqual(len(context.processes),4)
        self.assertTrue(all(process.terminated and process.joined for process in context.processes))


if __name__ == '__main__': main()
