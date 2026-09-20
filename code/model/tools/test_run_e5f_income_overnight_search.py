import tempfile
import unittest
from pathlib import Path

import run_e5f_income_overnight_search as search
from run_e5f_income_overnight_search import OBJECTIVE, PARAMETERS, load_seeds, proposals, run_search


class OvernightControllerTest(unittest.TestCase):
    def setUp(self):
        self.pilot = dict(zip(PARAMETERS, [.99, .3377, .3978, 1.05, 8.1, .081, .085, .265, 2.3]))
        self.plan = {"schema":"e5f_income_candidate_calibration_v1", "objective_canonical_sha256":OBJECTIVE,
                     "candidate_payload_fingerprint":"candidate-pin", "parameter_bounds":{n:[.01,20.] for n in PARAMETERS},
                     "pilot_parameters":self.pilot, "pilot_initial_psi":.15, "controller_sha256":"old-helper",
                     "prior_successful_cases":[{"label":f"prior_{i}","parameters":{**self.pilot,"theta0":.081+.01*(i+1)}} for i in range(4)]}
        self.plan["parameter_bounds"].update(beta_annual=[.94,.99],theta0=[0.,8.],first_birth_fixed_cost=[0.,8.],h_P=[.1,2.3])

    def receipt(self, loss, params):
        return {"schema":"e5f_initial_minimum_distance_result_v1","contract_sha256":OBJECTIVE,"loss":loss,"objective":loss,
                "normalization":{"target":2.1,"absolute_gap":0.,"completed_fertility":2.1,"psi_child":.15},"candidate_payload_fingerprint":"candidate-pin",
                "parameters":[{"parameter":n,"estimate":params[n],"structural_coordinate":True} for n in PARAMETERS]+[{"parameter":f"nuisance_{i}","estimate":1.,"structural_coordinate":False} for i in range(8)],
                "target_fit":[{"restriction_id":f"m{i}","target":1.,"model":1.,"gap":0.,"actual_weight":1.,"loss_contribution":loss/12 if i<12 else 0.,"scored":i<12} for i in range(13)],
                "_summary":{"status":"verified_scored_candidate","objective_canonical_sha256":OBJECTIVE,"loss":loss,"repetitions":2,"exact_loss_equality":True}}

    def test_seeds_and_exact96_unique_multivariate_proposals(self):
        seeds=load_seeds(self.plan); self.assertEqual(len(seeds),5)
        self.assertEqual({x["label"] for x in seeds},{"originalseed","prior_0","prior_1","prior_2","prior_3"})
        ps=proposals(self.plan); self.assertEqual(len(ps),96); self.assertEqual({p["scale"] for p in ps},{"near","wide"})
        self.assertTrue(all(p["proposal_seed"]==20260919 for p in ps)); self.assertEqual(len({search.fingerprint(p["parameters"]) for p in ps}),96)
        known={search.fingerprint(x["parameters"]) for x in seeds}; self.assertTrue(known.isdisjoint({search.fingerprint(p["parameters"]) for p in ps}))
        self.assertEqual({p["seed_label"] for p in ps},{"originalseed","prior_0","prior_1","prior_2","prior_3"})
        for label in {p["seed_label"] for p in ps}: self.assertEqual({p["scale"] for p in ps if p["seed_label"]==label},{"near","wide"})

    def test_fake_loop_failure_exact_repeat_and_17_plot_copy(self):
        incumbent=self.receipt(100.,self.pilot); calls=[]
        def fake(params,case,**kw):
            calls.append(kw); case.mkdir(parents=True,exist_ok=True); evaluation=case/"evaluation"; rep=2 if kw["repetitions"]==2 else 1
            diag=evaluation/"raw"/f"repetition_{rep:02d}"/"standard_diagnostics"; diag.mkdir(parents=True,exist_ok=True)
            for i in range(17): (diag/f"diagnostic_{i:02d}.png").write_bytes(b"PNG fixture")
            if kw["repetitions"]==2: out=self.receipt(88.,params); out["_evaluation"]=str(evaluation); return out
            if len(calls)==1: raise RuntimeError("numerical failure")
            out=self.receipt(90.-len(calls),params); out["_evaluation"]=str(evaluation); return out
        with tempfile.TemporaryDirectory() as d:
            out=Path(d)/"search"; result=run_search(self.plan,out,fake,incumbent=incumbent,max_proposals=2,workers=1,search_seconds=4000,verification_seconds=10)
            self.assertEqual(result["status"],"verified_selection"); self.assertTrue(any(x["status"]=="numerical_failure" for x in result["cases"]))
            self.assertEqual(calls[-1]["repetitions"],2); self.assertEqual(len(list((out/"selected_standard_diagnostics").glob("*.png"))),17)

    def test_incumbent_parameters_come_from_receipt(self):
        params=dict(self.pilot); params["theta0"]=.321; incumbent=self.receipt(100.,params); seen=[]
        def fake(p,case,**kw): seen.append(p); return self.receipt(100.,p)
        with tempfile.TemporaryDirectory() as d: result=run_search(self.plan,Path(d)/"s",fake,incumbent=incumbent,max_proposals=1,workers=1,search_seconds=4000,verification_seconds=10)
        self.assertEqual(result["selected"]["parameters"]["theta0"],.321); self.assertEqual(seen[-1]["theta0"],.321)

    def test_fatal_old_contract_error_sets_stop_and_propagates(self):
        old_error,old_stop,original=search.old_helper().ContractError,search.old_helper().STOP,search.validate_score_contract; old_stop.clear()
        try:
            search.validate_score_contract=lambda *args: (_ for _ in ()).throw(old_error("contract mismatch"))
            with tempfile.TemporaryDirectory() as d:
                with self.assertRaises(old_error): run_search(self.plan,Path(d)/"s",lambda *a,**k:None,incumbent=self.receipt(100.,self.pilot),max_proposals=1)
            self.assertTrue(old_stop.is_set())
        finally: search.validate_score_contract=original; old_stop.clear()

if __name__ == "__main__": unittest.main()
