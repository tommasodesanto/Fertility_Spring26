"""H16 rough path adapter over the unchanged, authenticated Estate-A native runtime."""
from __future__ import annotations
import copy
import json
import math
from pathlib import Path
import time
import numpy as np
import two_shock as original
import two_shock_runtime as native


def prune_map_packets(folder,reason):
    """Discard only newly written diagnostic bitmaps; retain numeric root evidence."""
    folder=Path(folder)
    record_path=folder/'native_record.json'
    if not record_path.is_file():return 0
    record=json.loads(record_path.read_text())
    packets=record.get('diagnostic_packets',[])
    removed=0
    for item in packets:
        path=Path(item['path']).resolve()
        original.require(path.is_relative_to(folder.resolve()) and path.name=='diagnostic_packet.pkl.gz',
                         'Only owned native diagnostic packets may be pruned')
        if path.is_file():path.unlink();removed+=1
    if packets:
        record['diagnostic_packets']=[]
        record['diagnostic_packets_pruned']=dict(count=len(packets),reason=reason)
        original.write(record_path,record)
        if (folder/'mapping.json').is_file():original.write(folder/'mapping.json',record)
    return removed


def no_arbitrage_prices(parameters,reference_price,endpoint_price,horizon):
    original.require(type(horizon) is int and horizon in (6,16),'Rough path horizon must be 6 or 16')
    A=float(parameters.R_gross)+float(parameters.delta)+float(parameters.tau_H)
    qT=float(endpoint_price)
    original.require(math.isfinite(A) and A>1 and math.isfinite(qT) and qT>0,'Invalid no-arbitrage inputs')
    rent=float(parameters.user_cost_rate)*float(reference_price)
    original.require(math.isfinite(rent) and rent>0,'Invalid selected-reference rent')
    steady=rent/(A-1)
    t=np.arange(horizon,dtype=float)
    q=steady+(qT-steady)*np.power(A,-(horizon-t))
    original.require(np.isfinite(q).all() and (q>0).all(),'Nonfinite no-arbitrage initialization')
    return q


class RoughStageAdapter(native.StageAdapter):
    def evaluate(self,*,psi,start_year,horizon,seed,gates,budget,endpoint_controls,path_controls,deadline,folder):
        from e5f_four_shock_acceleration import extend_measured_jacobian,solve_joint_with_acceleration
        folder=Path(folder);start=self.calls;native_start=self.rt.total_native_calls
        terminal,endpoint=self._endpoint(psi,folder/'endpoint',min(deadline,time.monotonic()+budget['endpoint_seconds']))
        p=path_controls;psi_path=np.full(horizon,psi)
        J=extend_measured_jacobian(seed,horizon)
        slopes=[float(np.median(np.abs(np.diag(J)[i*horizon:(i+1)*horizon]))) for i in range(2)]
        original.require(all(math.isfinite(s) and s>0 for s in slopes),'Native measured own slopes absent')
        warm,warm_kind=self.select_warm(horizon,psi)
        kind='no_arbitrage_backward_from_fresh_endpoint' if warm is None else warm_kind
        original.write(folder/'warm_start.json',dict(kind=kind,target_psi_hex=float(psi).hex(),horizon=horizon,
            source_psi_hex=None if warm is None else warm['psi_hex'],identity=self.identity(),
            source_folder=None if warm is None else warm['source_folder'],fresh_native_mapping_required=True,
            residuals_or_fertility_reused=False))
        q=no_arbitrage_prices(self.rt.P,self.q,endpoint['price'],horizon) if warm is None else warm['prices'].copy()
        b=np.full(horizon,terminal['parameters'].pension) if warm is None else warm['fiscal_values'].copy()
        latest={};attempted=0;completed=0
        path_deadline=min(deadline,time.monotonic()+budget['path_seconds'])
        def evaluate(q,b):
            nonlocal attempted,completed
            attempted+=1;number=attempted;map_folder=folder/f'map_{number:03d}'
            mapped,record=self._mapping(terminal,endpoint,q,b,psi_path,map_folder,path_deadline)
            previous=latest.get('number')
            if previous is not None:
                prune_map_packets(folder/f'map_{previous:03d}',reason='Superseded by a completed root mapping')
            completed+=1
            latest.update(native=mapped,record=record,number=number,pin=original.pin(map_folder/'native_record.json'))
            original.write(folder/'latest_completed.json',dict(map_number=number,record=latest['pin'],
                market_maximum_residual=max(map(abs,record['market_residual'])),
                fiscal_maximum_residual=max(map(abs,record['fiscal_residual']))))
            return dict(mapping_valid=all(record['gates'].values()),market_residual=record['market_residual'],
                        fiscal_residual=record['fiscal_residual'])
        root=solve_joint_with_acceleration(closure='fixed_tax',initial_prices=q,initial_fiscal_values=b,
            evaluate=evaluate,project_prices=lambda x:np.clip(x,self.q*p['price_bound_ratios'][0],
                                                       self.q*p['price_bound_ratios'][1]),
            fiscal_bounds=[self.pension*x for x in p['pension_bound_ratios']],
            market_tolerance=gates['market_tolerance'],fiscal_tolerance=gates['fiscal_tolerance'],
            market_slope=slopes[0],fiscal_slope=slopes[1],max_log_step=p['max_log_step'],
            damping=p['damping'],max_evaluations=p['max_evaluations'],deadline_monotonic=path_deadline,
            max_condition_number=self.plan['fit']['max_condition_number'],worsening_factor=self.plan['fit']['worsening_factor'],
            final_reproduction_tolerance=gates['final_reproduction_tolerance'],initial_jacobian=J if warm is None else warm['final_jacobian'].copy(),
            callback=lambda row:original.write(folder/'root_progress.json',row))
        native.retained.original_modules()[0].inner.write(folder/'root.json',root)
        actual_calls=self.rt.total_native_calls-native_start
        completed_calls=self.calls-start
        unfinished_calls=actual_calls-completed_calls
        budget_interrupted=bool(root.get('status')=='time_or_evaluation_budget' and root.get('converged') is False)
        original.require(unfinished_calls>=0 and (unfinished_calls==0 or budget_interrupted),
                         'Rough native-call ledger differs outside unfinished bounded budget')
        original.write(folder/'root_operation_receipt.json',dict(
            root_status=root.get('status'),root_converged=bool(root.get('converged')),
            path_evaluations_attempted=attempted,path_mappings_completed=completed,
            last_successful_map_number=latest.get('number'),
            last_successful_map_pin=latest.get('pin'),
            actual_native_calls=actual_calls,completed_operation_calls=completed_calls,
            unfinished_native_calls=unfinished_calls,budget_interrupted=budget_interrupted,
            final_reproduction_max_abs=root.get('final_reproduction_max_abs'),
            missing_final_replay=root.get('final_reproduction_max_abs') is None))
        original.require(bool(latest),'Path root produced no completed native mapping')
        record=latest['record']
        terminal_check=self.rt.terminal_checks(terminal,endpoint,latest['native'],psi_path,
            tolerance=1e-3,raw_queue_tolerance=1e-3)
        if root['converged']:
            self.retain_warm(horizon,psi,root,folder)
        return dict(identity=self.identity(),reference_manifest_sha256=self.rt.identity()['reference_sha256'],
            source_pins=self.rt.identity()['source_pins'],housing=self.rt.housing,
            shock_contract=dict(start_year=self.start_year,psi=psi,expectations='permanent_until_next_surprise'),
            psi=psi,horizon=horizon,accounting_valid=bool(record['accounting_valid'] and all(record['gates'].values())),
            policy_calls=actual_calls,completed_operation_calls=completed_calls,actual_native_calls=actual_calls,
            unfinished_native_calls=unfinished_calls,root_pass=bool(root['converged']),
            stationary_pass=endpoint['stationary_pass'],terminal_pass=bool(terminal_check['all_checks_pass']),
            replay_pass=bool(root['gates']['market_replay'] and root['gates']['fiscal_replay']) if 'gates' in root else False,
            stationary_renewal_gap=endpoint['stationary_renewal_gap'],
            market_maximum_residual=max(map(abs,record['market_residual'])),
            fiscal_maximum_residual=max(map(abs,record['fiscal_residual'])),
            replay_maximum_gap=root['final_reproduction_max_abs'] if root['final_reproduction_max_abs'] is not None else float('inf'),
            terminal=terminal_check,rows=record['rows'],fertility=record['fertility'],
            final_mapping_pin=latest['pin'],path_evaluations=attempted,path_mappings_completed=completed,
            latest_completed_map_number=latest['number'],native_reply=latest['native'],
            terminal_packet=terminal,endpoint=endpoint,psi_path=psi_path,root=root,
            original_horizon_gates_tested=False)


class RoughRuntime(native.NativeRuntime):
    def prune_candidate_packets(self,folder):
        return sum(prune_map_packets(path,reason='Superseded scalar candidate')
                   for path in sorted(Path(folder).glob('map_*')))

    def measure_seed(self,*,stage,start_year,inherited_state,folder,deadline,**unused):
        adapter=RoughStageAdapter(self.rt,self.plan,inherited_state,start_year,self.initialization if stage else None)
        self.adapters[stage]=adapter
        if stage==0:
            receipt=adapter.measure_seed(**self.plan['seed'],gates=self.plan['gates'],
                budget=self.plan['budget'],deadline=deadline,folder=folder)
            receipt.update(identity=adapter.identity())
        else:
            receipt=adapter.measure_inherited_seed(Path(folder),deadline)
        receipt['source_evidence']=[original.pin(path) for path in sorted(Path(folder).glob('map_*/native_record.json'))]
        self.seeds[stage]=copy.deepcopy(receipt)
        return receipt


def build_runtime(*,plan,output,smoke=False):
    original.require(plan['smoke'] is bool(smoke),'Rough runtime smoke mode differs')
    return RoughRuntime(plan,output)
