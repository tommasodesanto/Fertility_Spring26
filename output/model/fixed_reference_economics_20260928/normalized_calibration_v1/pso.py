"""Bounded eight-particle search; the caller owns the unchanged native objective."""
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np


def particle_swarm(objective, lower, spans, bounds, coordinates, config, chain, out, maxeval):
    design=next(s for s in config['pso_swarms'] if s['chain']==chain)
    lo=np.asarray([config['pso_initialization_ranges'][k][0] for k in coordinates])
    scale=np.asarray([config['pso_initialization_ranges'][k][1]-config['pso_initialization_ranges'][k][0] for k in coordinates])
    hardlo=(lower-lo)/scale; hardhi=(lower+spans-lo)/scale
    positions=(np.asarray([[p[k] for k in coordinates] for p in design['particles']])-lo)/scale
    velocity=np.zeros_like(positions); personal=positions.copy(); scores=np.full(8,np.inf)
    rng=np.random.default_rng(design['seed']); global_best=None; global_score=float('inf'); calls=0;generation=0
    settings=config['pso']
    def checkpoint(status,particle):
        state=dict(status=status,chain=chain,seed=design['seed'],generation=generation,particle=particle,objective_calls=calls,positions_initial_range_scaled=positions.tolist(),velocities_initial_range_scaled=velocity.tolist(),personal_best_positions=personal.tolist(),personal_best_losses=[float(x) if np.isfinite(x) else None for x in scores],global_best_position=None if global_best is None else global_best.tolist(),global_best_loss=global_score if np.isfinite(global_score) else None,rng_state=rng.bit_generator.state,original_hard_bounds=bounds,initialization_ranges=config['pso_initialization_ranges'],settings=settings,no_auto_resume=True)
        target=Path(out)/'pso_state.json';tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2,allow_nan=False)+'\n');tmp.replace(target)
    while calls<maxeval:
        for i in range(8):
            if calls>=maxeval:break
            checkpoint('before_objective',i)
            physical=lo+scale*positions[i]
            try:
                value=float(objective((physical-lower)/spans))
            except BaseException:
                checkpoint('stopped_or_failed_objective',i)
                raise
            calls+=1
            # The native numerical-rejection penalty must not establish a personal/global best.
            if value<1e12 and np.isfinite(value) and value<scores[i]:
                scores[i]=value;personal[i]=positions[i].copy()
                if value<global_score:global_score=value;global_best=positions[i].copy()
            checkpoint('completed_objective',i)
        if calls>=maxeval:break
        # A generation with no admissible particle explores by bounded random movement.
        attractor=positions if global_best is None else global_best[None,:]
        r1=rng.random(positions.shape);r2=rng.random(positions.shape)
        velocity=settings['inertia']*velocity+settings['cognitive']*r1*(personal-positions)+settings['social']*r2*(attractor-positions)
        if global_best is None:velocity=rng.uniform(-.25,.25,positions.shape)
        velocity=np.clip(velocity,-settings['velocity_cap_initial_span'],settings['velocity_cap_initial_span'])
        proposed=positions+velocity; repaired=np.clip(proposed,hardlo,hardhi)
        velocity[proposed!=repaired]=0.;positions=repaired;generation+=1
        checkpoint('velocity_update',None)
    checkpoint('objective_call_cap',None)
    return SimpleNamespace(success=False,message='particle swarm objective-call cap',nfev=calls)
