"""Torch-only pure tests for exact-zero menus and source overlay."""
import ast
import importlib.util
import unittest
from pathlib import Path

import numpy as np
import patch_solver
import worker
from types import SimpleNamespace


class ExactZeroChoice(unittest.TestCase):
    def choice(self, menu):
        namespace = {'np': np, 'DEAD_VALUE_CUTOFF': -1e90}
        exec(patch_solver.HELPER, namespace)
        return namespace['exact_zero_birth_choice'](np.asarray(menu, dtype=float))

    def test_rankings_tie_and_dead(self):
        menu = np.array([[[[3., 2.]]], [[[2., 3.]]], [[[2., 2.]]], [[[-1e100, -1e100]]]])
        value, probability = self.choice(menu)
        np.testing.assert_array_equal(value[:, 0, 0], [3., 3., 2., -1e100])
        np.testing.assert_array_equal(probability[:, 0, 0, :],
                                      [[1, 0], [0, 1], [1, 0], [0, 0]])

    def test_positive_source_is_preserved(self):
        root = Path(__file__).resolve().parents[4]
        source = (root/patch_solver.SOLVER).read_text()
        changed = patch_solver.patch_text(source)
        first = ('lf = Vfa / P.kappa_fert\n'
                 '                    ls, pr = logsumexp(lf, axis=3)\n'
                 '                    pr[np.max(Vfa, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0')
        later = ('l2, p2 = logsumexp(V2 / kf_cont, axis=3)\n'
                 '                            p2[np.max(V2, axis=3) <= DEAD_VALUE_CUTOFF, :] = 0.0')
        self.assertIn(first, source)
        self.assertIn(later, source)
        self.assertIn(first.replace('                    ', '                        '), changed)
        self.assertIn(later.replace('                            ', '                                '), changed)
        self.assertIn('choice_value, pr = exact_zero_birth_choice(Vfa)', changed)
        self.assertIn('continuation_value, p2 = exact_zero_birth_choice(V2)', changed)

    def test_isolated_binder_both_named_zero_restrictions(self):
        names=('kappa_fert','kappa_fert_continuation','first_birth_fixed_cost')
        bounds={'kappa_fert':(0.02,50.),'kappa_fert_continuation':(0.02,50.),
                'first_birth_fixed_cost':(0.,8.)}
        point=dict(kappa_fert=.12,kappa_fert_continuation=.22,first_birth_fixed_cost=.35)
        seen=[]
        def native_bind(value):
            seen.append(dict(value))
            for key in names:
                lo,hi=bounds[key]
                if not lo<=value[key]<=hi:raise ValueError('native bounds')
            return SimpleNamespace(**value,eps_fert=value['kappa_fert'],
                                   sequential_births=True,child_state_mode='independent_count',
                                   joint_nested_choice=False)
        for fixed,other in (('kappa_fert','kappa_fert_continuation'),
                            ('kappa_fert_continuation','kappa_fert')):
            actual=dict(point,**{fixed:0.})
            P=worker.bind_external_scale(native_bind,actual,names,bounds,fixed,0.)
            self.assertEqual(seen[-1][fixed],.02)
            self.assertEqual(getattr(P,fixed),0.)
            self.assertEqual(getattr(P,other),point[other])
            self.assertEqual(P.eps_fert,P.kappa_fert)
            bad=dict(actual,first_birth_fixed_cost=9.)
            with self.assertRaises(RuntimeError):
                worker.bind_external_scale(native_bind,bad,names,bounds,fixed,0.)
            with self.assertRaises(RuntimeError):
                worker.bind_external_scale(native_bind,dict(actual,extra=1),names,bounds,fixed,0.)

    def test_positive_binary_menu_values_and_probabilities_exact(self):
        root=Path(__file__).resolve().parents[4]
        source=(root/patch_solver.SOLVER).read_text()
        changed=patch_solver.patch_text(source)
        spec=importlib.util.spec_from_file_location('isolated_native_utils',
            root/'code/model/intergen_eqscale_seq_optimized/utils.py')
        native=importlib.util.module_from_spec(spec);spec.loader.exec_module(native)

        def statements(text,first_name,count):
            # Select the actual Bellman block by its assignment, retaining its
            # own RHS, native logsumexp call, and dead-row masking statement.
            matches=[]
            for node in ast.walk(ast.parse(text)):
                for _,body in ast.iter_fields(node):
                    if not isinstance(body,list):continue
                    for i,statement in enumerate(body):
                        target=statement.targets[0] if isinstance(statement,ast.Assign) and len(statement.targets)==1 else None
                        names=([target.id] if isinstance(target,ast.Name) else
                               [item.id for item in target.elts if isinstance(item,ast.Name)]
                               if isinstance(target,ast.Tuple) else [])
                        if first_name in names:
                            if first_name=='lf':
                                if count==3:
                                    following=body[i+count] if i+count<len(body) else None
                                    # Original active branch saves its two-action
                                    # menu; legacy nonsequential saves every slot.
                                    if (following is None or 'fert_probs' not in ast.unparse(following)
                                            or ':2]' not in ast.unparse(following)):
                                        continue
                                else:
                                    # In the patched active branch the save sits
                                    # outside the positive-scale else body.
                                    if (len(body[i:i+count])!=4 or
                                            'choice_value = P.kappa_fert * ls' not in
                                            ast.unparse(body[i+3])):
                                        continue
                            matches.append(body[i:i+count])
            self.assertEqual(len(matches),1)
            return matches[0]

        def execute(nodes,menu,kappa,branch):
            block=ast.fix_missing_locations(ast.Module(body=nodes,type_ignores=[]))
            namespace=dict(np=np,logsumexp=native.logsumexp,DEAD_VALUE_CUTOFF=-1e90,
                           P=SimpleNamespace(kappa_fert=kappa),kf_cont=kappa,
                           Vfa=menu.copy(),V2=menu.copy())
            exec(compile(block,'<actual-bellman-choice-block>','exec'),namespace)
            if branch=='first':
                return (namespace.get('choice_value',kappa*namespace['ls']),namespace['pr'])
            return (namespace.get('continuation_value',kappa*namespace['l2']),namespace['p2'])

        rng=np.random.default_rng(20260929)
        menu=rng.normal(size=(3,2,2,2))
        menu[0,0,0,:]=-1e100
        menu[1,0,0,:]=2.0
        for branch,kappa,first_name,old_count,new_count in (
                ('first',.10854911052270635,'lf',3,4),
                ('later',.22239802621622237,'l2',2,3)):
            old_nodes=statements(source,first_name,old_count)
            new_nodes=statements(changed,first_name,new_count)
            old_value,old_probability=execute(old_nodes,menu,kappa,branch)
            new_value,new_probability=execute(new_nodes,menu,kappa,branch)
            self.assertTrue(np.array_equal(old_value,new_value),branch)
            self.assertTrue(np.array_equal(old_probability,new_probability),branch)


if __name__ == '__main__':
    unittest.main()
