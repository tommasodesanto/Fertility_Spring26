"""Guard the observed birth calendar, aggregation, scales and diagnostic labels."""
import copy
import unittest

from e5f_matched_pf_birth_path import compare_birth_path


class BirthPathTests(unittest.TestCase):
    def setUp(self):
        self.model = [dict(calendar_year=y, birth_children=b*.8,
            birth_children_topcode_adjusted=b, adult_population=1+k*.1)
            for k, (y,b) in enumerate(zip((2007,2011,2015,2019),(10.,9.,8.,7.)))]
        self.data = [dict(decision_year=y, birth_year_start=y+1, birth_year_end=y+4,
            live_births_total=b*1000) for y,b in zip((2007,2011,2015,2019),(10.,9.,8.,7.))]

    def test_exact_shape_and_excluded_2023(self):
        rows = self.model + [dict(calendar_year=2023,birth_children=999,
            birth_children_topcode_adjusted=999,adult_population=1.)]
        result = compare_birth_path(rows,self.data)
        self.assertEqual(result['shape_mean_squared_gap'],0)
        self.assertEqual(result['rows'][-1]['birth_year_end'],2023)
        self.assertFalse(result['female_tfr_comparable'])
        self.assertEqual(result['informative_shape_blocks'],3)

    def test_calendar_shift_rejected(self):
        self.data[0]['birth_year_start']=2007
        with self.assertRaisesRegex(ValueError,'t\\+1'):
            compare_birth_path(self.model,self.data)

    def test_normalization_does_not_hide_anchor_movement(self):
        for row in self.model:
            row['birth_children']*=2
            row['birth_children_topcode_adjusted']*=2
        result=compare_birth_path(self.model,self.data,anchor_first_block_births=10.)
        self.assertEqual(result['shape_mean_squared_gap'],0)
        self.assertEqual(result['first_block_change_from_anchor'],1.)
        self.assertGreater(result['common_anchor_mean_squared_gap'],0)

    def test_duplicate_and_missing_blocks_fail(self):
        with self.assertRaises(ValueError):
            compare_birth_path(self.model,self.data+[self.data[0]])
        with self.assertRaises(ValueError):
            compare_birth_path(self.model[:-1],self.data)

    def test_bad_top_bin_and_nonfinite_fail(self):
        for value in (0.1,float('nan')):
            rows=copy.deepcopy(self.model)
            rows[1]['birth_children_topcode_adjusted']=value
            with self.assertRaises(ValueError):
                compare_birth_path(rows,self.data)


if __name__=='__main__':
    unittest.main()
