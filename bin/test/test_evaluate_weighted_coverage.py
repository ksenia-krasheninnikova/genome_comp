from unittest import TestCase

from get_coverage_stats import evaluate_weighted_coverage

class TestEvaluate_weighted_coverage(TestCase):
    def test_evaluate_weighted_coverage(self):
        header = [1]
        vals = [10]
        result= evaluate_weighted_coverage(header, vals)
        self.assertEquals(result, 1)

        header = [1,2,3]
        vals = [10,5,0]
        result= evaluate_weighted_coverage(header, vals)
        self.assertEquals(result, 1.5)


