from unittest import TestCase

from maf_extractor import process

class TestProcess(TestCase):

    def test_process(self):
        process('resources/test_process.bed','resources/test_process.maf')