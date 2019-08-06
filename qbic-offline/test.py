import unittest

import os

import qbic

class TestQbicOffline(unittest.TestCase):

    def test_inittbl(self):
        filename = "testing_resources/test.vcf"
        chromosome_version = "hg19"

        main.inittbl()

if __name__ == '__main__':
    unittest.main()
