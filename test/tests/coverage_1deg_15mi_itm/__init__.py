import subprocess
import unittest
import filecmp
import utils
import sys
import os


"""Test of a run of SPLAT's area plot around Summit Park, UT.
   This uses 1-sec arcs and the default distance."""

class TestObject1(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Set up our paths and the names of the local files we're interested in."""
        curdir = os.path.dirname(os.path.realpath(__file__))
        os.chdir(curdir)
        self.splat = utils.findsplat(curdir)
        self.srtms = os.path.join("..","..","srtm1")

        self.img = os.path.join(curdir, "coverage.png")
        self.img142 = os.path.join(curdir, "coverage_142.png")
        self.report = os.path.join(curdir, "tx-site_report.txt")
        self.report142 = os.path.join(curdir, "tx-site_report_142.txt")

    @classmethod
    def tearDownClass(self):
        #maybe don't remove these if there's an error?
        #os.remove(self.img)
        #os.remove(self.report)
        pass

    # tests start here

    def test_01_runsplat(self):
        """Run SPLAT to generate the files"""
        (imgname, suffix) = os.path.splitext(self.img)
        cmd="-hd -d %s -t tx -L 10.0 -R 15 -o %s" % (self.srtms, imgname)
        splatargs=cmd.split()
        splatargs.insert(0, self.splat)
        print(*splatargs)
        #cp = subprocess.run(splatargs, capture_output=True)
        cp = subprocess.run(splatargs)
        self.assertEqual(cp.returncode, 0)

    def test_02_compare_pngs(self):
        """Compare the area plot"""
        self.assertGreater(utils.pyssim(self.img, self.img142), 0.995, "Images don't match closely enough!")

    def test_03_compare_reports(self):
        """Compare the site_report"""
        utils.striplinefromfile(self.report, "Site Analysis Report For:")
        utils.striplinefromfile(self.report142, "Site Analysis Report For:")
        self.assertTrue(filecmp.cmp(self.report, self.report142), "tx-site_reports don't match")

