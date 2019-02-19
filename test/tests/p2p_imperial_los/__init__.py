import subprocess
import unittest
import filecmp
import utils
import sys
import os


"""Test of a normal single-threaded run of SPLAT's area plot around Summit Park, UT.
   This uses 3-sec arcs and the default distance."""

class TestObject1(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        """Set up our paths and the names of the local files we're interested in."""
        curdir = os.path.dirname(os.path.realpath(__file__))
        os.chdir(curdir)
        self.splat = utils.findsplat(curdir)
        self.srtms = os.path.join("..","..","srtm3")

        self.report = os.path.join(curdir, "tx-site_report.txt")
        self.report142 = os.path.join(curdir, "tx-site_report_142.txt")
        self.txrx = os.path.join(curdir, "tx-to-rx.txt")
        self.txrx142 = os.path.join(curdir, "tx-to-rx_142.txt")
        self.elev = os.path.join(curdir, "elevation_profile.png")
        self.elev142 = os.path.join(curdir, "elevation_profile_142.png")
        self.height = os.path.join(curdir, "height_profile.png")
        self.height142 = os.path.join(curdir, "height_profile_142.png")
        self.normheight = os.path.join(curdir, "normalized_height_profile.png")
        self.normheight142 = os.path.join(curdir, "normalized_height_profile_142.png")

    @classmethod
    def tearDownClass(self):
        #maybe don't remove these if there's an error?
        #os.remove(self.report)
        #os.remove(self.txrx)
        #os.remove(self.elev)
        #os.remove(self.height)
        #os.remove(self.normheight)
        pass

    # tests start here

    def test_01_runsplat(self):
        """Run SPLAT to generate the files"""
        cmd="-d %s -t tx -r rx -e %s -h %s -H %s" % (self.srtms, self.elev, self.height, self.normheight)
        splatargs=cmd.split()
        splatargs.insert(0, self.splat)
        print(*splatargs)
        #cp = subprocess.run(splatargs, capture_output=True)
        cp = subprocess.run(splatargs)
        self.assertEqual(cp.returncode, 0)

    def test_01_compare_reports(self):
        """Compare the tx site report"""
        utils.striplinefromfile(self.report, "SPLAT!")
        utils.striplinefromfile(self.report142, "SPLAT!")
        self.assertTrue(filecmp.cmp(self.report, self.report142), "site_reports don't match")

    def test_02_compare_txrx(self):
        """Compare the tx-to-rx"""
        utils.striplinefromfile(self.txrx, "SPLAT!")
        utils.striplinefromfile(self.txrx142, "SPLAT!")
        self.assertTrue(filecmp.cmp(self.txrx, self.txrx142), "tx-to-rx reports don't match")

    def test_04_compare_elevation_plots(self):
        """Compare the elevation plots"""
        self.assertGreater(utils.pyssim(self.elev, self.elev142), 0.94, "Images don't match closely enough!")

    def test_05_compare_height_plots(self):
        """Compare the height plots"""
        self.assertGreater(utils.pyssim(self.height, self.height142), 0.94, "Images don't match closely enough!")

    def test_06_compare_normalized_height_plots(self):
        """Compare the normalized height plots"""
        self.assertGreater(utils.pyssim(self.normheight, self.normheight142), 0.94, "Images don't match closely enough!")

