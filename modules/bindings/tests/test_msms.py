import unittest
from ost import *
from ost import settings
from ost.bindings import msms

class TestMSMSBindings(unittest.TestCase):
  
  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb")

    self.num_tri=111546
    self.num_vert=55762
    self.ases=20051.928
    self.asas=21574.324
    self.volume=60229.848

  def testCalculateSurface(self):
    surf=msms.CalculateSurface(self.protein,
                               msms_env='MSMSSERVER',
                               density=3,
                               radius=1.4,
                               attach_asa="asa",
                               attach_esa="esa")[0]
    self.assertEqual(self.num_vert, len(surf.GetVertexIDList()))
    self.assertEqual(self.num_tri, len(surf.GetTriIDList()))
 
  def testCalculateSurfaceArea(self):
    (msms_ases, msms_asas)=msms.CalculateSurfaceArea(self.protein, \
                                                     msms_env='MSMSSERVER',
                                                     density=3,
                                                     radius=1.4)
    self.assertAlmostEqual(self.ases, msms_ases[0], 1)
    self.assertAlmostEqual(self.asas, msms_asas[0], 1)

  def testCalculateSurfaceVolume(self):
    volume=msms.CalculateSurfaceVolume(self.protein, \
                                       msms_env='MSMSSERVER',
                                       density=3,
                                       radius=1.4)
    self.assertAlmostEqual(self.volume, volume, 1)

if __name__ == "__main__":
  # test if msms package is available on system, otherwise ignore tests
  VERSION_REQUIRED = '2.6.1'
  try:
    msms._GetExecutable(msms_exe=None, msms_env='MSMSSERVER')
  except(settings.FileNotFound):
    print("Could not find msms executable: ignoring unit tests")
    exit(0)
  version = msms.GetVersion(msms_exe=None, msms_env='MSMSSERVER')
  if version!=VERSION_REQUIRED:
    print("MSMS version (%s) does not match required version %s: ignoring unit tests"%(version, VERSION_REQUIRED))
  from ost import testutils
  testutils.RunTests()
