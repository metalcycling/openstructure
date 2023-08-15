import unittest
from ost import *
from ost import settings
from ost.bindings import cadscore
from ost import testutils

class TestCADBindings(unittest.TestCase):
  
  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb")


  def testCADClassic(self):

    try:
      # all of the following need to be present
      cad_calc_path = settings.Locate('CADscore_calc.bash')  
      cad_read_g_path = settings.Locate('CADscore_read_global_scores.bash')  
      cad_read_l_path = settings.Locate('CADscore_read_local_scores.bash')
      executable_path = settings.Locate('voroprot2')
    except settings.FileNotFound:
      self.skipTest("Could not find CAD score classic executables: ignoring unit tests")

    cad_result = cadscore.CADScore(self.protein, self.protein, 
                                   label="cad_classic")

    # model and reference are the same, we expect a global CAD score of 1
    self.assertEqual(cad_result.globalAA, 1.0)

    # one score per residue
    self.assertEqual(len(cad_result.localAA), len(self.protein.residues))

    # model and reference are the same, we expect local CAD scores of 0.0
    for score in cad_result.localAA.values():
      self.assertEqual(score, 0.0)

    # check whether this score is assigned to each residue as float property
    for r in self.protein.residues:
      self.assertTrue(r.HasProp("cad_classic"))
      self.assertEqual(r.GetFloatProp("cad_classic"), 0.0)


  def testCADVoronota(self):

    try:
      # all of the following need to be present
      voronota_cadscore_path = settings.Locate("voronota-cadscore")
      executable_path = settings.Locate("voronota")   
    except settings.FileNotFound:
      self.skipTest("Could not find CAD score voronota executables: ignoring unit tests")

    cad_result = cadscore.CADScore(self.protein, self.protein, mode="voronota",
                                   label="cad_voronota")

    # model and reference are the same, we expect a global CAD score of 1
    self.assertEqual(cad_result.globalAA, 1.0)

    # one score per residue
    self.assertEqual(len(cad_result.localAA), len(self.protein.residues))

    # model and reference are the same, we expect local CAD scores of 1.0
    for score in cad_result.localAA.values():
      self.assertEqual(score, 1.0)

    # check whether this score is assigned to each residue as float property
    for r in self.protein.residues:
      self.assertTrue(r.HasProp("cad_voronota"))
      self.assertEqual(r.GetFloatProp("cad_voronota"), 1.0)
 

if __name__ == "__main__":
  testutils.RunTests()
