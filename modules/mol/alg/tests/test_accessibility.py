from ost import io, mol, settings
import unittest
import os
from ost.bindings import naccess
from ost.bindings import dssp


class AccessibilityContainer:
  def __init__(self):
    self.overall_acc = None
    self.per_atom_acc = list()
    self.per_atom_acc_single_chain = list()


def AccessibilitiesRaw(ent, use_naccess = False):

  acc_container = AccessibilityContainer()

  if use_naccess:
    acc_container.overall_acc = naccess.CalculateSurfaceArea(ent)
  else:
    acc_container.overall_acc = mol.alg.Accessibility(ent) 

  for a in ent.atoms:
    acc_container.per_atom_acc.append(a.GetFloatProp("asaAtom"))

  for c in ent.chains:
    chain_selection = ent.Select("cname="+c.GetName())
    if use_naccess:
      naccess.CalculateSurfaceArea(chain_selection)
    else:
      mol.alg.Accessibility(chain_selection)

  for a in ent.atoms:
    acc_container.per_atom_acc_single_chain.append(a.GetFloatProp("asaAtom"))

  return acc_container


def AccessibilitiesOligo(ent):

  acc_container = AccessibilityContainer()

  acc_container.overall_acc = mol.alg.Accessibility(ent, oligo_mode=True)
  
  for a in ent.atoms:
    acc_container.per_atom_acc.append(a.GetFloatProp("asaAtom"))
    acc = a.GetFloatProp("asaAtom_single_chain")
    acc_container.per_atom_acc_single_chain.append(acc)

  return acc_container


def Compare(acc_one, acc_two):

  # 1 digits is sufficient for summed overall accessibility. 
  # In both cases, we should have the full accessibility of the full complex
  if abs(acc_one.overall_acc - acc_two.overall_acc) > 0.1:
    return False

  # the per atom values should be a bit more accurate => 2 digits
  # let's check the per atom accessibility of the full complex
  if len(acc_one.per_atom_acc) != len(acc_two.per_atom_acc):
    return False
  for a,b in zip(acc_one.per_atom_acc, acc_two.per_atom_acc):
    if abs(a - b) > 0.01:
      return False

  # let's check the per atom accessibility when only single chains 
  # are considered
  if len(acc_one.per_atom_acc_single_chain) != \
     len(acc_two.per_atom_acc_single_chain):
    return False
  for a,b in zip(acc_one.per_atom_acc_single_chain, \
                 acc_two.per_atom_acc_single_chain):
    if abs(a - b) > 0.01:
      return False

  return True


class TestAccessibility(unittest.TestCase):
  
  def testAccNACCESS(self):

    # tests oligo mode by comparing the results from doing the
    # corresponding calculations manually
    ent_one = io.LoadPDB(os.path.join("testfiles", "1a0s.pdb"))
    ent_two = io.LoadPDB(os.path.join("testfiles", "1a0s.pdb"))

    # we're only interested in peptide stuff...
    ent_one = ent_one.Select("peptide=true")
    ent_two = ent_two.Select("peptide=true")

    acc_classic = AccessibilitiesRaw(ent_one)
    acc_oligo = AccessibilitiesOligo(ent_two)

    self.assertTrue(Compare(acc_classic, acc_oligo))

    # if there is naccess around, we also check for equality with 
    # naccess results 
    try:
      naccess_path = settings.Locate("naccess")
      ent_three = io.LoadPDB(os.path.join("testfiles", "1a0s.pdb"))
      ent_three = ent_three.Select("peptide=true")
      acc_naccess = AccessibilitiesRaw(ent_three, use_naccess=True)
      self.assertTrue(Compare(acc_classic, acc_naccess))
    except:
      print("Could not find NACCESS, could not compare Accessiblity function...")


  def testAccDSSP(self):

    # only relevant if dssp there
    try:
      # same check used in dssp binding
      dssp_path = settings.Locate(['dsspcmbi', 'dssp', 'mkdssp'],
                                  env_name='DSSP_EXECUTABLE')
    except:
      print("Could not find DSSP, could not compare Accessibility function...")
      return

    # we assume oligo mode to be working as it is tested in 
    # testAccNACCESS. So we only test the single residue
    # accessibilitities
    ent_one = io.LoadPDB(os.path.join("testfiles", "1a0s.pdb"))
    ent_two = io.LoadPDB(os.path.join("testfiles", "1a0s.pdb"))
    ent_one = ent_one.Select("peptide=true")
    ent_two = ent_two.Select("peptide=true")

    dssp.AssignDSSP(ent_one, extract_burial_status=True, dssp_bin=dssp_path)
    mol.alg.Accessibility(ent_two, algorithm=mol.alg.AccessibilityAlgorithm.DSSP)

    for a,b in zip(ent_one.residues, ent_two.residues):

      # overall accessibility
      if a.HasProp("solvent_accessibility") and b.HasProp("asaAbs"):
        diff = abs(a.GetFloatProp("solvent_accessibility") -\
                   round(b.GetFloatProp("asaAbs")))
        self.assertTrue(diff < 0.01)

      # relative accessibility
      if a.HasProp("relative_solvent_accessibility") and b.HasProp("asaRel"):
        diff = abs(a.GetFloatProp("relative_solvent_accessibility") -\
                   b.GetFloatProp("asaRel"))
        self.assertTrue(diff < 0.01)


if __name__ == "__main__":

  from ost import testutils
  testutils.RunTests()
