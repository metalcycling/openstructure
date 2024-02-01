import sys
import unittest
from ost import *
from ost import settings
from ost import seq
from ost.bindings.clustalw import *
from ost.seq.alg import renumber

try:
  clustalw_path=settings.Locate(('clustalw', 'clustalw2'))
except(settings.FileNotFound):
  clustalw_path=None

class TestRenumber(unittest.TestCase):
  
  def setUp(self):
    self.target_seq = io.LoadSequence("testfiles/peptide.fasta")
    self.peptide_original = io.LoadEntity("testfiles/peptide_original.pdb")
    self.peptide_plus_5 = io.LoadEntity("testfiles/peptide_plus_5.pdb")
    self.peptide_random = io.LoadEntity("testfiles/peptide_random.pdb")
    self.peptide_del_1_2 = io.LoadEntity("testfiles/peptide_del_1_2.pdb")
    self.peptide_del_4 = io.LoadEntity("testfiles/peptide_del_4.pdb")
    self.peptide_mutation_3 = io.LoadEntity("testfiles/peptide_mutation_3.pdb")
    
  def testRenumbersChainsBasedOnSequence(self): 
    aln_seq = seq.CreateSequence('A', 'MP-T---NA')
    aln_seq.AttachView(self.peptide_original.Select(''))

    renumbered = renumber.Renumber(aln_seq)

    res = renumbered.residues
    self.assertEqual(len(res), 5)
    self.assertEqual(res[0].number, 1)
    self.assertEqual(res[1].number, 2)
    self.assertEqual(res[2].number, 4)
    self.assertEqual(res[3].number, 8)
    self.assertEqual(res[4].number, 9)

  def testRenumbersChainsBasedOnAlignment(self): 
    aln_seq = seq.CreateSequence('A', 'MP-T---NA')
    aln_seq.AttachView(self.peptide_original.Select(''))
    aln = seq.CreateAlignment(seq.CreateSequence('A', 'MP-T-XXNA'), aln_seq)

    renumbered = renumber.Renumber(aln)

    res = renumbered.residues
    self.assertEqual(len(res), 5)
    self.assertEqual(res[0].number, 1)
    self.assertEqual(res[1].number, 2)
    self.assertEqual(res[2].number, 4)
    self.assertEqual(res[3].number, 8)
    self.assertEqual(res[4].number, 9)

  def testRenumberPreservesBonds(self):
    aln_seq = seq.CreateSequence('A', 'MP-T---NA')
    aln_seq.AttachView(self.peptide_original.Select(''))

    renumbered = renumber.Renumber(aln_seq)

    self.assertTrue(mol.BondExists(renumbered.chains[0].FindAtom(1, 'N'), 
                                   renumbered.chains[0].FindAtom(1, 'CA')))

  def testClustalWPeptidePlusFive(self):
    """
    All residue numbers shifted by 5.
    Check whether internal atom order changes while renumbering
    (a new entity is generated in the edit_mode)
    TODO: add more basic tests: are all properties preserved?
    """
    if not clustalw_path:
      return
    model_seq=seq.SequenceFromChain(" ", self.peptide_plus_5.chains[0]) 
    model_seq.name="model"
    aln=ClustalW(self.target_seq,model_seq)
    aln.AttachView(1,self.peptide_plus_5.Select(""))
    ent_n=renumber.Renumber(aln)

    for r_counter in range(len(self.peptide_original.residues)):
      assert self.peptide_original.residues[r_counter].qualified_name == ent_n.residues[r_counter].qualified_name, \
              "Renumbering failed on residue level: restoring from ResNum+5"
      for atom_nr in range(len(self.peptide_original.residues[r_counter].atoms)):
        assert self.peptide_original.residues[r_counter].atoms[atom_nr].name==ent_n.residues[r_counter].atoms[atom_nr].name, \
                "Renumbering failed on atom level: restoring from ResNum+5"
  

  def testClustalWPeptideRandom(self):
    """
    Change residue names in random order
    """
    if not clustalw_path:
      return
    model_seq=seq.SequenceFromChain(" ", self.peptide_random.chains[0]) 
    model_seq.name="model"
    aln=ClustalW(self.target_seq,model_seq)
    aln.AttachView(1,self.peptide_random.Select(""))
    ent_n=renumber.Renumber(aln)

    for r_counter in range(len(self.peptide_original.residues)):
      assert self.peptide_original.residues[r_counter].qualified_name == ent_n.residues[r_counter].qualified_name, \
             "Renumbering failed on residue level: restoring from random residue numbers"
      for atom_nr in range(len(self.peptide_original.residues[r_counter].atoms)):
        assert self.peptide_original.residues[r_counter].atoms[atom_nr].name==ent_n.residues[r_counter].atoms[atom_nr].name, \
               "Renumbering failed on atom level: restoring from random residue numbers"


  def testClustalWPeptideDel_1_2(self):
    """
    First two residues were removed
    """
    if not clustalw_path:
      return
    model_seq=seq.SequenceFromChain(" ", self.peptide_del_1_2.chains[0]) 
    model_seq.name="model"
    aln=ClustalW(self.target_seq,model_seq)
    aln.AttachView(1,self.peptide_del_1_2.Select(""))
    ent_n=renumber.Renumber(aln)

    for r_counter in range(len(self.peptide_original.residues)):
      if r_counter==0 or r_counter==1:
        continue
      assert self.peptide_original.residues[r_counter].qualified_name == ent_n.residues[r_counter-2].qualified_name, \
             "Renumbering failed on residue level: restoring from random residue numbers"
      for atom_nr in range(len(self.peptide_original.residues[r_counter].atoms)):
        assert self.peptide_original.residues[r_counter].atoms[atom_nr].name==ent_n.residues[r_counter-2].atoms[atom_nr].name, \
               "Renumbering failed on atom level: restoring from random residue numbers"


  def testClustalWPeptideDel_4(self):
    """
    Residues in the middle (position 4) was removed
    """
    if not clustalw_path:
      return
    model_seq=seq.SequenceFromChain(" ", self.peptide_del_4.chains[0]) 
    model_seq.name="model"
    aln=ClustalW(self.target_seq,model_seq)
    aln.AttachView(1,self.peptide_del_4.Select(""))
    ent_n=renumber.Renumber(aln)

    flag=0
    for r_counter in range(len(self.peptide_original.residues)):
      if r_counter==3:
        flag=1
        continue
      assert self.peptide_original.residues[r_counter].qualified_name == ent_n.residues[r_counter-flag].qualified_name, \
             "Renumbering failed on residue level: restoring from random residue numbers"
      for atom_nr in range(len(self.peptide_original.residues[r_counter].atoms)):
        assert self.peptide_original.residues[r_counter].atoms[atom_nr].name==ent_n.residues[r_counter-flag].atoms[atom_nr].name, \
               "Renumbering failed on atom level: restoring from random residue numbers"


  def testClustalWPeptideMutation_3(self):
    """
    Mutation to GLY at postion 3
    """
    if not clustalw_path:
      return
    model_seq=seq.SequenceFromChain(" ", self.peptide_mutation_3.chains[0]) 
    model_seq.name="model"
    aln=ClustalW(self.target_seq,model_seq)
    aln.AttachView(1,self.peptide_mutation_3.Select(""))
    self.assertRaises(RuntimeError, renumber.Renumber, aln)



if __name__ == "__main__":
  # test renumbering
  # test if clustalw package is available on system, otherwise ignore tests
  if not clustalw_path:
    print("Could not find clustalw executable: ignoring some renumber unit tests")
  from ost import testutils
  testutils.RunTests()
