import unittest
from ost import *
from ost import settings
from ost import seq


def fixture():
  e=mol.CreateEntity()
  ede=e.EditXCS()
  chain=ede.InsertChain('A')
  for res in 'ABCDEFGH':
    r=ede.AppendResidue(chain, res)
    r.SetOneLetterCode(res)
    ede.InsertAtom(r, "XXX", geom.Vec3())
    
  return e
  
class TestSeq(unittest.TestCase):
  
  def setUp(self):
    self.ent=fixture()
    
  def testSeqSlice(self):
    seq_string='abcdefg'
    s=seq.CreateSequence('A', seq_string)
    self.assertEqual(s[1:5], seq_string[1:5])
    self.assertEqual(s[:-1], seq_string[:-1])
    self.assertEqual(s[-3:-2], seq_string[-3:-2])
    self.assertEqual(s[-3:], seq_string[-3:])
    self.assertEqual(s[3:4], seq_string[3:4])
    
  def testSeqListSlice(self):
    a=seq.CreateSequence('A','aaaa')
    b=seq.CreateSequence('B','bbbb')
    c=seq.CreateSequence('C','cccc')
    d=seq.CreateSequence('D','dddd')
    sl=seq.CreateSequenceList(a, b, c, d)
    sliced=sl[1:]
    self.assertEqual(len(sliced), 3)
    self.assertEqual(str(sliced[0]), str(b))
    self.assertEqual(str(sliced[1]), str(c))
    self.assertEqual(str(sliced[2]), str(d))
    sliced=sl[:-1]
    self.assertEqual(len(sliced), 3)
    self.assertEqual(str(sliced[0]), str(a))
    self.assertEqual(str(sliced[1]), str(b))
    self.assertEqual(str(sliced[2]), str(c))
    sliced=sl[-1:]
    self.assertEqual(len(sliced), 1)
    self.assertEqual(str(sliced[0]), str(d))
    
  def testAlnSlice(self):
    a=seq.CreateSequence('A','abcd')
    b=seq.CreateSequence('B','efgh')
    c=seq.CreateSequence('C','ijkl')
    d=seq.CreateSequence('D','mnop')
    aln=seq.CreateAlignment(a, b, c, d)
    sliced=aln[1:]
    self.assertEqual(len(sliced), 3)
    self.assertEqual(str(sliced[0]), 'bfjn')
    self.assertEqual(str(sliced[1]), 'cgko')
    self.assertEqual(str(sliced[2]), 'dhlp')
    sliced=aln[:-1]
    self.assertEqual(len(sliced), 3)
    self.assertEqual(str(sliced[0]), 'aeim')
    self.assertEqual(str(sliced[1]), 'bfjn')
    self.assertEqual(str(sliced[2]), 'cgko')
    sliced=aln[-1:]
    self.assertEqual(len(sliced), 1)
    self.assertEqual(str(sliced[0]), 'dhlp')

  def testViewsFromSequences_01(self):
    seq_a=seq.CreateSequence("A", "ABCD-FGH")
    seq_a.AttachView(self.ent.Select('rname=A,B,C,D,F,G,H'))
    seq_b=seq.CreateSequence("B", "ABCDEFGH")
    seq_b.AttachView(self.ent.Select(''))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'ABCDFGH')
    self.assertEqual(string_b, 'ABCDFGH')
  def testSeqIterBZDNG148(self):
    s=seq.CreateSequence('A', 'abcdef')
    for x in s: 
      pass
  def testViewsFromSequences_02(self):
    seq_a=seq.CreateSequence("A", "ABCD-FGH")
    seq_a.AttachView(self.ent.Select('rname=A,B,C,D,F,G,H'))
    seq_b=seq.CreateSequence("B", "ABCD-FGH")
    seq_b.AttachView(self.ent.Select('rname=A,B,C,D,F,G,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'ABCDFGH')
    self.assertEqual(string_b, 'ABCDFGH')
    
  def testViewsFromSequences_03(self):
    seq_a=seq.CreateSequence("A", "ABCD--GH")
    seq_a.AttachView(self.ent.Select('rname=A,B,C,D,G,H'))
    seq_b=seq.CreateSequence("B", "ABCD-FGH")
    seq_b.AttachView(self.ent.Select('rname=A,B,C,D,F,G,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'ABCDGH')
    self.assertEqual(string_b, 'ABCDGH')
 
  def testViewsFromSequences_04(self):
    seq_a=seq.CreateSequence("A", "ABCD-FGH")
    seq_a.AttachView(self.ent.Select('rname=A,B,C,D,F,G,H'))
    seq_b=seq.CreateSequence("B", "ABCD--GH")
    seq_b.AttachView(self.ent.Select('rname=A,B,C,D,G,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'ABCDGH')
    self.assertEqual(string_b, 'ABCDGH')
    
  def testViewsFromSequences_05(self):
    seq_a=seq.CreateSequence("A", "ABCD-F--")
    seq_a.AttachView(self.ent.Select('rname=A,B,C,D,F'))
    seq_b=seq.CreateSequence("B", "ABCDEFGH")
    seq_b.AttachView(self.ent.Select('rname=A,B,C,D,E,F,G,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'ABCDF')
    self.assertEqual(string_b, 'ABCDF')
    
  def testViewsFromSequences_06(self):
    seq_a=seq.CreateSequence("A", "--CD-FGH")
    seq_a.AttachView(self.ent.Select('rname=C,D,F,G,H'))
    seq_b=seq.CreateSequence("B", "ABCDEFGH")
    seq_b.AttachView(self.ent.Select('rname=A,B,C,D,E,F,G,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'CDFGH')
    self.assertEqual(string_b, 'CDFGH')


  def testViewsFromSequences_07(self):
    seq_a=seq.CreateSequence("A", "AB-D-FGH")
    seq_a.AttachView(self.ent.Select('rname=A,B,D,F,G,H'))
    seq_b=seq.CreateSequence("B", "AB-DEF-H")
    seq_b.AttachView(self.ent.Select('rname=A,B,D,E,F,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'ABDFH')
    self.assertEqual(string_b, 'ABDFH')


  def testViewsFromSequences_08(self):
    seq_a=seq.CreateSequence("A", "A-C-E-G")
    seq_a.AttachView(self.ent.Select('rname=A,C,E,G'))
    seq_b=seq.CreateSequence("B", "-B-D-H-")
    seq_b.AttachView(self.ent.Select('rname=B,D,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, '')
    self.assertEqual(string_b, '')

  def testViewsFromSequences_09(self):
    seq_a=seq.CreateSequence("A", "B-D-FGH")
    seq_a.AttachView(self.ent.Select('rname=A,B,D,F,G,H'))
    seq_a.offset=1
    seq_b=seq.CreateSequence("B", "B-DEF-H")
    seq_b.offset=1
    seq_b.AttachView(self.ent.Select('rname=A,B,D,E,F,H'))
    a, b=seq.ViewsFromSequences(seq_a, seq_b)
    string_a=''.join([r.one_letter_code for r in a.residues])
    string_b=''.join([r.one_letter_code for r in b.residues])
    self.assertEqual(string_a, 'BDFH')
    self.assertEqual(string_b, 'BDFH')

  def testNormalise(self):
    seq_a=seq.CreateSequence("A", "B-D-FGH")
    self.assertEqual("B-D-FGH", seq_a.GetString())
    seq_a.Normalise()
    self.assertEqual("BDFGH", seq_a.GetString())
    seq_a=seq.CreateSequence("A", "b.d-fgh")
    self.assertEqual("b.d-fgh", seq_a.GetString())
    seq_a.Normalise()
    self.assertEqual("BDFGH", seq_a.GetString())

if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()

