import unittest
from ost import conop, mol


class TestProcessor(unittest.TestCase):
  def testPyWrap(self):
    class MyProc(conop.Processor):
      def __init__(self):
        conop.Processor.__init__(self)
        self.count =0
      def DoProcess(self, diag, ent):
        self.count+=1
    p = MyProc()
    ent = mol.CreateEntity()
    p.Process(ent)
    self.assertEqual(p.count, 1)

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()


