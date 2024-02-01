import inspect
import sys
import unittest

from ost import xmlrunner
from ost.conop import GetDefaultLib


def RunTests():
  """
  This function behaves as a custom TestLoader for python unittests.

  With no system arguments, the default unittest TestRunner is used.

  If the first system argument (sys.argv[1]) is set to 'xml', a XMLTestRunner
  is used, which produces a JUnit compatible XML output file. Within the current
  module, each function is identified which is a subclass of unittest.TestCase
  and for each TestCase, a test suite is executed, producing an individual
  output file for each TestCase. The output file has the name,
  'PYTEST-<TestCaseName>.xml'.

  Example of a Python testcase:

  .. code-block:: python

    import unittest

    class TestRenumber(unittest.TestCase):

      def setUp(self):
        # prepare stuff"
        pass

      def testSomeFunction(self):
        # do some asserts
        pass

    if __name__ == "__main__":
      from ost import testutils
      testutils.RunTests()

  """
  try:
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      import __main__
      for name, obj in inspect.getmembers(__main__):
        if (isinstance(obj, type) and
                            issubclass(obj, unittest.TestCase)):
          suite = unittest.TestLoader().loadTestsFromTestCase(obj)
          stream = open('PYTEST-%s.xml' % name, 'w')
          xmlrunner.XMLTestRunner(stream).run(suite)
          stream.close()

    else:
      unittest.main()
  except Exception as e:
    print(e)


def DefaultCompoundLibIsSet():
  """
  This function checks if a default compound library is set.

  :return: True, if a compound library was found and set to be accessed with
           :func:`ost.conop.GetDefaultLib`. False otherwise.
  """
  # check if already there
  if GetDefaultLib():
    return True
  else:
    return False
