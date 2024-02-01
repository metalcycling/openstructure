import unittest
import ost

# helper log sink to campture messages
class _CapturingLogSink(ost.LogSink):
  def __init__(self):
    ost.LogSink.__init__(self)
  def LogMessage(self, message, severity):
    self.message = message
    self.severity = severity

# Altough the logging system might appear to be too simple to be worth writing a 
# specific test case for, it actually isn't. The python export is very fragile 
# and seemingly trivial changes can break the code in unexpected ways. So let's 
# check for some invariants
class TestLog(unittest.TestCase):

  def testGetLogSink(self):
    logsink = ost.GetCurrentLogSink()
    self.assertTrue(hasattr(logsink, 'LogMessage'))
    # Check if the return type of logsink is sane
    ost.PushLogSink(ost.GetCurrentLogSink())

  def testPushPopLogSink(self):
    class MyLogSink(ost.LogSink):
       def __init__(self):
         ost.LogSink.__init__(self)
    ls = MyLogSink()
    ost.PushLogSink(ls)
    self.assertEqual(ls, ost.GetCurrentLogSink())
    ost.PopLogSink()
    self.assertNotEqual(ls, ost.GetCurrentLogSink())

  def testLogMessage(self):
    ls = _CapturingLogSink()
    ost.PushVerbosityLevel(1)
    ost.PushLogSink(ls)
    ost.LogError('error message')
    self.assertEqual(ls.message, 'error message\n')
    self.assertEqual(ls.severity, 0)
    ost.LogWarning(1, 2, 3)
    self.assertEqual(ls.message, '1 2 3\n')
    self.assertEqual(ls.severity, 1)
    ost.PopLogSink()

  def testLogMultipleMessages(self):
    # observation: converting non-strings in logging can break following calls
    ls = _CapturingLogSink()
    ost.PushVerbosityLevel(1)
    ost.PushLogSink(ls)
    ost.LogWarning('foo')
    self.assertEqual(ls.message, 'foo\n')
    ost.LogWarning(1)
    self.assertEqual(ls.message, '1\n')
    ost.LogWarning('bar')
    self.assertEqual(ls.message, 'bar\n')
    ost.PopLogSink()

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()