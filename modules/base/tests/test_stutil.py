import unittest
from ost import stutil

class TestStUtils(unittest.TestCase):

  def setUp(self):
    self.data1 = [0,1,2,3,4,5,6,-5,-4,-3,-2,-1,1]
    self.data2 = [1,2,3,4,5,6,7,8,9,10,11,12,13]
    self.data3 = [0.1,0.5,0.7,2.4,0.5,4.1,0.9,-1.1,-0.5,0.7,-1.4,-7.5,8.5]
    self.data_short = [0,1]
    self.data_short2 = [0,1,2]

    self.mean = 0.5384615385
    self.mean2 = 0.6076923077
    self.median = 1
    self.median2 = 0.5
    self.median_short = 0.5
    self.median_short2 = 1.0
    self.stddev = 3.3192998478
    self.stddev2 = 3.4192805223
    self.minimum = -5
    self.minimum2 = -7.5
    self.maximum = 6
    self.maximum2 = 8.5
    self.correl = -0.39639313
    self.correl2 = -0.0619291504


  def isSimilar(self, value, reference, accuracy):
    if abs(value - reference) <= accuracy:
      return True
    return False

  def testMean(self):
    assert self.isSimilar(stutil.Mean(self.data1), self.mean, 0.001), \
           "Mean (%f) does not correspond to precalculated mean (%f)" % \
           (stutil.Mean(self.data1), self.mean)
    assert self.isSimilar(stutil.Mean(self.data3), self.mean2, 0.001), \
           "Mean (%f) does not correspond to precalculated mean (%f)" % \
           (stutil.Mean(self.data3), self.mean2)

  def testMedian(self):
    assert self.isSimilar(stutil.Median(self.data1), self.median, 0.001), \
           "Median (%f) does not correspond to precalculated median (%f)" % \
           (stutil.Median(self.data1), self.median)
    assert self.isSimilar(stutil.Median(self.data3), self.median2, 0.001), \
           "Median (%f) does not correspond to precalculated median (%f)" % \
           (stutil.Median(self.data3), self.median2)
    assert self.isSimilar(stutil.Median(self.data_short), self.median_short, 0.001), \
           "Median (%f) does not correspond to precalculated median (%f)" % \
           (stutil.Median(self.data_short), self.median_short)
    assert self.isSimilar(stutil.Median(self.data_short2), self.median_short2, 0.001), \
           "Median (%f) does not correspond to precalculated median (%f)" % \
           (stutil.Median(self.data_short2), self.median_short2)

  def testStddev(self):
    assert self.isSimilar(stutil.StdDev(self.data1), self.stddev, 0.001), \
           "StdDev (%f) does not correspond to precalculated StdDev (%f)" % \
           (stutil.StdDev(self.data1), self.stddev)
    assert self.isSimilar(stutil.StdDev(self.data3), self.stddev2, 0.001), \
           "StdDev (%f) does not correspond to precalculated StdDev (%f)" % \
           (stutil.StdDev(self.data3), self.stddev2)

  def testMinimum(self):
    assert self.isSimilar(stutil.Min(self.data1), self.minimum, 0.001), \
           "Minimum (%f) does not correspond to precalculated minimum (%f)" % \
           (stutil.Min(self.data1), self.minimum)
    assert self.isSimilar(stutil.Min(self.data3), self.minimum2, 0.001), \
           "Minimum (%f) does not correspond to precalculated minimum (%f)" % \
           (stutil.Min(self.data3), self.minimum2)

  def testMaximum(self):
    assert self.isSimilar(stutil.Max(self.data1), self.maximum, 0.001), \
           "Maximum (%f) does not correspond to precalculated maximum (%f)" % \
           (stutil.Max(self.data1), self.maximum)
    assert self.isSimilar(stutil.Max(self.data3), self.maximum2, 0.001), \
           "Maximum (%f) does not correspond to precalculated maximum (%f)" % \
           (stutil.Max(self.data3), self.maximum2)

  def testCorrel(self):
    assert self.isSimilar(stutil.Correl(self.data1, self.data2), self.correl, 0.001), \
           "Correl (%f) does not correspond to precalculated correl (%f)" % \
           (stutil.Correl(self.data1, self.data2), self.correl)
    assert self.isSimilar(stutil.Correl(self.data3, self.data2), self.correl2, 0.001), \
           "Correl (%f) does not correspond to precalculated correl (%f)" % \
           (stutil.Correl(self.data3, self.data2), self.correl2)

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()