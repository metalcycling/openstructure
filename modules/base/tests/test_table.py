'''
Unit tests for Table class

Author: Tobias Schmidt
'''

import os
import unittest
import glob
from ost.table import *
import ost

HAS_NUMPY=True
HAS_SCIPY_STATS=True
HAS_SCIPY_NDIMG=True
HAS_MPL=True
HAS_PIL=True
try:
  import numpy as np
except ImportError:
  HAS_NUMPY=False
  print("Could not find numpy: ignoring some table class unit tests")

try:
  import scipy.stats.mstats
except ImportError:
  HAS_SCIPY_STATS=False
  print("Could not find scipy.stats.mstats: ignoring some table class unit tests")
  
try:
  import scipy.ndimage
except ImportError:
  HAS_SCIPY_NDIMG=False
  print("Could not find scipy.ndimage: ignoring some table class unit tests")
  

try:
  import matplotlib
  matplotlib.use('Agg')
except ImportError:
  HAS_MPL=False
  print("Could not find matplotlib: ignoring some table class unit tests")
  
try:
  from PIL import Image
  from PIL import ImageChops
except ImportError:
  HAS_PIL=False
  print("Could not find python imagine library: ignoring some table class unit tests")

# setting up an OST LogSink to capture error messages
class _FetchLog(ost.LogSink):
  def __init__(self):
    ost.LogSink.__init__(self)
    self.messages = dict()

  def LogMessage(self, message, severity):
    levels=['ERROR', 'WARNING', 'INFO', 'VERBOSE', 'DEBUG', 'TRACE']
    level=levels[severity]
    if not level in list(self.messages.keys()):
      self.messages[level] = list()
    self.messages[level].append(message.strip())
  
class TestTable(unittest.TestCase):

  def tearDown(self):
    for filename in glob.glob('*_out.*'):
      os.remove(filename)

  def setUp(self):
    ost.PushVerbosityLevel(3)

  def CreateTestTable(self):
    '''
    creates a table with some test data
    
      first  second  third 
    ----------------------
     x            3     NA
     foo         NA  2.200
     NA           9  3.300

    '''
    tab = Table()
    tab.AddCol('first', 'string')
    tab.AddCol('second', 'int')
    tab.AddCol('third', 'float')
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 0)
    self.CompareColTypes(tab, ['first','second', 'third'], 'sif')
    tab.AddRow(['x',3, None], overwrite=None)
    tab.AddRow(['foo',None, 2.2], overwrite=None)
    tab.AddRow([None,9, 3.3], overwrite=None)
    return tab

  def CompareRowCount(self, t, row_count):
    '''
    Compare the number of rows
    '''
    self.assertEqual(len(t.rows),
                     row_count,
                     "row count (%i) different from expected value (%i)" \
                     %(len(t.rows), row_count))
  
  def CompareColCount(self, t, col_count):
    '''
    Compare the number of columns
    '''
    self.assertEqual(len(t.col_names),
                     col_count,
                     "column count (%i) different from expected value (%i)" \
                     %(len(t.col_names), col_count))
  
  def CompareColNames(self, t, col_names):
    '''
    Compare all column names of the table with a list of reference col names
    '''
    self.CompareColCount(t, len(col_names))
    for i, (col_name, ref_name) in enumerate(zip(t.col_names, col_names)):
        self.assertEqual(col_name,
                         ref_name,
                         "column name (%s) different from expected name (%s) at col %i" \
                         %(col_name, ref_name, i))
  
  def CompareDataFromDict(self, t, data_dict):
    '''
    Compare all values of a table with reference values given in the form of a
    dictionary containing a list of values for each column.
    '''
    self.CompareColCount(t, len(data_dict))
    for k, v in data_dict.items():
      self.CompareDataForCol(t, k, v)
      
  def CompareDataForCol(self, t, col_name, ref_data):
    '''
    Compare the values of each row of ONE column specified by its name with
    the reference values specified as a list of values for this column.
    '''
    self.CompareRowCount(t, len(ref_data))
    idx = t.GetColIndex(col_name)
    col_type = t.col_types[idx]
    for i, (row, ref) in enumerate(zip(t.rows, ref_data)):
      if (isinstance(ref, float) or isinstance(ref, int)) and (isinstance(row[idx], float) or isinstance(row[idx], int)):
        self.assertAlmostEqual(row[idx],
                               ref,
                               msg="data (%s) in col (%s), row (%i) different from expected value (%s)" \
                               %(row[idx], col_name, i, ref))
      else:
        self.assertEqual(row[idx],
                         ref,
                         "data (%s) in col (%s), row (%i) different from expected value (%s)" \
                         %(row[idx], col_name, i, ref))

  def CompareColTypes(self, t, col_names, ref_types):
    '''
    Compare the types of n columns specified by their names with reference
    values specified either as a string consisting of the short type names 
    (e.g 'sfb') or a list of strings consisting of the long type names
    (e.g. ['string','float','bool'])
    '''
    if type(ref_types)==str:
      trans = {'s' : 'string', 'i': 'int', 'b' : 'bool', 'f' : 'float'}
      ref_types = [trans[rt] for rt in ref_types]
    if type(col_names)==str:
      col_names = [col_names]
    self.assertEqual(len(col_names),
                     len(ref_types),
                     "number of col names (%i) different from number of reference col types (%i)" \
                     %(len(col_names), len(ref_types)))
    idxs = [t.GetColIndex(x) for x in col_names]
    for idx, ref_type in zip(idxs, ref_types):
      self.assertEqual(t.col_types[idx],
                       ref_type,
                       "column type (%s) at column %i, different from reference col type (%s)" \
                       %(t.col_types[idx], idx, ref_type))

  def CompareImages(self, img1, img2):
    '''
    Compares two images based on all pixel values. This function needs the
    python imaging library (PIL) package.
    '''
    if not HAS_PIL:
      return
    diff = ImageChops.difference(img1, img2)
    self.assertEqual(diff.getbbox(),None)

  def testAllowsToSearchColNames(self):
    tab = self.CreateTestTable()
    self.assertEqual(tab.SearchColNames('d$'), ['second', 'third'])
    self.assertEqual(tab.SearchColNames('(first|third)'), ['first','third'])

  def testProvidesDirectAccessToColumns(self):
    tab = Table(['x', 'two'], 'ii')
    tab.AddRow([1,2])
    self.assertEqual([1], list(tab.x))
    self.assertEqual([2], list(tab.two))
    self.assertRaises(AttributeError, tab.__getattr__, 'one')
  def testZip(self):
    tab=Table(['col1', 'col2', 'col3', 'col4'], 'sssi')
    tab.AddRow(['a', 'b', 'c', 1])
    tab.AddRow(['x', 'y', 'z', 2])
    z=tab.Zip('col1', 'col3')
    self.assertEqual(len(z), 2)
    self.assertEqual(z[0][0], 'a')
    self.assertEqual(z[0][1], 'c')
    self.assertEqual(z[1][0], 'x')
    self.assertEqual(z[1][1], 'z')
    z=tab.Zip('col3', 'col1')
    self.assertEqual(len(z), 2)
    self.assertEqual(z[0][0], 'c')
    self.assertEqual(z[0][1], 'a')
    self.assertEqual(z[1][0], 'z')
    self.assertEqual(z[1][1], 'x')
    z=tab.Zip('col1', 'col1')
    self.assertEqual(len(z), 2)
    self.assertEqual(z[0][0], 'a')
    self.assertEqual(z[0][1], 'a')
    self.assertEqual(z[1][0], 'x')
    self.assertEqual(z[1][1], 'x')
    z=tab.Zip('col1', 'col4')
    self.assertEqual(type(z[0][0]),str)
    self.assertEqual(type(z[1][0]),str)
    self.assertEqual(type(z[0][1]),int)
    self.assertEqual(type(z[1][1]),int)
    self.assertRaises(ValueError, tab.Zip, 'col5', 'col3')
  def testPercentiles(self):
    tab = Table(['nums'], 'i')
    self.assertEqual(tab.Percentiles('nums', [0,100]), [None, None])
    self.assertRaises(ValueError, tab.Percentiles, 'nums', [101])
    self.assertRaises(ValueError, tab.Percentiles, 'nums', [-1])
    for i in (35,15,50,40,20):
      tab.AddRow([i])
    self.assertEqual(tab.Percentiles('nums', [0,30,40,100]), [15,20,35,50])
  def testTableInitEmpty(self):
    '''
    empty table
    '''
    tab = Table()
    self.CompareColCount(tab, 0)
    self.CompareRowCount(tab, 0)
    self.assertRaises(ValueError, tab.GetColIndex, 'a')
  def testTableInitFromOldTable(self):
    old_tab = Table(['a','b'],'ff')
    new_tab = Table(old_tab.col_names,old_tab.col_types)
    new_tab.AddCol('c','f')
    self.assertEqual(old_tab.col_names, ['a','b'])
  def testGuessColumnType(self):
    self.assertEqual(GuessColumnType(['1', '1.3', '2']), 'float')
    self.assertEqual(GuessColumnType(['1', '1', '2']), 'int')
    self.assertEqual(GuessColumnType(['NONE', '1', '1', '2']), 'int')
    self.assertEqual(GuessColumnType(['NONE', '1', '1', '2']), 'int')
    self.assertEqual(GuessColumnType(['NONE', '1', '1', 'a']), 'string')
    self.assertEqual(GuessColumnType(['NONE', 'TRUE', 'False']), 'bool')
    self.assertEqual(GuessColumnType(['NONE']), 'string')
  def testTableInitSingleColEmpty(self):
    '''
    empty table with one float column:
    
     x 
    ---
    
    '''
    tab = Table(['x'], 'f')
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 0)
    self.CompareColNames(tab, ['x'])
    self.CompareColTypes(tab, 'x', 'f')
    
  def testTableInitMultiColEmpty(self):
    '''
    empty table with multiple column with different types:
    
     x  y  z  a 
    ------------
    
    '''
    tab = Table(['x','y','z','a'], 'sfbi')
    self.CompareColCount(tab, 4)
    self.CompareRowCount(tab, 0)
    self.CompareColNames(tab, ['x','y','z','a'])
    self.CompareColTypes(tab, ['x','y','z','a'], 'sfbi')
    self.CompareColTypes(tab, ['x','y','z','a'], ['string','float','bool','int'])
      
  def testTableInitSingleColSingleValueNonEmpty(self):
    '''
    table with one column and one row:
    
       x   
    -------
      5.000

    '''
    tab = Table(['x'], 'f', x=5)
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 1)
    self.CompareColNames(tab, ['x'])
    self.CompareColTypes(tab, 'x', 'f')
    
  def testTableInitMultiColSingleValueNonEmpty(self):
    '''
    table with three columns and one row:
    
       x      a      z   
    ---------------------
      5.000 False   1.425
    
    '''
    tab = Table(['x','a','z'], 'fbf', x=5, z=1.425, a=False)
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 1)
    self.CompareColNames(tab, ['x','a','z'])
    self.CompareColTypes(tab, ['z','x','a'], 'ffb')
    self.CompareDataFromDict(tab, {'x': [5], 'z': [1.425], 'a': [False]})
    
  def testTableInitMultiColSingleValueAndNoneNonEmpty(self):
    '''
    table with three columns and one row with two None values:
    
       x    a1  zzz 
    ----------------
      5.000 NA   NA 
    '''
    tab = Table(['x','a1','zzz'], 'fbf', x=5)
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 1)
    self.CompareColNames(tab, ['x','a1','zzz'])
    self.CompareColTypes(tab, ['zzz','x','a1'], 'ffb')
    self.CompareDataFromDict(tab, {'x': [5], 'zzz': [None], 'a1': [None]})
  
  def testTableInitSingleColMultiValueNonEmpty(self):
    '''
    table with one column and five row:
    
       x   
    -------
      0.000
      1.000
      2.000
      3.000
      4.000

    '''
    tab = Table(['x'], 'f', x=list(range(5)))
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 5)
    self.CompareColNames(tab, ['x'])
    self.CompareColTypes(tab, 'x', 'f')
    
  def testTableFilterColNamesTypes(self):
    """
    make sure the col_names and col_types are copied. 
    We don't want them to be referenced to the original table.
    This leads to strange surprises.
    """
    t = Table(['a', 'b'], 'ii')
    t.AddRow([1,2])
    t.AddRow([2,3])
    t.AddRow([2,3])
    t.AddRow([3,3])
    t.AddRow([4,3])
    t.AddRow([5,3])
    filt = t.Filter(a=2)
    filt.AddCol('c', 'i')
    self.assertEqual(len(t.col_names), 2)
    self.assertEqual(len(t.col_types), 2)

  def testTableInitMultiColMultiValueNonEmpty(self):
    '''
    table with two column and four rows:
    
      foo     bar
    ---------------
      i         10
      love      11
      unit      12
      tests     13

    '''
    
    tab = Table(['foo', 'bar'], 'si', bar=list(range(10,14)), foo=['i','love','unit','tests'])
    self.CompareColCount(tab, 2)
    self.CompareRowCount(tab, 4)
    self.CompareColNames(tab, ['foo','bar'])
    self.CompareColTypes(tab, ['foo', 'bar'], 'si')
    self.CompareDataFromDict(tab, {'bar': [10,11,12,13], 'foo': ['i','love','unit','tests']})
    
  def testTableInitMultiColMissingMultiValue(self):
    '''
    test if error is raised when creating rows with missing data
    '''
    
    self.assertRaises(ValueError, Table, ['foo', 'bar'], 'si',
                      bar=list(range(10,14)), foo=['i','love','tests'])
    
    
  def testTableInitMultiColMultiValueAndNoneNonEmpty(self):
    '''
    table with two column and four rows with None values:
    
      foo     bar
    ---------------
      i         NA
      love      NA
      unit      NA
      tests     NA

    '''
    tab = Table(['foo', 'bar'], 'si', foo=['i','love','unit','tests'])
    self.CompareColCount(tab, 2)
    self.CompareRowCount(tab, 4)
    self.CompareColNames(tab, ['foo','bar'])
    self.CompareColTypes(tab, ['foo', 'bar'], 'si')
    self.CompareDataFromDict(tab, {'bar': [None,None,None,None], 'foo': ['i','love','unit','tests']})
  
  def testTableAddSingleCol(self):
    '''
    init empty table, add one empty column:
    
     first 
    -------
    
    '''
    tab = Table()
    self.CompareColCount(tab, 0)
    self.CompareRowCount(tab, 0)
    tab.AddCol('first', 'string')
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 0)
    self.CompareColNames(tab, ['first'])
    self.CompareColTypes(tab, 'first', 's')
    
  def testTableAddSingleRow(self):
    '''
    init table with one col, add one row:
    
     first 
    -------
          2
    '''
    tab = Table(['first'],'i')
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 0)
    tab.AddRow([2], overwrite=None)
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 1)
    self.CompareColNames(tab, ['first'])
    self.CompareColTypes(tab, 'first', 'i')
    self.CompareDataFromDict(tab, {'first': [2]})
    
  def testTableAddSingleColSingleRow(self):
    '''
    init empty table, add one col, add one row:
    
     first 
    -------
          2
    '''
    tab = Table()
    tab.AddCol('first', 'int')
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 0)
    tab.AddRow([2], overwrite=None)
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 1)
    self.CompareColNames(tab, ['first'])
    self.CompareColTypes(tab, 'first', 'i')
    self.CompareDataFromDict(tab, {'first': [2]})
  
  def testTableAddSingleColWithRow(self):
    '''
    init table with two cols, add row with data, add third column:
    
     first  second  third 
    ----------------------
     x            3  3.141

    '''
    tab = Table(['first','second'],'si')
    self.CompareColCount(tab, 2)
    self.CompareRowCount(tab, 0)
    self.CompareColTypes(tab, ['first','second'], 'si')
    tab.AddRow(['x',3], overwrite=None)
    self.CompareColCount(tab, 2)
    self.CompareRowCount(tab, 1)
    tab.AddCol('third', 'float', 3.141)
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 1)
    self.CompareColTypes(tab, ['first','third','second'], 'sfi')
    self.CompareDataFromDict(tab, {'second': [3], 'first': ['x'], 'third': [3.141]})
    
  def testTableAddMultiColMultiRow(self):
    '''
    init empty table add three cols, add three rows with data:
    
      first  second  third 
    ----------------------
     x            3  1.000
     foo          6  2.200
     bar          9  3.300

    '''
    tab = Table()
    tab.AddCol('first', 'string')
    tab.AddCol('second', 'int')
    tab.AddCol('third', 'float')
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 0)
    self.CompareColTypes(tab, ['first','second', 'third'], 'sif')
    tab.AddRow(['x',3, 1.0], overwrite=None)
    tab.AddRow(['foo',6, 2.2], overwrite=None)
    tab.AddRow(['bar',9, 3.3], overwrite=None)
    tab.AddCol('fourth','bool',itertools.cycle([True,False]))
    self.CompareColCount(tab, 4)
    self.CompareRowCount(tab, 3)
    self.CompareDataFromDict(tab, {'second': [3,6,9], 'first': ['x','foo','bar'], 'third': [1,2.2,3.3], 'fourth': [True,False,True]})

  def testTableAddMultiColMultiRowFromDict(self):
    '''
    init empty table add three cols, add three rows with data:
    
      first  second  third 
    ----------------------
     x            3  1.000
     foo          6  2.200
     bar          9  3.300

    '''
    tab = Table()
    tab.AddCol('first', 'string')
    tab.AddCol('second', 'int')
    tab.AddCol('aaa', 'float')
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 0)
    self.CompareColTypes(tab, ['first','second', 'aaa'], 'sif')
    tab.AddRow({'first':'x','second':3, 'aaa':1.0}, overwrite=None)
    tab.AddRow({'aaa':2.2, 'second':6, 'first':'foo'}, overwrite=None)
    tab.AddRow({'second':9, 'aaa':3.3, 'first':'bar'}, overwrite=None)
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 3)
    self.CompareDataFromDict(tab, {'second': [3,6,9], 'first': ['x','foo','bar'], 'aaa': [1,2.2,3.3]})
    
  def testTableAddMultiRowMultiCol(self):
    '''
    init empty table add one col, add three rows with data,
    add one col without data, add one col with data:
    
      first  second  third 
    ----------------------
     x        NA     3.141
     foo      NA     3.141
     bar      NA     3.141

    '''
    tab = Table()
    tab.AddCol('first', 'string')
    self.CompareColCount(tab, 1)
    self.CompareRowCount(tab, 0)
    self.CompareColTypes(tab, ['first'], 's')
    tab.AddRow(['x'], overwrite=None)
    tab.AddRow(['foo'], overwrite=None)
    tab.AddRow(['bar'], overwrite=None)
    tab.AddCol('second', 'int')
    tab.AddCol('third', 'float', 3.141)
    self.CompareColCount(tab, 3)
    self.CompareRowCount(tab, 3)
    self.CompareDataFromDict(tab, {'second': [None,None,None],
                                   'first': ['x','foo','bar'],
                                   'third': [3.141, 3.141, 3.141]})

  def testAddMultiRowFromDict(self):
    tab = Table(['x','y','z'], 'fff')
    data = {'x': [1.2, 1.5], 'z': [1.6, 2.4]}
    tab.AddRow(data)
    self.CompareDataFromDict(tab, {'x': [1.2, 1.5],
                                   'y': [None, None],
                                   'z': [1.6, 2.4]})

    data = {'y': [5.1, 3.4, 1.5]}
    tab.AddRow(data)
    self.CompareDataFromDict(tab, {'x': [1.2, 1.5, None, None, None],
                                   'y': [None, None, 5.1, 3.4, 1.5],
                                   'z': [1.6, 2.4, None, None, None]})

    # must raise since length of data is not the same
    data = {'x': [1.2, 1.5], 'y': 1.2, 'z': [1.6, 2.4]}
    self.assertRaises(ValueError, tab.AddRow, data)

    # must raise since length of data is not the same
    data = {'x': [1.2, 1.5], 'y': [1.2], 'z': [1.6, 2.4]}
    self.assertRaises(ValueError, tab.AddRow, data)

    # overwrite certain rows
    data = {'x': [1.2, 1.9], 'z': [7.9, 3.5]}
    tab.AddRow(data, overwrite='x')
    self.CompareDataFromDict(tab, {'x': [1.2, 1.5, None, None, None, 1.9],
                                   'y': [None, None, 5.1, 3.4, 1.5, None],
                                   'z': [7.9, 2.4, None, None, None, 3.5]})

  def testAddRowFromDictWithOverwrite(self):
    '''
    add rows from dictionary with overwrite (i.e. overwrite third row with additional data)
    
      x     foo   bar 
    ------------------
     row1  True      1
     row2    NA      2
     row3  False     3
     
    '''
    tab = Table()
    tab.AddCol('x', 'string')
    tab.AddCol('foo', 'bool')
    tab.AddCol('bar', 'int')
    tab.AddRow(['row1',True, 1])
    tab.AddRow(['row2',None, 2])
    tab.AddRow(['row3',False, None])
    self.CompareDataFromDict(tab, {'x': ['row1', 'row2', 'row3'],
                                   'foo': [True, None, False],
                                   'bar': [1, 2, None]})
    tab.AddRow({'x':'row3', 'bar':3}, overwrite='x')
    self.CompareDataFromDict(tab, {'x': ['row1', 'row2', 'row3'],
                                   'foo': [True, None, False],
                                   'bar': [1, 2, 3]})
    
  def testAddRowFromListWithOverwrite(self):
    '''
    add rows from list with overwrite (i.e. overwrite third row with additional data)
    
      x     foo   bar 
    ------------------
     row1  True      1
     row2    NA      2
     row3  True      3
     
    '''
    
    tab = Table()
    tab.AddCol('x', 'string')
    tab.AddCol('foo', 'bool')
    tab.AddCol('bar', 'int')
    tab.AddRow(['row1',True, 1])
    tab.AddRow(['row2',None, 2])
    tab.AddRow(['row3',False, None])
    self.CompareDataFromDict(tab, {'x': ['row1', 'row2', 'row3'],
                                   'foo': [True, None, False],
                                   'bar': [1, 2, None]})
    tab.AddRow(['row3', True, 3], overwrite='x')
    self.CompareDataFromDict(tab, {'x': ['row1', 'row2', 'row3'],
                                   'foo': [True, None, True],
                                   'bar': [1, 2, 3]})

  def testRaiseErrorOnWrongDataLengthAddCol(self):
    tab = Table()
    tab.AddCol('a','f',[4.2,4.2,4.2])
    self.assertRaises(ValueError, tab.AddCol, 'b', 'f', [4.2,4.2])

  def testRaiseErrorColNameAlreadyExists(self):
    tab = Table()
    tab.AddCol('awesome','f')
    self.assertRaises(ValueError, tab.AddCol, 'awesome', 'f')

  def testRaiseErrorOnWrongColumnTypes(self):
    # wrong columns types in init
    self.assertRaises(ValueError, Table, ['bla','bli'], 'ab')
    
    # wrong column types in Coerce
    tab = Table()
    self.assertRaises(ValueError, tab._Coerce, 'bla', 'a')

    # wrong column types in AddCol
    self.assertRaises(ValueError, tab.AddCol, 'bla', 'a')
    
  def testParseColumnTypes(self):
    types = Table._ParseColTypes(['i','f','s','b'])
    self.assertEqual(types, ['int','float','string','bool'])
    
    types = Table._ParseColTypes(['int','float','string','bool'])
    self.assertEqual(types, ['int','float','string','bool'])
    
    types = Table._ParseColTypes(['i','float','s','bool'])
    self.assertEqual(types, ['int','float','string','bool'])

    types = Table._ParseColTypes(['i','fLOAT','S','bool'])
    self.assertEqual(types, ['int','float','string','bool'])
    
    types = Table._ParseColTypes('ifsb')
    self.assertEqual(types, ['int','float','string','bool'])
    
    types = Table._ParseColTypes('int,float,string,bool')
    self.assertEqual(types, ['int','float','string','bool'])
    
    types = Table._ParseColTypes('int,f,s,bool')
    self.assertEqual(types, ['int','float','string','bool'])
    
    types = Table._ParseColTypes('INT,F,s,bOOL')
    self.assertEqual(types, ['int','float','string','bool'])

    types = Table._ParseColTypes('boOl')
    self.assertEqual(types, ['bool'])
    
    types = Table._ParseColTypes('S')
    self.assertEqual(types, ['string'])
    
    types = Table._ParseColTypes(['i'])
    self.assertEqual(types, ['int'])
    
    types = Table._ParseColTypes(['FLOAT'])
    self.assertEqual(types, ['float'])

    self.assertRaises(ValueError, Table._ParseColTypes, 'bfstring')
    self.assertRaises(ValueError, Table._ParseColTypes, ['b,f,string'])
    self.assertRaises(ValueError, Table._ParseColTypes, 'bi2')
    self.assertRaises(ValueError, Table._ParseColTypes, ['b',2,'string'])
    self.assertRaises(ValueError, Table._ParseColTypes, [['b'],['f','string']])
    self.assertRaises(ValueError, Table._ParseColTypes, 'a')
    self.assertRaises(ValueError, Table._ParseColTypes, 'foo')
    self.assertRaises(ValueError, Table._ParseColTypes, ['a'])
    self.assertRaises(ValueError, Table._ParseColTypes, ['foo'])
  
  def testShortLongColumnTypes(self):
    tab = Table(['x','y','z','a'],['i','f','s','b'])
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],['int','float','string','bool'])
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],['i','float','s','bool'])
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],['i','fLOAT','S','bool'])
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],'ifsb')
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],'int,float,string,bool')
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],'int,f,s,bool')
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x','y','z','a'],'INT,F,s,bOOL')
    self.CompareColTypes(tab, ['x','y','z','a'], 'ifsb')
    
    tab = Table(['x'], 'boOl')
    self.CompareColTypes(tab, ['x'], 'b')
    tab = Table(['x'], 'B')
    self.CompareColTypes(tab, ['x'], 'b')
    tab = Table(['x'], ['b'])
    self.CompareColTypes(tab, ['x'], 'b')
    tab = Table(['x'], ['Bool'])
    self.CompareColTypes(tab, ['x'], 'b')
    
    self.assertRaises(ValueError, Table, ['x','y','z'], 'bfstring')
    self.assertRaises(ValueError, Table, ['x','y','z'], ['b,f,string'])
    self.assertRaises(ValueError, Table, ['x','y','z'], 'bi2')
    self.assertRaises(ValueError, Table, ['x','y','z'], ['b',2,'string'])
    self.assertRaises(ValueError, Table, ['x','y','z'], [['b'],['f','string']])
    self.assertRaises(ValueError, Table, ['x'], 'a')
    self.assertRaises(ValueError, Table, ['x'], 'foo')
    self.assertRaises(ValueError, Table, ['x'], ['a'])
    self.assertRaises(ValueError, Table, ['x'], ['foo'])
    
  def testCoerce(self):
    tab = Table()
    
    # None values
    self.assertEqual(tab._Coerce('NA', 'x'), None)
    self.assertEqual(tab._Coerce(None, 'x'), None)
    
    # int type
    self.assertTrue(isinstance(tab._Coerce(2 ,'int'), int))
    self.assertEqual(tab._Coerce(2 ,'int'), 2)
    self.assertTrue(isinstance(tab._Coerce(2.2 ,'int'), int))
    self.assertEqual(tab._Coerce(2.2 ,'int'), 2)
    self.assertEqual(tab._Coerce(True ,'int'), 1)
    self.assertEqual(tab._Coerce(False ,'int'), 0)
    self.assertRaises(ValueError, tab._Coerce, "foo" , 'int')
    
    # float type
    self.assertTrue(isinstance(tab._Coerce(2 ,'float'), float))
    self.assertEqual(tab._Coerce(2 ,'float'), 2.000)
    self.assertTrue(isinstance(tab._Coerce(3.141 ,'float'), float))
    self.assertEqual(tab._Coerce(3.141 ,'float'), 3.141)
    self.assertRaises(ValueError, tab._Coerce, "foo" , 'float')
    
    # string type
    self.assertTrue(isinstance(tab._Coerce('foo' ,'string'), str))
    self.assertTrue(isinstance(tab._Coerce('this is a longer string' ,'string'), str))
    self.assertTrue(isinstance(tab._Coerce(2.2 ,'string'), str))
    self.assertTrue(isinstance(tab._Coerce(2 ,'string'), str))
    self.assertTrue(isinstance(tab._Coerce(True ,'string'), str))
    self.assertTrue(isinstance(tab._Coerce(False ,'string'), str))
    
    # bool type
    self.assertEqual(tab._Coerce(True ,'bool'), True)
    self.assertEqual(tab._Coerce(False ,'bool'), False)
    self.assertEqual(tab._Coerce('falSE' ,'bool'), False)
    self.assertEqual(tab._Coerce('no' ,'bool'), False)
    self.assertEqual(tab._Coerce('not false and not no','bool'), True)
    self.assertEqual(tab._Coerce(0, 'bool'), False)
    self.assertEqual(tab._Coerce(1, 'bool'), True)
    
    # unknown type
    self.assertRaises(ValueError, tab._Coerce, 'bla', 'abc')
    
  def testRemoveCol(self):
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    tab.RemoveCol("second")
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'third': [None,2.2,3.3]})
    
    # raise error when column is unknown
    tab = self.CreateTestTable()
    self.assertRaises(ValueError, tab.RemoveCol, "unknown col")
    
  def testSortTable(self):
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    tab.Sort('first', '-')
    self.CompareDataFromDict(tab, {'first': [None,'foo','x'], 'second': [9,None,3], 'third': [3.3,2.2,None]})
    tab.Sort('first', '+')
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    tab.Sort('third', '+')
    self.CompareDataFromDict(tab, {'first': [None,'foo','x'], 'second': [9,None,3], 'third': [3.3,2.2,None]})

  def testGuessFormat(self):
    self.assertEqual(Table._GuessFormat('table_test.csv'), 'csv')
    self.assertEqual(Table._GuessFormat('table_test.pickle'), 'pickle')
    self.assertEqual(Table._GuessFormat('table_test.tab'), 'ost')
    self.assertEqual(Table._GuessFormat('table_test.ost'), 'ost')
    self.assertEqual(Table._GuessFormat('table_test.xyz'), 'ost')
    
  def testSaveLoadTableAutoFormat(self):   
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})

    # write to disc
    tab.Save("saveloadtable_filename_out.csv", format='csv')
    tab.Save("saveloadtable_filename_out.tab", format='ost')
    tab.Save("saveloadtable_filename_out.pickle", format='pickle')
    
    # read from disc: csv
    in_stream_csv = open("saveloadtable_filename_out.csv", 'r')
    tab_loaded_stream_csv = Table.Load(in_stream_csv)
    in_stream_csv.close()
    tab_loaded_fname_csv = Table.Load('saveloadtable_filename_out.csv')

    # check content: csv
    self.CompareDataFromDict(tab_loaded_stream_csv, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    self.CompareDataFromDict(tab_loaded_fname_csv, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
  
    # read from disc: pickle
    in_stream_pickle = open("saveloadtable_filename_out.pickle", 'rb')
    tab_loaded_stream_pickle = Table.Load(in_stream_pickle)
    in_stream_pickle.close()
    tab_loaded_fname_pickle = Table.Load('saveloadtable_filename_out.pickle')

    # check content: pickle
    self.CompareDataFromDict(tab_loaded_stream_pickle, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    self.CompareDataFromDict(tab_loaded_fname_pickle, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})

    # read from disc: ost
    in_stream_ost = open("saveloadtable_filename_out.tab", 'r')
    tab_loaded_stream_ost = Table.Load(in_stream_ost)
    in_stream_ost.close()
    tab_loaded_fname_ost = Table.Load('saveloadtable_filename_out.tab')

    # check content: ost
    self.CompareDataFromDict(tab_loaded_stream_ost, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    self.CompareDataFromDict(tab_loaded_fname_ost, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
  

  def testLoadTableOSTUnknownType(self):
    self.assertRaises(ValueError, Table.Load, os.path.join('testfiles','ost-table-unknown-type.tab'))

  def testLoadTableOSTNoType(self):
    tab = Table.Load(os.path.join('testfiles','ost-table-notype.tab'))
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    
  def testLoadOSTDifficultHeaders(self):
    tab = Table.Load(os.path.join('testfiles','ost-table-difficult-headers.tab'))
    self.assertEqual(tab.col_types, ['float','float','float','float','float'])

  def testSaveLoadTableOST(self):
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    
    # write to disc
    tab.Save("saveloadtable_filename_out.tab")
    out_stream = open("saveloadtable_stream_out.tab", 'w')
    tab.Save(out_stream)
    out_stream.close()
    
    # read from disc
    in_stream = open("saveloadtable_stream_out.tab", 'r')
    tab_loaded_stream = Table.Load(in_stream)
    in_stream.close()
    tab_loaded_fname = Table.Load('saveloadtable_filename_out.tab')
    
    # check content
    self.CompareDataFromDict(tab_loaded_stream, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    self.CompareDataFromDict(tab_loaded_fname, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    
    # check Errors for empty/non existing files
    self.assertRaises(IOError, Table.Load, 'nonexisting.file')
    self.assertRaises(IOError, Table.Load, os.path.join('testfiles','emptytable.tab'))
    in_stream = open(os.path.join('testfiles','emptytable.csv'), 'r')
    self.assertRaises(IOError, Table.Load, in_stream)
    in_stream.close()
    
  def testSaveLoadTableOSTWithSpaces(self):
    tab = self.CreateTestTable()
    tab.AddRow(['hello spaces',10, 10.1], overwrite=None)
    self.CompareDataFromDict(tab, {'first': ['x','foo',None,'hello spaces'], 'second': [3,None,9,10], 'third': [None,2.2,3.3,10.1]})

    # write to disc
    tab.Save("saveloadtable_withspaces_filename_out.tab")

    # read from disc
    tab_loaded_fname = Table.Load('saveloadtable_withspaces_filename_out.tab')
    self.CompareDataFromDict(tab_loaded_fname, {'first': ['x','foo',None,'hello spaces'], 'second': [3,None,9,10], 'third': [None,2.2,3.3,10.1]})
  def testSaveTableHTML(self):
    import io
    tab = self.CreateTestTable()
    stream = io.StringIO()
    tab.Save(stream, format='html')
    self.assertEqual(stream.getvalue(), '<table><tr><th>first</th><th>second</th><th>third</th></tr><tr><td>x</td><td>3</td><td></td></tr><tr><td>foo</td><td></td><td>2.200</td></tr><tr><td></td><td>9</td><td>3.300</td></tr></table>')
  def testSaveTableContext(self):
    import io
    tab = self.CreateTestTable()
    stream = io.StringIO()
    tab.Save(stream, format='context')
    self.assertEqual(stream.getvalue(), 
                     '\\starttable[l|r|i3r|]\n\\HL\n\\NC \\bf first\\NC \\bf second\\NC \\bf third \\AR\\HL\n\\NC x\\NC 3\\NC --- \\AR\n\\NC foo\\NC ---\\NC 2.200 \\AR\n\\NC ---\\NC 9\\NC 3.300 \\AR\n\\HL\n\\stoptable')

  def testSaveLoadTableCSV(self):
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})

    # write to disc
    tab.Save("saveloadtable_filename_out.csv", format='csv')
    out_stream = open("saveloadtable_stream_out.csv", 'w')
    tab.Save(out_stream, format='csv')
    out_stream.close()
    
    # read from disc
    in_stream = open("saveloadtable_stream_out.csv", 'r')
    tab_loaded_stream = Table.Load(in_stream, format='csv')
    in_stream.close()
    tab_loaded_fname = Table.Load('saveloadtable_filename_out.csv', format='csv')

    # check content
    self.CompareDataFromDict(tab_loaded_stream, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    self.CompareDataFromDict(tab_loaded_fname, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
  
  def testSaveLoadTablePickle(self):
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    # write to disc
    tab.Save("saveloadtable_filename_out.pickle", format='pickle')
    out_stream = open("saveloadtable_stream_out.pickle", 'wb')
    tab.Save(out_stream, format='pickle')
    out_stream.close()

    # read from disc
    in_stream = open("saveloadtable_stream_out.pickle", 'rb')
    tab_loaded_stream = Table.Load(in_stream, format='pickle')
    in_stream.close()
    tab_loaded_fname = Table.Load('saveloadtable_filename_out.pickle', format='pickle')

    # check content
    self.CompareDataFromDict(tab_loaded_stream, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})
    self.CompareDataFromDict(tab_loaded_fname, {'first': ['x','foo',None], 'second': [3,None,9], 'third': [None,2.2,3.3]})

  def testMergeTable(self):
    '''
    Merge the following two tables:
    
    x     y           x     u   
    -------           -------
    1 |  10           1 | 100
    2 |  15           3 | 200
    3 |  20           4 | 400  
    
    to get (only_matching=False):
    
    x      y     u
    ---------------
    1 |   10 |  100
    2 |   15 |   NA
    3 |   20 |  200
    4 |   NA |  400
    
    or (only_matching=True):
    
    x      y     u
    ---------------
    1 |   10 |  100
    3 |   20 |  200
    
    '''
    tab1 = Table(['x','y'],['int','int'])
    tab1.AddRow([1,10])
    tab1.AddRow([2,15])
    tab1.AddRow([3,20])
    
    tab2 = Table(['x','u'],['int','int'])
    tab2.AddRow([1,100])
    tab2.AddRow([3,200])
    tab2.AddRow([4,400])
    
    tab_merged = Merge(tab1, tab2, 'x', only_matching=False)
    tab_merged.Sort('x', order='-')
    self.CompareDataFromDict(tab_merged, {'x': [1,2,3,4], 'y': [10,15,20,None], 'u': [100,None,200,400]})
    
    tab_merged = Merge(tab1, tab2, 'x', only_matching=True)
    tab_merged.Sort('x', order='-')
    self.CompareDataFromDict(tab_merged, {'x': [1,3], 'y': [10,20], 'u': [100,200]})
    
  def testFilterTable(self):
    tab = self.CreateTestTable()
    tab.AddRow(['foo',1,5.15])
    tab.AddRow(['foo',0,1])
    tab.AddRow(['foo',1,12])
    
    # filter on one column
    tab_filtered = tab.Filter(first='foo')
    self.CompareDataFromDict(tab_filtered, {'first':['foo','foo','foo','foo'],
                                            'second':[None,1,0,1],
                                            'third':[2.2,5.15,1.0,12.0]})
    
    # filter on multiple columns
    tab_filtered = tab.Filter(first='foo',second=1)
    self.CompareDataFromDict(tab_filtered, {'first':['foo','foo'],
                                            'second':[1,1],
                                            'third':[5.15,12.0]})
    
    # raise Error when using non existing column name for filtering
    self.assertRaises(ValueError,tab.Filter,first='foo',nonexisting=1)
    
  def testMinTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[True,True,False])

    self.assertEqual(tab.Min('first'),'foo')
    self.assertEqual(tab.Min('second'),3)
    self.assertAlmostEqual(tab.Min('third'),2.2)
    self.assertEqual(tab.Min('fourth'),False)
    self.assertRaises(ValueError,tab.Min,'fifth')
    
    self.assertEqual(tab.MinIdx('first'),1)
    self.assertEqual(tab.MinIdx('second'),0)
    self.assertAlmostEqual(tab.MinIdx('third'),1)
    self.assertEqual(tab.MinIdx('fourth'),2)
    self.assertRaises(ValueError,tab.MinIdx,'fifth')
    
    self.assertEqual(tab.MinRow('first'),['foo', None, 2.20, True])
    self.assertEqual(tab.MinRow('second'),['x', 3, None, True])
    self.assertEqual(tab.MinRow('third'),['foo', None, 2.20, True])
    self.assertEqual(tab.MinRow('fourth'),[None, 9, 3.3, False])
    self.assertRaises(ValueError,tab.MinRow,'fifth')
    
  def testMaxTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[False,True,True])
    
    self.assertEqual(tab.Max('first'),'x')
    self.assertEqual(tab.Max('second'),9)
    self.assertAlmostEqual(tab.Max('third'),3.3)
    self.assertEqual(tab.Max('fourth'),True)
    self.assertRaises(ValueError,tab.Max,'fifth')
    
    self.assertEqual(tab.MaxIdx('first'),0)
    self.assertEqual(tab.MaxIdx('second'),2)
    self.assertAlmostEqual(tab.MaxIdx('third'),2)
    self.assertEqual(tab.MaxIdx('fourth'),1)
    self.assertRaises(ValueError,tab.MaxIdx,'fifth')
    
    self.assertEqual(tab.MaxRow('first'),['x', 3, None, False])
    self.assertEqual(tab.MaxRow('second'),[None, 9, 3.3, True])
    self.assertEqual(tab.MaxRow('third'),[None, 9, 3.3, True])
    self.assertEqual(tab.MaxRow('fourth'),['foo', None, 2.2, True])
    self.assertRaises(ValueError,tab.MaxRow,'fifth')
    
  def testSumTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[False,True,False])
    tab.AddCol('fifth','string',['foo','bar',None])
    
    self.assertRaises(TypeError,tab.Sum,'first')
    self.assertEqual(tab.Sum('second'),12)
    self.assertAlmostEqual(tab.Sum('third'),5.5)
    self.assertEqual(tab.Sum('fourth'),1)
    self.assertRaises(TypeError,tab.Sum,'fifth')
    self.assertRaises(ValueError,tab.Sum,'sixth')
    
  def testMedianTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[False,True,False])
    tab.AddCol('fifth','string',['foo','bar',None])
    
    self.assertRaises(TypeError,tab.Median,'first')
    self.assertEqual(tab.Median('second'),6.0)
    self.assertAlmostEqual(tab.Median('third'),2.75)
    self.assertEqual(tab.Median('fourth'),False)
    self.assertRaises(TypeError,tab.Median,'fifth')
    self.assertRaises(ValueError,tab.Median,'sixth')
    
  def testMeanTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[False,True,False])
    tab.AddCol('fifth','string',['foo','bar',None])
    
    self.assertRaises(TypeError,tab.Mean,'first')
    self.assertAlmostEqual(tab.Mean('second'),6.0)
    self.assertAlmostEqual(tab.Mean('third'),2.75)
    self.assertAlmostEqual(tab.Mean('fourth'),0.33333333)
    self.assertRaises(TypeError,tab.Mean,'fifth')
    self.assertRaises(ValueError,tab.Mean,'sixth')
    
  def testRowMeanTable(self):
    '''
      first  second  third fourth
    -----------------------------
     x            3     NA      1
     foo         NA  2.200      2
     NA           9  3.300      3
     NA          NA     NA     NA
    '''
    tab = self.CreateTestTable()
    tab.AddCol('fourth','float',[1,2,3])
    tab.AddRow([None, None, None, None])
    
    self.assertRaises(TypeError, tab.RowMean, 'mean', ['first', 'second'])
    tab.RowMean('mean', ['third', 'second', 'fourth'])
    self.CompareDataFromDict(tab, {'mean': [2,2.1,5.1,None],
                                   'first': ['x','foo',None,None],
                                   'second': [3,None,9,None],
                                   'third': [None,2.2,3.3,None],
                                   'fourth': [1,2,3,None]})
    
    
  def testStdDevTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[False,True,False])
    tab.AddCol('fifth','string',['foo','bar',None])
    
    self.assertRaises(TypeError,tab.StdDev,'first')
    self.assertAlmostEqual(tab.StdDev('second'),3.0)
    self.assertAlmostEqual(tab.StdDev('third'),0.55)
    self.assertAlmostEqual(tab.StdDev('fourth'),0.47140452079)
    self.assertRaises(TypeError,tab.StdDev,'fifth')
    self.assertRaises(ValueError,tab.StdDev,'sixth')
    
  def testCountTable(self):
    tab = self.CreateTestTable()
    tab.AddCol('fourth','bool',[False,True,False])
    
    self.assertEqual(tab.Count('first'),2)
    self.assertEqual(tab.Count('second'),2)
    self.assertEqual(tab.Count('third'),2)
    self.assertEqual(tab.Count('fourth'),3)
    self.assertEqual(tab.Count('first', ignore_nan=False),3)
    self.assertEqual(tab.Count('second', ignore_nan=False),3)
    self.assertEqual(tab.Count('third', ignore_nan=False),3)
    self.assertEqual(tab.Count('fourth', ignore_nan=False),3)
    self.assertRaises(ValueError,tab.Count,'fifth')
    
  def testCalcEnrichment(self):
    enrx_ref = [0.0, 0.041666666666666664, 0.083333333333333329, 0.125, 0.16666666666666666, 0.20833333333333334, 0.25, 0.29166666666666669, 0.33333333333333331, 0.375, 0.41666666666666669, 0.45833333333333331, 0.5, 0.54166666666666663, 0.58333333333333337, 0.625, 0.66666666666666663, 0.70833333333333337, 0.75, 0.79166666666666663, 0.83333333333333337, 0.875, 0.91666666666666663, 0.95833333333333337, 1.0]
    enry_ref = [0.0, 0.16666666666666666, 0.33333333333333331, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.66666666666666663, 0.66666666666666663, 0.83333333333333337, 0.83333333333333337, 0.83333333333333337, 0.83333333333333337, 0.83333333333333337, 0.83333333333333337, 0.83333333333333337, 1.0, 1.0, 1.0, 1.0]
    
    tab = Table(['score', 'rmsd', 'classific'], 'ffb',
                score=[2.64,1.11,2.17,0.45,0.15,0.85,1.13,2.90,0.50,1.03,1.46,2.83,1.15,2.04,0.67,1.27,2.22,1.90,0.68,0.36,1.04,2.46,0.91,0.60],
                rmsd=[9.58,1.61,7.48,0.29,1.68,3.52,3.34,8.17,4.31,2.85,6.28,8.78,0.41,6.29,4.89,7.30,4.26,3.51,3.38,0.04,2.21,0.24,7.58,8.40],
                classific=[False,True,False,True,True,False,False,False,False,False,False,False,True,False,False,False,False,False,False,True,False,True,False,False])
    
    enrx,enry = tab.ComputeEnrichment(score_col='score', score_dir='-',
                                      class_col='rmsd', class_cutoff=2.0,
                                      class_dir='-')
    
    for x,y,refx,refy in zip(enrx,enry,enrx_ref,enry_ref):
      self.assertAlmostEqual(x,refx)
      self.assertAlmostEqual(y,refy)
    
    enrx,enry = tab.ComputeEnrichment(score_col='score', score_dir='-',
                                      class_col='classific')
    
    for x,y,refx,refy in zip(enrx,enry,enrx_ref,enry_ref):
      self.assertAlmostEqual(x,refx)
      self.assertAlmostEqual(y,refy)
    
    tab.AddCol('bad','string','col')
    
    self.assertRaises(TypeError, tab.ComputeEnrichment, score_col='classific',
                      score_dir='-', class_col='rmsd', class_cutoff=2.0,
                      class_dir='-')
    
    self.assertRaises(TypeError, tab.ComputeEnrichment, score_col='bad',
                      score_dir='-', class_col='rmsd', class_cutoff=2.0,
                      class_dir='-')
    
    self.assertRaises(TypeError, tab.ComputeEnrichment, score_col='score',
                      score_dir='-', class_col='bad', class_cutoff=2.0,
                      class_dir='-')
    
    self.assertRaises(ValueError, tab.ComputeEnrichment, score_col='score',
                      score_dir='x', class_col='rmsd', class_cutoff=2.0,
                      class_dir='-')
    
    self.assertRaises(ValueError, tab.ComputeEnrichment, score_col='score',
                      score_dir='+', class_col='rmsd', class_cutoff=2.0,
                      class_dir='y')
    
  def testPlot(self):
    if not HAS_MPL or not HAS_NUMPY:
      return
    tab = self.CreateTestTable()
    self.assertRaises(ValueError, tab.Plot, 'second', x_range=1)
    self.assertRaises(ValueError, tab.Plot, x='second', y='third', y_range=[1,2,3])
    self.assertRaises(ValueError, tab.Plot, x='second', y='third', z_range='st')

  def testHexbin(self):
    if not HAS_MPL or not HAS_NUMPY:
      return
    tab = self.CreateTestTable()
    self.assertRaises(ValueError, tab.PlotHexbin, x='second', y='third', x_range=1)
    self.assertRaises(ValueError, tab.PlotHexbin, x='second', y='third', x_range=[1,2,3])

  def testPlotEnrichment(self):
    if not HAS_MPL or not HAS_PIL:
      return
    tab = Table(['score', 'rmsd', 'classific'], 'ffb',
                score=[2.64,1.11,2.17,0.45,0.15,0.85,1.13,2.90,0.50,1.03,1.46,2.83,1.15,2.04,0.67,1.27,2.22,1.90,0.68,0.36,1.04,2.46,0.91,0.60],
                rmsd=[9.58,1.61,7.48,0.29,1.68,3.52,3.34,8.17,4.31,2.85,6.28,8.78,0.41,6.29,4.89,7.30,4.26,3.51,3.38,0.04,2.21,0.24,7.58,8.40],
                classific=[False,True,False,True,True,False,False,False,False,False,False,False,True,False,False,False,False,False,False,True,False,True,False,False])
 
    pl = tab.PlotEnrichment(score_col='score', score_dir='-',
                            class_col='rmsd', class_cutoff=2.0,
                            class_dir='-',
                            save=os.path.join("testfiles","enrichment-out.png"))
    img1 = Image.open(os.path.join("testfiles","enrichment-out.png"))
    #img2 = Image.open(os.path.join("testfiles","enrichment.png"))
    #self.CompareImages(img1, img2)
    #pl.show()
    
  def testCalcEnrichmentAUCwithNone(self):
    if not HAS_NUMPY:
      return
    tab = Table(['pred_bfactors','ref_distances'], 'ff',
                ref_distances=[None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 2.445, 2.405, 2.361, 2.124, 1.957, 1.897, 1.422, 1.348, 1.247, 1.165, 1.153, 1.011, 0.992, 0.885, 0.852, 0.775, 0.757, 0.755, 0.735, 0.71, 0.656, 0.636, 0.609, 0.607, 0.604, 0.595, 0.572, 0.549, 0.458, 0.438, 0.41, 0.345, 0.304, 0.254, 0.241, 0.227, 2.68, 1.856, 1.312, 0.453],
                pred_bfactors=[1.85000,  2.01000,  2.12000,  2.14000,  2.15000,  2.18000,  2.20000,  2.26000,  2.28000,  2.31000,  2.37000,  2.38000,  2.39000,  2.39000,  2.43000,  2.43000,  2.49000,  2.51000,  2.56000,  2.58000,  2.65000,  2.67000,  2.72000,  2.75000,  2.77000,  2.81000,  2.91000,  2.95000,  3.09000,  3.12000,  3.25000,  3.30000,  3.33000,  3.38000,  3.39000,  3.41000,  3.41000,  3.45000,  3.57000,  3.59000,  3.64000,  3.76000,  3.76000,  3.92000,  3.95000,  3.95000,  4.05000,  4.06000,  4.07000,  4.14000,  4.14000,  4.18000,  4.24000,  4.28000,  4.40000,  4.43000,  4.43000,  4.48000,  4.50000,  4.51000,  4.54000,  4.63000,  4.64000,  4.79000,  4.93000,  5.07000,  5.12000,  5.20000,  5.41000,  5.42000,  5.44000,  5.52000,  5.68000,  5.78000,  5.80000,  5.93000,  6.11000,  6.31000,  6.50000,  6.53000,  6.55000,  6.60000,  6.73000,  6.79000,  6.81000,  7.44000,  8.45000,  8.81000,  9.04000,  9.29000,  9.30000, 10.99000, 11.42000, 12.55000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 99.99000, 99.99000, 99.99000, 99.99000])

    auc = tab.ComputeEnrichmentAUC(score_col='pred_bfactors', class_col='ref_distances')
    self.assertAlmostEqual(auc, 0.50714285714285)

    # when removing all None lines, no true positive is left
    auc = tab.ComputeEnrichmentAUC(score_col='ref_distances', class_col='pred_bfactors')
    self.assertEqual(auc, None)

    # when increasing the cutoff, we have again true positives
    auc = tab.ComputeEnrichmentAUC(score_col='ref_distances', class_col='pred_bfactors', class_cutoff=60)
    self.assertAlmostEqual(auc, 0.52013888888)

  def testCalcEnrichmentAUC(self):
    if not HAS_NUMPY:
      return
    auc_ref = 0.65277777777777779
    tab = Table(['score', 'rmsd', 'classific'], 'ffb',
                score=[2.64,1.11,2.17,0.45,0.15,0.85,1.13,2.90,0.50,1.03,1.46,2.83,1.15,2.04,0.67,1.27,2.22,1.90,0.68,0.36,1.04,2.46,0.91,0.60],
                rmsd=[9.58,1.61,7.48,0.29,1.68,3.52,3.34,8.17,4.31,2.85,6.28,8.78,0.41,6.29,4.89,7.30,4.26,3.51,3.38,0.04,2.21,0.24,7.58,8.40],
                classific=[False,True,False,True,True,False,False,False,False,False,False,False,True,False,False,False,False,False,False,True,False,True,False,False])
 
    auc = tab.ComputeEnrichmentAUC(score_col='score', score_dir='-',
                                   class_col='rmsd', class_cutoff=2.0,
                                   class_dir='-')
    
    self.assertAlmostEqual(auc, auc_ref)

  def testPlotROC(self):
    if not HAS_MPL or not HAS_PIL:
      return
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, True, False, True, True, True, False, False, True, False, True, False, True, False, False, False, True, False, True, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    pl = tab.PlotROC(score_col='score', score_dir='+',
                     class_col='classific',
                     save=os.path.join("testfiles","roc-out.png"))
    img1 = Image.open(os.path.join("testfiles","roc-out.png"))
    #img2 = Image.open(os.path.join("testfiles","roc.png"))
    #self.CompareImages(img1, img2)

    # no true positives
    tab = Table(['classific', 'score'], 'bf',
                classific=[False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    pl = tab.PlotROC(score_col='score', score_dir='+',
                     class_col='classific',
                     save=os.path.join("testfiles","roc-out.png"))
    self.assertEqual(pl, None)

  def testPlotLogROC(self):
    if not HAS_MPL or not HAS_PIL:
      return
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, True, False, True, True, True, False, False, True, False, True, False, True, False, False, False, True, False, True, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    pl = tab.PlotLogROC(score_col='score', score_dir='+',
                        class_col='classific',
                        save=os.path.join("testfiles","logroc-out.png"))
    img1 = Image.open(os.path.join("testfiles","logroc-out.png"))
    #img2 = Image.open(os.path.join("testfiles","roc.png"))
    #self.CompareImages(img1, img2)

    # no true positives
    tab = Table(['classific', 'score'], 'bf',
                classific=[False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    pl = tab.PlotLogROC(score_col='score', score_dir='+',
                        class_col='classific',
                        save=os.path.join("testfiles","logroc-out.png"))
    self.assertEqual(pl, None)

  def testPlotROCSameValues(self):
    if not HAS_MPL or not HAS_PIL:
      return
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, True, False, True, True, True, False, False, True, False, True, False, True, False, False, False, True, False, True, False],
                score=[0.9, 0.8, 0.7, 0.7, 0.7, 0.7, 0.53, 0.52, 0.51, 0.505, 0.4, 0.4, 0.4, 0.4, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    pl = tab.PlotROC(score_col='score', score_dir='+',
                     class_col='classific',
                     save=os.path.join("testfiles","roc-same-val-out.png"))
    img1 = Image.open(os.path.join("testfiles","roc-same-val-out.png"))
    #img2 = Image.open(os.path.join("testfiles","roc-same-val.png"))
    #self.CompareImages(img1, img2)
    #pl.show()

  def testCalcROCAUCwithNone(self):
    if not HAS_NUMPY:
      return
    tab = Table(['pred_bfactors','ref_distances'], 'ff',
                ref_distances=[None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 2.445, 2.405, 2.361, 2.124, 1.957, 1.897, 1.422, 1.348, 1.247, 1.165, 1.153, 1.011, 0.992, 0.885, 0.852, 0.775, 0.757, 0.755, 0.735, 0.71, 0.656, 0.636, 0.609, 0.607, 0.604, 0.595, 0.572, 0.549, 0.458, 0.438, 0.41, 0.345, 0.304, 0.254, 0.241, 0.227, 2.68, 1.856, 1.312, 0.453],
                pred_bfactors=[1.85000,  2.01000,  2.12000,  2.14000,  2.15000,  2.18000,  2.20000,  2.26000,  2.28000,  2.31000,  2.37000,  2.38000,  2.39000,  2.39000,  2.43000,  2.43000,  2.49000,  2.51000,  2.56000,  2.58000,  2.65000,  2.67000,  2.72000,  2.75000,  2.77000,  2.81000,  2.91000,  2.95000,  3.09000,  3.12000,  3.25000,  3.30000,  3.33000,  3.38000,  3.39000,  3.41000,  3.41000,  3.45000,  3.57000,  3.59000,  3.64000,  3.76000,  3.76000,  3.92000,  3.95000,  3.95000,  4.05000,  4.06000,  4.07000,  4.14000,  4.14000,  4.18000,  4.24000,  4.28000,  4.40000,  4.43000,  4.43000,  4.48000,  4.50000,  4.51000,  4.54000,  4.63000,  4.64000,  4.79000,  4.93000,  5.07000,  5.12000,  5.20000,  5.41000,  5.42000,  5.44000,  5.52000,  5.68000,  5.78000,  5.80000,  5.93000,  6.11000,  6.31000,  6.50000,  6.53000,  6.55000,  6.60000,  6.73000,  6.79000,  6.81000,  7.44000,  8.45000,  8.81000,  9.04000,  9.29000,  9.30000, 10.99000, 11.42000, 12.55000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 50.00000, 99.99000, 99.99000, 99.99000, 99.99000])

    auc = tab.ComputeROCAUC(score_col='pred_bfactors', class_col='ref_distances')
    self.assertAlmostEqual(auc, 0.55714285714285705)

    # when removing all None lines, no true positive is left
    auc = tab.ComputeROCAUC(score_col='ref_distances', class_col='pred_bfactors')
    self.assertEqual(auc, None)

    # when increasing the cutoff, we have again true positives
    auc = tab.ComputeROCAUC(score_col='ref_distances', class_col='pred_bfactors', class_cutoff=60)
    self.assertAlmostEqual(auc, 0.701388888888888)

  def testCalcROCAUC(self):
    if not HAS_NUMPY:
      return
    auc_ref = 0.68
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, True, False, True, True, True, False, False, True, False, True, False, True, False, False, False, True, False, True, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    auc = tab.ComputeROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertAlmostEqual(auc, auc_ref)

    # no true positives
    tab = Table(['classific', 'score'], 'bf',
                classific=[False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    auc = tab.ComputeROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertEqual(auc, None)
    
  def testLogROCAUCforPerfectCurve(self):
    if not HAS_NUMPY:
      return
    auc_ref = 1.0
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, True, True, True, True, True, False, False, False, False, False, False],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35])
    
    # test logAUC
    auc = tab.ComputeLogROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertAlmostEqual(auc, auc_ref)
    
    # test linear AUC
    auc = tab.ComputeROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertAlmostEqual(auc, auc_ref)
      
  def testCalcLogROCAUCRandomCurve(self):
    if not HAS_NUMPY:
      return
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, False, True, False, True, False, True, False, True, False, True, False],
                score=[0.9, 0.9, 0.7, 0.7, 0.55, 0.55, 0.4, 0.4, 0.3, 0.3, 0.3, 0.3])
    
    # test logAUC
    auc_ref = 0.1440197405305
    auc = tab.ComputeLogROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertAlmostEqual(auc, auc_ref)
    
    # test linear AUC    
    auc_ref = 0.5
    auc = tab.ComputeROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertAlmostEqual(auc, auc_ref)


  def testCalcROCAUCWithCutoff(self):
    if not HAS_NUMPY:
      return
    tab = Table(['classific', 'score'], 'ff',
                classific=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1],
                score=[0.9, 0.8, 0.7, 0.6, 0.55, 0.54, 0.53, 0.52, 0.51, 0.505, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    auc = tab.ComputeROCAUC(score_col='score', class_col='classific', class_cutoff=0.5)
    self.assertEqual(auc, 1.0)

    # no true positives
    auc = tab.ComputeROCAUC(score_col='score', class_col='classific', class_cutoff=1.0)
    self.assertEqual(auc, None)

  def testCalcROCFromFile(self):
    if not HAS_NUMPY:
      return
    tab = Table.Load(os.path.join('testfiles','roc_table.dat'))
    auc = tab.ComputeROCAUC(score_col='prediction', class_col='reference', class_cutoff=0.4)
    self.assertEqual(auc, 1.0)
      

  def testCalcROCAUCSameValues(self):
    if not HAS_NUMPY:
      return
    auc_ref = 0.685
    tab = Table(['classific', 'score'], 'bf',
                classific=[True, True, False, True, True, True, False, False, True, False, True, False, True, False, False, False, True, False, True, False],
                score=[0.9, 0.8, 0.7, 0.7, 0.7, 0.7, 0.53, 0.52, 0.51, 0.505, 0.4, 0.4, 0.4, 0.4, 0.36, 0.35, 0.34, 0.33, 0.30, 0.1])
    auc = tab.ComputeROCAUC(score_col='score', score_dir='+', class_col='classific')
    self.assertAlmostEqual(auc, auc_ref)

  def testCalcMCC(self):
    log = _FetchLog()
    ost.PushLogSink(log)

    tab = Table(['score', 'rmsd', 'class_rmsd', 'class_score', 'class_wrong'], 'ffbbb',
                score=      [2.64, 1.11, 2.17, 0.45,0.15,0.85, 1.13, 2.90, 0.50, 1.03, 1.46, 2.83, 1.15, 2.04, 0.67, 1.27, 2.22, 1.90, 0.68, 0.36,1.04, 2.46, 0.91,0.60],
                rmsd=[9.58,1.61,7.48,0.29,1.68,3.52,3.34,8.17,4.31,2.85,6.28,8.78,0.41,6.29,4.89,7.30,4.26,3.51,3.38,0.04,2.21,0.24,7.58,8.40],
                class_rmsd= [False,True, False,True,True,False,False,False,False,False,False,False,True, False,False,False,False,False,False,True,False,True,False,False],
                class_score=[False,False,False,True,True,True, False,False,True, False,False,False,False,False,True, False,False,False,True, True,False,False,True,True],
                class_wrong=[False,False,False,False,False,False, False,False,False, False,False,False,False,False,False, False,False,False,False, False,False,False,False,False])
    
    mcc = tab.ComputeMCC(score_col='score', score_dir='-', class_col='rmsd', class_dir='-', score_cutoff=1.0, class_cutoff=2.0)
    self.assertAlmostEqual(mcc, 0.1490711984)
    mcc = tab.ComputeMCC(score_col='class_score', class_col='class_rmsd')
    self.assertAlmostEqual(mcc, 0.1490711984)
    mcc = tab.ComputeMCC(score_col='score', score_dir='+', class_col='rmsd', class_dir='+', score_cutoff=1.0, class_cutoff=2.0)
    self.assertAlmostEqual(mcc, 0.1490711984)
    mcc = tab.ComputeMCC(score_col='score', score_dir='-', class_col='rmsd', class_dir='+', score_cutoff=1.0, class_cutoff=2.0)
    self.assertAlmostEqual(mcc, -0.1490711984)
    mcc = tab.ComputeMCC(score_col='score', score_dir='+', class_col='rmsd', class_dir='-', score_cutoff=1.0, class_cutoff=2.0)
    self.assertAlmostEqual(mcc, -0.1490711984)
    mcc = tab.ComputeMCC(score_col='class_wrong', class_col='class_rmsd')
    self.assertEqual(mcc,None)
    self.assertEqual(log.messages['WARNING'][0],
    'Could not compute MCC: MCC is not defined since factor (tp + fp) is zero')

  def testCalcMCCPreclassified(self):
    tab = Table(['reference', 'prediction1', 'prediction2'],'bbb',
                reference=  [False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True,  False, False, True,  False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, True, False, False, True,  False, False, True,  False, False, False, False, False, False, False, False, False, False, True, False, False, True, False, False, False, False, False, False, False, False],
                prediction1=[False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True,  False, True, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, True, False, True,  False, False, False, False, False, False],
                prediction2=[False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, True,  False, False, True,  False, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, True, False, False, True,  False, True,  True,  False, False, False, False, False, False, False, False, False, False, True, False, False, True, False, False, False, True,  False, False, False, False])
    mcc = tab.ComputeMCC(score_col='prediction1', class_col='reference')
    self.assertAlmostEqual(mcc, 0.538389277)
    mcc = tab.ComputeMCC(score_col='prediction2', class_col='reference')
    self.assertAlmostEqual(mcc, 0.882089673321)

  def testTableAsNumpyMatrixAsArray(self):
    if not HAS_NUMPY:
      return

    '''
    checks numpy matrix 
    
      first  second  third  fourth
    -------------------------------
     x            3     NA  True
     foo         NA  2.200  False
     NA           9  3.300  False
    '''
    
    tab = self.CreateTestTable()
    tab.AddCol('fourth','b',[True, False, False])
    m = tab.GetNumpyMatrixAsArray('second')
    mc = np.array([[3],[None],[9]])
    self.assertTrue(np.all(m==mc))
    mc = np.array([[3],[None],[10]])
    self.assertFalse(np.all(m==mc))
    m = tab.GetNumpyMatrixAsArray('third')
    mc = np.array([[None],[2.200],[3.300]])
    self.assertTrue(np.all(m==mc))
    m = tab.GetNumpyMatrixAsArray('second','third')
    mc = np.array([[3, None],[None, 2.200],[9, 3.300]])
    self.assertTrue(np.all(m==mc))
    m = tab.GetNumpyMatrixAsArray('third','second')
    mc = np.array([[None, 3],[2.200, None],[3.300, 9]])
    self.assertTrue(np.all(m==mc))

    self.assertRaises(TypeError, tab.GetNumpyMatrixAsArray, 'fourth')
    self.assertRaises(TypeError, tab.GetNumpyMatrixAsArray, 'first')
    self.assertRaises(RuntimeError, tab.GetNumpyMatrixAsArray)
    
  def testOptimalPrefactors(self):
    if not HAS_NUMPY:
      return
    tab = Table(['a','b','c','d','e','f'],
                'ffffff',
                a=[1,2,3,4,5,6,7,8,9],
                b=[2,3,4,5,6,7,8,9,10],
                c=[1,3,2,4,5,6,8,7,9],
                d=[0.1,0.1,0.1,0.2,0.3,0.3,0.4,0.5,0.8],
                e=[1,1,1,1,1,1,1,1,1],
                f=[9,9,9,9,9,9,9,9,9])
    
    pref = tab.GetOptimalPrefactors('c','a','b')
    self.assertAlmostEqual(pref[0],0.799999999)
    self.assertAlmostEqual(pref[1],0.166666666666)
    
    pref = tab.GetOptimalPrefactors('c','b','a')
    self.assertAlmostEqual(pref[0],0.166666666666)
    self.assertAlmostEqual(pref[1],0.799999999)
    
    pref = tab.GetOptimalPrefactors('c','b','a',weights='e')
    self.assertAlmostEqual(pref[0],0.166666666666)
    self.assertAlmostEqual(pref[1],0.799999999)
    
    pref = tab.GetOptimalPrefactors('c','b','a',weights='f')
    self.assertAlmostEqual(pref[0],0.166666666666)
    self.assertAlmostEqual(pref[1],0.799999999)
    
    pref = tab.GetOptimalPrefactors('c','a','b',weights='d')
    self.assertAlmostEqual(pref[0],0.6078825445851)
    self.assertAlmostEqual(pref[1],0.3394613806088)
    
    self.assertRaises(RuntimeError, tab.GetOptimalPrefactors, 'c','a','b',weight='d')
    self.assertRaises(RuntimeError, tab.GetOptimalPrefactors, 'c',weights='d')

  def testGaussianSmooth(self):
    if not HAS_SCIPY_NDIMG:
      return
    tab = Table(['a','b','c','d','e','f'],'fffffi',
                a=[0.5,1.0,2.0,3.0,2.5,1.0,0.5,2.3,1.0],
                b=[0.5,1.0,2.0,3.0,2.5,1.0,0.5,2.3,1.0],
                c=[0.5,1.0,2.0,3.0,2.5,1.0,0.5,2.3,1.0],
                d=[0.5,1.0,2.0,3.0,2.5,1.0,0.5,2.3,1.0],
                e=[0.5,1.0,2.0,3.0,2.5,1.0,0.5,2.3,1.0],
                f=[2,6,5,3,8,7,4,4,4])

    tab.GaussianSmooth('a')
    tab.GaussianSmooth('b', std=2.0)
    tab.GaussianSmooth('c', padding='wrap')
    tab.GaussianSmooth('d', padding='constant')
    tab.GaussianSmooth('e', padding='constant',c=3.0)
    tab.GaussianSmooth('f')

    ref_list=[]

    ref_list.append([0.74729766,1.20875404,1.93459464,2.39849076,2.11504816,
                     1.42457403,1.20524937,1.41025075,1.3557406])
    ref_list.append([1.23447249,1.41295267,1.65198705,1.79959835,1.78131778,
                     1.64501718,1.49728102,1.40589715,1.37147629])
    ref_list.append([0.9315564,1.24131027,1.93698455,2.39855767,2.11504816,
                     1.42450711,1.20285946,1.37769451,1.17148186])
    ref_list.append([0.5630556,1.17705895,1.93224488,2.39842384,2.11504816,
                     1.4244402,1.2005097,1.34599942,0.9872398 ])
    ref_list.append([1.46464039,1.35272941,1.94594196,2.39882533,2.11504816,
                     1.42484169,1.21420677,1.52166988,1.88882459])
    ref_list.append([3,4,4,5,6,6,4,4,4])

    tab_list=[[],[],[],[],[],[]]
    
    for row in tab.rows:
      for i,v in enumerate(row):
        tab_list[i].append(v)

    for i in range(len(ref_list[0])):
      self.assertAlmostEqual(tab_list[0][i],ref_list[0][i])
      self.assertAlmostEqual(tab_list[1][i],ref_list[1][i])
      self.assertAlmostEqual(tab_list[2][i],ref_list[2][i])
      self.assertAlmostEqual(tab_list[3][i],ref_list[3][i])
      self.assertAlmostEqual(tab_list[4][i],ref_list[4][i])
      self.assertAlmostEqual(tab_list[5][i],ref_list[5][i])
     

  def testIsEmpty(self):
    tab = Table()
    self.assertTrue(tab.IsEmpty())
    self.assertTrue(tab.IsEmpty(ignore_nan=False))
    self.assertRaises(ValueError, tab.IsEmpty, 'a')
    
    # empty table
    tab = Table(['a','b','c'], 'fff')
    self.assertTrue(tab.IsEmpty())
    self.assertTrue(tab.IsEmpty('a'))
    self.assertTrue(tab.IsEmpty('b'))
    self.assertTrue(tab.IsEmpty('c'))
    self.assertTrue(tab.IsEmpty(ignore_nan=False))
    self.assertTrue(tab.IsEmpty('a', ignore_nan=False))
    self.assertTrue(tab.IsEmpty('b', ignore_nan=False))
    self.assertTrue(tab.IsEmpty('c', ignore_nan=False))
    self.assertRaises(ValueError, tab.IsEmpty, 'd')
    
    # fill row with NAN values
    tab.AddRow([None,None,None])
    self.assertTrue(tab.IsEmpty())
    self.assertTrue(tab.IsEmpty('a'))
    self.assertTrue(tab.IsEmpty('b'))
    self.assertTrue(tab.IsEmpty('c'))
    self.assertFalse(tab.IsEmpty(ignore_nan=False))
    self.assertFalse(tab.IsEmpty('a', ignore_nan=False))
    self.assertFalse(tab.IsEmpty('b', ignore_nan=False))
    self.assertFalse(tab.IsEmpty('c', ignore_nan=False))
    
    # fill some values into column 'c' only
    tab.AddRow([None,None,1.0])
    self.assertFalse(tab.IsEmpty())
    self.assertTrue(tab.IsEmpty('a'))
    self.assertTrue(tab.IsEmpty('b'))
    self.assertFalse(tab.IsEmpty('c'))
    self.assertFalse(tab.IsEmpty('a', ignore_nan=False))
    self.assertFalse(tab.IsEmpty('b', ignore_nan=False))
    self.assertFalse(tab.IsEmpty('c', ignore_nan=False))
    
    # fill some values into all columns
    tab.AddRow([2.0,3.0,1.0])
    self.assertFalse(tab.IsEmpty())
    self.assertFalse(tab.IsEmpty('a'))
    self.assertFalse(tab.IsEmpty('b'))
    self.assertFalse(tab.IsEmpty('c'))
    self.assertFalse(tab.IsEmpty('a', ignore_nan=False))
    self.assertFalse(tab.IsEmpty('b', ignore_nan=False))
    self.assertFalse(tab.IsEmpty('c', ignore_nan=False))
    
  def testUnique(self):
    tab = self.CreateTestTable()
    tab.AddRow(['foo',4, 3.3])
    tab.AddRow([None,5, 6.3])
    self.assertEqual(tab.GetUnique('first'), ['x','foo'])
    self.assertEqual(tab.GetUnique('first', ignore_nan=False), ['x','foo', None])
    self.assertEqual(tab.GetUnique('second'), [3,9,4,5])
    self.assertEqual(tab.GetUnique('second', ignore_nan=False), [3,None,9,4,5])
    self.assertEqual(tab.GetUnique('third'), [2.2, 3.3, 6.3])
    self.assertEqual(tab.GetUnique('third', ignore_nan=False), [None, 2.2, 3.3, 6.3])
    
  def testCorrel(self):
    tab = self.CreateTestTable()
    self.assertEqual(tab.Correl('second','third'), None)
    tab.AddRow(['foo',4, 3.3])
    tab.AddRow([None,5, 6.3])
    tab.AddRow([None,8, 2])
    self.assertAlmostEqual(tab.Correl('second','third'), -0.4954982578)
    
  def testSpearmanCorrel(self):
    if not HAS_SCIPY_STATS:
      return
    tab = self.CreateTestTable()
    self.assertEqual(tab.SpearmanCorrel('second','third'), None)
    tab.AddRow(['foo',4, 3.3])
    tab.AddRow([None,5, 6.3])
    tab.AddRow([None,8, 2])
    self.assertAlmostEqual(tab.SpearmanCorrel('second','third'), -0.316227766)
    
  def testExtend(self):
    '''
     first  second  third 
    ----------------------
     x            3     NA
     foo         NA  2.200
     NA           9  3.300
    '''
    
    # simple extend of the same table
    tab = self.CreateTestTable()
    self.CompareDataFromDict(tab, {'first': ['x','foo',None],
                                   'second': [3,None,9],
                                   'third': [None,2.2,3.3]})
    
    tab.Extend(tab)
    self.CompareDataFromDict(tab, {'first': ['x','foo',None,'x','foo',None],
                                   'second': [3,None,9,3,None,9],
                                   'third': [None,2.2,3.3,None,2.2,3.3]})
    
    # simple extend of different tables with the same data
    tab = self.CreateTestTable()
    tab2 = self.CreateTestTable()
    tab.Extend(tab2)
    self.CompareDataFromDict(tab, {'first': ['x','foo',None,'x','foo',None],
                                   'second': [3,None,9,3,None,9],
                                   'third': [None,2.2,3.3,None,2.2,3.3]})
    self.CompareDataFromDict(tab2, {'first': ['x','foo',None],
                                    'second': [3,None,9],
                                    'third': [None,2.2,3.3]})
    
    # add additional columns to current table
    tab = self.CreateTestTable()
    tab2 = self.CreateTestTable()
    tab2.AddCol('foo','i',[1,2,3])
    tab.Extend(tab2)
    self.CompareDataFromDict(tab, {'first': ['x','foo',None,'x','foo',None],
                                   'second': [3,None,9,3,None,9],
                                   'third': [None,2.2,3.3,None,2.2,3.3],
                                   'foo': [None,None,None,1,2,3]})     
    
    # different order of the data
    tab = self.CreateTestTable()
    tab2 = Table(['third','second','first'],
                  'fis',
                  third=[None,2.2,3.3],
                  first=['x','foo',None],
                  second=[3, None, 9])
    self.CompareDataFromDict(tab2, {'first': ['x','foo',None],
                                    'second': [3,None,9],
                                    'third': [None,2.2,3.3]})
    tab.Extend(tab2)
    self.CompareDataFromDict(tab, {'first': ['x','foo',None,'x','foo',None],
                                   'second': [3,None,9,3,None,9],
                                   'third': [None,2.2,3.3,None,2.2,3.3]})
    
    # with overwrite (additional column)
    tab = self.CreateTestTable()
    tab2 = self.CreateTestTable()
    tab2.AddCol('foo','i',[1,2,3])
    tab.Extend(tab2, overwrite='first')
    self.CompareDataFromDict(tab, {'first': ['x','foo',None],
                                   'second': [3,None,9],
                                   'third': [None,2.2,3.3],
                                   'foo': [1,2,3]})
    
    # with overwrite (no matching value)
    tab = self.CreateTestTable()
    tab2 = Table(['third','second','first'],
                  'fis',
                  third=[None,2.2,3.3],
                  first=['a','bar','bla'],
                  second=[3, None, 9])
    tab.Extend(tab2, overwrite='first')
    self.CompareDataFromDict(tab, {'first': ['x','foo',None,'a','bar','bla'],
                                   'second': [3,None,9,3,None,9],
                                   'third': [None,2.2,3.3,None,2.2,3.3]})
    
    # with overwrite (with matching values)
    tab = self.CreateTestTable()
    tab2 = Table(['third','second','first'],
                  'fis',
                  third=[None,2.2,3.4],
                  first=['a','bar','bla'],
                  second=[3, None, 9])
    tab.Extend(tab2, overwrite='third')
    self.CompareDataFromDict(tab, {'first': ['a','bar',None,'bla'],
                                   'second': [3,None,9,9],
                                   'third': [None,2.2,3.3,3.4]})
    
    # cannot extend if types are different
    tab = Table('aaa','s',a=['a','b'])
    tab2 = Table('aaa','i',a=[1,2])
    self.assertRaises(TypeError, tab.Extend, tab2)

  def testSelection(self):

    #create table with some data
    tab = Table(['a','b','c','d'],'ifbs',
                a=[0,1,2,3,4,5,6,7,8,9],
                b=[0.5,1.0,1.1,1.2,1.3,1.4,1.5,2.3,2.5,3.0],
                c=[True,True,True,True,False,False,False,False,None,None],
                d=['a','b','c','d','e','f','g','h','i','j'])


    #check wether error gets raised in problematic cases like: a=1=b
    self.assertRaises(ValueError, tab.Select, 'a=1=b and c=True')

    #check wether errors get raised in stupid cases like: =1,2 and b=1.0
    self.assertRaises(ValueError, tab.Select, '=1,2 and c=True')
    self.assertRaises(ValueError, tab.Select, '1,2=')
    self.assertRaises(ValueError, tab.Select, '=1:2 and c=True')
    self.assertRaises(ValueError, tab.Select, '1:2=')

    #check wether errors get raised in parenthesis mismatches
    self.assertRaises(ValueError,tab.Select,'a=1,2 and ((b=1.0 or c=True)')
    self.assertRaises(ValueError,tab.Select,'a=1,2 and (b=1.0 or c=True))')

    #check wether error gets raised when two operands are not separated by an operator
    self.assertRaises(ValueError,tab.Select,'a=1 b=1.0')

    #check whether error gets raised when inconsistent types are compared
    self.assertRaises(ValueError,tab.Select,'c>d')
    self.assertRaises(ValueError,tab.Select,'c+d=5')

    from ost.table_selector import TableSelector

    selector=TableSelector(tab.col_types,tab.col_names,'')

    #check some examples for dijkstras shunting yard algorithm
    query_one='a=1:4 and ((b>5.0 or b<2.0) and c=True)'
    split_exp_one=selector._ExpressionLexer(query_one)
    parsed_exp_one=selector._ParseExpression(split_exp_one)
    rpn_one=[1, 'a', '<=', 'a', 4, '<=', 'and', 'b', 5.0, '>', 'b', 2.0, '<', 'or', 'c', True, '=', 'and', 'and']
    query_two='(a=1,2) or (b>5.0+a or b<2.0)'
    split_exp_two=selector._ExpressionLexer(query_two)
    parsed_exp_two=selector._ParseExpression(split_exp_two)
    rpn_two=['a', 1, '=', 'a', 2, '=', 'or', 'b', 5.0, 'a', '+', '>', 'b', 2.0, '<', 'or', 'or']
    self.assertEqual(selector._ShuntingYard(parsed_exp_one),rpn_one)
    self.assertEqual(selector._ShuntingYard(parsed_exp_two),rpn_two)

    #check operator evaluations
    self.assertTrue(selector._EvaluateOperator('=',False,False))
    self.assertFalse(selector._EvaluateOperator('=',False,True))
    self.assertTrue(selector._EvaluateOperator('and',True,True))
    self.assertFalse(selector._EvaluateOperator('and',True,False))
    self.assertTrue(selector._EvaluateOperator('or',True,False))
    self.assertTrue(selector._EvaluateOperator('or',False,True))
    self.assertFalse(selector._EvaluateOperator('or',False,False))
    self.assertTrue(selector._EvaluateOperator('!=',False,True))
    self.assertFalse(selector._EvaluateOperator('!=',True,True))
    self.assertTrue(selector._EvaluateOperator('<=',1.0,2.0))
    self.assertTrue(selector._EvaluateOperator('<=',1.0,1.0))
    self.assertFalse(selector._EvaluateOperator('<=',2.0,1.0))
    self.assertFalse(selector._EvaluateOperator('>=',1.0,2.0))
    self.assertTrue(selector._EvaluateOperator('>=',2.0,1.0))
    self.assertTrue(selector._EvaluateOperator('>=',1.0,1.0))
    self.assertTrue(selector._EvaluateOperator('<',1.0,2.0))
    self.assertFalse(selector._EvaluateOperator('<',2.0,1.0))
    self.assertFalse(selector._EvaluateOperator('<',1.0,1.0))
    self.assertFalse(selector._EvaluateOperator('>',1.0,2.0))
    self.assertTrue(selector._EvaluateOperator('>',2.0,1.0))
    self.assertFalse(selector._EvaluateOperator('>',1.0,1.0))
    self.assertEqual(selector._EvaluateOperator('+',1,1),2)
    self.assertEqual(selector._EvaluateOperator('-',1,1),0)
    self.assertEqual(selector._EvaluateOperator('*',2,2),4)
    self.assertEqual(selector._EvaluateOperator('/',2,2),1)
    self.assertEqual(selector._EvaluateOperator('+',None,5),None)
    self.assertEqual(selector._EvaluateOperator('-',2,None),None)
    self.assertEqual(selector._EvaluateOperator('/',2,None),None)
    self.assertEqual(selector._EvaluateOperator('*',None,3),None)

    #check a few rpn evaluation examples

    rpn_one=[1, 2, '<=', 2, 4, '<=', 'and', 6.0, 5.0, '>', 2.0, 2.0, '<', 'or', True, True, '=', 'and', 'and']
    rpn_two=[1, 0, '<=', 2, 4, '<=', 'and', 6.0, 5.0, '>', 2.0, 2.0, '<', 'or', True, True, '=', 'and', 'and']
    rpn_three=[1, 0, '<=', 2, 4, '<=', 'or', 6.0, 5.0, '>', 2.0, 2.0, '<', 'or', True, True, '=', 'and', 'and']

    self.assertTrue(selector._EvaluateRPN(rpn_one))
    self.assertFalse(selector._EvaluateRPN(rpn_two))
    self.assertTrue(selector._EvaluateRPN(rpn_three))

    #check a few selection examples

    query_one='c'
    query_two='a=2:9 and c'
    query_three='d=e,f,j and a!=4'
    query_four='(b=2.0:3.0 and c=False) or (b<1.0 and c)'
    query_five='b=0.0:1.2,2.5:8.0'
    query_six='c=None'
    query_seven='c!=None'
    query_eight='b/0.5=2'

    self.assertEqual([0,1,2,3], list(r[0] for r in tab.Select(query_one).rows))
    self.assertEqual([2,3], list(r[0] for r in tab.Select(query_two).rows))
    self.assertEqual([5,9], list(r[0] for r in tab.Select(query_three).rows))
    self.assertEqual([0,7], list(r[0] for r in tab.Select(query_four).rows))
    self.assertEqual([0,1,2,3,8,9], list(r[0] for r in tab.Select(query_five).rows))
    self.assertEqual([8,9], list(r[0] for r in tab.Select(query_six).rows))
    self.assertEqual([0,1,2,3,4,5,6,7], list(r[0] for r in tab.Select(query_seven).rows))
    self.assertEqual([1], list(r[0] for r in tab.Select(query_eight).rows))


    
if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
