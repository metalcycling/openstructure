:mod:`~ost.table` - Working with tabular data
================================================================================

.. module:: ost.table
  :synopsis: Working with tabular data
  
This module defines the table class that provides convenient functionality to
work with tabular data. It features functions to calculate statistical moments,
e.g. mean, standard deviations as well as functionality to plot the data using
matplotlib.

Basic Usage
--------------------------------------------------------------------------------

Populate table with data and plot the data:

.. code-block:: python

  from ost.table import *

  # create table with two columns, x and y both of float type
  tab = Table(['x', 'y'], 'ff')
  for x in range(1000):
    tab.AddRow([x, x**2])

  # create a plot
  plt = tab.Plot('x', 'y')

  # save resulting plot to png file
  plt.savefig('x-vs-y.png')

Iterating over table items:

.. code-block:: python

  # load table from file
  tab = Table.Load(...)

  # get column index for col 'foo'
  idx = tab.GetColIndex('foo')

  # iterate over all rows
  for row in tab.rows:
    # print complete row
    print(row)

    # print value for column 'foo'
    print(row[idx])

  # iterate over all rows of selected columns
  for foo, bar in tab.Zip('foo','bar'):
    print(foo, bar)

Doing element wise mathematical operations on entire colums:

.. code-block:: python

  # create table with two columns, foo and bar both of int type
  # and fill with values
  tab = Table(['foo', 'bar'], 'ii', foo=[1,2,3,4], bar=[1,4,9,16])

  # add new column by doing an element wise
  # addition of column foo and column bar
  tab.AddCol('qux', 'f', tab['foo']+tab['bar'])

  print(tab)


Select part of the table based on a query:

.. code-block:: python

  # create table with two columns, foo and bar both of int type
  # and fill with values
  tab = Table(['foo', 'bar'], 'ii', foo=[1,2,3,4], bar=[1,4,9,16])

  # select all rows where foo>=2 and bar<10
  subtab = tab.Select('foo>=2 and bar<10')
  print(subtab)

  # select all rows where foo>3 or bar=1
  subtab = tab.Select('foo>3 or bar=1')
  print(subtab)


Functions You Might be Interested In
--------------------------------------------------------------------------------
============================================= ============================================
**Adding/Removing/Reordering data**
:meth:`~ost.table.Table.AddRow`               add a row to the table
:meth:`~ost.table.Table.AddCol`               add a column to the table
:meth:`~ost.table.Table.RemoveCol`            remove a column from the table
:meth:`~ost.table.Table.RenameCol`            rename a column
:meth:`~ost.table.Table.Extend`               append a table to the end of another table
:meth:`~ost.table.Merge`                      merge two tables together
:meth:`~ost.table.Table.Sort`                 sort table by column
:meth:`~ost.table.Table.Filter`               filter table by values
:meth:`~ost.table.Table.Select`               select subtable based on query
:meth:`~ost.table.Table.Zip`                  extract multiple columns at once
:meth:`~ost.table.Table.SearchColNames`       search for matching column names

**Input/Output**
:meth:`~ost.table.Table.Save`                 save a table to a file
:meth:`~ost.table.Table.Load`                 load a table from a file
:meth:`~ost.table.Table.ToString`             convert a table to a string for printing

**Simple Math**
:meth:`~ost.table.Table.Min`                  compute the minimum of a column
:meth:`~ost.table.Table.Max`                  compute the maximum of a column
:meth:`~ost.table.Table.Sum`                  compute the sum of a column
:meth:`~ost.table.Table.Mean`                 compute the mean of a column
:meth:`~ost.table.Table.RowMean`              compute the mean for each row
:meth:`~ost.table.Table.Median`               compute the median of a column
:meth:`~ost.table.Table.StdDev`               compute the standard deviation of a column
:meth:`~ost.table.Table.Count`                compute the number of items in a column

**More Sophisticated Math**
:meth:`~ost.table.Table.Correl`               compute Pearson's correlation coefficient
:meth:`~ost.table.Table.SpearmanCorrel`       compute Spearman's rank correlation coefficient
:meth:`~ost.table.Table.ComputeMCC`           compute Matthew's correlation coefficient
:meth:`~ost.table.Table.ComputeROC`           compute receiver operating characteristics (ROC)
:meth:`~ost.table.Table.ComputeEnrichment`    compute enrichment
:meth:`~ost.table.Table.GetOptimalPrefactors` compute optimal coefficients for linear combination of columns

**Plot**
:meth:`~ost.table.Table.Plot`                 Plot data in 1, 2 or 3 dimensions
:meth:`~ost.table.Table.PlotHistogram`        Plot data as histogram
:meth:`~ost.table.Table.PlotROC`              Plot receiver operating characteristics (ROC)
:meth:`~ost.table.Table.PlotEnrichment`       Plot enrichment
:meth:`~ost.table.Table.PlotHexbin`           Hexagonal density plot
:meth:`~ost.table.Table.PlotBar`              Bar plot


============================================= ============================================

The Table class
--------------------------------------------------------------------------------


.. autoclass:: ost.table.Table
  :members:
  :undoc-members: SUPPORTED_TYPES

.. autofunction:: ost.table.Merge
