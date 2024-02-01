..  Note on large code blocks: keep max. width to 100 or it will look bad
                               on webpage!

=========================
Molecular Checker (Molck)
=========================

--------------------------------------
Where can I find the Molck executable? 
--------------------------------------

The Molck executable can be found at <YOUR-OST-STAGE-DIR>/bin

-----------
Basic Usage 
-----------

To check one PDB file (struc1.pdb) with Molck, use the following command:

.. code-block:: bash

    molck --complib <PATH TO COMPOUND LIB> struc1.pdb

The checked and cleaned file will be saved by default ad struc1-molcked.pdb.

Similarly it is possible to check a list of PDB files:

.. code-block:: bash

    molck --complib <PATH TO COMPOUND LIB> struc1.pdb struc2.pdb struc3.pdb


-----------
All Options 
-----------

The molck executable supports several other command line options,
please find them following:

.. code-block:: bash

    Usage: molck [options] file1.pdb [file2.pdb [...]]
    Options
      --complib=path       Location of the compound library file. If not provided,
                           the following locations are searched in this order:
                           1. Working directory,
                           2. OpenStructure standard library location (if the
                              executable is part of a standard OpenStructure
                              installation)
      --rm=<a>,<b>         Remove atoms and residues matching some criteria:
                           - zeroocc - Remove atoms with zero occupancy
                           - hyd - Remove hydrogen atoms
                           - oxt - Remove terminal oxygens
                           - nonstd - Remove all residues not one of the
                                      20 standard amino acids
                           - unk - Remove unknown atoms not following
                                   the nomenclature
                           Default: hyd
      --fix-ele            Clean up element column
      --stdout             Write cleaned file(s) to stdout
      --out=filename       Write cleaned file(s) to disk. % characters in the
                           filename are replaced with the basename of the input
                           file without extension. Default: %-molcked.pdb
      --color=auto|on|off  Whether output should be colored. Delault: auto
      --map-nonstd         Maps modified residues back to the parent amino acid,
                           for example: MSE -> MET, SEP -> SER.

================
Molck Python API
================

Within OST, one can also call the :func:`~ost.mol.alg.Molck` function directly
on entities to get the same effect as with the binary.
