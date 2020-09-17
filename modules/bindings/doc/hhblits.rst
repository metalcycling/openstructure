:mod:`~ost.bindings.hhblits` - Search related sequences in databases
================================================================================

Introduction
--------------------------------------------------------------------------------

HHblits is a sequence search tool like BLAST but able to find more distant
homologs. This is achieved by aligning hidden Markov models (HMM)
in the search process as opposed to `sequence-sequence` searches in BLAST.
HHblits works on a HMM database, usually that one is provided, queried with 
a HMM representing your target sequence. The latter one needs to be calculated 
before the actual search. The software suite needed for HHblits can be found on
`github <https://github.com/soedinglab/hh-suite>`_.
Alternatively, the deprecated HHblits 2.x suite can be found here:
`here <http://wwwuser.gwdg.de/~compbiol/data/hhsuite/releases/all/>`_.

On HHblits Versions
--------------------------------------------------------------------------------

The binding for HHblits 3 has internally been forked from the HHblits 2 binding. 
The binding for HHblits 2 is considered deprecated and doesn't receive bugfixes
anymore. Also the documentation refers to the HHblits 3 binding. The different
bindings can be imported explicitely:

.. code-block:: python

  from ost.bindings import hhblits2  
  from ost.bindings import hhblits3

Alternatively you can let OpenStructure figure out the HHblits version you're
using and import the appropriate binding for you under the base name hhblits. 
This assumes the HHblits binary (hhblits) to be in your path and raises an error 
otherwise.

.. code-block:: python

  from ost.bindings import hhblits


Examples
--------------------------------------------------------------------------------

A typical search: Get an instance of the binding, build the search profile out
of the query sequence, run the search and iterate results. 

.. code-block:: python

  from ost.bindings import hhblits3

  # Create a SequenceHandle, alternatively you can load any sequence in 
  # FASTA format using ost.io.LoadSequence(<PATH_TO_SEQUENCE>)
  query_seq = seq.CreateSequence('Query',
                                 'TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN')

  # set up the search environment
  # lets assume a default installation with hhblits binary at
  # <PATH_TO_HHBLITS_INSTALL>/bin/hhblits
  hh = hhblits3.HHblits(query_seq, '<PATH_TO_HHBLITS_INSTALL>')

  # now create a search profile for the query sequence against uniclust30 
  # which you can get with instructions in the hh-suite user guide (github)
  # <PATH_TO_DB>/uniclust30_2018_08 is just the prefix common to
  # all db files, so `ls <PATH_TO_DB>/uniclust30_2018_08*` would list all 
  # of them
  a3m_file = hh.BuildQueryMSA(nrdb='<PATH_TO_DB>/uniclust30_2018_08')

  # lets load the data in the a3m_file and display the generated
  # multiple sequence alignment note that ParseA3M is not a class method 
  # but a module function
  a3m_data = hhblits3.ParseA3M(open(a3m_file))
  print(a3m_data['msa'])

  # search time! we just search against uniclust again, but every HHblits db is
  # working here, e.g. one build from all the sequences in PDB
  hit_file = hh.Search(a3m_file, '<PATH_TO_DB>/uniclust30_2018_08')

  # lets have a look at the resuls
  with open(hit_file) as hit_fh:
      header, hits = hhblits3.ParseHHblitsOutput(hit_fh)
  for hit in hits:
      print(hit.aln)

  # cleanup
  hh.Cleanup()



Binding API
--------------------------------------------------------------------------------

.. automodule:: ost.bindings.hhblits3
   :synopsis: Search related sequences in databases
   :members:

..  LocalWords:  HHblits homologs
