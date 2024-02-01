IO Profiles for entity importer
================================================================================

.. currentmodule:: ost.io

As of version 1.1, OpenStructure introduces IO profiles to fine-tune the
behaviour of the molecule importers. A profile aggregates flags and methods
that affect the import of molecular structures and influence both the behaviour
of :mod:`~ost.conop` and :mod:`~ost.io`.

Basic usage of IO profiles
--------------------------------------------------------------------------------

You are most certainly reading this document because you were having trouble
loading PDB files. In that case, as a first step you will want to set the
profile parameter of  :func:`LoadPDB`. The profile parameter can either be the
name of a profile or an instance of :class:`IOProfile`. Both of the following
two examples are equivalent:

.. code-block:: python

  ent = io.LoadPDB('weird.pdb', profile=io.profiles['SLOPPY'])
  ent = io.LoadPDB('weird.pdb', profile='SLOPPY')

Profiles is a dictionary-like object containing all the profiles known to
OpenStructure. You can add new ones by inserting them into the dictionary.
If you are loading a lot of structures, you may want to set the default profile
to avoid having to pass the profile every time you load a structure.
This is done by assigning a different profile to ``DEFAULT``:

.. code-block:: python

  io.profiles['DEFAULT']='SLOPPY'
  ent = io.LoadPDB('weird.pdb')

Again, you can either assign the name of the profile, or the profile itself.
If none of the profiles available by default suits your needs, feel free to
create one to your liking.

Available default profiles
--------------------------------------------------------------------------------

The following profiles are available by default. For a detailed description of
what the different parameters mean, consult the documentation of
:class:`IOProfile`.

STRICT

  This profile is the default (also available as DEFAULT) and is known to
  work very well with PDB files coming from the official PDB website. It
  is equivalent to the following profile:

  .. code-block:: python

    IOProfile(dialect='PDB', fault_tolerant=False, quack_mode=False,
              processor=conop.RuleBasedProcessor(conop.GetDefaultLib()))

SLOPPY:

  This profile loads essentially everything

  .. code-block:: python

    IOProfile(dialect='PDB', fault_tolerant=True, quack_mode=True,
              processor=conop.RuleBasedProcessor(conop.GetDefaultLib()))

CHARMM:

  This format is the default when importing CHARMM trajectories and turns on the
  CHARMM specific compound dictionary.

  .. code-block:: python

    IOProfile(dialect='CHARMM', fault_tolerant=True, quack_mode=False,
              processor=conop.RuleBasedProcessor(conop.GetDefaultLib()))

.. note:: 

  The profiles are setup at the first import of the io module, i.e. something
  like  ``from ost import io`` or ``from ost.io import LoadPDB``. The processor
  parameter is set as stated above IF :func:`ost.conop.GetDefaultLib()` returns
  a valid compound library at that point in time. If not, the processor is set
  to :class:`ost.conop.HeuristicProcessor()`. Calling
  :func:`ost.conop.SetDefaultLib()` has thus no immediate effect on the default
  profiles! Two exceptions: :func:`ost.io.LoadPDB` and
  :class:`ost.io.LoadMMCIF()` have a logic in place to override the processor of
  the default profiles with :func:`ost.conop.GetDefaultLib`, using
  :class:`HeuristicProcessor` respectively. This logic does not apply to user
  defined profiles. 


The IOProfile Class
--------------------------------------------------------------------------------

.. class:: IOProfile(dialect='PDB', quack_mode=False, fault_tolerant=False,\
                     join_spread_atom_records=False, no_hetatms=False,\
                     calpha_only=False, read_conect=False, processor=None)

  Aggregates flags that control the import of molecular structures.

  .. attribute:: dialect

    :type: str

    The dialect to be used for PDB files. At the moment, this is either CHARMM
    or PDB. More will most likely come in the future. By setting the dialect to
    CHARMM, the loading is optimized for CHARMM PDB files. This turns on
    support for chain names with length up to 4 characters (column 72-76) and
    increase the size of the residue name to 4 residues.

  .. attribute:: quack_mode

    :type: bool

    Read/write property. When quack_mode is enabled, the chemical class for
    unknown residues is guessed based on its atoms and connectivity. Turn this
    on if you are working with non-standard conforming PDB files and are
    experiencing problems with the rendering of the backbone trace and/or see
    peptidic residues with unknown chemical classes.

  .. attribute:: fault_tolerant

    :type: bool

    If true, the import will succeed, even if the PDB contains faulty records.
    The faulty records will be ignored and import continues as if the records
    are not present.

  .. attribute::   join_spread_atom_records

    :type: bool
  
    If set to true, atom records belonging to the same residue are joined, even 
    if they do not appear sequentially in the PDB file.

  .. attribute:: no_hetatms

    :type: bool
  
    If set to true, HETATM records are ignored during import.

  .. attribute:: calpha_only

    :type: bool

    When set to true, forces the importer to only load atoms named CA. This is
    most useful in combination with protein-only PDB files to speed up
    subsequent processing and importing.

  .. attribute:: read_conect

    :type: bool

    Only relevant when reading files in PDB format. When set to true, reads CONECT
    statements and applies them in the PDB reader. This can result in
    hydrogen bonds, salt bridges etc. to be connected. Check the PDB format
    definition for more info. This may cause issues in subsequent processing,
    such as bonds being overriden, or extra, inconsistent bonds, as the
    processor suddenly has two separate sources of connectivity.
    For the use case where the input PDB file contains valid CONECT
    statements for all hetatms, you may want to disable processing of bonds
    between them in :attr:`ost.conop.Processor.connect_hetatm`

  .. attribute:: processor

    :type: :class:`ost.conop.HeuristicProcessor`/:class:`ost.conop.RuleBasedProcessor`

    Controls connectivity processing of loaded :class:`ost.mol.EntityHandle`.
    Even though its a keyword argument, processing will fail if not given.
