:mod:`~ost.io` - Input and Output of Sequences, Structures and Maps
================================================================================

.. toctree::
  :hidden:
  
  structure_formats
  sequence_formats
  sequence_profile_formats
  image_formats
  mmcif
  profile

.. module:: ost.io
  :synopsis: Input and output of sequences, alignments, structures, images and density maps.

The io module deals with the input and output of :class:`entities 
<ost.mol.EntityHandle>`, :class:`alignments <ost.seq.AlignmentHandle>`, 
:class:`sequences <ost.seq.SequenceHandle>`, :class:`images 
<ost.img.ImageHandle>`. Importers for common file formats containing molecules 
such as PDB, SDF and CHARMM trajectory files are available. Sequence and 
alignment file formats such as FASTA and CLUSTALW are supported as well as  
various image data (e.g. png, dm3) and density map files (e.g. CCP4, MRC). 

.. _mol-io:

Molecular Structures
--------------------------------------------------------------------------------

Loading Molecular Structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :mod:`~ost.io` modules offer several ways to load molecular structures 
depending on your requirements. The most general way is offered by 
:func:`~ost.io.LoadEntity`, which will automatically detect the file format based 
on the file extension.

.. function:: LoadEntity(filename, format='auto')

  Load entity from disk. If format is set to 'auto', the function guesses the 
  filetype based on the extension of the file e.g. files ending in '.pdb', '.ent', 
  '.ent.gz' and '.pdb.gz' will automatically be loaded as PDB files. For files 
  without or exotic extensions, the format can be set explicitly as the second 
  parameter. 
  
  .. code-block:: python

    # Recognizes SDF file by file extension
    ent = io.LoadEntity('file.sdf')

    # In this case, there is no file extensions, so you have to say it's a 
    # SDF file explicitly
    ent = io.LoadEntity('file', 'sdf')

  For a list of file formats supported by :func:`LoadEntity`, see :doc:`structure_formats`.
  
  :raises: :exc:`~ost.io.IOUnknownFormatException` if the format string supplied 
      is not recognized or the file format can not be detected based on the 
      file extension.
      
      :exc:`~ost.io.IOException` if the import fails due to an erroneous or 
      inexistent file.

Some of the formats have a dedicated function that allows you to tweak many
parameters that affect the import. PDB files can be loaded with  :func:`LoadPDB`
and mmCIF files with :func:`LoadMMCIF` (this also gives you access to the
:class:`MMCifInfo` class). It offers tighter control over the exact loading
behaviour.


.. autofunction:: ost.io.LoadPDB

.. function:: PDBStrToEntity(pdb_string, profile=IOProfile(), process=False)

  Load entity from a string in PDB format. By default the entity is loaded with
  an empty IO Profile and is not processed with the :class:`~ost.conop.Processor`,
  even if one is available in the IO Profile.

  :param pdb_string: A PDB file as a string.

  :param profile: The IO Profile to read the entity with. For more information
      on the IO Profiles available, see :doc:`profile`.

  :param process: If set to True, run the :class:`~ost.conop.Processor`
      contained in the IO Profile.

  :rtype: :class:`~ost.mol.EntityHandle`.

  To get an entity equivalent to one loaded with :func:`LoadPDB`, set the
  `profile` and `process` arguments as follows:

  .. code-block:: python

    with open('protein.pdb') as pdb_fd:
        pdb_str = pdb.read()
        ent = io.PDBStrToEntity(pdb_str, ost.io.profiles['DEFAULT'], True)

Loading Molecular Structures From Remote Repositories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`LoadPDB` already provides access to selected structural databases in 
the internet when enabling the *remote* flag. Predefined 
:class:`ost.io.remote.RemoteRepository` objects are available as

.. code-block:: python

  from ost.io import remote
  repo_name = 'smtl'
  repo = remote.REMOTE_REPOSITORIES.get(repo_name)

  # url for entry with id 1ake.1
  print(repo.URLForID('1ake.1'))

where *repo_name* can be one of ['pdb', 'cif', 'pdb_redo', 'smtl'].
Instead of explicit access, you can directly fetch data using:

.. autofunction:: ost.io.remote.RemoteGet
 
.. autofunction:: ost.io.remote.RemoteLoad

.. autoclass:: ost.io.remote.RemoteRepository
  :members:

Saving Molecular Structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Saving a complete entity or a view is a matter of calling 
:func:`~ost.io.SaveEntity`.

.. code-block:: python
  
  ent = io.LoadEntity('protein.pdb')
  # save full entity
  io.SaveEntity(ent, 'full.pdb')
  # only save C-alpha atoms
  io.SaveEntity(ent.Select('aname=CA and peptide=true'), 'calpha.pdb')
  
:func:`~ost.io.SavePDB` provides a simple way to save several entities into one 
file:


.. code-block:: python
  
  ent = io.LoadEntity('protein.pdb')
  # Save complete entity
  io.SavePDB(ent, 'full.pdb')
  # Save chain A and chain B separately
  io.SavePDB([ent.Select('cname=A'), ent.Select('cname=B')], 'split.pdb')


.. function:: SaveEntity(ent, filename, format='auto')
  
  Save entity to disk. If format is set to 'auto', the function guesses the 
  filetype based on the file extension, otherwise the supplied format is checked 
  against the available export plugins.
  
  :param ent: The entity to be saved
  :type  ent: :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`
  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string
  
  :raises: :exc:`~ost.io.IOUnknownFormatException` if the format string supplied 
      is not recognized or the file format can not be detected based on the 
      file extension.
      
.. autofunction:: ost.io.SavePDB

.. function:: EntityToPDBStr(ent, profile=IOProfile())

  Return entity as a string in PDB format. 

  :param entity: The :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`

  :param profile: The IO Profile to read the entity with. For more information
      on the IO Profiles available, see :doc:`profile`.
  :raises: IOException if the restrictions of the PDB format are not satisfied
           (see :meth:`ost.io.SavePDB`)

  :rtype: string.

.. _seq-io:

Sequences and Alignments
--------------------------------------------------------------------------------

Loading sequence or alignment files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: LoadSequence(filename, format='auto')

  Load sequence data from disk. If format is set to 'auto', the function guesses 
  the filetype based on the extension of the file. Files ending in '.fasta', 
  '.aln' will automatically be loaded.

  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string
  :rtype: :class:`~ost.seq.SequenceHandle`
  
  For files with non-standard extensions, the format can be set explicitly 
  specifying the `format` parameter. 
  
  .. code-block:: python

    # recognizes FASTA file by file extension
    myseq = io.LoadSequence('seq.fasta')
    # for obtaining a SequenceList
    seqlist = io.LoadSequenceList('seqs.fasta')
    # or for multiple alignments (here from CLUSTAL)
    aln = io.LoadAlignment('algnm.aln', format="clustal")
    
  For a list of file formats supported by :func:`LoadSequence` see
  :doc:`sequence_formats`.
  
  :raises: :exc:`~ost.io.IOUnknownFormatException` if the format string supplied 
      is not recognized or the file format can not be detected based on the 
      file extension.
      
      :exc:`~ost.io.IOException` if the import fails due to an erroneous or 
      inexistent file.

.. function:: LoadSequenceList(filename, format='auto')

  For a description of how to use :func:`LoadSequenceList` please refer to 
  :func:`LoadSequence`. For a list of file formats supported by
  :func:`LoadSequenceList` see :doc:`sequence_formats`.

  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string
  :rtype: :class:`~ost.seq.SequenceList`

.. function:: LoadAlignment(filename, format='auto')

  For a description of how to use :func:`LoadAlignment` please refer to 
  :func:`LoadSequence`. For a list of file formats supported by 
  :func:`LoadAlignment` see :doc:`sequence_formats`.

  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string
  :rtype: :class:`~ost.seq.AlignmentHandle`

.. function:: LoadSequenceProfile(filename, format='auto')

  Load sequence profile data from disk. If format is set to 'auto', the function
  guesses the filetype based on the extension of the file. Files ending in
  '.hhm' (output of HHblits) and '.pssm' (ASCII Table (PSSM) output of PSI-BLAST
  as generated with blastpgp and flag -Q) will automatically be loaded.
  
  For files with non-standard extensions, the format can be set explicitly 
  specifying the `format` parameter. 
  
  .. code-block:: python

    # recognizes hhm file by file extension
    myprof = io.LoadSequenceProfile('myhmm.hhm')
    # recognizes pssm file by file extension
    myprof = io.LoadSequenceProfile('myprof.pssm')

    # to override format
    myprof = io.LoadSequenceProfile('myfile', format='hhm')
    myprof = io.LoadSequenceProfile('myfile', format='pssm')
    
  For a list of file formats supported by :func:`LoadSequenceProfile` see
  :doc:`sequence_profile_formats`.

  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string
  
  :rtype: :class:`~ost.seq.ProfileHandle`

  :raises: :exc:`~ost.io.IOUnknownFormatException` if the format string supplied 
      is not recognized or the file format can not be detected based on the 
      file extension.
      
      :exc:`~ost.io.IOException` if the import fails due to an erroneous or 
      inexistent file.

.. function:: AlignmentFromString(data, format)

  Load alignment from string.

  The format argument is mandatory. For a list of supported formats,
  see :doc:`sequence_formats`.

  :param data: The alignment
  :type  data: string
  :param format: Name of the format
  :type  format: string
  :rtype: :class:`~ost.seq.AlignmentHandle`

.. function:: SequenceFromString(data, format)

  Load sequence from string.

  The format argument is mandatory. For a list of supported formats,
  see :doc:`sequence_formats`.

  :param data: The sequence
  :type  data: string
  :param format: Name of the format
  :type  format: string
  :rtype: :class:`~ost.seq.SequenceHandle`

.. function:: SequenceListFromString(data, format)

  Load a list of sequences from string.

  The format argument is mandatory. For a list of supported formats,
  see :doc:`sequence_formats`.

  :param data: The list of sequences
  :type  data: string
  :param format: Name of the format
  :type  format: string
  :rtype: :class:`~ost.seq.SequenceList`

.. function:: SequenceProfileFromString(data, format)

  Load sequence profile from string.

  The format argument is mandatory.

  :param data: String containing the data you would read from a file with
               specified format.
  :type  data: string
  :param format: Name of the format. Can either be "hhm" or "pssm".
  :type  format: string
  :rtype: :class:`~ost.seq.ProfileHandle`

Saving Sequence Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SaveSequence(sequence, filename, format='auto')

  Saving sequence data is performed by calling :func:`SaveSequence`.
  For files with non-standard extensions, the format can be set explicitly 
  specifying the 'format' parameter. 
  
  .. code-block:: python

    # recognizes FASTA file by file extension
    io.SaveSequence(myseq, 'seq.fasta')
    # for saving a SequenceList
    io.SaveSequenceList(seqlist, 'seqlist.fasta')
    # or for multiple alignments (here in FASTA format)
    io.SaveAlignment(aln, 'aln.fasta')
    
  For a list of file formats supported by :func:`SaveSequence` see
  :doc:`sequence_formats`.

  :param sequence: The sequence
  :type  sequence: :class:`~ost.seq.SequenceHandle`
  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string

  :raises: :exc:`~ost.io.IOUnknownFormatException` if the format string supplied 
      is not recognized or the file format can not be detected based on the 
      file extension.
      
      :exc:`~ost.io.IOException` if the import fails due to an erroneous or 
      inexistent file.

.. function:: SaveSequenceList(seq_list, filename, format='auto')

  For a desription of how to use :func:`SaveSequenceList` please refer to 
  :func:`SaveSequence`. For a list of file formats supported by 
  :func:`SaveSequenceList` see :doc:`sequence_formats`.

  :param seq_list: The sequence list
  :type  seq_list: :class:`~ost.seq.SequenceList`
  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string

.. function:: SaveAlignment(aln, filename, format='auto')

  For a desription of how to use :func:`SaveAlignment` please refer to 
  :func:`SaveSequence`.
  
  For a list of file formats supported by :func:`SaveAlignment` see 
  :doc:`sequence_formats`.

  :param aln: The alignment
  :type  aln: :class:`~ost.seq.AlignmentHandle`
  :param filename: The filename
  :type  filename: string
  :param format: Name of the format
  :type  format: string

.. function:: AlignmentToString(ali, format)

  Return alignment as a string.

  The format argument is mandatory. For a list of supported formats
  see :doc:`sequence_formats`.

  :param ali: The alignment
  :type  ali: :class:`~ost.seq.AlignmentHandle`
  :param format: Name of the format
  :type  format: string
  :rtype: string

.. function:: SequenceToString(seq, format)

  Return sequence as a string.

  The format argument is mandatory. For a list of supported formats
  see :doc:`sequence_formats`.

  :param seq: The sequence
  :type  seq: :class:`~ost.seq.SequenceHandle`
  :param format: Name of the format
  :type  format: string
  :rtype: string

.. function:: SequenceListToString(seq_list, format)

  Return sequence list as a string.

  The format argument is mandatory. For a list of supported formats
  see :doc:`sequence_formats`.

  :param seq_list: The sequence list
  :type  seq: :class:`~ost.seq.SequenceList`
  :param format: Name of the format
  :type  format: string
  :rtype: string

.. _img-io: 

Density Maps
--------------------------------------------------------------------------------

Loading Density Maps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: LoadImage(filename)

  Load density map from disk with the extension being guessed by the function.

  :param filename: The filename
  :type  filename: string

.. function:: LoadImage(filename, format)

  Load density map from disk. If no format is given, the function guesses the 
  filetype based on the extension of the file. If the extension is unknown or not present the
  filetype will be guessed based on the content of the file if possible.
  
  :param filename: The filename
  :type  filename: string
  :param format: The file format

  :raises: :exc:`~ost.io.IOUnknownFormatException` if the format supplied 
      is not recognized or the file format can not be detected based on the 
      file extension and content.
      
      :exc:`~ost.io.IOException` if the import fails due to an erroneous or 
      inexistent file.

  .. code-block:: python

    # recognizes mrc file by file extension
    ent = io.LoadImage('file.mrc')

    # it is always possible to explicitly set the image format
    # DAT file explicitly
    ent = io.LoadImage('file', Dat())

  For a list of file formats supported by :func:`LoadImage`, see :doc:`image_formats`.
  



Saving Density Maps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: SaveImage(image,filename)

  Save density map to disk with the function guessing the 
  filetype based on the file extension.

.. function:: SaveImage(image,filename, format)

  Save density map to disk. If no format is set, the function guesses the 
  filetype based on the file extension.
  
  :param image: The density map to be saved
  :type  image: :class:`~ost.img.IMageHandle`
  :param filename: The filename
  :type  filename: string
  :param format: The file format
  
  :raises: :exc:`~ost.io.IOUnknownFormatException` if the file format can not be detected based on the 
      file extension

  For a list of file formats supported by :func:`SaveImage`, see :doc:`image_formats`.
  
  .. code-block:: python

    # load density map
    image = io.LoadImage('density_map.ccp4')
    # save density map
    io.SaveImage(image, 'new_map.map', CCP4())


Stereochemical Parameters
--------------------------------------------------------------------------------

In order to check the structure for some stereo-chemical and steric clashes
before computing the lDDT scores it is required to pass parameter file based on
Engh and Huber parameters, and on the atomic radii as defined in the Cambridge
Structural Database. OpenStructure ships with default file called
`stereo_chemical_props.txt` located in `$OST_ROOT/share/openstructure`
directory. A function :func:`~ost.io.ReadStereoChemicalPropsFile` is used to
read this file.



.. function:: ReadStereoChemicalPropsFile(filename="", check=True)

  Read stereochemical parameters - if not provided a local version will be used.

  :param filename: The path to the parameter file that will be used. If set
                   to "", it reads the default file shipped with OpenStructure.
  :type filename: :class:`str`
  :param check: Raise an error when any of the resulting tables are empty.
  :type check:  :class:`bool`
  :return: Object containing stereochemical parameters
  :rtype: :class:`~ost.mol.alg.StereoChemicalProps`

.. function:: GetStereoChemicalPropsFile()

  Get the default path to the stereochemical paramteres file.
