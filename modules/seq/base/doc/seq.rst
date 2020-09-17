:mod:`~ost.seq` -- Sequences and Alignments
================================================================================

.. module:: ost.seq
   :synopsis: Contains classes and functions to deal with sequences and 
              alignments

The :mod:`seq` module helps you working with sequence data of various kinds. It 
has classes for :class:`single sequences <SequenceHandle>`, :class:`lists of 
sequences <SequenceList>` and :class:`alignments <AlignmentHandle>` of two or
more sequences. 


.. _attaching-views:

Attaching Structures to Sequences
--------------------------------------------------------------------------------


As OpenStructure is a computational structural biology framework, it is not
surprising that the sequence classes have been designed to work together with
structural data. Each sequence can have an attached :class:`~ost.mol.EntityView`
allowing for fast mapping between residues in the entity view and position in
the sequence.

.. _sequence-offset:

Sequence Offset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When using sequences and structures together, often the start of the structure 
and the beginning of the sequence do not fall together. In the following case, 
the alignment of sequences B and C only covers a subsequence of structure A::

  A acefghiklmnpqrstuvwy
  B     ghiklm
  C     123-45
  
We would now like to know which residue in protein A is aligned to which residue 
in sequence C. This is achieved by setting the sequence offset of sequence C to 
4. In essence, the sequence offset influences all the mapping operations from 
position in the sequence to residue index and vice versa. By default, the 
sequence offset is 0.

Loading and Saving Sequences and Alignments
--------------------------------------------------------------------------------

The :mod:`~ost.io` module supports input and output of common sequence formats.
Single  sequences can be loaded from disk with :func:`~ost.io.LoadSequence`,
alignments are loaded with :func:`~ost.io.LoadAlignment` and lists of sequences
are loaded with :func:`~ost.io.LoadSequenceList`. In addition to the file based
input  methods, sequences can also be loaded from a string:

.. code-block:: python

  seq_string = '''>sequence
  abcdefghiklmnop'''
  s = io.SequenceFromString(seq_string, 'fasta')
  print(s.name, s) # will print "sequence abcdefghiklmnop"
  
Note that, in that case specifying the format is mandatory.

The SequenceHandle
--------------------------------------------------------------------------------

.. function:: CreateSequence(name, sequence, role="UNKNOWN")

  Create a new :class:`SequenceHandle` with the given name and sequence. 

  :param name: Name of the sequence
  :type  name: :class:`str`
  :param sequence: String of characters representing the sequence. Only   
       'word' characters (no digits), '?', '-' and '.' are allowed. In an
       upcoming release, '?' and '.' will also be forbidden so its best to
       translate those to 'X' or '-'.
  :type sequence: :class:`str`
  :param role: Role of the sequence (optional)
  :type role:  :class:`str`
  :raises InvalidSequence: When the sequence string contains forbidden
       characters. In the future, '?' and '.' will also raise this exception.

.. class:: SequenceHandle
           ConstSequenceHandle

  Represents a sequence. New instances are created with :func:`CreateSequence`.
  
  .. method:: GetPos(residue_index)
  
    Get position of residue with index in sequence. This is best illustrated in 
    the following example:
    
    .. code-block:: python
      
      s=seq.CreateSequence("A", "abc---def")
      print(s.GetPos(1)) # prints 1
      print(s.GetPos(3)) # prints 6
    
    The reverse mapping, that is from position in the sequence to residue index 
    can be achieved with :meth:`GetResidueIndex`.
  
  .. method:: GetResidueIndex(pos)
     
    Get residue index of character at given position. This method is the
    inverse of :meth:`GetPos`. If the sequence contains a gap at that position,
    an :exc:`Error` is raised. Admires the
    :ref:`sequence offset <sequence-offset>`.
    
    .. code-block:: python
      
      s=seq.CreateSequence("A", "abc--def")
      print(s.GetResidueIndex(1)) # prints 1
      print(s.GetResidueIndex(6)) # prints 4
      # the following line raises an exception of type
      # Error with the message "requested position contains 
      # a gap"
      print(s.GetResidueIndex(3))

  .. method:: GetResidue(pos)
     
    As, :meth:`GetResidueIndex`, but directly returns the residue view. If no
    view is attached, or if the position is a gap, an invalid residue view
    is returned.
    
    :rtype: :class:`~ost.mol.ResidueView`
    
  .. method:: GetLastNonGap()
     
    Get position of last non-gap character in sequence. In case of an empty
    sequence, or, a sequence only consisting of hyphens, -1 is returned
     
  .. method:: GetFirstNonGap()
  
    Get position of first non-gap character in sequence. In case of an empty
    sequence, or, a sequence only consisting of hyphens, -1 is returned.

  .. method:: AttachView(view)
              AttachView(view, chain_name)
    
    Attach an :class:`~mol.EntityView` to sequence. The first signature requires
    that the view contains one chain. If not, an :exc:`IntegrityError` is
    raised. The second signature will select the chain with the given name. If 
    no such chain exists, an :exc:`IntegrityError` is raised.
    
  .. method:: HasAttachedView()
  
    Returns True when the sequence has a view attached, False if not.
    
  .. method:: GetAttachedView()
  
    Returns the attached :class:`~mol.EntityView`, or an invalid
    :class:`~mol.EntityView` if no view has been attached. Also available as 
    the property :attr:`attached_view`.
    
  .. method:: GetName()
  
    Returns the name of the sequence. Also available as the property
    :attr:`name`
  
  .. method:: SetOffset()
  
    Set the :ref:`sequence offset <sequence-offset>`. By default, the offset is
    0. Also available as the property :attr:`offset`.
    
  .. method:: GetOffset()
    
    Returns the :ref:`sequence offset <sequence-offset>`. Also available as
    :attr:`offset`.
    
    
  .. method:: GetGaplessString()
     
    Returns a string version of this sequence with all hyphens removed. Also
    available as the property :attr:`gapless_string`.
     
  .. method:: Normalise()
     
    Remove '-' and '.' as gaps from the sequence and make it all upper case.
    Works in place.
   
  .. method:: SetName()
  
    Set name of the sequence. Also available as the property :attr:`name`.

  .. method:: GetOneLetterCode(pos)
              __getitem__(pos)
              __getitem__(slice)

    :return: Character at position *pos* of sequence (also supports pythonic
             slicing with [] operator)
    :rtype:  :class:`str`
  
  .. attribute:: gapless_string
     
    Shorthand for :meth:`GetGaplessString()`
     
  .. attribute:: name
  
    Shorthand for :meth:`GetName`/:meth:`SetName`
  
  .. attribute:: attached_view
  
    Shorthand for :meth:`GetAttachedView`.

  .. attribute:: offset
  
    Shorthand for :meth:`GetOffset`/:meth:`SetOffset`

  .. attribute:: role

    Role of this sequence.

    :type: :class:`str`

  .. method:: __len__()
    
    :return: The length of the sequence (including insertions and deletions)
    
  .. method:: __str__()

    :return: The sequence as a string.

  .. method:: Copy()
    
    Create a deep copy of the sequence. The newly created sequence has the same
    attached view (not a deep copy of the view!).

.. function:: SequenceFromChain(name, chain)

  :return: Sequence extracted from one letter codes in given *chain* with a
           view to the chain attached to it
  :rtype:  :class:`SequenceHandle`
  :param name: Name of the sequence
  :type  name: :class:`str`
  :param chain: Chain from which to extract sequence
  :type chain:  :class:`~ost.mol.ChainHandle` / :class:`~ost.mol.ChainView`

.. function:: Match(s1, s2)

  :param s1: The first sequence
  :param s2: The second sequence
  :type s1: :class:`SequenceHandle`, or :class:`str`
  :type s2: :class:`SequenceHandle`, or :class:`str`

  Check whether the two sequences s1 and s2 match. This function performs are
  case-insensitive comparison of the two sequences. The character  'X' is
  interpreted as a wild card character that always matches the other sequence.

The SequenceList
--------------------------------------------------------------------------------

.. function:: CreateSequenceList()

  Creates and returns a new :class:`SequenceList` with no sequences.

.. class:: SequenceList
           ConstSequenceList

  Represents a list of sequences. The class provides a row-based interface.

  .. method:: GetCount()
              __len__()
  
    :return: Number of sequences in the list.
    :rtype:  :class:`int`

  .. method:: AddSequence(sequence)

    Append a sequence to the list.

  .. method:: GetMinLength()
              GetMaxLength()

    :return: Minimal / maximal length of the sequences in this list.
    :rtype:  :class:`int`

  .. method:: FindSequence(name)
  
    Find sequence with given *name*. If the alignment contains several sequences
    with the same name, the first sequence is returned.

  .. method:: SequencesHaveEqualLength()

    :return: True if all sequences have same length.

  .. method:: Take(n)

    :return: First *n* (or last *-n* if *n* negative) sequences.

  .. method:: Slice(first, n)

    :return: *n* sequences starting from *first*.

  .. method:: __getitem__(key)

    :return: Sequence(s) indexed by *key* (supports pythonic slicing).


The AlignmentHandle
--------------------------------------------------------------------------------

The :class:`AlignmentHandle` represents a list of aligned sequences. In
contrast to :class:`SequenceList`, an alignment requires all sequences to be of 
the same length. New instances of alignments are created with 
:func:`CreateAlignment` and :func:`AlignmentFromSequenceList`.

Typically sequence alignments are used column-based, i.e by looking at an  
aligned columns in the sequence alignment. To get a row-based (sequence) view
on the sequence list, use :meth:`~AlignmentHandle.GetSequences()`. 

All functions that operate on an alignment will again produce a valid alignment. 
This mean that it is not possible to change the length of one sequence, without  
adjusting the other sequences, too.

The following example shows how to iterate over the columns and sequences of
an alignment:

.. code-block:: python

  aln=io.LoadAlignment('aln.fasta')
  # iterate over the columns
  for col in aln:
    print(col)

  # iterate over the sequences
  for s in aln.sequences:
    print(s)

.. function:: CreateAlignment()

  Creates and returns a new :class:`AlignmentHandle` with no sequences.
  
.. function:: AlignmentFromSequenceList(sequences)
  
  Create a new alignment from the given list of sequences
  
  :param sequences: the list of sequences
  :type sequences: :class:`ConstSequenceList`
  
  :raises: :exc:`InvalidAlignment` if the sequences do not have the same length.

.. class:: AlignmentHandle
  
  .. method:: GetSequence(index)
  
    :return: Sequence at the given index, raising an IndexError when trying
             to access an inexistent sequence.
    :rtype:  :class:`ConstSequenceHandle`
    
  .. method:: GetSequences()
  
    :return: List of all sequence of the alignment. Also available as
             :attr:`sequences`.
    :rtype:  :class:`ConstSequenceList`
    
  .. method:: GetLength()
              __len__()
  
    :return: Length of the alignment.
    :rtype:  :class:`int`
    
  .. method:: GetCount()
  
    :return: Number of sequences in the alignment. Also available as
             :attr:`sequence_count`.
    :rtype:  :class:`int`
  
  .. method:: ToString(width=80)
  
    :return: Formatted string version of the alignment. The sequences are
             split into smaller parts to fit into the number columns specified.
    :rtype:  :class:`str`
    
    .. code-block:: python
    
      aln=seq.CreateAlignment()
      aln.AddSequence(seq.CreateSequence("A", "abcdefghik"))
      aln.AddSequence(seq.CreateSequence("B", "1234567890"))
      # The following command will print the output given below
      print(aln.ToString(7))
      # A abcde
      # B 12345
      #
      # A fghik
      # B 67890

  .. method:: FindSequence(name)
  
    :return: Sequence with given *name*. If the alignment contains several
             sequences with the same name, the first sequence is returned.
    
  .. method:: SetSequenceName(seq_index, name)
  
    Set the name of the sequence at index `seq_index` to `name`
    (see :attr:`SequenceHandle.name`).
    
  .. method:: Copy()
    
    Create a deep copy of the alignment by copying each contained sequence
    (see :meth:`SequenceHandle.Copy`)

  .. method:: GetPos(seq_index, res_index)
    
    :return: Position of residue with index equal to `res_index` in sequence at
             index `seq_index` (see :meth:`SequenceHandle.GetPos`)
    
  .. method:: GetResidueIndex(seq_index, pos)
  
    :return: Residue index of residue at position `pos` in sequence at index
             `seq_index` (see :meth:`SequenceHandle.GetResidueIndex`)
    
  .. method:: GetResidue(seq_index, pos)
  
    :return: Attached residue at position `pos` in sequence at index `seq_index`
             (see :meth:`SequenceHandle.GetResidue`).
  
  .. method:: AttachView(seq_index, view)
              AttachView(seq_index, view, chain_name)
    
    Attach the given view to the sequence at index `seq_index`
    (see :meth:`SequenceHandle.AttachView`).
    
  .. method:: Cut(start, end)
  
    Removes the columns in the half-closed interval `start`, `end` from the
    alignment. Note that this function does not update offsets!
    
    .. code-block:: python
    
      aln=seq.CreateAlignment()
      aln.AddSequence(seq.CreateSequence("A", "abcd---hik"))
      aln.AddSequence(seq.CreateSequence("B", "1234567890"))
      aln.Cut(4, 7)
      
      print(aln.ToString(80))
      # will print
      # A abcdhik
      # B 1234890
   
  .. method:: Replace(new_region, start, end)
  
    Replace the columns in the half-closed interval `start`, `end` with the
    columns in `new_region`.
    
    :param new_region: The region to be inserted
    :type new_region: :class:`AlignedRegion` or :class:`AlignmentHandle`
  
  .. method:: GetMatchingBackboneViews(index1=0, index2=1)
  
    Returns a tuple of entity views containing matching backbone atoms for the 
    two sequences at index1 and index2, respectively. For each aligned column in
    the alignment, backbone atoms are added to the view if both aligned residues 
    have them. It is guaranteed that the two views contain the same number of 
    atoms and that the order of the atoms in the two views is the same.
    
    The output of this function can be used to superpose two structures with
    :func:`~ost.mol.alg.SuperposeSVD`.
    
    
    :param index1: The index of the first sequence
    
    :param index2: The index of the second sequence.
    
    :raises: In case one of the two sequences doesn't have an attached view, a 
       :exc:`RuntimeError` is raised.

  .. method:: AddSequence(sequence)

    Append a sequence to the alignment. The sequence must have the same length as
    sequences already present in the alignment.

    :raises: :exc:`RuntimeError` if the sequence length does not match
    :param sequence: Sequence to be added
    :type sequence: :class:`ConstSequenceHandle`

  .. method:: GetSequenceOffset(index)
              SetSequenceOffset(index, offset)

    Get/set the offset for sequence at *index*
    (see :attr:`SequenceHandle.offset`).

    :param index: The index of the sequence
    :type index: :class:`int`
    :param offset: The new offset
    :type offset: :class:`int`
    :rtype: :class:`int`
  
  .. method:: GetSequenceRole(index)
              SetSequenceRole(index, role)

    Get/Set the sequence role for sequence at *index*
    (see :attr:`SequenceHandle.role`).

    :param index: The index of the sequence
    :type index: :class:`int`
    :param role: The new role
    :type role: :class:`str`
    :rtype: :class:`str`

  .. method:: GetCoverage(index)

    Get coverage of sequence at *index* to the first sequence.

    :param index: The index of the sequence
    :type index: :class:`int`
    :returns: Coverage as a number between 0 and 1.
  
  .. method:: RemoveSequence(index)

    Remove sequence at *index* from the alignment.
  
  .. attribute:: sequences
  
    Shorthand for :meth:`GetSequences`

  .. attribute:: sequence_count
  
    Shorthand for :meth:`GetCount`

  .. method:: __getitem__(pos)

    :return: Column at position *pos* of alignment.
    :rtype:  :class:`AlignedColumn`

  .. method:: __getitem__(slice)

    :return: Columns defined by by pythonic slicing.
    :rtype:  :class:`AlignedRegion`


.. class:: AlignedRegion

  Represents a slice of an :class:`AlignmentHandle`.

  .. method:: GetAlignmentHandle()

    :return: Alignment from which we slices.
    :rtype:  :class:`AlignmentHandle`

  .. method:: GetLength()
              __len__()

    :return: Number of columns in the slice.

  .. method:: __getitem__(pos)

    :return: Column at position *pos* within this slice.
    :rtype:  :class:`AlignedColumn`

  .. attribute:: start

    Starting position in alignment.

  .. attribute:: end

    One after end position in alignment.


.. class:: AlignedColumn

  .. method:: GetIndex()

    :return: Position in alignment.

  .. method:: GetRowCount()

    :return: Number of rows in the column.

  .. method:: GetResidue(row)
  
    :return: Attached residue for sequence at given *row* of this column
             (see :meth:`AlignmentHandle.GetResidue`).

  .. method:: __getitem__(row)

    :return: Character at given *row* of this column.
    :rtype:  :class:`str`

  .. method:: __str__()

    :return: String representation of column in alignment.


Extracting views from sequences
--------------------------------------------------------------------------------

.. function:: ViewsFromSequences(seq1, seq2)

  Returns a tuple of entity views containing only the atoms of the aligned
  residues. The order of residues in the two views is guaranteed to be the same
  but the order of atoms within each residue may differ. If the order of atoms
  is crucial (e.g. for :func:`~ost.mol.alg.SuperposeSVD`) either prefilter the
  attached views to include only one atom per residue or use the slower (approx.
  50% more runtime) :meth:`AlignmentHandle.GetMatchingBackboneViews`.

  :return: Pair of views including all the aligned residues of the two given
           sequences. An alignment is 
  :rtype:  :class:`tuple` with two :class:`~ost.mol.EntityView`

  :raises: :class:`Exception` if sequence lengths do not match or if any of the
           sequences is lacking an attached view.

.. function:: ViewsFromAlignment(aln, index1=0, index2=1)

  :return: Pair of views as in :meth:`ViewsFromSequences`.
  :rtype:  :class:`tuple` with two :class:`~ost.mol.EntityView`

  :param aln: Alignment from which to extract sequences.
  :type aln:  :class:`AlignmentHandle`
  :param index1: Index of first sequence in *aln* to use.
  :type index1:  :class:`int`
  :param index2: Index of second sequence in *aln* to use.
  :type index2:  :class:`int`


Handling Sequence Profiles
--------------------------------------------------------------------------------

The :class:`ProfileHandle` provides a simple container for profiles for each
residue. It mainly contains:

- *N* :class:`ProfileColumn` objects (*N* = number of residues in sequence)
  which each contains 20 amino acid frequencies
- a :attr:`~ProfileHandle.sequence` (:class:`str`) of length *N*
- a :attr:`~ProfileHandle.null_model` to use for this profile

Optionally, HMM-related information can be added. This is transition 
probabilities between Match, Insertion or Deletion states or neff values
(number of effective sequences, a measure of local sequence diversity).

.. class:: HMMTransition

  The possible HMM-transitions between Match(M), Insertion(I) and Deletion(D)
  states. Transitions between Deletion and Insertion are disallowed:

  HMM_M2M, HMM_M2I, HMM_M2D, HMM_I2M, HMM_I2I, HMM_D2M, HMM_D2D


.. class:: HMMData

  Data container for HMM-related information that can be assigned to profile 
  columns.

    .. attribute:: neff

      Local diversity of all sequences that have a residue at this column of 
      the full alignment.

    .. attribute:: neff_i

      Local diversity of all sequences that have an insertion at this column of 
      the full alignment.

    .. attribute:: neff_d

      Same for deletion.

    .. method:: GetProb(transition)

      Get HMM transition probability

      :param transition: The transition
      :type transition: :class:`HMMTransition`

    .. method:: SetProb(transition, prob)

      Set HMM transition probability

      :param transition: The transition
      :param prob: The probablity to be set

      :type transition: :class:`HMMTransition`
      :type prob: :class:`float`



.. class:: ProfileColumn

  .. method:: BLOSUMNullModel()

    Static method, that returns a new :class:`ProfileColumn` with amino acid
    frequencies given from the BLOSUM62 substitution matrix.

  .. method:: HHblitsNullModel()

    Static method, that returns a new :class:`ProfileColumn` with amino acid
    frequencies as set in HHblits output.

  .. method:: GetFreq(aa)

    :return: Frequency of *aa*
    :rtype:  :class:`float`
    :param aa: One letter code of standard amino acid
    :type aa:  :class:`str`

  .. method:: SetFreq(aa,freq)

    :param aa:  One letter code of standard amino acid
    :param freq:  The frequency of the given amino acid
    :type aa:  :class:`str`
    :type freq:  :class:`float`

  .. method:: GetScore(other, null_model)

    :return: Column score as in Soeding-2005 paper.
    :rtype:  :class:`float`
    :param other: Other column to compute score with.
    :type other:  :class:`ProfileColumn`
    :param null_model: Null model to use for weighting.
    :type null_model:  :class:`ProfileColumn`

  .. method:: SetHMMData(data)
    
    :param data: Data to be set
    :type data: :class:`HMMData`


  .. method:: GetHMMData()

    Returns previously set :class:`HMMData` object.

    :rtype: :class:`HMMData`

    :raises: :exc:`~exceptions.Error` if data has never been set.

  .. attribute:: entropy

    Shannon entropy based on the columns amino acid frequencies

    :type: :class:`float`


.. class:: ProfileHandle

  .. method:: __len__()
    
    Returns the length of the sequence for which we have profile.

    :rtype: :class:`int`

  .. method:: AddColumn(col, olc='X')

    Appends column in the internal column list.

    :param col: Column to add to :attr:`columns`
    :type col:  :class:`ProfileColumn`
    :param olc: One letter code to add to :attr:`sequence`
    :type col:  :class:`str`

  .. method:: Extract(from, to)

    :param from:  Col Idx to start from
    :param to:  End Idx, not included in sub-ProfileHandle

    :type from:  :class:`int`
    :type to:  :class:`int`

    :returns: sub-profile as defined by given indices
              (:attr:`null_model` and :attr:`neff` are copied,
              you might want to manually reset neff)
    :rtype: :class:`ProfileHandle`

    :raises: :exc:`~exceptions.Error` if *to* <= *from* or
              *to* > :meth:`__len__`.

  .. method:: GetAverageScore(other, offset=0)

    :return: Average column score between *other.columns[i]* and this object's
             *columns[i+offset]* for *i* in [*0, len(other)-1*] using this
             object's :attr:`null_model`. See :meth:`ProfileColumn.GetScore`.
    :rtype:  :class:`float`
    :param other: Other profile to compare with.
    :type other:  :class:`ProfileHandle`
    :param offset: Start comparison at column *offset* of this object.
    :type offset:  :class:`int`

    :raises: :exc:`~exceptions.Error` if any *columns[i+offset]* out of bounds.


  .. attribute:: sequence

    Sequence for which we have this profile. When setting a new value, the
    length and the number of profile columns must match (exception thrown
    otherwise).

    :type: :class:`str`

  .. attribute:: columns

    Iterable columns of the profile (read-only).

    :type: :class:`ProfileColumnList`

  .. attribute:: null_model

    Null model of the profile. By default this is set to
    :meth:`ProfileColumn.HHblitsNullModel`.

    :type: :class:`ProfileColumn`

  .. attribute:: avg_entropy

    Average entropy of all the columns (read-only).

    :type: :class:`float`

  .. attribute:: neff

    Measure for sequence diversity which is defined as the average of the 
    per-column neff values. However, this is just a convenience attribute
    which can be set to arbitrary values but there is no guarantee that 
    it's the actual average of the per-column values.


.. class:: ProfileDB

  A simple database to gather :class:`ProfileHandle` objects. It is possible
  to save them to disk in a compressed format with limited accuracy
  (4 digits for each frequency).

  .. method:: Save(filename)

    :param filename:  Name of file that will be generated on disk.
    :type filename:  :class:`str`

  .. method:: Load(filename)

    Static loading method

    :param filename:  Name of file from which the database should be loaded.
    :type filename:  :class:`str`
    :returns:  The loaded database

  .. method:: AddProfile(name, prof)

    :param name:  Name of profile to be added
    :param prof:  Profile to be added

    :type name:  :class:`str`
    :type prof:  :class:`ProfileHandle`
    :raises:  :class:`Exception` when filename is longer than 255 characters.

  .. method:: GetProfile(name)

    :param name:  Name of profile to be returned
    :type name:  :class:`str`
    :returns:  The requested :class:`ProfileHandle`
    :raises:  :class:`Exception` when no :class:`ProfileHandle` for **name** exists.

  .. method:: Size()

    :returns: Number of :class:`ProfileHandle` objects in the database

  .. method:: GetNames()

    :returns: A nonsorted list of the names of all :class:`ProfileHandle`
              objects in the database
