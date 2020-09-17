:mod:`seq.alg <ost.seq.alg>` -- Algorithms for Sequences
================================================================================

.. module:: ost.seq.alg
  :synopsis: Algorithms for sequences

Algorithms for Alignments
--------------------------------------------------------------------------------

.. function:: MergePairwiseAlignments(pairwise_alns, ref_seq)

  :param pairwise_alns: A list of pairwise alignments
  :type pairwise_alns: :class:`~ost.seq.AlignmentList`

  :param ref_seq: The reference sequence
  :type ref_seq: :class:`~ost.seq.SequenceHandle`

  :returns: The merged alignment
  :rtype: :class:`~ost.seq.AlignmentHandle`

  Merge a list of pairwise alignments into a multiple sequence alignments. This
  function uses the reference sequence as the anchor and inserts gaps where
  needed. This is also known as the *star method*.

  The resulting multiple sequence alignment provides a simple way to map between 
  residues of pairwise alignments, e.g. to compare distances in two structural 
  templates.
  
  There are a few things to keep in mind when using this function:
  
   - The reference sequence mustn't contain any gaps
   
   - The first sequence of each pairwise alignments corresponds to the reference 
     sequence. Apart from the presence of gaps, these two sequences must be 
     completely identical.

   - If the reference sequence has an offset, the first sequence of each pairwise alignment 
     must have the same offset. This offset is inherited by the first sequence of the final
     output alignment.
   
   - The resulting multiple sequence alignment is by no means optimal. For 
     better results, consider using a multiple-sequence alignment program such 
     as MUSCLE or ClustalW.
   
   - Residues in columns where the reference sequence has gaps should not be 
     considered as aligned. There is no information in the pairwise alignment to 
     guide the merging, the result is undefined.


     **Example:**

     .. code-block:: python

       ref_seq = ost.seq.CreateSequence('ref', 'acdefghiklmn')
       seq_a1 = seq.CreateSequence('A1', 'acdefghikl-mn')
       seq_a2 = seq.CreateSequence('A2', 'atd-fghikllmn')
       seq_b1 = seq.CreateSequence('B1', 'acdefg-hiklmn')
       seq_b2 = seq.CreateSequence('B2', 'acd---qhirlmn')

       aln_a = seq.CreateAlignment()
       aln_a.AddSequence(seq_a1)
       aln_a.AddSequence(seq_a2)
       print(aln_a)
       # >>> A1  acdefghikl-mn
       # >>> A2  atd-fghikllmn

       aln_b = seq.CreateAlignment()
       aln_b.AddSequence(seq_b1)
       aln_b.AddSequence(seq_b2)
       print(aln_b)
       # >>> B1  acdefg-hiklmn
       # >>> B2  acd---qhirlmn

       aln_list = ost.seq.AlignmentList()
       aln_list.append(aln_a)
       aln_list.append(aln_b)

       merged_aln = ost.seq.alg.MergePairwiseAlignments(aln_list, ref_seq)
       print(merged_aln)
       # >>> ref  acdefg-hikl-mn
       # >>> A2   atd-fg-hikllmn
       # >>> B2   acd---qhirl-mn


.. autofunction:: ValidateSEQRESAlignment

.. autofunction:: AlignToSEQRES

.. autofunction:: AlignmentFromChainView

.. function:: Conservation(aln, assign=true, prop_name="cons", ignore_gap=false)

  Calculates conservation scores for each column in the alignment, according to
  the ConSurf method (Armon et al., J. Mol. Biol. (2001) 307, 447-463).
  
  The conservation score is a value between 0 and 1. The bigger the number 
  the more conserved the aligned residues are. 
  
  :param aln: An alignment handle
  :type aln: :class:`~ost.seq.AlignmentHandle`
  :param assign: If true, the conservation scores are assigned to attached 
      residues. The name of the property can be changed with the *prop_name* 
      parameter. Useful when coloring entities based on sequence conservation.
  :param prop_name: The property name for assigning the conservation to 
      attached residues. Defaults to 'cons'.
  :param ignore_gap: If true, the dissimilarity between two gaps is increased to
      6.0 instead of 0.5 as defined in the original version. Without this, a
      stretch where in the alignment there is only one sequence which is
      aligned to only gaps, is considered highly conserved (depending on the
      number of gap sequences).

.. function:: LocalAlign(seq1, seq2, subst_weight, gap_open=-5, gap_ext=-2)

  Performs a Smith/Waterman local alignment of *seq1* and *seq2* and returns
  the best-scoring alignments as a list of pairwise alignments.
  
  **Example:**
  
  .. code-block:: python
  
    seq_a = seq.CreateSequence('A', 'acdefghiklmn')
    seq_b = seq.CreateSequence('B', 'acdhiklmn')
    alns = seq.alg.LocalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    print(alns[0].ToString(80))
    # >>> A acdefghiklmn
    # >>> B acd---hiklmn

  :param seq1: A valid sequence
  :type seq1: :class:`~ost.seq.ConstSequenceHandle`
  :param seq2: A valid sequence  
  :type seq2: :class:`~ost.seq.ConstSequenceHandle`
  :param subst_weigth: The substitution weights matrix
  :type subst_weight: :class:`SubstWeightMatrix`
  :param gap_open: The gap opening penalty. Must be a negative number
  :param gap_ext: The gap extension penalty. Must be a negative number
  :returns: A list of best-scoring, non-overlapping alignments of *seq1* and 
     *seq2*. Since alignments always start with a replacement, the start is
     stored in the sequence offset of the two sequences.


.. function:: GlobalAlign(seq1, seq2, subst_weight, gap_open=-5, gap_ext=-2)

  Performs a Needleman/Wunsch global alignment of *seq1* and *seq2* and returns
  the best-scoring alignment.
  
  **Example:**
  
  .. code-block:: python
  
    seq_a = seq.CreateSequence('A', 'acdefghiklmn')
    seq_b = seq.CreateSequence('B', 'acdhiklmn')
    alns = seq.alg.GlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    print(alns[0].ToString(80))
    # >>> A acdefghiklmn
    # >>> B acd---hiklmn

  :param seq1: A valid sequence
  :type seq1: :class:`~ost.seq.ConstSequenceHandle`
  :param seq2: A valid sequence  
  :type seq2: :class:`~ost.seq.ConstSequenceHandle`
  :param subst_weigth: The substitution weights matrix
  :type subst_weight: :class:`SubstWeightMatrix`
  :param gap_open: The gap opening penalty. Must be a negative number
  :param gap_ext: The gap extension penalty. Must be a negative number
  :returns: Best-scoring alignment of *seq1* and *seq2*.

.. function:: ShannonEntropy(aln, ignore_gaps=True)

  Returns the per-column Shannon entropies of the alignment. The entropy
  describes how conserved a certain column in the alignment is. The higher
  the entropy is, the less conserved the column. For a column with no amino 
  aids, the entropy value is set to NAN.

  :param aln: Multiple sequence alignment
  :type aln: :class:`~ost.seq.AlignmentHandle`
  :param ignore_gaps: Whether to ignore gaps in the column.
  :type ignore_gaps: bool

  :returns: List of column entropies

.. function:: SemiGlobalAlign(seq1, seq2, subst_weight, gap_open=-5, gap_ext=-2)

  Performs a semi-global alignment of *seq1* and *seq2* and returns the best-
  scoring alignment. The algorithm is Needleman/Wunsch same as GlobalAlign, but
  without any gap penalty for starting or ending gaps. This is prefereble 
  whenever one of the sequences is significantly shorted than the other.
  This make it also suitable for fragment assembly.
  
  **Example:**
  
  .. code-block:: python
  
    seq_a = seq.CreateSequence('A', 'abcdefghijklmnok')
    seq_b = seq.CreateSequence('B', 'cdehijk')
    alns = seq.alg.GlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    print(alns[0].ToString(80))
    # >>> A abcdefghijklmnok
    # >>> B --cde--hi-----jk
    alns = seq.alg.SemiGlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    print(alns[0].ToString(80))
    # >>> A abcdefghijklmnok
    # >>> B --cde--hijk-----

  :param seq1: A valid sequence
  :type seq1: :class:`~ost.seq.ConstSequenceHandle`
  :param seq2: A valid sequence  
  :type seq2: :class:`~ost.seq.ConstSequenceHandle`
  
  :param subst_weigth: The substitution weights matrix
  :type subst_weight: :class:`SubstWeightMatrix`
  :param gap_open: The gap opening penalty. Must be a negative number
  :param gap_ext: The gap extension penalty. Must be a negative number
  :returns: best-scoring alignment of *seq1* and *seq2*.

.. autofunction:: ost.seq.alg.renumber.Renumber

.. _contact-prediction:

Contact Prediction
--------------------------------------------------------------------------------

This is a set of functions for predicting pairwise contacts from a multiple
sequence alignment (MSA). The core method here is mutual information which uses 
coevolution to predict contacts. Mutual information is complemented by two other 
methods which score pairs of columns of a MSA from the likelyhood of certain
amino acid pairs to form contacts (statistical potential) and the likelyhood
of finding certain substitutions of aminio-acid pairs in columns of the MSA
corresponding to interacting residues.

.. class:: ContactPredictionScoreResult
  
  Object containing the results form a contact prediction. 

  .. attribute:: matrix

    An *NxN* :class:`~ost.FloatMatrix` where *N* is the length of the alignment.
    The element *i,j* corresponds to the score of the corresponding
    columns of the MSA. High scores correspond to high likelyhood of
    a contact.

  .. attribute:: sorted_indices

    List of all indices pairs *i,j*, containing (N*N-1)/2 elements,
    as the **matrix** is symmetrical and elements in the diagonal
    are ignored. The indices are sorted from the pair most likely to form
    a contact to the least likely one.

  .. method:: GetScore(i,j)

    returns **matrix(i,j)**

    :param i: First index
    :param j: Second index
    :type i:  :class:`int`
    :type j:  :class:`int`

  .. method:: SetScore(i,j,score)

    Sets **matrix(i,j)** to **score**

    :param i: First index
    :param j: Second index
    :param score: The score
    :type i:  :class:`int`
    :type j:  :class:`int`
    :type score:  :class:`float`

.. autofunction:: PredictContacts

.. function:: CalculateMutualInformation(aln, \
                weights=LoadConstantContactWeightMatrix(), \
                apc_correction=true, zpx_transformation=true, \
                small_number_correction=0.05)

    Calculates the mutual information (MI) from a multiple sequence alignemnt. Contributions of each pair of amino-acids are weighted using the matrix **weights** (weighted mutual information). The average product correction (**apc_correction**) correction and transformation into Z-scores (**zpx_transofrmation**) increase prediciton accuracy by reducing the effect of phylogeny and other noise sources. The small number correction reduces noise for alignments with small number of sequences of low diversity.

    :param aln: The multiple sequences alignment
    :type aln:  :class:`~ost.seq.AlignmentHandle`
    :param weights: The weight matrix
    :type weights:  :class`ContactWeightMatrix`
    :param apc_correction: Whether to use the APC correction
    :type apc_correction:  :class:`bool`
    :param zpx_transformation:  Whether to transform the scores into Z-scores
    :type zpx_transformation: :class:`bool`
    :param small_number_correction: initial values for the probabilities of having a given pair of amino acids *p(a,b)*.
    :type small_number_correction: :class:`float`

.. function:: CalculateContactScore(aln, \
                weights=LoadDefaultContactWeightMatrix())
  
  Calculates the Contact Score (*CoSc*) from a multiple sequence alignment. For each pair of residues *(i,j)* (pair of columns in the MSA), *CoSc(i,j)* is the average over the values of the **weights** corresponding to the amino acid pairs in the columns.

  :param aln: The multiple sequences alignment
  :type aln:  :class:`~ost.seq.AlignmentHandle`
  :param weights: The contact weight matrix
  :type weights:  :class`ContactWeightMatrix`

.. function:: CalculateContactSubstitutionScore(aln, ref_seq_index=0, \
                weights=LoadDefaultPairSubstWeightMatrix())

  Calculates the Contact Substitution Score (*CoEvoSc*) from a multiple sequence alignment. For each pair of residues *(i,j)* (pair of columns in the MSA), *CoEvoSc(i,j)* is the average over the values of the **weights** corresponding to substituting the amino acid pair in the reference sequence (given by **ref_seq_index**) with all other pairs in columns *(i,j)* of the **aln**.

  :param aln: The multiple sequences alignment
  :type aln:  :class:`~ost.seq.AlignmentHandle`
  :param weights: The pair substitution weight matrix
  :type weights:  :class`ContactWeightMatrix`

.. function:: LoadDefaultContactWeightMatrix()
  
  :returns: *CPE*, a :class:`ContactWeightMatrix` that was calculated from a large (>15000) set of
    high quality crystal structures as *CPE=log(CF(a,b)/NCF(a,b))* and then normalised so that all its elements are comprised between 0 and 1. *CF(a,b)* is the frequency of amino acids *a* and *b* for pairs of contacting residues and *NCF(a,b)* is the frequency of amino acids *a* and *b* for pairs of non-contacting residues. Apart from weights for the standard amino acids, this matrix gives a weight of 0 to all pairs for which at least one amino-acid is a gap.

.. function:: LoadConstantContactWeightMatrix()
  
  :returns: A :class:`ContactWeightMatrix`. This matrix gives a weight of one to all pairs of
   standard amino-acids and a weight of 0 to pairs for which at least one amino-acid is a gap.

.. function:: LoadDefaultPairSubstWeightMatrix()
  
  :returns: *CRPE*, a :class:`PairSubstWeightMatrix` that was calculated from a large (>15000) set of
    high quality crystal structures as *CRPE=log(CRF(ab->cd)/NCRF(ab->cd))* and then normalised so that all its elements are comprised between 0 and 1. *CRF(ab->cd)* is the frequency of replacement of a pair of amino acids  *a* and *b* by a pair *c* and *d* in columns of the MSA corresponding to contacting residues and *NCRF(ab->cd)* is the frequency of replacement of a pair of amino acids  *a* and *b* by a pair *c* and *d* in columns of the MSA corresponding to non-contacting residues. Apart from weights for the standard amino acids, this matrix gives a weight of 0 to all pair substitutions for which at least one amino-acid is a gap.


.. class:: PairSubstWeightMatrix(weights, aa_list)

  This class is used to associate a weight to any substitution from one amino-acid pair *(a,b)* to any other pair *(c,d)*.

  .. attribute:: weights

    A :class:`~ost.FloatMatrix4` of size *NxNxNxN*, where *N=len(aa_list)*

  .. attribute:: aa_list

    A :class:`CharList` of one letter codes of the amino acids for which weights are found in the **weights** matrix.

.. class:: ContactWeightMatrix(weights, aa_list)

  This class is used to associate a weight to any pair of amino-acids.

  .. attribute:: weights

    A :class:`~ost.FloatMatrix` of size *NxN*, where *N=len(aa_list)*

  .. attribute:: aa_list

    A :class:`CharList` of one letter codes of the amino acids for which weights are found in the **weights** matrix.

Get and analyze distance matrices from alignments
--------------------------------------------------------------------------------

Given a multiple sequence alignment between a reference sequence (first sequence
in alignment) and a list of structures (remaining sequences in alignment with an
attached view to the structure), this set of functions can be used to analyze
differences between the structures.

**Example:**

.. code-block:: python
  
  # SETUP: aln is multiple sequence alignment, where first sequence is the
  #        reference sequence and all others have a structure attached

  # clip alignment to only have parts with at least 3 sequences (incl. ref.)
  # -> aln will be cut and clip_start is 1st column of aln that was kept
  clip_start = seq.alg.ClipAlignment(aln, 3)
  
  # get variance measure and distance to mean for each residue pair
  d_map = seq.alg.CreateDistanceMap(aln)
  var_map = seq.alg.CreateVarianceMap(d_map)
  dist_to_mean = seq.alg.CreateDist2Mean(d_map)

  # report min. and max. variances
  print("MIN-MAX:", var_map.Min(), "-", var_map.Max())
  # get data and json-strings for further processing
  var_map_data = var_map.GetData()
  var_map_json = var_map.GetJsonString()
  dist_to_mean_data = dist_to_mean.GetData()
  dist_to_mean_json = dist_to_mean.GetJsonString()

.. function:: ClipAlignment(aln, n_seq_thresh=2, set_offset=true, \
                            remove_empty=true)

  Clips alignment so that first and last column have at least the desired number
  of structures.

  :param aln: Multiple sequence alignment. Will be cut!
  :type aln:  :class:`~ost.seq.AlignmentHandle`
  :param n_seq_thresh: Minimal number of sequences desired.
  :type n_seq_thresh:  :class:`int`
  :param set_offset: Shall we update offsets for attached views?
  :type set_offset:  :class:`bool`
  :param remove_empty: Shall we remove sequences with only gaps in cut aln?
  :type remove_empty:  :class:`bool`
  :returns: Starting column (0-indexed), where cut region starts (w.r.t.
            original aln). -1, if there is no region in the alignment with
            at least the desired number of structures.
  :rtype:   :class:`int`

.. function:: CreateDistanceMap(aln)

  Create distance map from a multiple sequence alignment.
  
  The algorithm requires that the sequence alignment consists of at least two
  sequences. The sequence at index 0 serves as a frame of reference. All the
  other sequences must have an attached view and a properly set sequence offset
  (see :meth:`~ost.seq.AlignmentHandle.SetSequenceOffset`).
  
  For each of the attached views, the C-alpha distance pairs are extracted and
  mapped onto the corresponding C-alpha distances in the reference sequence.

  :param aln: Multiple sequence alignment.
  :type aln:  :class:`~ost.seq.AlignmentHandle`
  :returns: Distance map.
  :rtype:   :class:`DistanceMap`
  :raises:  Exception if *aln* has less than 2 sequences or any sequence (apart
            from index 0) is lacking an attached view.

.. function:: CreateVarianceMap(d_map, sigma=25)

  :returns: Variance measure for each entry in *d_map*.
  :rtype:   :class:`VarianceMap`
  :param d_map: Distance map as created with :func:`CreateDistanceMap`.
  :type d_map:  :class:`DistanceMap`
  :param sigma: Used for weighting of variance measure
                (see :meth:`Distances.GetWeightedStdDev`)
  :type sigma:  :class:`float`
  :raises:  Exception if *d_map* has no entries.

.. function:: CreateDist2Mean(d_map)

  :returns: Distances to mean for each structure in *d_map*.
            Structures are in the same order as passed when creating *d_map*.
  :rtype:   :class:`Dist2Mean`
  :param d_map: Distance map as created with :func:`CreateDistanceMap`.
  :type d_map:  :class:`DistanceMap`
  :raises:  Exception if *d_map* has no entries.

.. class:: Distances
  
  Container used by :class:`DistanceMap` to store a pair wise distance for each
  structure. Each structure is identified by its index in the originally used
  alignment (see :func:`CreateDistanceMap`).

  .. method:: GetDataSize()

    :returns: Number of pairwise distances.
    :rtype:   :class:`int`

  .. method:: GetAverage()

    :returns: Average of all distances.
    :rtype:   :class:`float`
    :raises:  Exception if there are no distances.

  .. method:: GetMin()
              GetMax()

    :returns: Minimal/maximal distance.
    :rtype:   :class:`tuple` (distance (:class:`float`), index (:class:`int`))
    :raises:  Exception if there are no distances.

  .. method:: GetDataElement(index)

    :returns: Element at given *index*.
    :rtype:   :class:`tuple` (distance (:class:`float`), index (:class:`int`))
    :param index: Index within list of distances (must be < :meth:`GetDataSize`).
    :type index:  :class:`int`
    :raises:  Exception if there are no distances or *index* out of bounds.

  .. method:: GetStdDev()

    :returns: Standard deviation of all distances.
    :rtype:   :class:`float`
    :raises:  Exception if there are no distances.

  .. method:: GetWeightedStdDev(sigma)

    :returns: Standard deviation of all distances multiplied by
              exp( :meth:`GetAverage` / (-2*sigma) ).
    :rtype:   :class:`float`
    :param sigma: Defines weight.
    :type sigma:  :class:`float`
    :raises:  Exception if there are no distances.

  .. method:: GetNormStdDev()

    :returns: Standard deviation of all distances divided by :meth:`GetAverage`.
    :rtype:   :class:`float`
    :raises:  Exception if there are no distances.

.. class:: DistanceMap

  Container returned by :func:`CreateDistanceMap`.
  Essentially a symmetric :meth:`GetSize` x :meth:`GetSize` matrix containing
  up to :meth:`GetNumStructures` distances (list stored as :class:`Distances`).
  Indexing of residues starts at 0 and corresponds to the positions in the
  originally used alignment (see :func:`CreateDistanceMap`).

  .. method:: GetDistances(i_res1, i_res2)

    :returns: List of distances for given pair of residue indices.
    :rtype:   :class:`Distances`
    :param i_res1: Index of residue.
    :type i_res1:  :class:`int`
    :param i_res2: Index of residue.
    :type i_res2:  :class:`int`

  .. method:: GetSize()

    :returns: Number of residues in map.
    :rtype:   :class:`int`

  .. method:: GetNumStructures()

    :returns: Number of structures originally used when creating the map
              (see :func:`CreateDistanceMap`).
    :rtype:   :class:`int`

.. class:: VarianceMap

  Container returned by :func:`CreateVarianceMap`.
  Like :class:`DistanceMap`, it is a symmetric :meth:`GetSize` x :meth:`GetSize`
  matrix containing variance measures.
  Indexing of residues is as in :class:`DistanceMap`.

  .. method:: Get(i_res1, i_res2)

    :returns: Variance measure for given pair of residue indices.
    :rtype:   :class:`float`
    :param i_res1: Index of residue.
    :type i_res1:  :class:`int`
    :param i_res2: Index of residue.
    :type i_res2:  :class:`int`

  .. method:: GetSize()

    :returns: Number of residues in map.
    :rtype:   :class:`int`

  .. method:: Min()
              Max()

    :returns: Minimal/maximal variance in the map.
    :rtype:   :class:`float`

  .. method:: ExportDat(file_name)
              ExportCsv(file_name)
              ExportJson(file_name)

    Write all variance measures into a file. The possible formats are:

    - "dat" file: a list of "*i_res1+1* *i_res2+1* variance" lines
    - "csv" file: a list of ";" separated variances (one line for each *i_res1*)
    - "json" file: a JSON formatted file (see :meth:`GetJsonString`)

    :param file_name: Path to file to be created.
    :type file_name:  :class:`str`
    :raises:  Exception if the file cannot be opened for writing.

  .. method:: GetJsonString()

    :returns: A JSON formatted list of :meth:`GetSize` lists with
              :meth:`GetSize` variances
    :rtype:   :class:`str`

  .. method:: GetData()

    Gets all the data in this map at once. Note that this is much faster (10x
    speedup observed) than parsing :meth:`GetJsonString` or using :meth:`Get`
    on each element.

    :returns: A list of :meth:`GetSize` lists with :meth:`GetSize` variances.
    :rtype:   :class:`list` of :class:`list` of :class:`float`

.. class:: Dist2Mean

  Container returned by :func:`CreateDist2Mean`.
  Stores distances to mean for :meth:`GetNumResidues` residues of
  :meth:`GetNumStructures` structures.
  Indexing of residues is as in :class:`DistanceMap`.
  Indexing of structures goes from 0 to :meth:`GetNumStructures` - 1 and is in
  the same order as the structures in the originally used alignment.

  .. method:: Get(i_res, i_str)

    :returns: Distance to mean for given residue and structure indices.
    :rtype:   :class:`float`
    :param i_res: Index of residue.
    :type i_res:  :class:`int`
    :param i_str: Index of structure.
    :type i_str:  :class:`int`

  .. method:: GetNumResidues()

    :returns: Number of residues.
    :rtype:   :class:`int`

  .. method:: GetNumStructures()

    :returns: Number of structures.
    :rtype:   :class:`int`

  .. method:: ExportDat(file_name)
              ExportCsv(file_name)
              ExportJson(file_name)

    Write all distance measures into a file. The possible formats are:

    - "dat" file: a list of "*i_res+1* distances" lines (distances are space
      separated)
    - "csv" file: a list of ";" separated distances (one line for each *i_res*)
    - "json" file: a JSON formatted file (see :meth:`GetJsonString`)

    :param file_name: Path to file to be created.
    :type file_name:  :class:`str`
    :raises:  Exception if the file cannot be opened for writing.

  .. method:: GetJsonString()

    :returns: A JSON formatted list of :meth:`GetNumResidues` lists with
              :meth:`GetNumStructures` distances.
    :rtype:   :class:`str`

  .. method:: GetData()

    Gets all the data in this map at once. Note that this is much faster (10x
    speedup observed) than parsing :meth:`GetJsonString` or using :meth:`Get`
    on each element.

    :returns: A list of :meth:`GetNumResidues` lists with
              :meth:`GetNumStructures` distances.
    :rtype:   :class:`list` of :class:`list` of :class:`float`


HMM Algorithms
--------------------------------------------------------------------------------
Openstructure implements basic HMM-related functionality that aims at
calculating an HMM-HMM alignment score as described in
Soding, Bioinformatics (2005) 21(7), 951-60. This is the score which is
optimized in the Viterbi algorithm of the hhalign tool. 
As a prerequisite, OpenStructure also implements adding pseudo counts to 
:class:`ost.seq.ProfileHandle` in order to avoid zero probabilities for 
unobserved transitions/emissions. Given these requirements, all functions
in this section require HMM related data (transition probabilities, neff values,
etc.) to be set, which is the case if you load a file in hhm format.

.. method:: HMMScore(profile_0, profile_1, aln, s_0_idx, s_1_idx, \
                     match_score_offset=-0.03,correl_score_weight=0.1, \
                     del_start_penalty_factor=0.6, \
                     del_extend_penalty_factor=0.6, \
                     ins_start_penalty_factor=0.6, \
                     ins_extend_penalty_factor=0.6)

  Scores an HMM-HMM alignment given in *aln* between *profile_0* and 
  *profile_1*.
  The score is described in Soding, Bioinformatics (2005) 21(7), 951-60 and 
  consists of three components: 

    * sum of column alignment scores of all aligned columns, the 
      *match_score_offset* is applied to each of those scores
    * sum of transition probability scores, the prefactor of those scores can 
      be controlled with penalty factors (*del_start_penalty_factor* etc.)
    * correlation score which rewards conserved columns occuring in clusters,
      *correl_score_weight* controls its contribution to the total score

  You have to make sure that proper pseudo counts are already assigned before 
  calling this function. You can find a usage example in this documentation.
  This score is not necessarily consistent with the output generated with 
  hhalign, i.e. you take the hhalign output alignment and directly feed it
  into this function with the same profiles and expect an equal score. 
  The reason is that by default, hhalign performs a re-alignment step but the 
  output score actually relates to the initial alignment coming from the 
  Viterbi alignment. To get consistent results, run hhalign with the 
  -norealign flag.

  :param profile_0:     First profile to be scored
  :param profile_1:     Second profile to be scored
  :param aln:           Alignment connecting the two profiles
  :param s_0_idx:       Idx of sequence in *aln* that describes *profile_0*
  :param s_1_idx:       Idx of sequence in *aln* that describes *profile_1*
  :param match_score_offset: Offset which is applied to each column alignment 
                             score
  :param correl_score_weight: Prefactor to control contribution of correlation 
                              score to total score
  :param del_start_penalty_factor: Factor which is applied for each transition 
                                   score starting a deletion
  :param del_extend_penalty_factor: Factor which is applied for each transition 
                                    score extending a deletion
  :param ins_start_penalty_factor: Factor which is applied for each transition 
                                   score starting an insertion
  :param ins_extend_penalty_factor: Factor which is applied for each transition 
                                    score extending an insertion

  :type profile_0:      :class:`ost.seq.ProfileHandle`
  :type profile_1:      :class:`ost.seq.ProfileHandle`
  :type aln:            :class:`ost.seq.AlignmentHandle`
  :type s_0_idx:        :class:`int`
  :type s_1_idx:        :class:`int`
  :type match_score_offset: :class:`float`
  :type correl_score_weight: :class:`float`
  :type del_start_penalty_factor: :class:`float`
  :type del_extend_penalty_factor: :class:`float`
  :type ins_start_penalty_factor: :class:`float`
  :type ins_extend_penalty_factor: :class:`float`

  :raises:  Exception if profiles don't have HMM information assigned or 
            specified sequences in *aln* don't match with profile SEQRES. 
            Potentially set sequence offsets are taken into account.


**Example with pseudo count assignment:**

.. code-block:: python

  from ost import io, seq

  prof_query = io.LoadSequenceProfile("query.hhm")
  prof_tpl = io.LoadSequenceProfile("tpl.hhm")
  aln = io.LoadAlignment("aln.fasta")

  # assign pseudo counts to transition probabilities
  seq.alg.AddTransitionPseudoCounts(prof_query)
  seq.alg.AddTransitionPseudoCounts(prof_tpl)

  # hhblits/hhalign 3 assign different pseudo counts to 
  # query and template. The reason is computational efficiency.
  # The more expensive Angermueller et al. pseudo counts
  # are assigned to the query.
  path_to_crf = "/path/to/hh-suite/data/context_data.crf"
  lib = seq.alg.ContextProfileDB.FromCRF(path_to_crf)
  seq.alg.AddAAPseudoCounts(prof_query, lib)

  # templates are assigned the computationally cheaper pseudo
  # counts derived from a Gonnet substitution matrix
  seq.alg.AddAAPseudoCounts(prof_tpl)

  # assign null model pseudo counts
  # this should be done AFTER you assigned pseudo counts to emission
  # probabilities as this affects the result
  seq.alg.AddNullPseudoCounts(prof_query)
  seq.alg.AddNullPseudoCounts(prof_tpl)

  print("score:", seq.alg.HMMScore(prof_query, prof_tpl, aln, 0, 1))


.. method:: AddNullPseudoCounts(profile)

  Adds pseudo counts to null model in *profile* as implemented in hhalign.
  Conceptually we're mixing the original null model with the frequencies 
  observed in the columns of *profile*. The weight of the original null model
  depends on the neff value of *profile*. This function should be called
  AFTER you already assigned pseudo counts to the emission probabilities
  as this affects the result. 

  :param profile:       Profile to add pseudo counts
  :type profile:        :class:`ost.seq.ProfileHandle`

  :raises:  Exception if profile doesn't have HMM information assigned


.. method:: AddTransitionPseudoCounts(profile, gapb=1.0, gapd=0.15, gape=1.0)

  Adds pseudo counts to the transition probabilities in *profile* as implemented 
  in hhalign with equivalent parameter naming and default parameterization.
  The original transition probabilities are mixed with prior
  probabilities that are controlled by *gapd* and *gape*. Priors:

    * priorM2I = priorM2D = *gapd* * 0.0286
    * priorM2M = 1.0 - priorM2D - priorM2I
    * priorI2I = priorD2D = 1.0 * *gape* / (gape - 1.0 + 1.0/0.75)
    * priorI2M = priorD2M = 1.0 - priorI2I

  Transition probabilities of column i starting from a match state are 
  then estimated with pM2X = (neff[i] - 1) * pM2X + *gape* * priorM2X.
  Starting from an insertion/deletion state we have 
  pI2X = neff_ins[i] * pI2X + *gape* * priorI2X. In the end, all
  probabilities are normalized such that (pM2M, pM2I, pM2D) sum up to one,
  (pI2M, pI2I) sum up to one and (pD2I, pD2D) sum up to one.

  :param profile:       Profile to add pseudo counts
  :type profile:        :class:`ost.seq.ProfileHandle`

  :raises:  Exception if profile doesn't have HMM information assigned


.. method:: AddAAPseudoCounts(profile, a=1.0, b=1.5, c=1.0)

  Adds pseudo counts to the emission probabilities in *profile* by mixing in
  probabilities from the Gonnet matrix as implemented in hhalign with equivalent 
  parameter naming and default parameterization. We only implement the 
  diversity-dependent mode for the mixing factor tau (default in hhalign), which
  for column *i* depends on *neff[i]* , *a* , *b* and *c* .

  :param profile:       Profile to add pseudo counts
  :type profile:        :class:`ost.seq.ProfileHandle`

  :raises:  Exception if profile doesn't have HMM information assigned


.. class:: ContextProfileDB

  Database that contains context profiles which will be used to add pseudo 
  counts as described by 
  Angermueller et al., Bioinformatics (2012) 28, 3240-3247. 

  .. method:: FromCRF(filename)

    Static load function which reads a crf file provided in an hh-suite 
    installation. Default location: "path/to/hhsuite/data/context_data.crf"

    :param filename:    Filename of CRF file
    :type filename:     :class:`str`

  .. method:: Save(filename)

    Saves database in OST-internal binary format which can be loaded faster than
    a crf file.

    :param filename:    Filename to save db
    :type filename:     :class:`str`

  .. method:: Load(filename)

    Static load function that loads database in OST-internal binary format.

    :param filename:    Filename of db
    :type filename:     :class:`str`


.. method:: AddAAPseudoCounts(profile, db, a=0.9, b=4.0, c=1.0)

  Adds pseudo counts to the emission probabilities in *profile* by utilizing 
  context profiles as described in 
  Angermueller et al., Bioinformatics (2012) 28, 3240-3247.
  We only implement the 
  diversity-dependent mode for the mixing factor tau (default in hhalign), which
  for column *i* depends on *neff[i]* , *a* , *b* and *c* .

  :param profile:       Profile to add pseudo counts
  :type profile:        :class:`ost.seq.ProfileHandle`
  :param db:            Database of context profiles
  :type db:             :class:`ContextProfileDB`

  :raises:  Exception if profile doesn't have HMM information assigned
