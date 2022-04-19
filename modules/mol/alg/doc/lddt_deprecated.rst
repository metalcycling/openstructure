:orphan:

lDDT (deprecated)
================================================================================

.. function:: LocalDistDiffTest(model, distance_list, tolerance_list, \
                                sequence_separation=0, \
                                local_lddt_property_string="")
  
  This function counts the number of conserved local contacts between a model
  and a reference structure which is needed to compute the Local Distance
  Difference Test score.

  The Local Distance Difference Test score is a number between zero and one,
  which measures the agreement of  local contacts between a model and a
  reference structure. One means complete agreement, and zero means no agreement
  at all. The calculation of this score does not require any superposition
  between the model and the reference structures.
  
  All distances between atoms in the reference structure that are shorter than a
  certain predefined length (inclusion radius) are compared with the
  corresponding distances in the model structure. If the difference between a
  reference distance and the corresponding model distance is smaller than a
  threshold value (tolerance), that distance is considered conserved. The final
  lDDT score is the fraction of conserved distances. Missing atoms in the model
  structure lead to non-conserved distances (and thus lower the final lDDT
  score).
  
  This function takes as an input a list of distances to be checked for
  conservation. Any number of threshold values  can be specified when the
  function is called. All thresholds are then applied in sequence and the return
  counts are averaged over all threshold values. A sequence separation parameter
  can be passed to the function. If this happens, only distances between
  residues whose separation in sequence is higher than the provided parameter
  are considered when the score is computed.

  If a string is passed as the last parameter, residue-based counts and the
  value of the residue-based Local Distance Difference Test score are saved in
  each ResidueHandle as int and float properties. Specifically, the local
  residue-based lddt score is stored in a float property named as the provided
  string, while the residue-based number of conserved and total distances are
  saved in two int properties named <string>_conserved and <string>_total.

  :param model: the model structure
  :type model: :class:`~ost.mol.EntityView`
  :param distance_list: the list of distances to check for conservation
  :type distance_list: :class:`~ost.mol.alg.GlobalRDMap`
  :param tolerance_list: a list of thresholds used to determine distance
                         conservation
  :param sequence_separation: sequence separation parameter used when computing
                              the score
  :param local_lddt_property_string: the base name for the ResidueHandle
                                     properties that store the local scores

  :returns: a tuple containing the counts of the conserved distances in the
            model and of all the checked distances

.. function:: LocalDistDiffTest(model, reference_list, distance_list, settings)

  Wrapper around :func:`LocalDistDiffTest` above.

  :param model: the model structure
  :type model: :class:`~ost.mol.EntityView`
  :param reference_list: the list of reference structures from which distances were derived
  :type reference_list: :class:`list` of :class:`~ost.mol.EntityView`
  :param distance_list: A residue distance map prepared with :func:`PreparelDDTGlobalRDMap`
    with *reference_list* and *settings* as parameters.
  :type distance_list:  :class:`~ost.mol.alg.GlobalRDMap`
  :param settings: lDDT settings
  :type settings: :class:`~ost.mol.alg.lDDTSettings`

  :returns: the Local Distance Difference Test score (conserved distances
            divided by all the checked distances)
  :rtype:   :class:`float`

.. function:: LocalDistDiffTest(model, target, cutoff, max_dist, \
                                local_lddt_property_string="")

  Wrapper around :func:`LocalDistDiffTest` above using:
  *distance_list* = :func:`CreateDistanceList` with *target* and *max_dist* as
  parameters and *tolerance_list* = [*cutoff*].

  :param model: the model structure
  :type model: :class:`~ost.mol.EntityView`
  :param target: the target structure from which distances are derived
  :type target: :class:`~ost.mol.EntityView`
  :param cutoff: single distance threshold to determine distance conservation
  :type cutoff:  :class:`float`
  :param max_dist: the inclusion radius in Angstroms (to determine which
                   distances are checked for conservation)
  :type max_dist:  :class:`float`
  :param local_lddt_property_string: the base name for the ResidueHandle
                                     properties that store the local scores

  :returns: the Local Distance Difference Test score (conserved distances
            divided by all the checked distances)
  :rtype:   :class:`float`


.. function:: LocalDistDiffTest(alignment, tolerance, radius, ref_index=0, \
                                mdl_index=1)

  Calculates the Local Distance Difference Test score (see previous function)
  starting from an alignment between a reference structure and a model. The
  AlignmentHandle parameter used to provide the  alignment to the function needs
  to have the two structures attached to it. By default the first structure in
  the alignment is considered to be the reference structure, and the second
  structure is taken as the model. This can however be changed by passing the
  indexes of the two structures in the AlignmentHandle as parameters to the
  function.

  .. note::

    This function uses the old implementation of the Local Distance Difference
    Test algorithm and will give slightly different results from the new one.

  :param alignment: an alignment containing the sequences of the reference and
                    of the model structures, with the structures themselves
                    attached
  :type alignment:  :class:`~ost.seq.AlignmentHandle`
  :param tolerance: a list of thresholds used to determine distance conservation
  :param radius: the inclusion radius in Angstroms (to determine which distances
                 are checked for conservation)
  :param ref_index: index of the reference structure in the alignment 
  :param mdl_index: index of the model in the alignment

  :returns: the Local Distance Difference Test score


.. function:: LDDTHA(model, distance_list, sequence_separation=0)

  This function calculates the Local Distance Difference Test, using the same
  threshold values as the GDT-HA test (the default set of thresholds used for
  the lDDT score) (See previous functions). The thresholds are 0.5, 1, 2, and 4
  Angstroms.

  The function only compares the input distance list to the first chain of the
  model structure.

  The local residue-based lDDT score values are stored in the ResidueHandles of
  the model passed to the function in a float property called "locallddt".

  A sequence separation parameter can be passed to the function. If this
  happens, only distances between residues whose separation is higher than the
  provided parameter are considered when computing the score.

  :param model: the model structure
  :type model:  :class:`~ost.mol.EntityView`
  :param distance_list: the list of distances to check for conservation
  :type distance_list:  :class:`~ost.mol.alg.GlobalRDMap`
  :param sequence_separation: sequence separation parameter

  :returns: the Local Distance Difference Test score


.. function:: DistanceRMSDTest(model, distance_list, cap_difference, \
                               sequence_separation=0, \
                               local_drmsd_property_string="")
  
  This function performs a Distance RMSD Test on a provided model, and
  calculates the two values that are necessary to determine the Distance RMSD
  Score, namely the sum of squared distance deviations and the number of
  distances on which the sum was computed.

  The Distance RMSD Test (or DRMSD Test) computes the deviation in the length of
  local contacts between a model and a reference structure and expresses it in
  the form of a score value. The score has an an RMSD-like form, with the
  deviations in the RMSD formula computed as contact distance differences. The
  score is open-ended, with a value of zero meaning complete agreement of local
  contact distances, and a positive value revealing a disagreement of magnitude
  proportional to the score value itself. This score does not require any
  superposition between the model and the reference.
  
  This function processes a list of distances provided by the user, together
  with their length in the reference structure. For each distance that is found
  in the model, its difference with the reference length is computed and used as
  deviation term in the RMSD-like formula.When a distance is not present in the
  model because one or both the atoms are missing, a default deviation value
  provided by the user is used.

  The function only processes distances between atoms that do not belong to the
  same residue, and considers only standard residues in the first chain of the
  model. For residues with symmetric sidechains (GLU, ASP, ARG, VAL, PHE, TYR),
  the naming of the atoms is ambiguous. For these residues, the function
  computes the Distance RMSD Test score that each naming convention would
  generate when considering all non-ambiguous surrounding atoms. The solution
  that gives the lower score is then picked to compute the final Distance RMSD
  Score for the whole model.
  
  A sequence separation parameter can be passed to the function. If this
  happens, only distances between residues whose separation is higher than the
  provided parameter are considered when computing the score.

  If a string is passed as last parameter to the function, the function computes
  the Distance RMSD Score for each residue and saves it as a float property in
  the ResidueHandle, with the passed string as property name. Additionally, the
  actual sum of squared deviations and the number of distances on which it was
  computed are stored as properties in the ResidueHandle. The property names are
  respectively <passed string>_sum (a float property) and <passed string>_count
  (an integer property).

  :param model: the model structure
  :type model:  :class:`~ost.mol.EntityView`
  :param distance_list: the list of distances to check (here we only use the
                        first of the two distance values stored, the second
                        is ignored)
  :type distance_list:  :class:`~ost.mol.alg.GlobalRDMap`
  :param cap_difference: a default deviation value to be used when a distance is
                         not found in the model
  :param sequence_separation: sequence separation parameter
  :param local_ldt_property_string: the base name for the ResidueHandle
                                    properties that store the local scores

  :returns: a tuple containing the sum of squared distance deviations, and the
            number of distances on which it was computed.


.. function:: DRMSD(model, distance_list, cap_difference, sequence_separation=0)

  This function calculates the Distance RMSD Test score (see
  :func:`DistanceRMSDTest`).
  
  The function only considers distances between atoms not belonging to the same
  residue, and only compares the input distance list to the first chain of the
  model structure. It requires, in addition to the model and the list
  themselves, a default deviation value to be used in the DRMSD Test when a
  distance is not found in the model.

  The local Local Distance Difference Test score values are stored in the
  ResidueHandles of the model passed to the function in a float property called
  "localdrmsd".

  A sequence separation parameter can be passed to the function. If this
  happens, only distances between residues whose separation is higher than the
  provided parameter are considered when computing the score.

  :param model: the model structure
  :type model:  :class:`~ost.mol.EntityView`
  :param distance_list: the list of distances as in :func:`DistanceRMSDTest`
  :type distance_list: :class:`~ost.mol.alg.GlobalRDMap`
  :param cap_difference: a default deviation value to be used when a distance is
                         not found in the model
  :param sequence_separation: sequence separation parameter
  :returns: the Distance RMSD Test score


.. function:: CreateDistanceList(reference, radius)
              CreateDistanceListFromMultipleReferences(reference_list, \
                                                       tolerance_list, \
                                                       sequence_separation, \
                                                       radius)

  Both these functions create lists of distances to be checked during a Local
  Distance Difference Test (see description of the functions above).

  .. note::

    These functions process only standard residues present in the first chain of
    the reference structures.

  The only difference between the two functions is that one takes a single
  reference structure and the other a list of reference structures. The
  structures in the list have to be properly prepared before being passed to the
  function. Corresponding residues in the structures must have the same residue
  number, the same chain name, etc. Gaps are allowed and automatically dealt
  with: if information about a distance is present in at least one of the
  structures, it will be considered.

  If a distance between two atoms is shorter than the inclusion radius in all
  structures in which the two atoms are present, it is included in the list.
  However, if the distance is longer than the inclusion radius in at least one
  of the structures, it is not considered to be a local interaction and is
  excluded from the list.

  The multiple-reference function takes care of residues with ambiguous
  symmetric sidechains. To decide which naming convention to use, the function
  computes a Local Distance Difference Test score foreach reference against the
  first reference structure in the list, using only non ambiguously-named atoms.
  It picks then the naming convention that gives the highest score, guaranteeing
  that all references are processed with the correct atom names.

  The cutoff list that will later be used to compute the Local Distance
  Difference Test score and the sequence separation parameter must be passed to
  the multi-reference function. These parameters do not influence the output
  distance list, which always includes all distances within the provided radius
  (to make it consistent with the single-reference corresponding function).
  However, the parameters are used when dealing with the naming convention of
  residues with ambiguous nomenclature.

  :param reference: a reference structure from which distances are derived
  :type reference:  :class:`~ost.mol.EntityView`
  :param reference_list: a list of reference structures from which distances are
                         derived
  :type reference_list:  list of :class:`~ost.mol.EntityView`
  :param tolerance_list: a list of thresholds used to determine distance
                         conservation when computing the lDDT score
  :param sequence_separation: sequence separation parameter used when computing
                              the lDDT score
  :param radius: inclusion radius (in Angstroms) used to determine the distances
                 included in the list
  
  :returns: :class:`~ost.mol.alg.GlobalRDMap`


.. function:: PreparelDDTGlobalRDMap(reference_list, cutoff_list, sequence_separation, max_dist)

  A wrapper around :func:`CreateDistanceList` and
  :func:`CreateDistanceListFromMultipleReferences`. Depending on the length of
  the ``reference_list`` it calls one or the other.

  :param reference_list: a list of reference structures from which distances are
    derived
  :type reference_list:  list of :class:`~ost.mol.EntityView`
  :param max_dist: the inclusion radius in Angstroms (to determine which
                   distances are checked for conservation)
  :type max_dist:  :class:`float`
  :param sequence_separation: sequence separation parameter ie. maximum distance
                              between two sequences.
  :type sequence_separation: :class:`int`
  :returns: :class:`~ost.mol.alg.GlobalRDMap`


.. function:: CleanlDDTReferences(reference_list)

  Prepares references to be used in lDDT calculation. It checks if all references
  has the same chain name and selects this chain for for further calculations.

  .. warning::

    This function modifies the passed *reference_list* list.

  :param reference_list: A list of reference structures from which distances are
                         derived
  :type reference_list:  :class:`list` of :class:`~ost.mol.EntityView`

.. function:: GetlDDTPerResidueStats(model, distance_list, structural_checks, label)

  Get the per-residue statistics from the lDDT calculation.

  :param model: The model structure
  :type model: :class:`~ost.mol.EntityHandle`
  :param distance_list: The list of distances to check for conservation
  :type distance_list: :class:`~ost.mol.alg.GlobalRDMap`
  :param structural_checks: Were structural checks performed on the model?
  :type structural_checks: :class:`bool`
  :param label: Label used for ResidueHandle properties that store the local
                scores.
  :type label: :class:`str`
  :returns: Per-residue local lDDT scores
  :rtype: :class:`list` of :class:`~ost.mol.alg.lDDTLocalScore`


.. function:: PrintlDDTPerResidueStats(scores, structural_checks, cutoffs_length)

  Print per-residue statistics from lDDT calculation.

  :param scores: Local lDDT scores
  :type scores: :class:`list` of :class:`~ost.mol.alg.lDDTLocalScore`
  :param structural_checks: Where structural checks performed on the model?
  :type structural_checks: :class:`bool`
  :param cutoffs_length: Length of the cutoffs list used to calculate lDDT
  :type cutoffs_length: :class:`int`


.. class:: lDDTLocalScore(cname, rname, rnum, is_assessed, quality_problems, \
                          local_lddt, conserved_dist, total_dist)

  Object containing per-residue information about calculated lDDT.

  :param cname: Sets :attr:`cname`
  :param rname: Sets :attr:`rname`
  :param rnum: Sets :attr:`rnum`
  :param is_assessed: Sets :attr:`is_assessed`
  :param quality_problems: Sets :attr:`quality_problems`
  :param local_lddt: Sets :attr:`local_lddt`
  :param conserved_dist: Sets :attr:`conserved_dist`
  :param total_dist: Sets :attr:`total_dist`

  .. attribute:: cname

    Chain name.

    :type: :class:`str`

  .. attribute:: rname

    Residue name.

    :type: :class:`str`

  .. attribute:: rnum

    Residue number.

    :type: :class:`int`

  .. attribute:: is_assessed

    Is the residue taken into account? Yes or No.

    :type: :class:`str`

  .. attribute:: quality_problems

    Does the residue have quality problems?
    No if there are no problems, NA if the problems were not assessed, Yes if
    there are sidechain problems and Yes+ if there are backbone problems.

    :type: :class:`str`

  .. attribute:: local_lddt

    Local lDDT score for residue.

    :type: :class:`float`

  .. attribute:: conserved_dist

    Number of conserved distances.

    :type: :class:`int`

  .. attribute:: total_dist

    Total number of distances.

    :type: :class:`int`

  .. method:: ToString(structural_checks)

    :return: String representation of the lDDTLocalScore object.
    :rtype:  :class:`str`

    :param structural_checks: Where structural checks applied during calculations?
    :type structural_checks: bool

  .. method:: GetHeader(structural_checks, cutoffs_length)

    Get the names of the fields as printed by ToString method.

    :param structural_checks: Where structural checks applied during calculations?
    :type structural_checks: bool
    :param cutoffs_length: Length of the cutoffs list used for calculations
    :type cutoffs_length: int





.. class:: lDDTScorer(reference, model, settings)

  Object to compute lDDT scores using :func:`LocalDistDiffTest` as in
  `Mariani et al. <https://dx.doi.org/10.1093/bioinformatics/btt473>`_.
  
  Example usage.
  
  .. code:: python
  
    #! /bin/env python
    """Run lDDT from within script."""
    from ost.io import LoadPDB
    from ost.mol.alg import (CleanlDDTReferences,
                             lDDTSettings, lDDTScorer)

    ent_full = LoadPDB('3ia3', remote=True)
    model_view = ent_full.Select('cname=A')
    references = [ent_full.Select('cname=C')]

    #
    # Initialize settings with default parameters and print them
    settings = lDDTSettings()
    settings.PrintParameters()

    # Clean up references
    CleanlDDTReferences(references)
    #
    # Calculate lDDT
    scorer = lDDTScorer(references=references, model=model_view, settings=settings)
    print("Global score:", scorer.global_score)
    scorer.PrintPerResidueStats()
  
  :param references: Sets :attr:`references`
  :param model: Sets :attr:`model`
  :param settings: Sets :attr:`settings`
  
  .. attribute:: references
  
    A list of reference structures.
    
    :type: list(:class:`~ost.mol.EntityView`)
  
  .. attribute:: model
  
    A model structure. 
    
    :type: :class:`~ost.mol.EntityView`
    
  .. attribute:: settings
  
    Settings used to calculate lDDT.
    
    :type: :class:`~ost.mol.alg.lDDTSettings`
  
  .. attribute:: global_dist_list
  
    Global map of residue properties.
    
    :type: :class:`~ost.mol.alg.GlobalRDMap`

  .. attribute:: global_score
  
    Global lDDT score. It is calculated as :attr:`conserved_contacts` divided
    by :attr:`total_contacts`.
    
    :type: float

  .. attribute:: conserved_contacts
  
    Number of conserved distances.
  
    :type: int
  
  .. attribute:: total_contacts
  
    Number of total distances.
  
    :type:
  
  .. attribute:: local_scores
  
    Local scores. For each of the residue lDDT is it is calculated as residue
    conserved contacts divided by residue total contacts.
  
    :type: list(:class:`~ost.mol.alg.lDDTLocalScore`)
  
  .. attribute:: is_valid
  
    Is the calculated score valid?
  
    :type: bool
  
  .. method:: PrintPerResidueStats
    
    Print per-residue statistics.


.. class:: UniqueAtomIdentifier(chain, residue_number, residue_name, atom_name)

  Object containing enough information to uniquely identify an atom in a
  structure.

  :param chain: A string containing the name of the chain to which the atom
                belongs
  :param residue_number: The number of the residue to which the atom belongs
  :type residue_number:  :class:`~ost.mol.ResNum`
  :param residue_name: A string containing the name of the residue to which
                       the atom belongs
  :param atom_name: A string containing the name of the atom

  .. method:: GetChainName() 

    Returns the name of the chain to which the atom belongs, as a String  

  .. method:: GetResNum() 

    Returns the number of the residue the atom belongs to, as a
    :class:`~ost.mol.ResNum` object

  .. method:: GetResidueName()
    
     Returns the name of the residue to which the atom belongs, as a String
 
  .. method:: GetAtomName()

     Returns the name of the atom, as a String

  .. method:: GetQualifiedAtomName()

     Returns the qualified name of the atom (the chain name, followed by a
     unique residue identifier and the atom name. For example: "A.GLY2.CA")


.. class:: ResidueRDMap

  Dictionary-like object containing the list of interatomic distances that
  originate from a single residue to be checked during a run of the Local
  Distance Difference Test algorithm
  (key = pair of :class:`UniqueAtomIdentifier`, value = pair of floats
  representing min and max distance observed in the structures used to build
  the map).

.. class:: GlobalRDMap

  Dictionary-like object containing all the :class:`~ost.mol.alg.ResidueRDMap` objects related to all the residues
  (key = :class:`~ost.mol.ResNum`, value = :class:`ResidueRDMap`).

  
.. function:: PrintResidueRDMap(residue_distance_list)

  Prints to standard output all the distances contained in a
  :class:`~ost.mol.alg.ResidueRDMap` object.


.. function:: PrintGlobalRDMap(global_distance_list)

  Prints to standard output all the distances contained in each of the
  :class:`~ost.mol.alg.ResidueRDMap` objects that make up a
  :class:`~ost.mol.alg.GlobalRDMap` object.