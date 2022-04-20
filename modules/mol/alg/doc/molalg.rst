:mod:`mol.alg <ost.mol.alg>` -- Algorithms for Structures
================================================================================

.. module:: ost.mol.alg
   :synopsis: Algorithms operating on molecular structures

Local Distance Test scores (lDDT, DRMSD)
--------------------------------------------------------------------------------

.. warning::

  The code that comes with
  `Mariani et al. <https://dx.doi.org/10.1093/bioinformatics/btt473>`_ is
  considered deprecated. lDDT has been re-implemented with a focus on
  supporting quaternary structure and compounds beyond the 20 standard
  proteinogenic amino acids. The old code is still available and documented
  :doc:`here <lddt_deprecated>`.

.. note::

  :class:`lddt.lDDTScorer` provides the raw Python API to compute lDDT but
  stereochemistry checks as described in
  `Mariani et al. <https://dx.doi.org/10.1093/bioinformatics/btt473>`_
  must be done seperately. You may want to check out the
  ``compare-structures`` action (:ref:`ost-compare-structures`) to
  compute lDDT with pre-processing and support for quaternary structures.


.. autoclass:: ost.mol.alg.lddt.lDDTScorer
  :members:

.. autoclass:: ost.mol.alg.lddt.SymmetrySettings
  :members:

.. autofunction:: ost.mol.alg.lddt.GetDefaultSymmetrySettings

.. function:: CheckStructure(ent, \
                             bond_table, \
                             angle_table, \
                             nonbonded_table, \
                             bond_tolerance, \
                             angle_tolerance)

  Perform structural checks and filters the structure.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityView`
  :param bond_table: List of bond stereo chemical parameters obtained from
    :class:`~ost.io.StereoChemicalParamsReader` or :func:`FillStereoChemicalParams`
  :type bond_table: :class:`~ost.mol.alg.StereoChemicalParams`
  :param angle_table: List of angle stereo chemical parameters obtained from
    :class:`~ost.io.StereoChemicalParamsReader` or :func:`FillStereoChemicalParams`
  :type angle_table: :class:`~ost.mol.alg.StereoChemicalParams`
  :param nonbonded_table: Information about the clashing distances obtained from
    :class:`~ost.io.StereoChemicalParamsReader` or :func:`FillClashingDistances`
  :type nonbonded_table: :class:`~ost.mol.alg.ClashingDistances`
  :param bond_tolerance: Tolerance in stddev for bonds
  :type bond_tolerance: :class:`float`
  :param angle_tolerance: Tolerance in stddev for angles
  :type angle_tolerance: :class:`float`

    
.. class:: StereoChemicalProps(bond_table, angle_table, nonbonded_table)
  
  Object containing the stereo-chemical properties read form stereochmical_props.txt
  file.

  :param bond_table: Sets :attr:`bond_table`
  :param angle_table: Sets :attr:`angle_table`
  :param nonbonded_table: Sets :attr:`nonbonded_table`

  .. attribute:: bond_table
  
    Object containing bond parameters
    
    :type: :class:`~ost.mol.alg.StereoChemicalParams`

  .. attribute:: angle_table
    
    Object containing angle parameters
    
    :type: :class:`~ost.mol.alg.StereoChemicalParams`

  .. attribute:: nonbonded_table
    
    Object containing clashing distances parameters
    
    :type: :class:`~ost.mol.alg.ClashingDistances`


.. class:: lDDTSettings(radius=15, \
                        sequence_separation=0, \
                        cutoffs=(0.5, 1.0, 2.0, 4.0), \
                        label="locallddt")

  Object containing the settings used for lDDT calculations.

  :param radius: Sets :attr:`radius`.
  :param sequence_separation: Sets :attr:`sequence_separation`.
  :param cutoffs: Sets :attr:`cutoffs`.
  :param label: Sets :attr:`label`.

  .. attribute:: radius

    Distance inclusion radius.

    :type: :class:`float`

  .. attribute:: sequence_separation

    Sequence separation.

    :type: :class:`int`

  .. attribute:: cutoffs

    List of thresholds used to determine distance conservation.

    :type: :class:`list` of :class:`float`

  .. attribute:: label

    The base name for the ResidueHandle properties that store the local scores.

    :type: :class:`str`

  .. method:: PrintParameters()

    Print settings.

  .. method:: ToString()

    :return: String representation of the lDDTSettings object.
    :rtype:  :class:`str`




:mod:`qsscoring <ost.mol.alg.qsscoring>` -- Quaternary Structure (QS) scores
--------------------------------------------------------------------------------

.. automodule:: ost.mol.alg.qsscoring
   :members:
   :synopsis: Scoring of quaternary structures

.. currentmodule:: ost.mol.alg


.. _steric-clashes:

Steric Clashes
--------------------------------------------------------------------------------

The following function detects steric clashes in atomic structures. Two atoms are clashing if their euclidian distance is smaller than a threshold value (minus a tolerance offset). 


.. function:: FilterClashes(entity, clashing_distances, always_remove_bb=False)

  This function filters out residues with non-bonded clashing atoms. If the
  clashing atom is a backbone atom, the complete residue is removed from the
  structure, if the atom is part of the sidechain, only the sidechain atoms are
  removed. This behavior is changed  by the *always_remove_bb* flag: when the
  flag is set to True the whole residue is removed even if a clash is just
  detected in the side-chain.

  The function returns a view containing all elements (residues, atoms) that
  have not been removed from the input structure, plus a
  :class:`~ost.mol.alg.ClashingInfo` object containing information about the
  detected clashes.
  
  Two atoms are defined as clashing if their distance is shorter than the
  reference distance minus a tolerance threshold. The information about the
  clashing distances and the tolerance thresholds for all possible pairs of
  atoms is passed to the function as a parameter.

  Hydrogen and deuterium atoms are ignored by this function.
  
  :param entity: The input entity
  :type entity: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param clashing_distances: information about the clashing distances
  :type clashing_distances: :class:`~ost.mol.alg.ClashingDistances`
  :param always_remove_bb: if set to True, the whole residue is removed even if
                           the clash happens in the side-chain
  :type always_remove_bb:  :class:`bool`

  :returns: A tuple of two elements: The filtered :class:`~ost.mol.EntityView`,
            and a :class:`~ost.mol.alg.ClashingInfo` object


.. function:: CheckStereoChemistry(entity, bond_stats, angle_stats, \
                                   bond_tolerance, angle_tolerance, \
                                   always_remove_bb=False)

  This function filters out residues with severe stereo-chemical violations. If
  the violation involves a backbone atom, the complete residue is removed from
  the structure, if it involves an atom that is part of the sidechain, only the
  sidechain is removed. This behavior is changed  by the *always_remove_bb*
  flag: when the flag is set to True the whole residue is removed even if a
  violation is just detected in the side-chain.

  The function returns a view containing all elements (residues, atoms) that
  have not been removed from the input structure, plus a
  :class:`~ost.mol.alg.StereoChemistryInfo` object containing information about
  the detected stereo-chemical violations.
    
  A violation is defined as a bond length that lies outside of the range:
  [mean_length-std_dev*bond_tolerance, mean_length+std_dev*bond_tolerance] or an
  angle width outside of the range [mean_width-std_dev*angle_tolerance,
  mean_width+std_dev*angle_tolerance ]. The information about the mean lengths
  and widths and the corresponding standard deviations is passed to the function
  using two parameters.

  Hydrogen and deuterium atoms are ignored by this function.

  :param entity: The input entity
  :type entity: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param bond_stats: statistics about bond lengths
  :type bond_stats: :class:`~ost.mol.alg.StereoChemicalParams`
  :param angle_stats: statistics about angle widths
  :type angle_stats: :class:`~ost.mol.alg.StereoChemicalParams`
  :param bond_tolerance: tolerance for bond lengths (in standard deviations)
  :type bond_tolerance:  :class:`float`
  :param angle_tolerance: tolerance for angle widths (in standard deviations)
  :type angle_tolerance:  :class:`float`
  :param always_remove_bb: if set to True, the whole residue is removed even if
                           the clash happens in the side-chain
  :type always_remove_bb:  :class:`bool`

  :returns: A tuple of two elements: The filtered :class:`~ost.mol.EntityView`, and a :class:`~ost.mol.alg.StereoChemistryInfo` object


.. class:: ClashingInfo

  This object is returned by the :func:`FilterClashes` function, and contains
  information about the clashes detected by the function.

  .. method:: GetClashCount()

    :return: number of clashes between non-bonded atoms detected in the
             input structure

  .. method:: GetAverageOffset()

    :return: a value in Angstroms representing the average offset by which
             clashing atoms lie closer than the minimum acceptable distance
             (which of course differs for each possible pair of elements)

  .. method:: GetClashList()

    :return: list of detected inter-atomic clashes
    :rtype:  :class:`list` of :class:`ClashEvent`


.. class:: ClashEvent

  This object contains all the information relative to a single clash detected
  by the :func:`FilterClashes` function

  .. method:: GetFirstAtom()
              GetSecondAtom()

    :return: atoms which clash
    :rtype:  :class:`~ost.mol.alg.UniqueAtomIdentifier`

  .. method:: GetModelDistance()

    :return: distance (in Angstroms) between the two clashing atoms as observed
             in the model

  .. method:: GetAdjustedReferenceDistance()

    :return: minimum acceptable distance (in Angstroms) between the two atoms
             involved in the clash, as defined in :class:`ClashingDistances`


.. class:: StereoChemistryInfo

  This object is returned by the :func:`CheckStereoChemistry` function, and
  contains information about bond lengths and planar angle widths in the
  structure that diverge from the parameters tabulated by Engh and Huber in the
  International Tables of Crystallography. Only elements that diverge from the
  tabulated value by a minimumnumber of standard deviations (defined when the
  CheckStereoChemistry function is called) are reported.

  .. method:: GetBadBondCount()

    :return: number of bonds where a serious violation was detected

  .. method:: GetBondCount()

    :return: total number of bonds in the structure checked by the
             CheckStereoChemistry function

  .. method:: GetAvgZscoreBonds()

    :return: average z-score of all the bond lengths in the structure, computed
             using Engh and Huber's mean and standard deviation values

  .. method:: GetBadAngleCount()

    :return: number of planar angles where a serious violation was detected

  .. method:: GetAngleCount()

    :return: total number of planar angles in the structure checked by the
             CheckStereoChemistry function

  .. method:: GetAvgZscoreAngles()

    :return: average z-score of all the planar angle widths, computed using Engh
             and Huber's mean and standard deviation values.

  .. method:: GetBondViolationList()

     :return: list of bond length violations detected in the structure
     :rtype:  :class:`list` of :class:`~ost.mol.alg.StereoChemicalBondViolation`

  .. method:: GetAngleViolationList()

     :return: list of angle width violations detected in the structure
     :rtype: :class:`list` of :class:`~ost.mol.alg.StereoChemicalAngleViolation`


.. class:: StereoChemicalBondViolation

  This object contains all the information relative to a single detected violation of stereo-chemical parameters in a bond length

  .. method:: GetFirstAtom()
              GetSecondAtom()

    :return: first / second atom of the bond
    :rtype:  :class:`~ost.mol.alg.UniqueAtomIdentifier`

  .. method:: GetBondLength()

    :return: length of the bond (in Angstroms) as observed in the model

  .. method:: GetAllowedRange()

    :return: allowed range of bond lengths (in Angstroms), according to Engh and
             Huber's tabulated parameters and the tolerance threshold used when
             the :func:`CheckStereoChemistry` function was called
    :rtype:  :class:`tuple` (minimum and maximum)


.. class:: StereoChemicalAngleViolation

  This object contains all the information relative to a single detected violation of stereo-chemical parameters in a planar angle width

  .. method:: GetFirstAtom()
              GetSecondAtom()
              GetThirdAtom()

    :return: first / second (vertex) / third atom that defines the planar angle
    :rtype:  :class:`UniqueAtomIdentifier`

  .. method:: GetAngleWidth()

    :return: width of the planar angle (in degrees) as observed in the model

  .. method:: GetAllowedRange()

    :return: allowed range of angle widths (in degrees), according to Engh and
             Huber's tabulated parameters and the tolerance threshold used when
             the :func:`CheckStereoChemistry` function was called
    :rtype:  :class:`tuple` (minimum and maximum)


.. class:: ClashingDistances

  Object containing information about clashing distances between non-bonded atoms

  .. method:: ClashingDistances()

    Creates an empty distance list

  .. method:: SetClashingDistance(ele1, ele2, clash_distance, tolerance)

    Adds or replaces an entry in the list

    :param ele1: string containing the first element's name 
    :param ele2: string containing the second element's name 
    :param clash_distance: minimum clashing distance (in Angstroms)
    :param tolerance: tolerance threshold (in Angstroms)

  .. method:: GetClashingDistance(ele1, ele2)

    :return: reference distance and a tolerance threshold (both in Angstroms)
             for two elements
    :rtype:  :class:`tuple` (minimum clashing distance, tolerance threshold)
    :param ele1: string containing the first element's name 
    :param ele2: string containing the second element's name 

  .. method:: GetAdjustedClashingDistance(ele1, ele2)

    :return: reference distance (in Angstroms) for two elements, already
             adjusted by the tolerance threshold
    :param ele1: string containing the first element's name
    :param ele2: string containing the second element's name

  .. method:: GetMaxAdjustedDistance()
 
    :return: longest clashing distance (in Angstroms) in the list, after
             adjustment with tolerance threshold

  .. method:: IsEmpty()

    :return: True if the list is empty (i.e. in an invalid, useless state)
 
  .. method:: PrintAllDistances()

    Prints all distances in the list to standard output


.. class:: StereoChemicalParams

  Object containing stereo-chemical information about bonds and angles. For each
  item (bond or angle in a specific residue), stores the mean and standard
  deviation

  .. method:: StereoChemicalParams()

    Creates an empty parameter list

  .. method:: SetParam(item, residue, mean, standard_dev)

    Adds or replaces an entry in the list

    :param item: string defining a bond (format: X-Y) or an angle (format:
                 X-Y-Z), where X,Y an Z are atom names
    :param residue: string containing the residue type for this entry
    :param mean: mean bond length (in Angstroms) or angle width (in degrees)
    :param standard_dev: standard deviation of the bond length (in Angstroms)
                         or of the angle width (in degrees)

  .. method GetParam(item, residue)

    :return: entry from the list as set in :meth:`SetParam`
    :rtype:  :class:`tuple` (mean, standard deviation)
    :param item: string as used in :meth:`SetParam`
    :param residue: string as used in :meth:`SetParam`

  .. method ContainsParam(item, residue)

    :return: True if a specific entry is present in the list, False if not
    :param item: string as used in :meth:`SetParam`
    :param residue: string as used in :meth:`SetParam`

  .. method:: IsEmpty()

    :return: True if the list is empty (i.e. in an invalid, useless state)
 
  .. method:: PrintAllParameters()

    Prints all entries in the list to standard output  


.. function:: FillClashingDistances(file_content)
              FillBondStereoChemicalParams(file_content)
              FillAngleStereoChemicalParams(file_content)

  These three functions fill a list of reference clashing distances, a list of
  stereo-chemical parameters for bonds and a list of stereo-chemical parameters
  for angles, respectively, starting from the content of a parameter file.

  :param file_content: list of lines from the parameter file
  :type file_content:  :class:`list` of :class:`str`

  :rtype: :class:`~ost.mol.alg.ClashingDistances` or
          :class:`~ost.mol.alg.StereoChemicalParams`


.. function:: FillClashingDistancesFromFile(filename)
              FillBondStereoChemicalParamsFromFile(filename)
              FillAngleStereoChemicalParamsFromFile(filename)

  These three functions fill a list of reference clashing distances, a list of
  stereo-chemical parameters for bonds and a list of stereo-chemical parameters
  for angles, respectively, starting from a file path.

  :param filename: path to parameter file
  :type filename:  :class:`str`

  :rtype: :class:`~ost.mol.alg.ClashingDistances` or
          :class:`~ost.mol.alg.StereoChemicalParams`


.. function:: DefaultClashingDistances()
              DefaultBondStereoChemicalParams()
              DefaultAngleStereoChemicalParams()

  These three functions fill a list of reference clashing distances, a list of
  stereo-chemical parameters for bonds and a list of stereo-chemical parameters
  for angles, respectively, using the default parameter files distributed with
  OpenStructure.

  :rtype: :class:`~ost.mol.alg.ClashingDistances` or
          :class:`~ost.mol.alg.StereoChemicalParams`


.. function:: ResidueNamesMatch(probe, reference)

  The function requires a reference structure and a probe structure. The
  function checks that all the residues in the reference structure that appear
  in the probe structure (i.e., that have the same ResNum) are of the same
  residue type. Chains are compared by order, not by chain name (i.e.: the first
  chain of the reference will be compared with the first chain of the probe
  structure, etc.)

  :param probe: the structure to test
  :type probe:  :class:`~ost.mol.EntityView`
  :param reference: the reference structure
  :type reference:  :class:`~ost.mol.EntityView`

  :return: True if the residue names are the same, False otherwise


Superposing structures
--------------------------------------------------------------------------------

.. autofunction:: Superpose

.. autofunction:: ParseAtomNames

.. autofunction:: MatchResidueByNum

.. autofunction:: MatchResidueByIdx

.. autofunction:: MatchResidueByLocalAln

.. autofunction:: MatchResidueByGlobalAln

.. class:: SuperpositionResult

  .. attribute:: rmsd

    RMSD of the superposed entities.

  .. attribute:: view1
                 view2

    Two :class:`~ost.mol.EntityView` used in superposition (not set if methods
    with :class:`~ost.geom.Vec3List` used).

  .. attribute:: transformation

    Transformation (:class:`~ost.geom.Mat4`) used to map :attr:`view1` onto
    :attr:`view2`.

  .. attribute:: fraction_superposed
                 rmsd_superposed_atoms
                 ncycles

    For iterative superposition (:func:`IterativeSuperposeSVD`): fraction and
    RMSD of atoms that were superposed with a distance below the given
    threshold and the number of iteration cycles performed.

.. method:: SuperposeSVD(view1, view2, apply_transform=True)
            SuperposeSVD(list1, list2)

  Superposition of two sets of atoms minimizing RMSD using a classic SVD based
  algorithm.

  Note that the atom positions in the view are taken blindly in the order in
  which the atoms appear.

  :param view1: View on the model entity
  :type view1:  :class:`~ost.mol.EntityView`
  :param view2: View on the reference entity
  :type view2:  :class:`~ost.mol.EntityView`
  :param list1: List of atom positions for model entity
  :type list1:  :class:`~ost.geom.Vec3List`
  :param list2: List of atom positions for reference entity
  :type list2:  :class:`~ost.geom.Vec3List`
  :param apply_transform: If True, the superposition transform is applied to
                          the (full!) entity handle linked to *view1*.
  :type apply_transform:  :class:`bool`

  :return: An instance of :class:`SuperpositionResult`.

.. method:: IterativeSuperposeSVD(view1, view2, max_iterations=5, \
                                  distance_threshold=3.0, apply_transform=True)
            IterativeSuperposeSVD(list1, list2, max_iterations=5, \
                                  distance_threshold=3.0)

  Iterative superposition of two sets of atoms. In each iteration cycle, we
  keep a fraction of atoms with distances below *distance_threshold* and get
  the superposition considering only those atoms.

  Note that the atom positions in the view are taken blindly in the order in
  which the atoms appear.

  :param view1: View on the model entity
  :type view1:  :class:`~ost.mol.EntityView`
  :param view2: View on the reference entity
  :type view2:  :class:`~ost.mol.EntityView`
  :param list1: List of atom positions for model entity
  :type list1:  :class:`~ost.geom.Vec3List`
  :param list2: List of atom positions for reference entity
  :type list2:  :class:`~ost.geom.Vec3List`
  :param max_iterations: Max. number of iterations to be performed
  :type max_iterations:  :class:`int`
  :param distance_threshold: Distance threshold defining superposed atoms
  :type distance_threshold:  :class:`float`
  :param apply_transform: If True, the superposition transform is applied to
                          the (full!) entity handle linked to *view1*.
  :type apply_transform:  :class:`bool`

  :return: An instance of :class:`SuperpositionResult`.

  :raises: Exception if atom counts do not match or if less than 3 atoms.

.. method:: CalculateRMSD(view1, view2, transformation=geom.Mat4())

  :return: RMSD of atom positions (taken blindly in the order in which the
           atoms appear) in the two given views.
  :rtype:  :class:`float`

  :param view1: View on the model entity
  :type view1:  :class:`~ost.mol.EntityView`
  :param view2: View on the reference entity
  :type view2:  :class:`~ost.mol.EntityView`
  :param transformation: Optional transformation to apply on each atom position
                         of *view1*.
  :type transformation:  :class:`~ost.geom.Mat4`


Algorithms on Structures
--------------------------------------------------------------------------------

.. method:: Accessibility(ent, probe_radius=1.4, include_hydrogens=False,\
                          include_hetatm=False, include_water=False,\
                          oligo_mode=False, selection="", asa_abs="asaAbs",\
                          asa_rel="asaRel", asa_atom="asaAtom", \
                          algorithm = NACCESS)
            
  Calculates the accesssible surface area for ever atom in *ent*. The algorithm
  mimics the behaviour of the bindings available for the NACCESS and DSSP tools 
  and has been tested to reproduce the numbers accordingly.

  :param ent:           Entity on which to calculate the surface
  :type ent:            :class:`~ost.mol.EntityView` /
                        :class:`~ost.mol.EntityHandle`

  :param probe_radius:  Radius of probe to determine accessible surface area
  :type probe_radius:  :class:`float`

  :param include_hydrogens: Whether to include hydrogens in the solvent
                            accessibility calculations. By default, 
                            every atom with ele=H,D is simply neglected.
  :type include_hydrogens:  :class:`bool`

  :param include_hetatms: Whether to include atoms flagged as hetatms 
                          , i.e. ligands, in the solvent
                          accessibility calculations. They are neglected 
                          by default.
  :type include_hetatms:  :class:`bool`

  :param include_water: Whether to include water in the solvent
                        accessibility calculations. By default, 
                        every residue with name "HOH" is neglected.
  :type include_water:  :class:`bool`

  :param oligo_mode:    A typical used case of accessibility calculations is to 
                        determine the solvent accessibility of a full complex
                        and then the accessibility of each chain individually.
                        Lots of calculations can be cached because only the 
                        values of the atoms close to an interface change.
                        This is exactly what happens when you activate the 
                        oligo mode. It returns exactly the same value but adds,
                        additionally to the values estimated in full complex, 
                        the values from each individual chain as float 
                        properties on every residue and atom. Example for atom 
                        accessible surface if the according property name is set 
                        to "asaAtom": Accessibility in the full complex is 
                        stored as "asaAtom", the accessibility when only 
                        considering that particular chain is stored as 
                        "asaAtom_single_chain".
                        The other properties follow the same naming scheme.
  :type oligo_mode:     :class:`bool`

  :param selection:     Selection statement, that gets applied on *ent* before 
                        doing anything. Everything that is not selected is 
                        neglected. The default value of "" results in no
                        selection at all.
  :type selection:      :class:`str`

  :param asa_abs:       Float property name to assign the summed solvent
                        accessible surface from each atom to a residue.
  :type asa_abs:        :class:`str`

  :param asa_rel:       Float property name to assign the relative solvent 
                        accessibility to a residue. This is the absolute 
                        accessibility divided by the maximum solvent 
                        accessibility of that particular residue. 
                        This maximum solvent accessibility is dependent on
                        the chosen :class:`AccessibilityAlgorithm`.
                        Only residues of the 20 standarad amino acids 
                        can be handled.

                        * In case of the **NACCESS** algorithm you can expect
                          a value in range [0.0, 100.0] and a value of -99.9
                          for non standard residues.

                        * In case of the **DSSP** algorithm you can expect a
                          value in range [0.0, 1.0], no float property is assigned
                          in case of a non standard residue.

  :type asa_rel:       :class:`str`

  :param asa_atom:      Float property name to assign the solvent accessible 
                        area to each atom.
  :type asa_atom:       :class:`str`   

  :param algorithm:     Specifies the used algorithm for solvent accessibility
                        calculations

  :type algorithm:      :class:`AccessibilityAlgorithm`   

  

  :return: The summed solvent accessibilty of each atom in *ent*.

.. class:: AccessibilityAlgorithm

  The accessibility algorithm enum specifies the algorithm used by
  the respective tools. Available are:

    NACCESS, DSSP



.. method:: AssignSecStruct(ent)

  Assigns secondary structures to all residues based on hydrogen bond patterns
  as described by DSSP.

  :param ent:           Entity on which to assign secondary structures
  :type ent:            :class:`~ost.mol.EntityView` /
                        :class:`~ost.mol.EntityHandle`



.. class:: FindMemParam

  Result object for the membrane detection algorithm described below

  .. attribute:: axis

    initial search axis from which optimal membrane slab could be found

  .. attribute:: tilt_axis

    Axis around which we tilt the membrane starting from the initial axis

  .. attribute:: tilt

    Angle to tilt around tilt axis

  .. attribute:: angle

    After the tilt operation we perform a rotation around the initial axis
    with this angle to get the final membrane axis

  .. attribute:: membrane_axis

    The result of applying the tilt and rotation procedure described above.
    The membrane_axis is orthogonal to the membrane plane and has unit length.

  .. attribute:: pos

    Real number that describes the membrane center point. To get the actual
    position you can do: pos * membrane_axis

  .. attribute:: width

    Total width of the membrane in A

  .. attribute:: energy

    Pseudo energy of the implicit solvation model

  .. attribute:: membrane_asa

    Membrane accessible surface area

  .. attribute:: membrane_representation

    Dummy atoms that represent the membrane. This entity is only valid if
    the according flag has been set to True when calling FindMembrane.


.. method:: FindMembrane(ent, assign_membrane_representation=True, fast=False)

  Estimates the optimal membrane position of a protein by using an implicit 
  solvation model. The original algorithm and the used energy function are 
  described in: Lomize AL, Pogozheva ID, Lomize MA, Mosberg HI (2006) 
  Positioning of proteins in membranes: A computational approach.

  There are some modifications in this implementation and the procedure is
  as follows:

  * Initial axis are constructed that build the starting point for initial 
    parameter grid searches.

  * For every axis, the protein is rotated so that the axis builds the z-axis

    * In order to exclude internal hydrophilic pores, only the outermost atoms
      with respect the the z-axis enter an initial grid search
    * The width and position of the membrane is optimized for different 
      combinations of tilt and rotation angles (further described in
      :class:`FindMemParam`). The top 20 parametrizations 
      (only top parametrization if *fast* is True) are stored for further 
      processing.

  * The 20 best membrane parametrizations from the initial grid search 
    (only the best if *fast* is set to True) enter a final 
    minimization step using a Levenberg-Marquardt minimizer. 


  :param ent:           Entity of a transmembrane protein, you'll get weird 
                        results if this is not the case. The energy term
                        of the result is typically a good indicator whether
                        *ent* is an actual transmembrane protein. The 
                        following float properties will be set on the atoms:

                        * 'asaAtom' on all atoms that are selected with 
                          ent.Select('peptide=true and ele!=H') as a result
                          of envoking :meth:`Accessibility`.
                        * 'membrane_e' the contribution of the potentially 
                          membrane facing atoms to the energy function. 
                          
  :type ent:            :class:`ost.mol.EntityHandle` / :class:`ost.mol.EntityView`

  :param assign_membrane_representation: Whether to construct a membrane 
                                         representation using dummy atoms

  :type assign_membrane_representation: :class:`bool`

  :param fast:          If set to false, the 20 best results of the initial grid
                        search undergo a Levenberg-Marquardt minimization and
                        the parametrization with optimal minimized energy is 
                        returned. 
                        If set to yes, only the best result of the initial grid
                        search is selected and returned after 
                        Levenberg-Marquardt minimization.

  :returns:             The results object
  :rtype:               :class:`ost.mol.alg.FindMemParam`

.. _traj-analysis:

Trajectory Analysis
--------------------------------------------------------------------------------

This is a set of functions used for basic trajectory analysis such as extracting
positions, distances, angles and RMSDs. The organization is such that most
functions have their counterpart at the individual :class:`frame level
<ost.mol.CoordFrame>` so that they can also be called on one frame instead of
the whole trajectory.

All these functions have a "stride" argument that defaults to stride=1, which is
used to skip frames in the analysis.


.. function:: SuperposeFrames(frames, sel, from=0, to=-1, ref=-1)

  This function superposes the frames of the given coord group and returns them
  as a new coord group.
  
  :param frames: The source coord group.
  :type frames: :class:`~ost.mol.CoordGroupHandle`
  :param sel: An entity view containing the selection of atoms to be used for     
    superposition. If set to an invalid view, all atoms in the coord group are 
    used.
  :type sel: :class:`ost.mol.EntityView`
  :param from: index of the first frame
  :param to: index of the last frame plus one. If set to -1, the value is set to 
     the number of frames in the coord group
  :param ref: The index of the reference frame to use for superposition. If set 
     to -1, the each frame is superposed to the previous frame.
     
  :returns: A newly created coord group containing the superposed frames.

.. function:: SuperposeFrames(frames, sel, ref_view, from=0, to=-1)

  Same as SuperposeFrames above, but the superposition is done on a reference
  view and not on another frame of the trajectory.
  
  :param frames: The source coord group.
  :type frames: :class:`~ost.mol.CoordGroupHandle`
  :param sel: An entity view containing the selection of atoms of the frames to be used for     
    superposition.
  :type sel: :class:`ost.mol.EntityView`
  :param ref_view: The reference view on which the frames will be superposed. The number
    of atoms in this reference view should be equal to the number of atoms in sel.
  :type ref_view: :class:`ost.mol.EntityView`
  :param from: index of the first frame
  :param to: index of the last frame plus one. If set to -1, the value is set to 
     the number of frames in the coord group     
  
  :returns: A newly created coord group containing the superposed frames.


.. function:: AnalyzeAtomPos(traj, atom1, stride=1)

  This function extracts the position of an atom from a trajectory. It returns 
  a vector containing the position of the atom for each analyzed frame.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.
  
  
.. function:: AnalyzeCenterOfMassPos(traj, sele, stride=1)

  This function extracts the position of the center-of-mass of a selection 
  (:class:`~ost.mol.EntityView`) from a trajectory and returns it as a vector.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param sele: The selection from which the center of mass is computed
  :type sele: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.
  

.. function:: AnalyzeDistanceBetwAtoms(traj, atom1, atom2, stride=1)

  This function extracts the distance between two atoms from a trajectory 
  and returns it as a vector.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The first :class:`~ost.mol.AtomHandle`.
  :param atom2: The second :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.  


.. function:: AnalyzeAngle(traj, atom1, atom2, atom3, stride=1)

  This function extracts the angle between three atoms from a trajectory 
  and returns it as a vector. The second atom is taken as being the central
  atom, so that the angle is between the vectors (atom1.pos-atom2.pos) and
  (atom3.pos-atom2.pos).
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The first :class:`~ost.mol.AtomHandle`.
  :param atom2: The second :class:`~ost.mol.AtomHandle`.
  :param atom3: The third :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.


.. function:: AnalyzeDihedralAngle(traj, atom1, atom2, atom3, atom4, stride=1)

  This function extracts the dihedral angle between four atoms from a trajectory 
  and returns it as a vector. The angle is between the planes containing the  
  first three and the last three atoms.
  
  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param atom1: The first :class:`~ost.mol.AtomHandle`.
  :param atom2: The second :class:`~ost.mol.AtomHandle`.
  :param atom3: The third :class:`~ost.mol.AtomHandle`.
  :param atom4: The fourth :class:`~ost.mol.AtomHandle`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.

.. function:: AnalyzeDistanceBetwCenterOfMass(traj, sele1, sele2, stride=1)

  This function extracts the distance between the center-of-mass of two 
  selections (:class:`~ost.mol.EntityView`) from a trajectory and returns it as 
  a vector.

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param sele1: The selection from which the first center of mass is computed
  :type sele1: :class:`~ost.mol.EntityView`.
  :param sele2: The selection from which the second center of mass is computed
  :type sele2: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.


.. function:: AnalyzeRMSD(traj, reference_view, sele_view, stride=1)

  This function extracts the rmsd between two :class:`~ost.mol.EntityView` and 
  returns it as a vector. The views don't have to be from the same entity. The 
  reference positions are taken directly from the reference_view, evaluated only 
  once. The positions from the sele_view are evaluated for each frame.
  If you want to compare to frame i of the trajectory t, first use 
  t.CopyFrame(i) for example:
 
  .. code-block:: python
     
    eh = io.LoadPDB(...)
    t = io.LoadCHARMMTraj(eh, ...)
    sele = eh.Select(...)
    t.CopyFrame(0)
    mol.alg.AnalyzeRMSD(t, sele, sele)

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param reference_view: The selection used as reference structure
  :type reference_view: :class:`~ost.mol.EntityView`.
  :param sele_view: The selection compared to the reference_view
  :type sele_view: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
      consecutive frames analyzed.

.. function:: AnalyzeMinDistance(traj, view1, view2, stride=1)

  This function extracts the minimal distance between two sets of atoms 
  (view1 and view2) for each frame in a trajectory and returns it as a vector.

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param view1: The first group of atoms
  :type view1: :class:`~ost.mol.EntityView`.
  :param view2: The second group of atoms
  :type view2: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.
     
.. function:: AnalyzeMinDistanceBetwCenterOfMassAndView(traj, view_cm, view_atoms, stride=1)

  This function extracts the minimal distance between a set of atoms 
  (view_atoms) and the center of mass of a second set of atoms (view_cm) 
  for each frame in a trajectory and returns it as a vector.

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param view_cm: The group of atoms from which the center of mass is taken
  :type view_cm: :class:`~ost.mol.EntityView`.
  :param view_atoms: The second group of atoms
  :type view_atoms: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.

.. function:: AnalyzeAromaticRingInteraction(traj, view_ring1, view_ring2, stride=1)

  This function is a crude analysis of aromatic ring interactions. For each frame in a trajectory, it calculates
  the minimal distance between the atoms in one view and the center of mass of the other
  and vice versa, and returns the minimum between these two minimal distances.
  Concretely, if the two views are the heavy atoms of two rings, then it returns the minimal
  center of mass - heavy atom distance betweent he two rings

  :param traj: The trajectory to be analyzed.
  :type traj: :class:`~ost.mol.CoordGroupHandle`
  :param view_ring1: First group of atoms
  :type view_ring1: :class:`~ost.mol.EntityView`.
  :param view_ring2: Second group of atoms
  :type view_ring2: :class:`~ost.mol.EntityView`.
  :param stride: Size of the increment of the frame's index between two 
     consecutive frames analyzed.  


:mod:`helix_kinks <ost.mol.alg.helix_kinks>` -- Algorithms to calculate Helix Kinks
---------------------------------------------------------------------------------------------------------------

.. automodule:: ost.mol.alg.helix_kinks
   :members:

:mod:`trajectory_analysis <ost.mol.alg.trajectory_analysis>` -- DRMSD, pairwise distances and more
---------------------------------------------------------------------------------------------------------------

.. automodule:: ost.mol.alg.trajectory_analysis
   :members:

:mod:`structure_analysis <ost.mol.alg.structure_analysis>` -- Functions to analyze structures
---------------------------------------------------------------------------------------------------------------

.. automodule:: ost.mol.alg.structure_analysis
   :members:


.. _mapping-functions:

Mapping functions
--------------------------------------------------------------------------------

.. currentmodule:: ost.mol.alg

The following functions help to convert one residue into another by reusing as
much as possible from the present atoms. They are mainly meant to map from
standard amino acid to other standard amino acids or from modified amino acids
to standard amino acids.

.. function:: CopyResidue(src_res, dst_res, editor)

  Copies the atoms of ``src_res`` to ``dst_res`` using the residue names
  as guide to decide which of the atoms should be copied. If ``src_res`` and
  ``dst_res`` have the same name, or ``src_res`` is a modified version of
  ``dst_res`` (i.e. have the same single letter code), :func:`CopyConserved`
  will be called, otherwise :func:`CopyNonConserved`.

  If a CBeta atom wasn't already copied from ``src_res``, a new one at a
  reconstructed position will be added to ``dst_res`` if it is not ``GLY`` and
  all backbone positions are available to do it.

  :param src_res: The source residue
  :type src_res: :class:`~ost.mol.ResidueHandle`
  :param dst_res: The destination residue (expected to be a standard amino acid)
  :type dst_res: :class:`~ost.mol.ResidueHandle`
  :param editor: Editor used to modify *dst_res*.
  :type editor: :class:`~ost.mol.XCSEditor`

  :returns: True if the residue could be copied as a conserved residue,
            False if it had to fallback to :func:`CopyNonConserved`.

.. function:: CopyConserved(src_res, dst_res, editor)

  Copies the atoms of ``src_res`` to ``dst_res`` assuming that the parent
  amino acid of ``src_res`` (or ``src_res`` itself) are identical to ``dst_res``.

  If ``src_res`` and ``dst_res`` are identical, all heavy atoms are copied
  to ``dst_res``. If ``src_res`` is a modified version of ``dst_res`` and the
  modification is a pure addition (e.g. the phosphate group of phosphoserine),
  the modification is stripped off and all other heavy atoms are copied to
  ``dst_res``. If the modification is not a pure addition, it falls back to
  :func:`CopyNonConserved`.

  Additionally, the selenium atom of ``MSE`` is converted to sulphur to map
  ``MSE`` to ``MET``.

  :param src_res: The source residue
  :type src_res: :class:`~ost.mol.ResidueHandle`
  :param dst_res: The destination residue (expected to be a standard amino acid)
  :type dst_res: :class:`~ost.mol.ResidueHandle`
  :param editor: Editor used to modify *dst_res*.
  :type editor: :class:`~ost.mol.XCSEditor`

  :returns: A tuple of bools stating whether the residue could be copied without
            falling back to :func:`CopyNonConserved` and whether the CBeta atom
            was copied from ``src_res`` to ``dst_res``.

.. function:: CopyNonConserved(src_res, dst_res, editor)

  Copies the heavy backbone atoms and CBeta (except for ``GLY``) of ``src_res``
  to ``dst_res``.

  :param src_res: The source residue
  :type src_res: :class:`~ost.mol.ResidueHandle`
  :param dst_res: The destination residue (expected to be a standard amino acid)
  :type dst_res: :class:`~ost.mol.ResidueHandle`
  :param editor: Editor used to modify *dst_res*.
  :type editor: :class:`~ost.mol.XCSEditor`

  :returns: A tuple of bools as in :func:`CopyConserved` with the first bool
            always being False.


Molecular Checker (Molck)
--------------------------------------------------------------------------------

Programmatic usage
##################

Molecular Checker (Molck) could be called directly from the code using Molck
function:

.. code-block:: python

  #! /bin/env python

  """Run Molck with Python API.


  This is an exemplary procedure on how to run Molck using Python API which is
  equivalent to the command line:

  molck <PDB PATH> --rm=hyd,oxt,nonstd,unk \
                   --fix-ele --out=<OUTPUT PATH> \
                   --complib=<PATH TO compounds.chemlib>
  """

  from ost.io import LoadPDB, SavePDB
  from ost.mol.alg import MolckSettings, Molck
                         
  from ost.conop import CompoundLib


  pdbid = "<PDB PATH>"
  lib = CompoundLib.Load("<PATH TO compounds.chemlib>")

  # Using Molck function
  ent = LoadPDB(pdbid)
  ms = MolckSettings(rm_unk_atoms=True,
                     rm_non_std=True,
                     rm_hyd_atoms=True,
                     rm_oxt_atoms=True,
                     rm_zero_occ_atoms=False,
                     colored=False,
                     map_nonstd_res=False,
                     assign_elem=True)
  Molck(ent, lib, ms)
  SavePDB(ent, "<OUTPUT PATH>")

It can also be split into subsequent commands for greater controll:

.. code-block:: python

  #! /bin/env python

  """Run Molck with Python API.


  This is an exemplary procedure on how to run Molck using Python API which is
  equivalent to the command line:

  molck <PDB PATH> --rm=hyd,oxt,nonstd,unk \
                   --fix-ele --out=<OUTPUT PATH> \
                   --complib=<PATH TO compounds.chemlib>
  """

  from ost.io import LoadPDB, SavePDB
  from ost.mol.alg import (RemoveAtoms, MapNonStandardResidues,
                           CleanUpElementColumn)
  from ost.conop import CompoundLib


  pdbid = "<PDB PATH>"
  lib = CompoundLib.Load("<PATH TO compounds.chemlib>")
  map_nonstd = False

  # Using function chain
  ent = LoadPDB(pdbid)
  if map_nonstd:
      MapNonStandardResidues(lib=lib, ent=ent)

  RemoveAtoms(lib=lib,
              ent=ent,
              rm_unk_atoms=True,
              rm_non_std=True,
              rm_hyd_atoms=True,
              rm_oxt_atoms=True,
              rm_zero_occ_atoms=False,
              colored=False)

  CleanUpElementColumn(lib=lib, ent=ent)
  SavePDB(ent, "<OUTPUT PATH>")

API
###

.. class:: MolckSettings(rm_unk_atoms=False, rm_non_std=False, \
                         rm_hyd_atoms=True, rm_oxt_atoms=False, \
                         rm_zero_occ_atoms=False, colored=False, \
                         map_nonstd_res=True, assign_elem=True)

  Stores settings used for Molecular Checker.

  :param rm_unk_atoms: Sets :attr:`rm_unk_atoms`.
  :param rm_non_std: Sets :attr:`rm_non_std`.
  :param rm_hyd_atoms: Sets :attr:`rm_hyd_atoms`.
  :param rm_oxt_atoms: Sets :attr:`rm_oxt_atoms`.
  :param rm_zero_occ_atoms: Sets :attr:`rm_zero_occ_atoms`.
  :param colored: Sets :attr:`colored`.
  :param map_nonstd_res: Sets :attr:`map_nonstd_res`.
  :param assign_elem: Sets :attr:`assign_elem`.

  .. attribute:: rm_unk_atoms

    Remove unknown and atoms not following the nomenclature.
    
    :type: :class:`bool`

  .. attribute:: rm_non_std

    Remove all residues not one of the 20 standard amino acids
    
    :type: :class:`bool`

  .. attribute:: rm_hyd_atoms

    Remove hydrogen atoms
    
    :type: :class:`bool`

  .. attribute:: rm_oxt_atoms

    Remove terminal oxygens
    
    :type: :class:`bool`

  .. attribute:: rm_zero_occ_atoms

    Remove atoms with zero occupancy
    
    :type: :class:`bool`

  .. attribute:: colored

    Whether output should be colored
    
    :type: :class:`bool`

  .. attribute:: map_nonstd_res

    Maps modified residues back to the parent amino acid, for example
    MSE -> MET, SEP -> SER
    
    :type: :class:`bool`

  .. attribute:: assign_elem

    Clean up element column
    
    :type: :class:`bool`


  .. method:: ToString()

    :return: String representation of the MolckSettings.
    :rtype:  :class:`str`

.. warning::

  The API here is set such that the functions modify the passed structure *ent*
  in-place. If this is not ok, please work on a copy of the structure.

.. function:: Molck(ent, lib, settings, [prune=True])

  Runs Molck on provided entity.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`
  :param settings: Molck settings
  :type settings: :class:`MolckSettings`
  :param prune: Whether to remove residues/chains that don't contain atoms 
                anymore after Molck cleanup
  :type prune: :class:`bool` 



.. function:: MapNonStandardResidues(ent, lib)

  Maps modified residues back to the parent amino acid, for example MSE -> MET.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`

.. function:: RemoveAtoms(ent, lib, rm_unk_atoms=False, rm_non_std=False, \
                          rm_hyd_atoms=True, rm_oxt_atoms=False, \
                          rm_zero_occ_atoms=False, colored=False)

  Removes atoms and residues according to some criteria.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`
  :param rm_unk_atoms: See :attr:`MolckSettings.rm_unk_atoms`
  :param rm_non_std: See :attr:`MolckSettings.rm_non_std`
  :param rm_hyd_atoms: See :attr:`MolckSettings.rm_hyd_atoms`
  :param rm_oxt_atoms: See :attr:`MolckSettings.rm_oxt_atoms`
  :param rm_zero_occ_atoms: See :attr:`MolckSettings.rm_zero_occ_atoms`
  :param colored: See :attr:`MolckSettings.colored`

.. function:: CleanUpElementColumn(ent, lib)

  Clean up element column.

  :param ent: Structure to check
  :type ent: :class:`~ost.mol.EntityHandle`
  :param lib: Compound library
  :type lib: :class:`~ost.conop.CompoundLib`
