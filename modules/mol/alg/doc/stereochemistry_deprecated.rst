:orphan:

Stereochemistry (deprecated)
================================================================================


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