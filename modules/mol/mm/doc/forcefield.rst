Forcefields
================================================================================

.. currentmodule:: ost.mol.mm

The forcefields are a dump for interactions with their parameters, but also
for atom specific information or residue definitions in the form of a 
:class:`BuildingBlock`. Objects for modifying residues can be set in form of 
:class:`BlockModifier` or :class:`HydrogenConstructor`.
They're also involved in dealing with the naming mess we're observing in the molecular mechanics
community and contain definable renaming rules that can be applied on an
:class:`EntityHandle` for renaming from e.g. PDB standard to the forcefield
specific standard. The standard forcefields in OpenStructure are loaded from
the files provided by Gromacs and the "standard" naming is therefore the same.
This has implications for controlling the protonation states for histidine.
If you e.g. want to enforce a d-protonated histidine you have to name
it HISD. Further reading can be found in the 
`Gromacs Manual <http://www.gromacs.org/Documentation/Manual>`_ 

Loading the standard forcefields provided by OpenStructure
--------------------------------------------------------------------------------

.. function:: LoadCHARMMForcefield()

   Loads the CHARMM27 forcefield read from Gromacs
   
   :returns: The loaded :class:`Forcefield`


.. function:: LoadAMBERForcefield()

   Loads the AMBER03 forcefield read from Gromacs
   
   :returns: The loaded :class:`Forcefield`


Reading forcefields
--------------------------------------------------------------------------------

.. class:: FFReader(base_dir)

  The :class:`FFReader` builds up a :class:`Forcefield`, that gets updated with
  every call to the read functions. If the read files contain preprocessor 
  statements as they are used in Gromacs, they will be applied to all
  subsequent lines read in. Parsed preprocessor statements are:
  #include, #define, #ifdef, #ifndef, #else and #endif

  Note that this class is rather experimental. It has nevertheless been 
  thoroughly tested for loading the CHARMM and AMBER forcefields in the
  Gromacs format. The reader is capable of resolving the preprocessor statements
  as they are used in Gromacs.

  :param base_dir:      Base path of the reader.
                        All loaded files must be defined relative to this base 
                        path.

  :type base_dir:       :class:`str`

  .. method:: ReadGromacsForcefield()

    Searches and reads the forcefield.itp and atomtypes.atp files 
    in the **base_dir** given at initialization. All atom specific 
    informations and bonded as well as nonbonded forces are read 
    this way.

  .. method:: ReadResidueDatabase(basename)

    Searches and reads all files belonging the the residue database
    defined by **basename**. With *basename=aminoacids* this function
    searches and reads all files in the **base_dir** matching *aminoacids.x*
    where *x* is *.rtp .arn .hdb .n.tdb .c.tdb .vsd .r2b*.
    Only the rtp file is mandatory, all others are neglected if not present.

    :param basename:    Basename of residue database to be loaded

    :type basename:     :class:`str`

  .. method:: ReadITP(basename)

    Searches and reads the itp file in the **base_dir**. *basename=amazing_ion*
    would therefore load the file *amazing_ion.itp*

    :param basename:    Basename of itp file to be loaded

    :type basename:     :class:`str`

  .. method:: SetForcefield(forcefield)

    Resets reader internal forcefield. Everything read so far is lost,
    except the already read preprocessor statements.

    :param forcefield:  Forcefield to be set

    :type forcefield:   :class:`Forcefield`

  .. method:: GetForcefield()

    Get the forcefield with everything read so far.

    :returns: The reader internal :class:`Forcefield` 


  .. code-block:: python
    
    path = "path_to_gromacs/share/top/charmm27.ff"
    reader = FFReader(path)

    #read in the data given in forcefield.itp and atomtypes.atp
    reader.ReadGromacsForcefield()

    #we also want to read several residue databases
    reader.ReadResidueDatabase("aminoacids")
    reader.ReadResidueDatabase("rna")
    reader.ReadResidueDatabase("dna")

    #ions and water are also nice to have, they're stored in itp files
    reader.ReadITP("tip3p")
    reader.ReadITP("ions")

    #let's finally get the reader internal forcefield out
    ff = reader.GetForcefield()

    #there is also an amazing ion definition in some other directory
    new_reader = FFReader("path/to/directory/with/itp/files")

    #we want to modify the previously read forcefield
    new_reader.SetForcefield(ff)

    #and read the amazing ion definition from an itp file
    #note, that any previously defined preprocessor statements
    #from the previous reader are lost 
    new_reader.ReadITP("amazing_ion")

    #the new forcefield finally contains everything we need, lets
    #extract it and save it down
    ff = new_reader.GetForcefield()
    ff.Save("charmm_forcefield.dat")

Generating forcefields with Antechamber
--------------------------------------------------------------------------------

The antechamber submodule of mol.mm defines functions to use Antechamber (from
AmberTools) to automatically generate force field parameters and load the
results into :class:`~ost.mol.mm.Forcefield` objects.

**Example usage**:

.. code-block:: python

  from ost.mol import mm

  # create parameters for RVP using PDB's component dictionary
  mm.antechamber.RunAntechamber('RVP', 'components.cif', base_out_dir='ligands')

  # create force field
  ff = mm.Forcefield()
  ff = mm.antechamber.AddFromPath(ff, 'ligands/RVP')
  # equivalent: ff = mm.antechamber.AddFromFiles(ff, 'ligands/RVP/frcmod',
  #                                              'ligands/RVP/out.mpdb')
  # since Antechamber cannot deal with ions, you can do it manually
  ff = mm.antechamber.AddIon(ff, 'CL', 'CL', 35.45, -1.0, 0.4401, 0.4184)
  # save it
  ff.Save('ligands/ff.dat')

**Functions**:

.. automodule:: ost.mol.mm.antechamber
  :members:

The Forcefield Class
--------------------------------------------------------------------------------

.. currentmodule:: ost.mol.mm

.. class:: Forcefield


  .. method:: Save(filename)

    Dumps forcefield into a binary file on disk

    :param filename:    Filename of the saved forcefield

    :type filename:     :class:`str` 



  .. staticmethod:: Load(filename)

    reads in binary forcefield file

    :param filename:    Filename of the forcefield to be loaded

    :type filename:     :class:`str`

    :returns:           loaded :class:`Forcefield`

    :raises:            :class:`RuntimeError` when **filename** can't be found



  .. method:: AddBond(bond)

    :param bond:        Bond to be added

    :type bond:         :class:`Interaction`

    :raises:            :class:`RuntimeError` when given interaction has
                                              no bond specific FuncType



  .. method:: AddAngle(angle)

    :param angle:       Angle to be added

    :type angle:        :class:`Interaction`

    :raises:            :class:`RuntimeError` when given interaction has
                                              no angle specific FuncType


  .. method:: AddDihedral(dihedral)

    :param dihedral:    Dihedral to be added

    :type dihedral:     :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no dihedral specific FuncType


  .. method:: AddImproper(improper)

    :param improper:    Improper to be added

    :type improper:     :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no improper specific FuncType


  .. method:: AddCMap(cmap)

    :param cmap:        CMap to be added

    :type cmap:         :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no cmap specific FuncType


  .. method:: AddImplicitGenborn(gb)

    :param gb:          GB to be added

    :type gb:           :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no gb specific FuncType


  .. method:: AddLJ(lj)

    :param lj:          LJ to be added

    :type lj:           :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no lj specific FuncType


  .. method:: AddLJPair(lj_pair)

    :param lj_pair:     LJPair to be added

    :type lj_pair:      :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no lj_pair specific FuncType


  .. method:: AddConstraint(constraint)

    :param constraint:  Constraint to be added

    :type constraint:   :class:`Interaction`     

    :raises:            :class:`RuntimeError` when given interaction has
                                              no constraint specific FuncType


  .. method:: AddMass(type, mass)

    :param type:        Type of atom
    :param mass:        Its mass

    :type type:         :class:`str`
    :type mass:         :class:`float`

  .. method:: SetFudgeLJ(factor)

    :param factor:      Factor with which the 1,4 Lennard Jones term
                        should be dampened

    :type factor:       :class:`float`


  .. method:: SetFudgeQQ(factor)

    :param factor:      Factor with which the 1,4 electrostatic term
                        should be dampened

    :type factor:       :class:`float`


  .. method:: SetGenPairs(gen_pairs)

    :param gen_pairs:   If set to false, all 1,4 interactions must be set
                        with AddLJPair. The Lorentz-Berthelot rule gets
                        used otherwise. 

    :type gen_pairs:    :class:`bool`


  .. method:: AddResidueRenamingRule(name, ff_main_name, ff_n_ter_name, ff_c_ter_name, ff_two_ter_name)

    :param name:        Original name of the residue 
                        (e.g. PDB/Gromacs standard)
    :param ff_main_name: Forcefield specific residue name
    :param ff_n_ter_name: Forcefield specific name if the residue
                          is N-Terminal
    :param ff_c_ter_name: Forcefield specific name if the residue
                                       is C-Terminal
    :param ff_two_ter_name: Forcefield specific name if the residue
                            is N- and C-Terminal

    :type name:            :class:`str`
    :type ff_main_name:    :class:`str`
    :type ff_n_ter_name:   :class:`str`
    :type ff_c_ter_name:   :class:`str`
    :type ff_two_ter_name: :class:`str`



  .. method:: AddAtomRenamingRule(res_name, old_atom_name, new_atom_name)

    :param res_name:    Forcefield specific name of the residue the
                                     atom belongs to

    :param old_atom_name: Atom name in PDB/Gromacs standard

    :param new_atom_name: FF specific atom name

    :type res_name:      :class:`str`
    :type old_atom_name: :class:`str`
    :type new_atom_name: :class:`str`


  .. method:: AddBuildingBlock(name, block)

    :param name:        Name of residue this :class:`BuildingBlock` 
                        is supposed to be related to

    :param block:       BuildingBlock to be added

    :type block:        :class:`BuildingBlock`
    :type name:         :class:`str`


  .. method:: AddHydrogenConstructor(name, h_constructor)

    :param name:        Name of residue this 
                        :class:`HydrogenConstructor` 
                        is supposed to be related to

    :param h_constructor: HydrogenConstructor to be added

    :type name:          :class:`str`
    :type h_constructor: :class:`HydrogenConstructor`


  .. method:: AddBlockModifier(name, modifier)

    :param name:        Name of residue this 
                        :class:`BlockModifier` 
                        is supposed to be related to

    :param modifier:    BlockModifier to be added

    :type name:         :class:`str`
    :type modifier:     :class:`BlockModifier`


  .. method:: SetStandardCTer(res_name, ter_name)

    Setting a standard CTer influences the behaviour of the GetCTerModifier 
    function. If no specific block modifier is defined there, this is the
    one that gets returned.

    :param res_name:    Forcefield specific residue name this block 
                        modifier is supposed to be related to

    :param ter_name:    Name of the default c-terminal block 
                        modifier for this residue

    :type res_name:     :class:`str`
    :type ter_name:     :class:`str`


  .. method:: SetStandardNTer(res_name, ter_name)

    Setting a standard NTer incluences the behaviour of the GetNTerModifier 
    function. If no specific block modifier is defined there, this is the
    one that gets returned.

    :param res_name:    Forcefield specific residue name this block 
                        modifier is supposed to be related to

    :param ter_name:    Name of the default n-terminal block 
                        modifier for this residue

    :type res_name:     :class:`str`
    :type ter_name:     :class:`str`
    


  .. method:: AssignFFSpecificNames(ent,[,reverse = False])

    This function does the forcefield specific renaming magic. It takes
    the given :class:`EntityHandle` and applies the rules set in
    AddResidueRenamingRule and AddAtomRenamingRule.

    :param ent:         Entity to be renamed

    :param reverse:     If False, the function does the renaming
                        from PDB/Gromacs naming to the forcefield
                        specific naming.
                        If True, the opposite happens.

    :type ent:          :class:`EntityHandle`
    :type reverse:      :class:`bool`


  .. method:: GetBond(type1, type2)

    :param type1:       Type of interacting particle 1
    :param type2:       Type of interacting particle 2

    :type type1:        :class:`str`
    :type type2:        :class:`str`

    :returns: an :class:`Interaction` with a bond FuncType

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found


  .. method:: GetAngle(type1, type2, type3)

    :param type1:       Type of interacting particle 1
    :param type2:       Type of interacting particle 2
    :param type3:       Type of interacting particle 3

    :type type1:        :class:`str`
    :type type2:        :class:`str`
    :type type3:        :class:`str`

    :returns: an :class:`Interaction` with a angle FuncType

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found


  .. method:: GetDihedrals(type1, type2, type3, type4)

    Several dihedral definitions can be merged to one dihedral function.
    This function therefore returns a list. 
    In a first step all dihedrals matching the given types are gathered
    and returned.
    If no dihedrals can be found, the search continues by including
    wildcard characters in the atom types (X). All found dihedrals
    matching with all possible combinations of wildcards are then gathered
    and returned.

    :param type1:       Type of interacting particle 1
    :param type2:       Type of interacting particle 2
    :param type3:       Type of interacting particle 3
    :param type4:       Type of interacting particle 4

    :type type1:        :class:`str`
    :type type2:        :class:`str`
    :type type3:        :class:`str`
    :type type4:        :class:`str`

    :returns: a :class:`list` of :class:`Interaction` objects with dihedral 
              FuncType matching given types

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found



  .. method:: GetImpropers(type1, type2, type3, type4)

    The same search strategy as in GetDihedrals is used to extract 
    the impropers.

    :param type1:       Type of interacting particle 1
    :param type2:       Type of interacting particle 2
    :param type3:       Type of interacting particle 3
    :param type4:       Type of interacting particle 4

    :type type1:        :class:`str`
    :type type2:        :class:`str`
    :type type3:        :class:`str`
    :type type4:        :class:`str`

    :returns: a :class:`list` of :class:`Interaction` objects with improper
              FuncType matching given types

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found


  .. method:: GetCMap(type1, type2, type3, type4, type5)

    :param type1:       Type of interacting particle 1
    :param type2:       Type of interacting particle 2
    :param type3:       Type of interacting particle 3
    :param type4:       Type of interacting particle 4
    :param type5:       Type of interacting particle 5

    :type type1:        :class:`str`
    :type type2:        :class:`str`
    :type type3:        :class:`str`
    :type type4:        :class:`str`
    :type type5:        :class:`str`

    :returns: an :class:`Interaction` with a cmap FuncType

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found


  .. method:: GetImplicitGenborn(type)

    :param type:        Type of particle

    :type type:         :class:`str`

    :returns: an :class:`Interaction` with a gb FuncType

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given type can be found


  .. method:: GetLJ(type)

    :param type:        Type of particle

    :type type:         :class:`str`

    :returns: an :class:`Interaction` with a lj FuncType

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given type can be found


  .. method:: GetLJ(type1, type2,[,pair=False])
    :noindex:

    :param type1:        Type of interacting particle 1
    :param type2:        Type of interacting particle 2
    :param pair:         If set to true, the interaction is
                         assumed to be a 1,4-interaction and
                         the set lj_pairs are first searched
                         for matches. In case of no success,
                         the function uses the Lorentz-Berthelot
                         rule to combine the sigma and epsilon 
                         parameters.
                         If set to false, the Lorentz-Berthelot
                         rule is applied directly.

    :type type1:        :class:`str`
    :type type2:        :class:`str`
    :type pair:         :class:`bool`


    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found
                                              or when pair is true and no 
                                              appropriate lj_pair is set 
                                              despite gen_pair flag being false.


  .. method:: GetConstraint(type1, type2)

    :param type1:       Type of interacting particle 1
    :param type2:       Type of interacting particle 2

    :type type1:        :class:`str`
    :type type2:        :class:`str`

    :returns: an :class:`Interaction` with a constraint FuncType

    :raises:            :class:`RuntimeError` when no :class:`Interaction`
                                              matching given types can be found


  .. method:: GetMass(type)

    :param type:        Type of particle

    :type type:         :class:`str`

    :returns:           :class:`float` - the mass

    :raises:            :class:`RuntimeError` if no mass has been set for this 
                        atom type


  .. method:: GetFudgeLJ()

    :returns:  :class:`float` - Factor with which the 1,4 Lennard Jones 
                term should be dampened

  .. method:: GetFudgeQQ()

    :returns:  :class:`float` - Factor with which the 1,4 
                electrostatic term should be dampened


  .. method:: GetAtomType(res_name, atom_name)

    :param res_name:    Forcefield specific residue name

    :param atom_name:   Forcefield specific atom name belonging
                        to that residue

    :type res_name:     :class:`str`
    :type atom_name:    :class:`str`

    :returns:           :class:`str` - atom type

    :raises:            :class:`RuntimeError` if forcefield has no such
                        :class:`BuildingBlock` or when atom is not present 
                        in that :class:`BuildingBlock`   


  .. method:: GetHydrogenConstructor(res_name)

    :param res_name:    Name of residue
    :type res_name:     :class:`str`

    :returns: :class:`HydrogenConstructor` for this name, invalid if it can't
              be found


  .. method:: GetBuildingBlock(res_name)

    :param res_name:    Name of residue
    :type res_name:     :class:`str`

    :returns:  :class:`BuildingBlock` for this name, invalid if it can't be 
               found

  .. method:: GetBuildingBlockNames()

    :returns:  :class:`list` of all building block names present in that 
               forcefield


  .. method:: GetBlockModifier(res_name)

    :param res_name:    Name of residue
    :type res_name:     :class:`str`

    :returns: :class:`BlockModifier` for this name, invalid if it can't
              be found


  .. method:: GetNTerModifier(res_name,[,ter_name=""])

    :param res_name:    Name of residue

    :param ter_name:    If not set, the ter_name
                        defined by SetStandardNTer gets used

    :type res_name:     :class:`str`
    :type ter_name:     :class:`str`


    :returns: :class:`BlockModifier` for this name, invalid if it can't
              be found


  .. method:: GetCTerModifier(name,[,ter_name=""])

    :param res_name:    Name of residue

    :param ter_name:    If not set, the ter_name
                        defined by SetStandardCTer gets used

    :type res_name:     :class:`str`
    :type ter_name:     :class:`str`

    :returns: :class:`BlockModifier` for this name, invalid if it can't
              be found




    
