.. currentmodule:: ost.conop

The compound library
================================================================================

Compound libraries contain information on chemical compounds, such as their 
connectivity, chemical class and one-letter-code. The compound library has 
several uses, but the most important one is to provide the connectivy 
information for the :class:`rule-based processor <RuleBasedBuilder>`. 

The compound definitions for standard PDB files are taken from the 
components.cif dictionary provided by the PDB. The dictionary is updated with 
every PDB release and augmented with the compound definitions of newly 
crystallized compounds. 

If you downloaded the bundle, a recent version of the compound library is 
already included. If you are compiling from source or want to incorporate the 
latest compound definitions, follow :ref:`these instructions <mmcif-convert>` to 
build the compound library manually.


.. function:: GetDefaultLib()

  :return: Default compound library set by :func:`SetDefaultLib`. If you got
           OpenStructure as a bundle or you :ref:`compiled <cmake-flags>`  it
           with a specified ``COMPOUND_LIB`` flag, this will return a compound
           library when executing scripts with ``ost``.
  :rtype:  :class:`CompoundLib` or None if no library set

.. function:: SetDefaultLib(lib)

  :param lib: Library to be set as default compound library.
  :type lib:  :class:`CompoundLib`


.. class:: CompoundLib

  .. staticmethod:: Load(database, readonly=True)
  
    Load the compound lib from database with the given name.
    
    :param readonly: Whether the library should be opened in read-only mode. It 
      is important to note that only one program at the time has write access to 
      compound library. If multiple programs try to open the compound library in 
      write mode, the programs can deadlock.
    :type readonly: :class:`bool`
    
    :returns: The loaded compound lib or None if it failed.
    
  .. staticmethod:: Create(database)
    
    Create a new compound library
    
  .. method:: FindCompound(tlc, dialect='PDB')
  
    Lookup compound by its three-letter-code, e.g ALA. If no compound with that 
    name exists, the function returns None. Compounds are cached after they have 
    been loaded with FindCompound. To delete the compound cache, use 
    :meth:`ClearCache`.
    
    :returns: The found compound
    :rtype: :class:`Compound`
  
  .. method:: Copy(dst_filename)
  
    Copy database to dst_filename. The new library will be an exact copy of the 
    database. The special name `:memory:` will create an in-memory version of 
    the database. At the expense of memory, database lookups will become much 
    faster.
    
    :returns: The copied compound library
    
    :rtype: :class:`CompoundLib`

  .. method:: ClearCache()
  
    Clear the compound cache.

  .. method:: SetChemLibInfo()

     When creating the new library the current date and the Version of OST used
     are stored into the table chemlib_info.

  .. method:: GetOSTVersionUsed()

     :return: OST version (ost_version_used from the table chemlib_info)
     :rtype:  :class:`str`

  .. method:: GetCreationDate()

     :return: creation date (creation_date from the table chemlib_info)
     :rtype:  :class:`str`


.. class:: Compound

  Holds the description of a chemical compound, such as three-letter-code, and
  chemical class.

  .. attribute:: id
  
    Alias for :attr:`three_letter_code`
    
  .. attribute:: three_letter_code
  
    Three-letter code of the residue, e.g. ALA for alanine. The three-letter 
    code is unique for each compound, always in uppercase letters and is between  
    1 and 3 characters long.
    
    code is always uppercase.
    
  .. attribute:: one_letter_code
  
    The one letter code of the residue, e.g. 'G' for glycine. If undefined, the 
    one letter code of the residue is set to '?'

  .. attribute:: formula
  
    The chemical composition, e.g. 'H2 O' for water. The elements are listed in 
    alphabetical order.
    
  .. attribute:: dialect
  
    The dialect of the compound.
    
  .. attribute:: atom_specs

    The atom definitions of this compound. Read-only.

    :type: list of :class:`AtomSpec`
          
  .. attribute:: bond_specs
  
    The bond definitions of this compound. Read-only.
    
    :type: list of :class:`BondSpec`
    
  .. attribute:: chem_class
  
    The :class:`~ost.mol.ChemClass` of this compound. Read-only.
    
    :type: :class:`str`
    
  .. attribute:: chem_type
  
    The :class:`~ost.mol.ChemType` of this compound. Read-only.
    
    :type: :class:`str`
    
  .. attribute:: inchi
  
    The InChI code of this compound, without the 'InChI=' part, e.g 
    '1S/H2O/h1H2' for water. Read-only.
    
    :type: :class:`str`
    
  .. attribute:: inchi
  
    The InChIKey of this compound without the 'InChIKey=' part, e.g.
    'XLYOFNOQVPJJNP-UHFFFAOYSA-N' for water. Read-only.
    
    :type: :class:`str`
    

.. class:: AtomSpec

  Definition of an atom
  
  .. attribute:: element
  
    The element of the atom
    
  .. attribute:: name
  
    The primary name of the atom
    
  .. attribute:: alt_name
  
    Alternative atom name. If the atom has only one name, this is identical to 
    :attr:`name`
    
  .. attribute:: is_leaving
  
    Whether this atom is required for a residue to be complete. The best example 
    of a leaving atom is the *OXT* atom of amino acids that gets lost when a 
    peptide bond is formed.

.. class:: BondSpec

  Definition of a bond
  
  .. attribute:: atom_one
    
    The first atom of the bond, encoded as index into the 
    :attr:`Compound.atom_specs` array.
    
  .. attribute:: atom_two
  
    The second atom of the bond, encoded as index into the 
    :attr:`Compound.atom_specs` array.
    
  .. attribute:: order
  
    The bond order, 1 for single bonds, 2 for double-bonds and 3 for 
    triple-bonds
    

Example: Translating SEQRES entries
--------------------------------------------------------------------------------

In this example we will translate the three-letter-codes given in the SEQRES record to one-letter-codes. Note that this automatically takes care of modified amino acids such as selenium-methionine.


.. code-block:: python

  compound_lib=conop.CompoundLib.Load('compounds.chemlib')
  seqres='ALA GLY MSE VAL PHE'
  sequence=''
  for tlc in seqres.split():
    compound=compound_lib.FindCompound(tlc)
    if compound:
       sequence+=compound.one_letter_code
  print(sequence) # prints 'AGMVF'

.. _mmcif-convert:

Creating a compound library
--------------------------------------------------------------------------------

The simplest way to create compound library is to use the :program:`chemdict_tool`. The programs allows you to import the chemical 
description of the compounds from a mmCIF dictionary, e.g. the components.cif dictionary provided by the PDB. The latest dictionary for can be downloaded from the `wwPDB site <http://www.wwpdb.org/ccd.html>`_. The files are rather large, it is therefore recommended to download the gzipped version.

After downloading the file use :program:`chemdict_tool` to convert the MMCIF  dictionary into our internal format.  

.. code-block:: bash
  
  chemdict_tool create <components.cif> <compounds.chemlib>

Note that the :program:`chemdict_tool` only understands `.cif` and `.cif.gz` 
files. If you have would like to use other sources for the compound definitions, consider writing a script by using the :doc:`compound library <compoundlib>` API.

If you are working with CHARMM trajectory files, you will also have to add the 
definitions for CHARMM. Assuming your are in the top-level source directory of 
OpenStructure, this can be achieved by:

.. code-block:: bash

  chemdict_tool update modules/conop/data/charmm.cif <compounds.chemlib> charmm


Once your library has been created, you need to tell cmake where to find it and 
make sure it gets staged.


.. code-block:: bash
  
  cmake -DCOMPOUND_LIB=compounds.chemlib
  make
