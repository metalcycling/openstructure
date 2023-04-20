Queries
================================================================================

.. currentmodule:: ost.mol

OpenStructure includes a powerful query system that allows you to perform custom 
selections on a molecular entity in a convenient way.


The Basics
--------------------------------------------------------------------------------

It is often convenient to highlight or focus certain parts of the structure. 
OpenStructure includes a powerful query system that allows you to perform custom 
selections in a convenient way. Selections are carried out mainly by calling the Select method made available by EntityHandle and :class:`EntityView` objects while providing a query string. Queries are written using a dedicated mini-language. For example, to select all arginine residues of a given structure, one would write:

.. code-block:: python

  arginines = model.Select('rname=ARG')

A simple selection query (called a predicate) consists of a property (here,
`rname`), a comparison operator (here, `=`) and an argument (here, `ARG`).  The 
return value of a call to the :meth:`EntityHandle.Select` method is always an 
:class:`EntityView`. The :class:`EntityView` always contains a full hierarchy of 
elements, never standalone separated elements. In the above example, the 
:class:`EntityView` called `arginines` will contain all chains from the 
structure called `model` that have at least one arginine. In turn these chains 
will contain all residues that have been identified as arginines. The residues 
themselves will contain references to all of their atoms. Of course, queries are 
not limited to selecting residues based on their type, it is also possible to 
select atom by name:

.. code-block:: python

  c_betas = model.Select('aname=CB') 

As before, `c_betas` is an instance of an :class:`EntityView` object and contains
a full hierarchy. The main difference to the previous example is that the 
selected residues do not contain a list of all of their atoms but only the 
C-beta. These examples clarify why the name 'view' was chosen for this result of 
a :meth:`~EntityHandle.Select` statement. It represents a reduced, restrained 
way of looking at the original structure.

Both the selection statements that have been used so far take strings as their arguments. However, selection properties such as `rnum` (residue number), take numeric arguments. With numeric arguments it is possible to use identity operators ( `!=` and `=`). It is also possible to compare them using the `>`, `<`, `>=` and `<=` operators. For example, the 20 N-terminal residues of a protein can be selected with:

.. code-block:: python

  n_term = model.Select('rnum<=20')

If you want to supply arguments with special characters they need to be put in
quotation marks (' or "). For instance, this is needed for any chain name
containing spaces as in:

.. code-block:: python

  model.Select('cname=" "')

Almost any name can be quoted with :meth:`QueryQuoteName`.


Combining predicates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Selection predicates can be combined with boolean operators. For example , you might want to select all C atoms with crystallographic occupancy higher than 50. These atoms must match the predicate `ele=C` in addition to the predicate `occ>50`. In the query language this can be written as:

.. code-block:: python

  model.Select('ele=C and occ>50')
  
Compact forms are available for several selection statements. For example, to select all arginines and aspargines, one could use a statement like:

.. code-block:: python
  
  arg_and_asn = model.Select('rname=ARG or rname=ASN')

However, this is rather cumbersome as it requires the word `rname` to be typed  twice. Since the only difference between the two parts of the selection is the argument that follows  the word `rname`, the statement can also be written in an abbreviated form:

.. code-block:: python
  
  arg_and_asn = model.Select('rname=ARG,ASN')

Another example: to select residues with numbers in the range 130 to 200, one could use the following statement

.. code-block:: python
  
  center = model.Select('rnum>=130 and rnum<=200')

or alternatively use the much nicer syntax:

.. code-block:: python

  center = model.Select('rnum=130:200')

This last statement is completely equivalent to the previous one. This syntax  
can be used when the selection statement requires a range of integer values 
within a closed interval.

Distance Queries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The query

.. code-block:: python

  around_center = model.Select('5 <> {0,0,0}')
  
selects all chains, residues and atoms that lie with 5 Å to the origin of the
reference system ({0,0,0}). The `<>` operator is called the 'within' operator.
Instead of a point, the within statements can also be used to return a view
containing all chains, residues and atoms within a radius of another selection
statement applied to the same entity. Square brackets are used to delimit the
inner query statement.

.. code-block:: python

  around_hem = model.Select('5 <> [rname=HEM]')
  model.Select('5 <> [rname=HEM and ele=C] and rname!=HEM')

Bonds and Queries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When an :class:`EntityView` is generated by a selection, it includes by default
only bonds for which both connected atoms  satisfy the query statement. This can
be changed by passing the parameters `EXCLUSIVE_BONDS` or `NO_BONDS` when
calling the Select method. `EXCLUSIVE_BONDS` adds bonds to the
:class:`EntityView` when at least one of the two atoms falls within the boundary
of the selection. `NO_BONDS` suppresses the bond inclusion step completely.

Whole Residue Queries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the parameter `MATCH_RESIDUES` is passed when the Select method is called,
the resulting :class:`EntityView` will include whole residues for which at least
one atom satisfies the query. This means that if at least one atom in the
residue  falls within the boundaries of the selection, all atoms of the residue
will be included in the View.

More Query Usage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The high level interface for queries are the Select methods of the 
EntityHandle and :class:`EntityView` classes. By passing in a query string, a view
consisting of a subset of the elements is returned. 

Queries also offer a second interface: `IsAtomSelected()`, 
`IsResidueSelected()` and `IsChainSelected()` take an atom, residue or 
chain as their argument and return true or false, depending on whether the 
element fulfills the predicates.

.. _genprop-in-queries:

Generic Properties in Queries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The query language can also be used for numeric generic properties (i.e.
float and int), but the syntax is slightly different. To access any generic
properties, it needs to be specified that they are generic and at which level 
they are defined. Therefore, all generic properties start with a `g`, followed by an `a`, `r` or `c` for atom, residue or chain level respectively. 

.. code-block:: python
   
   # set generic properties for atom, residue, chain
   atom_handle.SetFloatProp("testpropatom", 5.2)
   resid_handle.SetFloatProp("testpropres", 1.1)
   chain_handle.SetIntProp("testpropchain", 10)
   
   # query statements
   sel_a = e.Select("gatestpropatom<=10.0")
   sel_r = e.Select("grtestpropres=1.0")
   sel_c = e.Select("gctestpropchain>5")

Since generic properties do not need to be defined for all parts of an entity
(e.g. it could be specified for one single :class:`AtomHandle`), the query
statement will throw an error unless you specify a default value in the query
statement which can be done using a ':' character:

.. code-block:: python

  # if one or more atoms have no generic properties
  
  sel = e.Select("gatestprop=5")
  # this will throw an error
  
  # you can specify a default value:
  sel = e.Select("gatestprop:1.0=5")
  # this will run through smoothly and use 1.0 as
  # the default value for all atoms that do not
  # have the generic property 'testprop'

Using this method, you will be warned if a generic property is not set for all
atoms, residues or chains unless you specify a default value. So, be careful
when you do.

Available Properties
--------------------------------------------------------------------------------

The following properties may be used in predicates. The type is given for each property.

Properties of Chains
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**cname/chain** (str) :attr:`Chain name<ChainHandle.name>` 
  
Properties of Residues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**rname** (str): :attr:`Residue name<ResidueHandle.name>`
 
**rnum** (int): :attr:`Residue number<ResidueHandle.number>`. Currently only the
numeric part is honored.
 
**rtype** (str): Residue type as given by the DSSP code (e.g. "H" for alpha
helix, "E" for extended), "helix" for all helix types, "ext" or "strand" for
all beta sheets or "coil" for any type of coil (see :class:`SecStructure`).
  
**rindex** (int): :attr:`Index<ResidueHandle.index>` of residue handle in chain.
This index is the same for views and handles.
  
**peptide** (bool): Whether the residue is :attr:`peptide linking
<ResidueHandle.peptide_linking>`.

**nucleotide** (bool): Whether the residue is :attr:`nucleotide linking
<ResidueHandle.nucleotide_linking>`.

**protein** (bool): Whether the residue is considered to be
:attr:`part of a connected protein <ResidueHandle.is_protein>`.

**rbfac** (float): average B (temperature) factor of residue

**ligand** (bool) Whether the residue is a :meth:`ligand
<ResidueHandle.is_ligand>`.

**water** (bool) Whether the residue is water.

Properties of Atoms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**aname** (str): :attr:`Atom name<AtomHandle.name>`
  
**ele** (str): :attr:`Atom element<AtomHandle.element>`

**occ** (float): :attr:`Atom occupancy<AtomHandle.occupancy>`

**abfac** (float): :attr:`Atom B-factor<AtomHandle.b_factor>`

**x** (float): :attr:`X<AtomHandle.pos>` coordinate of atom.

**y** (float): :attr:`Y<AtomHandle.pos>` coordinate of atom.

**z** (float): :attr:`Z<AtomHandle.pos>` coordinate of atom.

**aindex** (int): :attr:`Atom index<AtomHandle.index>`

**ishetatm** (bool): Whether the atom is a :attr:`heterogenous
atom<AtomHandle.is_hetatom>`.

**acharge** (float): :attr:`Atom charge<AtomHandle.charge>`


Query API documentation
--------------------------------------------------------------------------------

In the following, the interface of the query class is described. In general, you will not have to use this interface but will pass the query as string directly. 

.. class:: Query(string='')
  
  Create a new query from the given string. The constructor does not throw any
  error in case the query contains syntax errors. Use :attr:`valid` to check
  whether the query was valid.


  .. attribute:: string
  
    The string used to create the query.
  
    :type: str
    
  .. attribute:: valid
  
    True, when the query could be compiled without syntax errors.
    
    :type: bool

  .. attribute:: error
    
    If :attr:`valid` is false, this attribute contains the error message.
    Otherwise it is set to an empty string
    
    :type: str
  .. method:: IsAtomSelected(atom)
  
    Returns true, when the given atom handle fulfills the predicates, false if
    not.

  .. method::  IsChainSelected(chain)
  
    Return true if at least one of the atomso of the chain matches the
    predicates.
    
  .. method:: IsResidueSelected(residue)
  
    Returns true, when at least one atom of the residue matches the predicates.

.. class:: QueryFlag

  Defines flags to change default behaviour of Select queries. Possible values:

  * ``EXCLUSIVE_BONDS`` - adds bonds to the :class:`EntityView` when at least
    one of the two bonded atoms was selected (by default both must be selected)
  * ``NO_BONDS`` - do not include any bonds (by default bonds are included)
  * ``MATCH_RESIDUES`` - include all atoms of a residue if any of its atoms is
    selected (by default only selected atoms are included)

.. method:: QueryQuoteName(name)

  Adds appropriate quotation marks to use *name* in a :class:`Query`. For
  instance the following code snippet would generate a query string selecting
  all chains from a list of chain names:

  .. code-block:: python

    query = "cname=" + ','.join(mol.QueryQuoteName(name) for name in names)

  Note that there is some limited support of wild card symbols (* and ?) which
  may have undesired effects in a query such as the code above.

  :param name: Name to put in quotation marks
  :type name:  :class:`str`
  :rtype: :class:`str`
  :raises: :exc:`~exceptions.Exception` if *name* cannot be used in queries.
           This happens if *name* includes both ' and " or if it ends with \\.
