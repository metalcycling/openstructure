Storing Custom Data
================================================================================

Introduction
--------------------------------------------------------------------------------

.. currentmodule:: ost

It is often very convenient to store any arbitrary data inside an Entity. A few 
examples are: 

  * calculated properties of atoms
  * sequence conservation of a residue
  * interaction energy of a substructure with its surrounding
  * fit of a fragment inside an electron density map

In OpenStructure this is supported by the use of generic properties. Most 
building blocks are derived from :class:`~ost.GenericPropContainer`, meaning
that arbitrary key-value pairs can be stored in them. In essence, the following 
classes support generic properties:

  * :class:`~ost.mol.EntityHandle` and :class:`~ost.mol.EntityView`
  * :class:`~ost.mol.ChainHandle` and :class:`~ost.mol.ChainView`
  * :class:`~ost.mol.ResidueHandle` and :class:`~ost.mol.ResidueView`
  * :class:`~ost.mol.AtomHandle` and :class:`~ost.mol.AtomView`
  * :class:`~ost.mol.BondHandle`
  * :class:`~ost.seq.SequenceHandle`

The view variants will reflect the generic properties of the handle variants.

A generic property key is always a string, and a value can be one of string, 
float, int or bool. For each of these data types, methods to retrieve and store 
values are available both in Python and C++.

Storing and Accessing Data
--------------------------------------------------------------------------------

All OpenStructure building blocks that are :class:`~ost.GenericPropContainer`\s,
have four different methods to store generic data, depending on the data type 
(i.e. string, float, int or bool).

To store a float value with the key 'myfloatprop' in all atoms of an entity:

.. code-block:: python
  
  import math
  for atom in entity.GetAtomList(): 
    val=5*math.sin(0.4*atom.GetPos().x)
    atom.SetFloatProp("myfloatprop", val)
  
If a GenericProp at a given level (i.e. atom, bond, residue, chain or entity) 
already exists, it will be overwritten. To check if it exists, use:

.. code-block:: python
  
  exists=atom.HasProp("myfloatprop")
  print(exists)
    
To access the value of a generic property, we first check if the property exists
and then access it, using the method suitable for the data type of the property. 
For the previously set property `myfloatprop` of the data type real, at the atom 
level:

.. code-block:: python
  
  for atom in entity.GetAtomList(): 
    if atom.HasProp("myfloatprop"):
      print(atom.GetFloatProp("myfloatprop"))
        
When trying to access a property that has not been set, or one that has been 
set, but at a different level, an error is thrown. The same is true when trying 
to access a property of a different data type, e.g.:

.. code-block:: python

  # all of the following lines will throw errors
  # error because the property does not exist 
  print(atom.GetFloatProp("unknownprop"))
  
  # error because the property was set at another level
  print(entity.GetFloatProp("myfloatprop"))
  
  # error because the data type of the property is different
  print(atom.GetStringProp("myfloatprop"))
      

Use of Generic Properties in Queries
--------------------------------------------------------------------------------

The :doc:`../mol/base/query` can also be used for numeric generic properties 
(i.e. bool, int, float), but the syntax is slightly different. To access any 
generic properties, it needs to be specified that they are generic and at which 
level (chain, residue, atom) they are defined. Therefore, all generic properties 
start with a 'g', followed by an 'a', 'r' or 'c' for atom, residue or chain 
level respectively. For more details see :doc:`../mol/base/query`. 


API documentation
--------------------------------------------------------------------------------

.. class:: GenericPropContainer

  .. method:: HasProp(key)
  
    checks existence of property. Returns true, if the the class contains a
    property with the given name, false if not.
  
  .. method:: GetPropAsString(key)
  
    Returns the string representation of a property, or the empty String if 
    the property addressed by key does not exist. Note that this is not the 
    same as trying to get a generic float/int/bool property as a string type; 
    the latter will result in a boost:get exception. Use this method to obtain 
    a representation suitable for output.
    
  .. method:: GetStringProp(key)
              GetStringProp(key, default_value)
  
    Get string property. The first signature raises a GenericPropError error if
    the property does not exist, the second returns the default value.
  
  
  .. method:: GetFloatProp(key)
              GetFloatProp(key, default_value)
   
    Get float property. The first signature raises a GenericPropError error if
    the property does not exist, the second returns the default value.


  .. method:: GetIntProp(key)
              GetIntProp(key, default_value)

    Get int property. The first signature raises a GenericPropError error if
    the property does not exist, the second returns the default value.

  .. method:: GetBoolProp(key)
              GetBoolProp(key, default_value)

    Get bool property. The first signature raises a GenericPropError error if
    the property does not exist, the second returns the default value.
        
  .. method:: ClearProps()
    
    Remove all generic properties
  
 
  .. method:: SetStringProp(key, value)
    
    Set string property, overriding an existing property with the same name
    
  .. method:: SetFloatProp(key, value)
    
    Set float property, overriding an existing property with the same name

  .. method:: SetIntProp(key, value)
  
    Set int property, overriding an existing property with the same name
  
  .. method:: SetBoolProp(key, value)
  
    Set bool property, overriding a property with the same name
    
  .. method:: RemoveProp(key)
  
    Removes the property with given key, regardless of its type. If the property 
    does not exist, the method has no effect.
    

