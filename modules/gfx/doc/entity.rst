The Graphical Entity
================================================================================

.. currentmodule:: ost.gfx
  
.. note::

  This document describes the API of the :class:`Entity` and related classes. 
  For an overview of the graphical entity, see :doc:`../intro-03`.

.. _render-mode:

Render Modes
--------------------------------------------------------------------------------

The entity supports different render modes. They are depicted in the following 
images together with the options they support. To change the render mode of the 
full entity or part of the entity, use the :meth:`Entity.SetRenderMode` method.

Render Modes at the Connectivity Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. data:: SIMPLE

  .. image:: simple.png 
  
  Renders bonds as lines and single atoms as crosses.

.. data:: CPK
  
  .. image:: cpk.png

  
  Renders atoms as spheres with radii proportional to their van der Waals 
  radius.

.. data:: CUSTOM
  
  .. image:: custom.png

Render Modes at the Trace Level
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. data:: TUBE

  .. image:: tube.png
  
  Renders a smooth backbone trace with circular profiles.

.. data:: SLINE

  .. image:: sline.png
  
  Renders a smooth backbone trace with lines.
  
.. data:: TRACE

  .. image:: trace.png
  
  Renders the Calpha atoms connected by tubes.
  
.. data:: LINE_TRACE

  .. image:: line_trace.png

  Renders the Calpha connected by lines.
  
.. data:: HSC

  .. image:: hsc.png

  Renders a helix/strand/coil cartoon.


.. _render-options:

Render Options
--------------------------------------------------------------------------------


**line_width**

  The line width. Defaults to 1.0
  
  :type: float

**aa_lines**

  Controls whether antialiased lines are enabled. Only available when
  compiled with shaders. Defaults to false.
  
  :type: bool
  
**draw_bond_order**

  Controls whether double/triple bonds are drawn as multiple lines. This 
  requires that the bond order information has been assigned. As of now, 
  only the  :class:`~ost.conop.RuleBasedBuilder` assigns this information.
  
  Defaults to off.
  
  :type: bool
  



**bond_order_distance**

  The distance between the lines depicting higher order bonds.


Loading and Displaying an Entity  
--------------------------------------------------------------------------------

This code examples shows how to load a structure and display it on the screen. It also shows how to programmatically combine different render modes in one graphical entity.

.. literalinclude:: /../../examples/demos/load_and_display.py


Entity API Reference
--------------------------------------------------------------------------------

.. class:: Entity(name, entity)
           Entity(name, render_mode, entity)

  Graphical representation of entities. 
  
  :param name: Name of the entity.
  :type  name: str
  :param render_mode: One of the :ref:`render-mode`.
     Defaults to :obj:`SIMPLE`.
  :param entity: The entity. Only atoms and bonds that part part of the entity
     view (or handle) are rendered.
  :type  entity: :class:`~ost.mol.EntityView`, or :class:`~ost.mol.EntityHandle`

  .. attribute:: name
  
    The name of the entity. Read-only. Also available as :meth:`GetName`
    
    :type: str
    
  .. attribute:: selection
  
    The selection of the entity is rendered with a green halo. As a shortcut, 
    the property an also be set as a string. The following two statements are 
    identical:
    
    .. code-block:: python
    
      go.selection = query
      go.selection = go.view.Select(query)
    
    :type: :class:`~ost.mol.EntityView`
  
  .. attribute:: center
  
    The geometric center of the graphical entity. Read-only.
    
    :type: :class:`~ost.geom.Vec3`
  
  
  .. method:: SetRenderMode(render_mode)
              SetRenderMode(render_mode, view, keep=False)
    
    Changes the render mode of the graphical entity. The first signature changes 
    the render mode of all atoms and bonds, whereas the second only changes the 
    display of atoms/bonds that are part of the view. If `keep` is set to false,
    the atoms and bonds that are part of the view will only be displayed in the 
    new render mode. If keep is set to true, the atoms and bonds will be 
    rendered with the new render mode in addition to whatever render mode they 
    were rendered previously.
    
    :param render_mode: One of the :ref:`render-mode`.
       Defaults to :obj:`SIMPLE`.
    :param view: A valid entity view.
    :type  view: :class:`~ost.mol.EntityView`

  .. method:: SetColor(color, selection='')
  
    Sets the atoms matching the selection predicates to `color`.
    
    :param color: The new color
    :type  color: :class:`Color`
    :param selection: A valid :doc:`query <../mol/base/query>`
    :type  selection: str
    
  .. method:: ColorByElement()
  
    Apply `CPK coloring <http://en.wikipedia.org/wiki/CPK_coloring>`_ coloring
    to atoms. This will color the carbon atoms in white, nitrogen in blue, 
    oxygen in red, phosphorous in pink and sulfur in yellow.
    
  
  .. method:: ColorBy(prop, color1, color2, min, max, hint='U')
              ColorBy(prop, color1, color2, hint='U')
              ColorBy(prop, gradient, min, max, hint='U')              
              ColorBy(prop, gradient, hint='U')
    
    Colors the entity by mapping a numeric property to a color. This property 
    can be a built-in a property such as atomic b-factor, charge, residue number 
    or be a :doc:`custom property <../base/generic>` assigned by the user or an 
    algorithm and be defined at any level (chain, residue, atom). The naming of 
    the built-in properties is identical to the properties available in 
    :doc:`query language <../mol/base/query>`.
    
    In case a generic property is undefined for a certain element, the property 
    value of min is used for that element.
    
    :param prop: The name of the numeric property
    :type  prop: str
    :param color1: Color corresponding to min value
    :type color1: :class:`Color`
    :param color2: Color corresponding to max value
    :type color2: :class:`Color`    
    :param gradient: The gradient to be used
    :type gradient: :class:`Gradient`
    :type min: Value mapped to `color1`, or the color at 0.0 in the gradient
    :type max: Value mapped to `color2` or the color at 1.0 in the gradient
    :param hint: When `prop` refers to a generic property, hint must be set to
       the hierarchy level at which the property is defined. 
    :type hint: str
    
  .. method:: UpdatePositions()
  
    Tells the graphical entity to update its vertex arrays. Use this function 
    when you modified coordinates of the underlying 
    :class:`~ost.mol.EntityHandle` and would like to see the changes on the 
    screen.
    
    :see: :meth:`UpdateViews`
    
  .. method:: UpdateView()
  
    Tells the graphical object to recalculate all views and vertex arrays. Use
    this function when you added new chains, residues, atoms or bonds to the 
    entity and would like to see the changes on the screen. 
    
    :see: :meth:`UpdatePositions`
