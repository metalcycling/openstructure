Introduction to the :mod:`~ost.gfx` Module
================================================================================

.. currentmodule:: ost.gfx

For the course of this tutorial, we assume that you have :ref:`DNG up and 
running<start-dng>`.

As already mentioned in the :doc:`intro-01`, the basic representation of 
molecular structures, sequences, images are completely independent from the 
graphical rendering. Thus, to display it on the screen, we first have to create 
a graphical representation of the object. This is illustrated for the entity 
below, but follows a similar scheme for density maps, surfaces and other types 
of data.

.. code-block:: python
  
  pdb = io.LoadPDB('/path/to/molecule.pdb')
  go = gfx.Entity('PROTEIN', pdb)
  scene.Add(go)
  scene.CenterOn(go)

The last command is not mandatory, but very convenient to center the camera on the object. 


The Scene
--------------------------------------------------------------------------------


The scene is the central registry for graphical objects and manages rendering 
parameters. Among other parameters, it is used to setup the lighting, fog, 
background color and the camera. The scene is a singleton, meaning that there is 
only one scene available. The instance can be accessed via :func:`gfx.Scene`. 
Because the scene is so important and commonly used, the scene is also available 
as the `scene` variable in the interactive python shell as well as from scripts. 
In fact, this is the preferred way to use functionality of the scene. As an easy example, let's change the background color:

.. code-block:: python

  scene.SetBackground(gfx.WHITE)

Now, the color of the screen has changed to white. To set it back to black, simply assign :obj:`gfx.BLACK` to the background property.

Objects are added to the scene by calling :meth:`Scene.Add`. Once added to the 
scene, the objects will be rendered automatically. Since all the graphical 
objects have a name that identifies them, they can be retrieved by name and 
manipulated:

.. code-block:: python

  # retrieving the previously added protein
  obj = scene['PROTEIN']
  # set color of the protein to red
  obj.SetColor(gfx.RED)


The last command assigns a uniform red color to all atoms.


The Graphical Entity
--------------------------------------------------------------------------------

Now let's look at the graphical entity in more detail. The graphical entity is 
responsible for the rendering of :class:`~ost.mol.EntityHandle` and 
:class:`~ost.mol.EntityView` instances.

Render Modes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Molecules are rendered in :ref:`different representations <render-mode>`. Some 
representations show the structure at atomic detail (connectivity level display 
styles), while others simplify the structure significantly (trace-based display 
styles). 

To change the display mode of an entity programmatically, use the :meth:`Entity.SetRenderMode` method:

.. code-block:: python

  # render the molecule with an interpolated line trace
  obj.SetRenderMode(gfx.SLINE)
  # switch back to simple line rendering
  obj.SetRenderMode(gfx.SIMPLE)
  # switch to helix-strand-coil cartoon
  obj.SetRenderMode(gfx.HSC)

It is also possible to change the render mode of only a part of the structure. 
To display the ligand with balls and sticks and everything else in cartoon mode,
we make use of the query language to select everything that is not part of the 
protein. Remember from above, that we stored our molecule in the `pdb` variable.

.. code-block:: python
  
  not_protein = pdb.Select('peptide=false')
  obj.SetRenderMode(gfx.HSC)
  obj.SetRenderMode(gfx.CUSTOM, not_protein)


The Basics of Coloring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The color of each atom can be set individually. The simplest coloring scheme is 
to set all atoms to the same color:

.. code-block:: python

  obj.SetColor(gfx.YELLOW) 

Assigning colors based on selections is also possible. To reproduce the 
well-known `CPK coloring <http://en.wikipedia.org/wiki/CPK_coloring>`_, 4 
consecutive coloring operations are sufficient:

.. code-block:: python

  obj.SetColor(gfx.WHITE, 'ele=C')
  obj.SetColor(gfx.BLUE, 'ele=N')
  obj.SetColor(gfx.RED, 'ele=O')
  obj.SetColor(gfx.YELLOW, 'ele=S')


But because this is such a common task, the entity provides a built-in method to color the atoms by element:

.. code-block:: python

  obj.ColorByElement()

As already seen above, the coloring operations can be chained to achieve the desired effect. For example, sometimes one wants to apply conventional CPK coloring to the protein and then change the color of C-alpha atoms to pink:

.. code-block:: python

  obj.ColorByElement()
  obj.SetColor(gfx.Color(1, 0, 1), 'aname=CA and ele=C')
  


Internally, the chain of calls to the coloring methods is stored as a list of 
coloring operations that is reapplied sequentially. This makes sure that the 
coloring remains when changing render modes or updating the geometry of the 
molecule. This storage scheme has one drawback. Consider the following code:

.. code-block:: python

  for i in range(100): 
    obj.SetColor(gfx.RED)

While the code obviously achieves the desired effect, it is far from optimal 
because all of the 100 color operations are applied every time the graphical 
entity needs an update. When you see a slow down in rendering after playing 
around with the coloring operations, consider calling 
:meth:`Entity.ClearColorOps` and reapply only the set of coloring operations 
that are required to achieve your rendering.

Coloring the Entity By Property
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most complex but also most powerful coloring method is 
:meth:`Entity.ColorBy`, which allows to color atoms by a numeric property. This 
property can either be a built-in property such as atomic b-factor, charge, 
residue number or be :doc:`custom properties <../base/generic>` assigned by the 
user or an algorithm and be defined at any level (chain, residue, atom). The 
naming of the built-in properties is identical to the properties available in 
:doc:`query language <../mol/base/query>`. For example, to achieve a heat-map 
coloring for b-factors:

.. code-block:: python  

  # will color the atoms by the atomic b-factor, setting atoms with 
  # a temperature factor of 0 to blue and atoms with a temperature factor of 
  # 100 and higher to red. Values in between are linearly interpolated.
  obj.ColorBy('abfac', gfx.BLUE, gfx.RED, 0, 100)

Fine-tuning the Coloring
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :obj:`HSC` render mode uses two colors to color the secondary structure 
elements. The main color affects the top and bottom of extended  and the outside 
of helical elements. The detail color is used for the inner side of helices and 
the rim of extended elements. This color is changed with 
:meth:`Entity.SetDetailColor`.
 
