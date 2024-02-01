Surfaces
================================================================================

.. currentmodule:: ost.mol

.. class:: SurfaceHandle

  TODO

  .. method:: Attach(entity, cutoff)
    
    Attach an :class:`EntityHandle` or :class:`EntityView` to a molecular
    surface. For each vertex of the surface, the atom which is closest to this
    vertex and within the cutoff distance is attached. Like this, the atom and
    its properties can be retrieved efficiently, which is e.g. used for coloring
    a surface according to the underlying entity (see
    :meth:`gfx.Surface.ColorBy`).
    
    Be aware that increasing the cutoff distance will significantly increase the
    time needed to attach an entity to the surface.
    
    :param entity: entity to be attached to the surface
    :type  entity: :class:`EntityHandle`, :class:`EntityView`
    :param cutoff: distance cutoff
    :type  cutoff: float
    :rtype: None

  .. method:: FindWithin(position, cutoff)

    Returns all :class:`surface vertices<SurfaceVertex>` within a certain cutoff 
    distance of a specified position.
     
    :param position: position
    :type  position: :class:`~ost.geom.Vec3`
    :param cutoff:   distance cutoff
    :type  cutoff:   float
    :rtype:          :class:`SurfaceVertexList`

  .. method:: GetTriIDList()
  
    Returns a list containing all IDs of all triangles of a surface.

    :rtype: list of ints

  .. method:: GetVertex(id)

    Returns a :class:`SurfaceVertex` with the specified ID.

    :param id: vertex id
    :type  id: int
    :rtype: :class:`SurfaceVertex`

  .. method:: GetVertexIDList()

    Returns a list containing all IDs of all vertices of a surface.
    
    :rtype: list of ints

  .. method:: Invert()

    Does an in place inversion of all vertex normal vectors.
     
    :rtype: None

.. class:: SurfaceVertex

  TODO
  
  .. method:: SurfaceVertex()
              SurfaceVertex(pos, norm, type, atom)
              
    Creates a new SurfaceVertex.
              
    :param pos:  vertex position
    :type  pos:  :class:`geom.Vec3`
    :param norm: vertex normal vector
    :type  norm: :class:`geom.Vec3`
    :param type: vertex type
    :type  type: int
    :param atom: attached atom
    :type  atom: :class:`AtomHandle` 
              
  .. attribute:: Atom
    
    :class:`AtomHandle` attached to the surface vertex.
    
    The AtomHandle is invalid, if no :class:`~ost.mol.EntityHandle` has been
    attached or if no atom was found within the cutoff distance.
    
    See :meth:`SurfaceHandle.Attach` for further details.
    
    :type: :class:`AtomHandle`
  
  .. attribute:: Normal
    
    Normal vector of the surface vertex.
  
    :type: :class:`~ost.geom.Vec3`
  
  .. attribute:: Position
  
    Spatial position of the surface vertex.
  
    :type: :class:`~ost.geom.Vec3`

.. class:: SurfaceVertexList

  .. method:: append(vertex)

    Add a single surface vertex.

    :param vertex: surface vertex
    :type  vertex: :class:`SurfaceVertex`
    :rtype: None

  .. method:: extend(list)
  
    Add a list of surface vertices.

    :param list: list of surface vertices
    :type  list: :class:`SurfaceVertexList`
    :rtype: None

.. class:: SurfaceTriIDList

  .. method:: append(id)
  
    Add a single triangle ID.

    :param id: trinangle vertex id
    :type  id: int
    :rtype: None

  .. method:: extend(list)

    Add a list of triangle IDs.

    :param list: list of triangle vertex ids
    :type  list: list of ints
    :rtype: None




