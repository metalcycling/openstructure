Geometric Objects
================================================================================

.. currentmodule:: ost.geom

Geometrical Objects in Two Dimensions
--------------------------------------------------------------------------------

.. class:: Line2()
           Line2(from, to)
  
  Parametric line in two dimensions as defined by an origin and a normalized
  direction vector. The first constructor creates a line with origin (0,0) and 
  direction along the x axis. The second signature creates a line originating 
  from `from` and pointing towards `to`.
  
  .. method:: At(t)
    
    Returns the point on the line at (signed) distance t from origin.
    
    :param t: free parameter
    :type  t: float
    :rtype:   :class:`Vec2`
  

  .. method:: GetOrigin()
  
    Returns the origin of the line: Also available as :attr:`origin`.
    
    :rtype: :class:`Vec2`
    
  .. method:: GetDirection()
  
    Returns the normalized direction vector. Also available as
    :attr:`direction`.
    
    :rtype: :class:`Vec2`
    
  .. attribute:: direction
  
  .. attribute:: origin

    
.. class:: Rectangle2()
           Rectangle2(top_left, bottom_right)
  
  Axis aligned rectangle. The first signature creates a rectangle with top-left 
  corner (-1, -1) and bottom-right corner (1, 1), wheras the second method 
  allows to set the top-left and bottom-right corners directly.
  
  :type top_left: :class:`Vec2`
  :param top_left: The top-left corner
  :param bottom_right: The bottom-right corner  
  :type bottom_right: :class:`Vec2`  

  .. method:: GetWidth()
    
    Returns the width of the rectangle. Also available as :attr:`width`.  
  
  .. method:: GetHeight()
  
    Returns the height of the rectangle. Also available as :attr:`height`.  
    
  .. attribute:: width
  
    :type: float
    
  .. attribute:: height

    :type: float  
    
  .. method:: GetStart()
  
    Get top-left corner
    
    :rtype: :class:`Vec2`
    
  .. method:: GetEnd()

    Get bottom-right corner

    :rtype: :class:`Vec2`    
  
  .. method:: SetStart(top_left)
  
    Set top-left corner, leaving the bottom-right corner untouched.
    
  .. method:: SetEnd(bottom_right)
  
    Set the bottom-right corner, leaving the top-left corner untouched.

.. class:: Circle2()
           Circle2(circle)
           Circle2(center, radius)

  The first signature creates a circle centered at (0, 0) and radius 1.0. The 
  second signature creates a circle with the same paramters as `circle`. The 
  third signature creates a new circle with given center and radius.
     
  .. method:: SetCenter(center)
  
    Set center of circle
    
    :type center: :class:`Vec2`
    :param center: The new center
    
  .. method:: SetRadius(radius)
  
    Set radius of circle
  
    :type radius: float
    :param center: The new radius
    
  .. method:: GetCenter()
  
    Returns the center of the circle
    
  .. method:: GetRadius()
  
    Returns the radius of the circle
    
  .. method:: GetArea()
  
    Returns the area of the circle
    
  .. method:: GetCircumference()
  
    Returns the circumference of the circle


.. class:: Ellipse2()
           Ellipse2(center, a, b, gamma)

  An ellipse is defined by a center, two principal axis and gamma that 
  defines the angle between the first principal axis an the x-axis.
  
  .. method:: At(t)
  
    ?
  
  .. method:: AtAngle(angle)
  
    ?
  
  .. method:: GetBoundingBox()
  
    Returns the bounding rectangle (axis-aligned) of the ellipse
    
    :rtype: :class:`Rectangle2`
    

  .. method:: GetA()
  
    Returns the first principal-axis
  
  .. method:: GetB()
  
    Returns the second principal-axis
    
  .. method:: GetGamma()
    
    Returns the angle of the first principal axis to the x-axis
    
  .. method:: GetArea()
  
    Returns the area of the ellipse
    
  .. method:: GetOrigin()
  
    Returns the center of the ellipse
    
  .. method:: SetA(a)
  
    Set the length of the first principal axis
    
  .. method:: SetB(b)
  
    Set the length of the second principal axis
    
  .. method:: SetGamma(gamma)
  
    Set the angle of the first principal axis to the x-axis
    
  .. method:: SetOrigin(ori)
  
    Set the center of the ellipse
    
Geometrical Objects in Three Dimensions
-------------------------------------------------------------------------------

.. class:: Line3()
           Line3(from, to)
  
  Parametric line in three dimensions as defined by an origin and a normalized
  direction vector. The first constructor creates a line with origin (0,0) and 
  direction along the x axis. The second signature creates a line originating 
  from `from` and pointing towards `to`.
  
  .. method:: At(t)
    
    Returns the point on the line at (signed) distance t from origin.
    
    :param t: free parameter
    :type  t: float
    :rtype:   :class:`Vec3`
  

  .. method:: GetOrigin()
  
    Returns the origin of the line: Also available as :attr:`origin`.
    
    :rtype: :class:`Vec3`
    
  .. method:: GetDirection()
  
    Returns the normalized direction vector. Also available as
    :attr:`direction`.
    
    :rtype: :class:`Vec3`
    
  .. attribute:: direction
  
    :type: :class:`Vec3`
  
  .. attribute:: origin
    
    :type: :class:`Vec3`

.. class:: Plane()
           Plane(p1, p2, p3)
           Plane(x, y, z, p)
           Plane(line, point)
           Plane(point, normal)
  
  A plane in 3d-space. The plane can be constructed by either passing in 3 
  points (p1, p2, p3), a normal and a point, the four parameters that define the 
  implicit plane equation (`x`, `y`, `z`, `p`) or a line and a point.
  
  .. method:: GetNormal()
  
    Returns the normal of the plane. Also available as :attr:`normal`
    
    :rtype: :class:`Vec3`
    
  .. method:: GetP()
  
    Returns the plane offset, i.e. the projection of any point on the plane onto
    the normal. Also available as :attr:`p`.
    
    :rtype: float
    
  .. method:: GetOrigin()
  
    Get the origin of the plane. Also available as :attr:`origin`.
    
    :rtype: :class:`Vec3`

  .. attribute:: origin
  
    :type: :class:`Vec3`
  .. attribute:: normal
  
    :type: :class:`Vec3`
      
  .. attribute:: p
  
    :type: float
 

.. class:: Sphere()
           Sphere(center, radius)
           
  Represents a sphere in 3d space. The first constructor creates a sphere with 
  radius 1, centered at (0, 0, 0), the second allows you to set the radius and 
  center directly.
  
  :param center: The center
  :type  center: :class:`Vec3`
  :param radius: The radius
  :type  radius: float
  
  .. attribute:: radius
    
    The radius of the sphere. Read-write. Also available as :meth:`GetRadius`, 
    :meth:`SetRadius`.
    
    :type: float
    
  .. attribute:: origin
  
    The center of the sphere. Read-write. Also available as :meth:`GetOrigin`,
    :meth:`SetOrigin`.
    
    :type: :class:`Vec3`

  .. method:: GetOrigin()  
  
    See :attr:`origin`
    
  .. method:: SetOrigin(origin)
    
    See :attr:`origin`
    
  .. method:: GetRadius()
  
    See :attr:`radius`
     
  .. method:: SetRadius(radius)
    
    See :attr:`radius`

.. class:: AlignedCuboid(min, max)
  
  Axis aligned cuboid is a cuboid whose axes are aligned to the x-, y-, and z- 
  axes of the coordinate system. For arbitrarily oriented bounding cuboid 
  class, see :class:`Cuboid`.
  
  .. method:: GetMin()
  
    Get minimum coordinate, i.e. the lower bound of x-, y-, and z for 
    any point in the cuboid
    
    :rtype: :class:`Vec3`
    
  .. method:: GetMax()
  
    Get maximum coordinate, i.e. the upper bound of x-, y-, and z for
    any point in the cuboid.
    
    :rtype: :class:`Vec3`

.. class:: CuboidAxis()
           CuboidAxis(dir, half_extent)
  
  A cuboid axis is defined by a half-extent, and a direction vector. This class 
  is used in together with the :class:`Cuboid` class.
  
  :param dir: Direction vector, will be normalized
  :type  dir: :class:`Vec3`
  :param half_extent: The half extent
  :type half_extent: float
  
  .. attribute:: vector
  
    The normalized direction vector of the cuboid axis. Also available as 
    :meth:`GetVector`
    
    :type: :class:`Vec3`
  
  .. attribute:: half_extent
  
    The half extent of the cuboid axis is the magnitude of the cuboid 
    axis measured from the center to the corner. Also available as
    :meth:`GetHalfExtent`
    
    :type: float
  
  .. attribute:: extent
  
    The extent of the cuboid axis. This value is always twice the 
    :attr:`half_extent`. Read-only. Also available as 
    :meth:`GetExtent`.
    
    :type: float
    
  .. method:: GetHalfExtent()
  
    See :attr:`half_extent`
  .. method:: GetExtent()
  
    See :attr:`extent`
    
  .. method:: GetVector()
  
    See :attr:`vector`
    
.. class:: Cuboid(center, axis_a, axis_b, axis_c)
  
  An arbitrarily oriented cuboid defined by a center and 3 axis. The 3 cuboid 
  axis are stored in the order they are passed to the constructor. This means, 
  that there is no guarantee that the 3 axes form a right-handed coordinate 
  system. If a right-handed coordinate system is a requirement, you have to 
  ensure this on your own:
  
  .. code-block:: python
  
    center=...
    axis_a=geom.CuboidAxis(...)
    axis_b=geom.CuboidAxis(...)
    axis_c=geom.CuboidAxis(geom.Cross(axis_a.vector, axis_b.vector), ...)
    
    cuboid=geom.Cuboid(center, axis_a, axis_b, axis_c)
  
  :param center: The center
  :type  center: :class:`Vec3`
  :param axis_a: The first axis
  :type  axis_a: :class:`CuboidAxis`
  :param axis_b: The second axis
  :type  axis_b: :class:`CuboidAxis`
  :param axis_c: The third axis
  :type  axis_c: :class:`CuboidAxis`
  
  .. attribute:: center
  
    The center of the cuboid.
    
    :type: :class:`Vec3`
    
  .. attribute:: axis_a
  
    The first cuboid axis
    
    :type: :class:`CuboidAxis`
    
  .. attribute:: axis_b

    The second cuboid axis

    :type: :class:`CuboidAxis`
    
  .. attribute:: axis_c

    The third cuboid axis

    :type: :class:`CuboidAxis`    

Operations on Geometrical Objects
--------------------------------------------------------------------------------

.. function:: Angle(lhs, rhs)

  Calculate the angle (in radians) between `lhs` and `rhs`. 
  
  :param lhs: First object
  :type  lhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :param rhs: Second object
  :type  rhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :rtype: float

.. function:: IntersectionPoint(lhs, rhs)

  Calculates and returns the intersection point between `lhs` and `rhs`
  
  :param lhs: First object
  :type  lhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :param rhs: Second object
  :type  rhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :raises: :exc:`~exceptions.RuntimeError` when the two objects do not intersect
  :rtype: :class:`Vec3` (:class:`Vec2` in case of :class:`Line2`)

.. function:: IntersectionLine(plane2, plane2)

  Returns the intersection line between `plane1` and `plane2`.
  
  :param plane1: The first plane
  :type  plane1: :class:`Plane`
  :param plane2: The second plane
  :type  plane2: :class:`Plane`  
  
  :raises: :exc:`~exceptions.RuntimeError` if the two planes are parallel.
  
.. function:: Distance(lhs, rhs)

  Returns the minimal distance between `lhs` and `rhs`. 
  
  :param lhs: First object
  :type  lhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :param rhs: Second object
  :type  rhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :rtype: float

.. function:: IsOnLine(line, point, epsilon=geom.EPSILON)

  Check whether `point` lies on `line` and returns true if point i no further 
  away than `epsilon`.
  
  :rtype: bool

.. function:: IsInPlane(plane, object, epsilon=geom.EPSILON)
  
  Check whether `object` lies in `plane` and returns true if the difference is 
  no bigger than `epsilon`.
  
  :param plane: The plane
  :type  plane: :class:`Plane`
  :param object: The object
  :type  object: :class:`Vec3` or :class:`Line3`
  
  :rtype: bool

.. function:: AreParallel(lhs, rhs, epsilon=geom.EPSILON)

  Check whether `lhs` and `rhs` are parallel and returns true, if the difference 
  is below the given treshold `epsilon`.
  
  :param lhs: First object
  :type  lhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :param rhs: Second object
  :type  rhs: :class:`Line2`, :class:`Line3`, :class:`Plane`
  
  :rtype: bool

.. function:: AreIntersecting(line1, line2, epsilon=geom.EPSILON)

  Check whether `line1` and `line2` are intersecting and returns true, if they 
  intersect below the given threshold `epsilon`.
  
  :param lhs: First line
  :type  lhs: :class:`Line2`, :class:`Line3`
  
  :param rhs: Second line
  :type  rhs: :class:`Line2`, :class:`Line3`
  
  :rtype: bool
  

  
.. function:: IsInSphere(sphere, point)

  Check whether the `sphere` contains `point`.
  
  :rtype: bool
  

.. function:: MinDistance(list1, list2)

  Calculate the minimal distance between two sets of points
  
  :param list1: the first set of points
  :type list1: :class:`~ost.geom.Vec3List`
  
  :param list2: the second set of points
  :type list2: :class:`~ost.geom.Vec3List`
  
  