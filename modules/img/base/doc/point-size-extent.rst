Point, Size, Extent
================================================================================

.. currentmodule:: ost.img


This document describes some of the classes used in the :mod`~ost.img` module to describe image properties such as size, pixel coordinates and extents.


Point
--------------------------------------------------------------------------------

.. class:: Point(x=0, y=0, z=0)
           Point(vec)
           Point(size)
           Point(point)

  :param x: integral x coordinate
  :type  x: int
  :param y: integral y coordinate
  :type  y: int
  :param z: integral z coordinate
  :type  z: int 
  :param vec: The vectors coordinates will be rounded to the next 
     integral value. In case of a :class:`~ost.geom.Vec4`, the coordinates will
     be divided by the homogenous coordinate.
  :type vec: :class:`~ost.geom.Vec3`, :class:`~ost.geom.Vec4`
  :param point:
  :type point: Initialize from point

  This class represents an image pixel. It is defined using three integer 
  numbers, corresponding to the pixel’s indexes along the width, height and 
  depth dimensions respectively. 
  
  .. method:: __getitem__(index)
    
    Read-write access to the coordinates.
    
    ..code-block:: python
    
      p=img.Point(1,2,3)
      print(p)
      p[1]=5

  .. method:: ToVec3()
    
    Converts the point to a :class:`~ost.geom.Vec3`
  
  .. method:: ToVec2()
  
    Converts the point to a :class:`~ost.geom.Vec2`. The z coordinate of the
    point is ommitted.
 
  .. method:: ToVec4()
  
    Converts the point to a :class:`~ost.geom.Vec4`, with the homogenous
    coordinate set to 1.
  
  


Size
--------------------------------------------------------------------------------

.. class:: Size(width, height, depth)
           Size(width, height)

  This property describes the size of an image. It is defined using three 
  integer numbers: When depth is not specified, the size is assumed to refer to 
  a 2D image two-dimensional (depth=1).

  The :class:`Size` class can be passed to the :func:`CreateImage` function to
  define the size of the image being created.

  .. code-block:: python 

    # uses the size class to create an image
    s=img.Size(40,20,30)
    i=img.CreateImage(s)
    # uses a temporary instance of the Size class
    i=img.CreateImage(img.Size(40,20,30))

Extent
--------------------------------------------------------------------------------

.. class:: Extent(first_point, last_point)
           Extent(first_point, size)
           Extent(size, center_point)

  This class describes the extent of an image in pixels. The image extent is 
  a rectangle in 2d and a cuboid in 3d, and can be defined in a number of ways: 
  giving the first and the last point of the extent, specifying the first point 
  and the size of the extent, or listing the size of the extent and its central 
  point.


  **Some examples:**

  .. code-block:: python 

    # defines a rectangular extent that goes from pixel (2,4) to pixel (5,10)  
    s=img.Extent(img.Point(2,4),img.Point(5,10))
    # defines the same extent using its first point and its size
    s=img.Extent(img.Point(2,4),img.Size(4,7))
    # again the same extent using its size and central pixel
    s=img.Extent(img.Size(4,7),img.Point(3,7))

  The :class:`Extent` class can be passed to the :func:`CreateImage` function to
  create an image with a specific extent.

  .. code-block:: python 

    # uses the Extent class to create an image
    e=img.Extent(img.Point(2,4),img.Size(4,7))
    i=img.CreateImage(e)
    # uses a temporary instance of the Extent class
    i=img.CreateImage(img.Extent(img.Point(2,4),img.Size(4,7)))

  Given an :class:`Extent`, is it possible to recover its full size, and also
  the length of each of its dimensions separately.

  **Examples:**

  .. code-block:: python 

    # gets the size of an extent (e is a 3D extent)
    s=e.GetSize()
    # gets the three dimensions separately
    w=e.GetWidth()
    h=e.GetHeight()
    d=e.GetDepth()

  When one needs to visit all the image pixels belonging to an :class:`Extent`
  in sequence, one can use an :class:`ExtentIterator`:

  .. code-block:: python 

    # defines an iterator over an extent e
    ei=ExtentIterator(e)
    # visits all the pixels in the extent and
    # prints out their values
    for pixel in ei:
      print(i.GetReal(pixel))
