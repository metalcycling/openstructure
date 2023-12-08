Vectors
================================================================================
.. currentmodule:: ost.geom

The :class:`Vec2`, :class:`Vec3`, :class:`Vec4` classes implement vectors in 2, 
3 and four dimensions. They support basic arithmetic via overloaded operators. 
Essentially, the following basic operations are available:

 * adding and subtracting two vectors
 * negation
 * multiplying and dividing by scalar value

This is shown in the following example:

.. code-block:: python
   
   vec_a=geom.Vec2(1, 0)
   vec_b=geom.Vec2(0, 1)
   print(vec_a, vec_b)
   print(vec_a+vec_b)
   print(vec_a*3-vec_b)

The standard vector operations are implemented as :ref:`free standing functions 
<vector-functions>`:


.. code-block:: python
   
  vec_a=geom.Vec3(1, 0, 0)
  vec_b=geom.Vec3(0, 1, 0)
  
  print(geom.Dot(vec_a, vec_b))
  print(geom.Cross(vec_a, vec_b))
  
  print(geom.Normalize(geom.Vec3(1, 1, 0)))
  
  print(geom.Length(geom.Vec3(1, 1, 1)))


Vector Classes
--------------------------------------------------------------------------------
.. class:: Vec2(x=0.0, y=0.0, z=0.0)
           Vec2(vec)
           
   Real-valued vector in 2 dimensions.
   
   :param   x: x coordinate
   :type    x: float or int   
   :param   y: y coordinate
   :type    y: float or int
   :param vec: the coordinates are set to the coordinates of vec. If vec is a
               :class:`Vec2`, the coordinates are copied directly, If vec is a
               :class:`Vec3`, the x and y coordinates are set to the
               coordinates of vec and z is silently swallowed. If vec is of
               type :class:`Vec4`, x and y are divided by the homogenous
               coordinate w, raising a DivideByZeroException when w is zero.
   :type  vec: Vec2, Vec3 or Vec4
   
   .. attribute:: x

    The x-coordinate of the vector.

    :type: float

   .. attribute:: y

     The y-coordinate of the vector.
.. class:: Vec3(x=0.0, y=0.0, z=0.0)
           Vec3(vec)
    
    Real-valued vector in 3 dimensions.
        
    :param   x: x coordinate
    :param   y: y coordinate
    :param   z: z coordinate
    :param vec: the coordinates are set to the coordinates of vec. If vec is a
                :class:`Vec3`, the coordinates are copied directly, If vec is a
                :class:`Vec2`, the x and y coordinates are set to the
                coordinates of vec and z is initialized to zero. If vec is of
                type :class:`Vec4`, x, y and z are divided by homogenous
                coordinate w, raising a DivideByZeroException when w is zero.
    :type  vec: Vec2, Vec3 or Vec4
    
    .. attribute:: x
    
     The x-coordinate of the vector.
    
     :type: float or int
     
    .. attribute:: y
      
      The y-coordinate of the vector.
    
      :type: float or int
      
    .. attribute:: z
    
      The z-coordinate of the vector.
      
      :type: float or int

.. class:: Vec4(x=0.0, y=0.0, z=0.0, w=1.0)
           Vec4(vec)

    Real-valued vector in 4 dimensions.

    :param   x: x coordinate
    :type    x: float or int
    :param   y: y coordinate
    :type    y: float or int
    :param   z: z coordinate
    :type    z: float or int       
    :param   w: w (homogenous) coordinate
    :type    w: float or int       
    :param vec: the coordinates are set to the coordinates of vec. If vec is a
                :class:`Vec4`, the coordinates are copied directly, If vec is a
                :class:`Vec2`, the x and y coordinates are set to the
                coordinates of vec and z and w are initialized to 0 and 1,
                respectively. If vec is of type :class:`Vec4`, x, y and z are
                divided by homogenous coordinate w, raising a
                DivideByZeroException when w is zero.

    .. attribute:: x

     The x-coordinate of the vector.

     :type: float or int

    .. attribute:: y

      The y-coordinate of the vector.

      :type: float or int

    .. attribute:: z

      The z-coordinate of the vector.

      :type: float or int
  
    .. attribute:: w

      The homogenous coordinate.

      :type: float or int


.. _vector-functions:

Functions Operating on Vectors     
--------------------------------------------------------------------------------

.. function:: Cross(vec_a, vec_b)

   Cross product of `vec_a` and `vec_b`

   :type vec_a: Vec3
   :type vec_b: Vec3

.. function:: Dot(vec_a, vec_b)

   Dot (scalar) product of `vec_a` and `vec_b`

   :param vec_a: first vector
   :type  vec_a: Vec3
   :param vec_b: second vector
   :type  vec_b: Vec3

.. function:: Length(vec)

   Length of vector

   :param vec:
   :type  vec: Vec2, Vec3 or Vec4

.. function:: Length2(vec)

  Returns the squared length of `vec`

  :param vec:
  :type  vec: Vec2, Vec3 or Vec4

.. function:: Normalize(vec)

   Returns a normalized version of `vec`

   :param vec:  Vector to be normalized
   :type  vec: Vec2, Vec3 or Vec4

