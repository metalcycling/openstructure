Matrices
================================================================================

.. currentmodule:: ost.geom

The :mod:`~ost.geom` module defines matrices in two, three and four dimensions. 
All matrices store the values in row-major order, meaning that, the matrix ((1, 
2), (3,4)) stores the values as (1, 2, 3, 4). This is illustrated in 
the following code examples:

.. code-block:: python

  m=geom.Mat2(1, 2, 3, 4)
  print(m) # will print {{1,2},{3,4}}
  print(m[(0,0)], m[(0,1)], m[(1,0)], m[(1,1)]) # will print 1, 2, 3, 4

Matrices support arithmetic via overloaded operators. The following operations are 
supported:

  * adding and subtracting two matrices
  * negation
  * multiplication of matrices
  * multiplying and dividing by scalar value

The Matrix Classes
--------------------------------------------------------------------------------

.. class:: Mat2()
           Mat2(d00, d01, d10, d11)

  2x2 real-valued matrix. The first signature creates a new identity matrix. The 
  second signature initializes the matrix in row-major order. 

  .. staticmethod:: Identity()
  
    Returns the 2x2 identity matrix
  
  
.. class:: Mat3()
           Mat3(d00, d01, d02, d10, d11, d12, d20, d21, d22)

  3x3 real-valued matrix. The first signature creates a new identity matrix. The 
  second signature initializes the matrix in row-major order.
  
  .. staticmethod:: Identity()
  
    Returns the 3x3 identity matrix
    
.. class:: Mat4()
           Mat4(d00, d01, d02, d03, d10, d11, d12, d13, d20, d21, d22, d23, d30, d31, d32, d33)

  4x4 real-valued matrix. The first signature creates a new identity matrix. The 
  second signature initializes the matrix in row-major order.
  
  .. method:: ExtractRotation()
    
    Returns the 3x3 submatrix
    
  .. method:: PasteRotation(mat)
  
    Set the 3x3 submatrix of the top-left corner to `mat`
    
  .. method:: ExtractTranslation()
  
    Extract translation component from matrix. Only meaningful when matrix 
    is a combination of rotation and translation matrices, otherwise the result 
    is undefined.
    
  .. PasteTranslation(trans)
  
    Set the translation component of the matrix to `trans`
    
    :param trans: The translation
    :type  trans: :class:`Vec3`
    
  .. staticmethod:: Identity()
  
    Returns the 4x4 identity matrix

Functions Operating on Matrices
--------------------------------------------------------------------------------
.. function:: Equal(lhs, rhs, epsilon=geom.EPSILON)

  Compares the two matrices `lhs` and `rhs` and returns True, if all 
  of the element-wise differences are smaller than epsilon. `lhs` 
  and `rhs` must be matrices of the same dimension.
  
  :param lhs: First matrix
  :type  lhs: :class:`Mat2`, :class:`Mat3` or :class:`Mat4`
  :param rhs: Second matrix
  :type  rhs: :class:`Mat2`, :class:`Mat3` or :class:`Mat4`

.. function:: Transpose(mat)

  Returns the transpose of `mat`
  
  :param mat: The matrix to be transposed
  :type  lhs: :class:`Mat2`, :class:`Mat3` or :class:`Mat4`

.. function:: Invert(mat)

  Returns the inverse of `mat`
  
  :param mat: The matrix to be inverted
  :type  mat: :class:`Mat2`, :class:`Mat3` or :class:`Mat4`
  
  What happens when determinant is 0?
  
.. function:: CompMultiply(lhs, rhs)

  Returns the component-wise product of `lhs` and `rhs`. `lhs` and 
  `rhs` must be vectors of the same dimension.
  
  :param lhs: The lefthand-side vector
  :type  lhs: :class:`~Vec2`, :class:`~Vec3` or 
              :class:`~Vec4`
  :param rhs: The righthand-side vector
  :type  rhs: :class:`~Vec2`, :class:`~Vec3` or 
              :class:`~Vec4`
              
.. function:: CompDivide(lhs, rhs)

  Returns the component-wise quotient of `lhs` divided by `rhs`. `lhs` 
  and `rhs` must be vectors of the same dimension.

  :param lhs: The lefthand-side vector
  :type  lhs: :class:`~Vec2`, :class:`~Vec3` or 
              :class:`~Vec4`
  :param rhs: The righthand-side vector
  :type  rhs: :class:`~Vec2`, :class:`~Vec3` or 
              :class:`~Vec4`

.. function:: Det(mat)

  Returns the determinant of `mat`
  :param mat: A matrix
  :type  mat: :class:`~Mat2`, :class:`~Mat3` or :class:`~Mat4`

.. function:: Minor(mat, i, j)
  
  Returns the determinant of the 2x2 matrix generated from `mat`  by 
  removing the ith row and jth column.

.. function:: EulerTransformation(phi, theta, xi)
  
  Returns a rotation matrix for the 3 euler angles `phi`, `theta`, and
  `xi`. The 3 angles are given in radians.
  
  
.. function:: AxisRotation(axis, angle)

  Returns a rotation matrix that represents a rotation of `angle` 
  around the `axis`. 
  
  :param axis: The rotation axis. Will be normalized
  :type  axis: :class:`Vec3`
  :param angle: Rotation angle (radians) in clockwise direction when 
      looking down the axis.
  
.. function:: OrthogonalVector(vec)

  Get arbitrary vector orthogonal to `vec`. The returned vector is of length 
  1, except when `vec` is a zero vector. In that case, the returned vector is 
  (0, 0, 0).
  
  :param vec: A vector of arbitrary length
  :type vec: :class:`Vec3`
