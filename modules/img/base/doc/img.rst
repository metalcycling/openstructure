:mod:`~ost.img` - Images and Density Maps
================================================================================

.. toctree::
  :hidden:

  point-size-extent

.. module:: ost.img
   :synopsis: Images and density maps
              
Introduction : The ImageHandle
--------------------------------------------------------------------------------

OpenStructure offers extensive processing capabilities for planar 2d images and
3d maps using the :mod:`img` module. The main datastructure for images is the 
:class:`ImageHandle` class.

ImageHandles provide a clean and efficient interface to interact with images and 
maps. An :class:`ImageHandle` can store an image in either real (:obj:`SPATIAL`) 
or Fourier (:obj:`FREQUENCY`) space and always keep track of the currently 
active  domain. This means,for example that one can apply a Fourier 
transformation to an :class:`ImageHandle` containing a :obj:`SPATIAL` image and 
the image will correctly  identify the new active domain as  :obj:`FREQUENCY`. 
The :class:`ImageHandle` also understands, for example, that applying a 
Fourier Transform to a centrosymmetric :obj:`FREQUENCY` image results in a real 
:obj:`SPATIAL` image, but applying it to a non-centrosymmetric one results in a 
complex :obj:`SPATIAL` image.

Furthermore, the image will make sure that real and Fourier space 
information about the image are always in sync. If, for example, the pixel 
sampling is changed while the current active domain is real space, the pixel 
sampling in Fourier space will be adjusted accordingly, and vice versa.

Creating and visualizing ImageHandles
--------------------------------------------------------------------------------

As a first step, enter the following lines in the OpenStructure python console:

.. code-block:: python
  
  im=img.CreateImage(img.Size(200,200))

This will create an empty, 2D image, with a height and width of 200 pixels, 
whose origin (ie the pixel with the coordinates <0,0>) is in the top-left 
corner.

.. code-block:: python
  
  v=gui.CreateDataViewer(im)

A viewer window will pop up (see below), showing a white frame on a black 
background. The inner area of the white frame is the image, which is empty.

Reading and writing into an image
-------------------------------------------------------------------------------

Data can be read and written from and into an image using the following 
commands:

.. code-block:: python

  # writes the real value 23.4 into pixel 10,10
  im.SetReal(img.Point(10,10),23.4)
  # reads the value in pixel 10,10
  val=im.GetReal(img.Point(10,10))
    
The complex equivalents are also available    

.. code-block:: python

  # writes the complex value value 2+3j into pixel 10,10
  im.SetComplex(img.Point(10,10),2+3j)
  # reads the value in pixel 10,10
  val=im.GetComplex(img.Point(10,10))

The image knows in which domain it is, and will adjust the type of data being 
written accordingly. For example, if one writes a complex value in a 
:obj:`SPATIAL` image, the value will be automatically converted to a real one by 
taking the amplitude of the complex number. On the other hand, if one writes a 
real value in a :obj:`FREQUENCY` image, the value is automatically
converted to complex by setting the imaginary part to 0.

Image properties
-----------------

.. _data-domain:

The data domain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The data domain of an image specifies wether the image contains data in the 
spatial or frequency domain. A :obj:`HALF_FREQUENCY` domain also exists, 
representing centrosymmetric frequency data (such as the data coming from the 
Fourier transform of an image from the real spatial domain)

.. data:: SPATIAL

  Real-valued spatial images

.. data:: COMPLEX_SPATIAL

  Complex-valued spatial images, i.e. resulting from a Fourier transform of
  the :obj:`FREQUENCY` domain.
  
.. data:: FREQUENCY
  
  Complex frequeny domain.
  
.. data:: HALF_FREQUENCY

  Centrosymmetric frequency images
    


.. _spatial-origin:

The spatial origin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The spatial origin of an image is the first pixel of its extent. Specifically, 
this is the top left pixel for 2D images and the top-front-left  corner for 3D
images.

.. _absolute-origin:

The absolute origin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The absolute origin of an image describes the coordinates, in the absolute 
reference system used by OpenStructure, of the pixel in with all indexes equal 
to 0. Please note that the pixel does not necessarily need to belong to the 
extent of the image.

.. _pixel-sampling:

Pixel sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The pixel sampling property of an image describes the size of its pixels. For
the same image, the size of pixels in the :obj:`SPATIAL` and in the
:obj:`FREQUENCY` :ref:`data-domain` are obviously interdependent.
OpenStructure takes care of the transformation and allows access to both pixel 
sampling irrespective of the current image domain.

The ImageHandle class
--------------------------------------------------------------------------------

The public interface of the ImageHandle class provides many ways to manipulate 
image properties. What follows is a brief description of the most important 
methods and attributes of the ImageHandle class.

.. class:: ImageHandle

  .. attribute:: size
  
    The size of the image. Read-only.
    
    :type: :class:`Size`
    
  .. attribute:: extent
  
    The extent of the image. Read-only.
    
    :type: :class:`Extent`
  
  .. attribute:: type
  
    The DataType of an image represents the nature of the data it contains. An
    image can contain :obj`REAL` or :obj:`COMPLEX` values.
  
  .. attribute:: absolute_origin
  
    The absolute origin of an image describes the coordinates, in the absolute 
    reference system used by OpenStructure, of the pixel in with all indexes 
    equal to 0. Please note that the pixel does not necessarily need to belong 
    to the extent of the image. Read-write.
    
    :type: :class:`~ost.geom.Vec3`
    

  .. attribute:: spatial_origin
  
    The spatial origin of an image is the first pixel of its extent. 
    Specifically, this is the top left pixel for 2D images and the 
    top-front-left  corner for 3Dimages.
    
    :type: :class:`Point`

  .. attribute:: domain
  
    The current domain of the image. See :ref:`data-domain`. Read-only.
    

  .. method:: Apply(algorithm)

     Applies an algorithm on an image and returns a new ImageHandle containing 
     the modified image. See :doc:`../alg/alg` 

     :param algorithm: algorithm
     :type  algorithm: Instance of an algorithm class - See :doc:`../alg/alg`.
     :rtype: :class:`ImageHandle`

  .. method:: ApplyIP(algorithm)

     Applies an algorithm on an image, overwriting the current content. See 
     `:doc:`../alg/alg`

     :param algorithm: algorithm
     :type  algorithm: Instance of an algorithm class - See: doc:`../alg/alg`.

  .. method:: CenterSpatialOrigin()

     Sets the :ref:`spatial-origin` of an image in such a way that the central 
     pixel has all 0 indexes, e.g. (0,0) or (0,0,0) for 3d images.
    
  .. method:: CoordToIndex(coord)

     Returns the indexes of an image corresponding to the specified absolute 
     coordinates. (See :ref:`absolute-origin`). A given set of absolute 
     coordinates will almost never fall exactly at the center of a pixel, so 
     this method return fractional indexes. 

     :param coord: Absolute coordinates
     :type  coord: :class:`~ost.geom.Vec3`
     :rtype: Vec3

  .. method:: Copy()
     
     Creates a and returns a copy of the ImageHandle. The new handle does not 
     point to the same underlying data as the old one, but to a complete and 
     separate copy of the original data.

     :rtype: :class:`ImageHandle`

  .. method:: Extract(extent)
     
     Creates and returns a new image that contains a copy of a portion of the 
     original image. The extracted image keeps the same :ref:`data-domain` of 
     the original image, but extractions from images in the 'FREQUENCY' or 
     'HALF FREQUENCY' domains result in :obj:`COMPLEX ` :obj:`SPATIAL` images. 
     This transformation is necessary, since the there is no guarantee that the 
     extracted :obj:`FREQUENCY` sub-image is centered around the origin and 
     hence back-transformable to :obj:`SPATIAL`.

     :param extent: Portion of the image to extract
     :type  extent: :class:`Extent`
     :rtype: :class:`ImageHandle`

  .. method:: FractionalIndexToCoord(frac_pixel)

     Same as :meth:`IndexToCoord`, but introduces subpixel precision by 
     accepting fractional numbers for pixel indexes.
     
     :param frac_pixel: Fractional pixel indexes
     :type  frac_pixel: :class:`~ost.geom.Vec3`
     :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetAbsoluteOrigin()
    
     Returns the :ref:`absolute-origin` of an image

     :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetComplex(pixel)

     Returns the complex value of the specified image pixel. If the image holds 
     data of the 'REAL' :ref:`data-domain`, the method return s complex value 
     with the pixel content as real part and a null imaginary part.

     :param pixel: Image pixel
     :type  pixel: :class:`Point`
     :rtype: complex

  .. method:: GetDomain()
    
     See :attr:`domain`

  .. method:: GetExtent()

    See :attr:`extent`


  .. method:: GetFrequencySampling()
     
     Returns the :ref:`pixel-sampling` of the image in the :obj:`FREQUENCY` 
     :ref:`data-domain`
       
     :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetIntpolComplex(frac_pixel)

     Returns the interpolated complex value of the virtual pixel corresponding 
     to the specified fractional indexes. This is computed by calculating a 
     weighted vector sum of the values of the surrounding pixels. If the image 
     holds data of the 'REAL' :ref:`data-domain`, the method computes the 
     interpolated value using bilinear interpolation (trilinear for 3D images), 
     then returns a complex value with the interpolated value as real part and a 
     null imaginary part.

     :param frac_pixel: Fractional pixel indexes
     :type  frac_pixel: :class:`~ost.geom.Vec3`
     :rtype: complex

  .. method:: GetIntpolReal(frac_pixel)

     Returns the interpolated value of the virtual pixel corresponding to the 
     specified fractional indexes. This is computed by using bilinear 
     interpolation (trilinear for 3D images). If the image holds data of the 
     :obj:`COMPLEX ` :ref:`data-domain`, the method computes the interpolated 
     value as a weighted vector sum of the values of the surrounding pixels, 
     then returns the amplitude of the interpolated value. 

     :param frac_pixel: Fractional pixel indexes
     :type  frac_pixel: :class:`~ost.geom.Vec3` 
     :rtype: float

  .. method:: GetPixelSampling()
   
     Returns the :ref:`pixel-sampling` of the image in the current 
     :ref:`data-domain`.

     :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetReal(pixel)

     Returns the value of the specified image pixel. If the image holds data of 
     the :obj:`COMPLEX ` :ref:`data-domain`, the method return the amplitude of 
     the pixel content.

     :param pixel: Image pixel
     :type  pixel: :class:`Point`
     :rtype: float

  .. method:: GetSize()

     Returns the :class:`Size` of the image.

     :rtype: :class:`Size`

  .. method:: GetSpatialOrigin()

     Returns the :ref:`spatial-origin` of the image.

     :rtype: :class:`Point`

  .. method:: GetSpatialSampling()

     Return the :ref:`pixel-sampling` of the image in the :obj:`SPATIAL`
     :ref:`data-domain`.

     :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetType()
     
     Returns the :ref:`data-domain` of the image (:obj:`REAL` or :obj:`COMPLEX`)
  
     :rtype: DataType  ???????????

  .. method:: IndexToCoord(pixel)
 
     Returns the absolute coordinates (See :ref:`absolute-origin`) corresponding 
     to the pixel with the specified indexes. Please note this method always
     returns the coordinates corresponding to the center of the pixel. When 
     sub-pixel precision is needed, the method :meth:`FractionalIndexToCoord` 
     can be used.
    
     :param pixel:
     :type  pixel: :class:`Point`
     :rtype: :class:`~ost.geom.Vec3`

  .. method:: IsFrequency()

     Returns true if the :ref:`data-domain` of the image is `FREQUENCY` or 
     `HALF-FREQUENCY`     

     :rtype: bool

  .. method:: IsSpatial()

     Returns true if the :ref:`data-domain` of the image is :obj:`SPATIAL`.  

     :rtype: bool

  .. method:: IsValid()

     Returns true, when the image handle is valid, false if not. It is an error 
     to call any method on the image handle other than `IsValid` when this 
     method returns false.

     :rtype: bool

  .. method:: Paste(source_image)
   
     Copies the content of a different image into the current one, overwriting 
     pre-existing data . The values of pixels with given indexes in the source 
     image are copied into pixels with the same indexes in the target image. If 
     the two images have different extents (see :class:`Extent`), pixels that do 
     not exist in both are ignored. Please notice that this method only copies 
     the pixel content: all other properties of the image are left untouched. 

     :param source_image: Source image that gets pasted into the current one
     :type  source_image: :class:`ImageHandle`

  .. method:: SetAbsoluteOrigin(coord)
  
     Sets the :ref:`absolute-origin` of the image to the specified coordinates  
  
     :param coord: Absolute coordinates
     :type  coord: :class:`~ost.geom.Vec3`
 
  .. method:: SetComplex(pixel, value)

     Sets the content of the specified pixel to the provided value.  If the 
     image holds data of the 'REAL' :ref:`data-domain`, the method sets the pixel 
     to the amplitude of the provided.
     value.

     :param pixel: Image pixel
     :type  pixel: :class:`Point`
     :param value: Value 
     :type  value: complex
 
  .. method:: SetPixelSampling(sampling)

     Sets the :ref:`pixel-sampling` of the image to the provided values in the 
     current :ref:`data-domain`. 

     :param sampling: Size of a pixel in x,y and z direction respectively
     :type  sampling: :class:`~ost.geom.Vec3`

  .. method:: SetReal(pixel, value)

     Sets the content of the specified pixel to the provided value.  If the 
     image holds data of the :obj:`COMPLEX` :ref:`data-domain`, the method sets 
     the pixel to a value has a real part equal to the provided value and a null 
     complex part.
 
     :param pixel: Image pixel
     :type  pixel: :class:`Point`
     :param value: Value
     :type  value: float

  .. method:: SetSpatialOrigin(pixel_indexes)

     Sets the :ref:`spatial-origin` of the image to the provided indexes.

     :param pixel_indexes: Indexes of the first pixel of the image
     :type  pixel_indexes: :class:`Point`
     :rtype: None

  .. method:: SetSpatialSampling(sampl)

     Sets the :ref:`pixel-sampling` if the image to the provided values in the 
     spatial :ref:`data-domain`

     :param sampl: Size of a pixel in x,y and z direction respectively
     :type  sampl: :class:`~ost.geom.Vec3`
 
  
