Introduction to the :mod:`~ost.img` Module
================================================================================

For the course of this tutorial, we assume that you have :ref:`DNG up and running<start-dng>`.


Loading Images and Density Maps
--------------------------------------------------------------------------------

OpenStructure features an :mod:`~ost.img` module that is dedicated to the 
manipulation of images/density maps. The images or density maps can either be 
one-, two- or three-dimensional. The most common formats used in X-ray/ 
electron crystallography and atomic force microscopy are supported in addition 
to several general purpose image formats. See
:doc:`supported file formats <io/image_formats>` for details. The :mod:`~ost.img`
module was originally developed as part of the Image Processing Library &
Toolbox IPLT. More documentation and examples can also be found on the
`IPLT website <http://www.iplt.org>`.

To load a density map, type

.. code-block:: python

  map = io.LoadImage('/path/to/examples/code_fragments/map/1ppt.map')

This will load the fragment density map from the specified file 'fragment.map' 
and store the result in map. 

Now let's inspect what we just loaded:

.. code-block:: python

  print(map.GetPixelSampling(), map.GetSize())
    
We can see that the sampling is set to 1.0 Angstrom in all three dimensions. The loaded map is an instance of :class:`~ost.img.ImageHandle`, a class to represent images in one, two and three dimensions.

Manipulating Images and Density Maps
--------------------------------------------------------------------------------

The algorithms used for manipulation of an image are found in the 
:mod:`img.alg <ost.img.alg>` module. Therefore before using an algorithm we 
first have to import the :mod:`img.alg <ost.img.alg>` module.

.. code-block:: python

  from ost.img import alg


The :mod:`img.alg <ost.img.alg>` module provides a wide range of algorithms to 
manipulate image data. Here for example we use a LowPassFilter to restrict the 
resolution of the density map to frequencies lower than a certain threshold.

.. code-block:: python

   map_filtered = map.Apply(alg.LowPassFilter(3.0))

The filtered map is stored in a new variable called `map_filtered`. A complete list of algorithms is available on the  :doc:`img/alg/alg` page.


Displaying Images and Density Maps
--------------------------------------------------------------------------------

Now that we have a filtered map it's time to have a look at it. There are 
fundamentally two ways to visualize 3-dimensional density maps. One is by 
drawing isocontour surfaces. These are conceptually similar to contour lines 
used in cartography: every point on an isosurface has the same density value. 
Isosurfaces are easy to create in OpenStructure:

.. code-block:: python

   go = gfx.MapIso("filtered", map_filtered,0.5)
   scene.Add(go)

The other way to visualize a 3-dimensional map is by showing one 2-dimensional 
density slice at a time, allowing the user to move through the slices. In 
OpenStructure this is achieved using a DataViewer). 
A DataViewer showing the filtered map is created using the following command:

.. code-block:: python

  gui.CreateDataViewer(map_filtered)

This command displays a panel showing one slice of the density map lying on a 
particular (x,y) plane in the coordinate reference system. The `z` and `x`
keys can be used to move to slices lying at a lower or higher coordinate along 
the `z` axis, allowing the examination of the full 3-dimensional volume.

A more detailed explanation of the :mod:`~ost.img` module can be found in the 
tutorial section for :mod:`~ost.img`.


Continue with :doc:`part three <intro-03>` of the tutorial.
