:mod:`img.alg <ost.img.alg>` - Image Processing Algorithms
================================================================================

.. module:: ost.img.alg
  :synopsis: Image processing algorithms
  
Applying algorithms
-------------------

While image properties are usually manipulated using method of the 
:class:`~ost.img.ImageHandle` class, their data content is manipulated using 
image algorithms. Image algorithms are objects. Each of them is a class, and its 
methods are used to handle the algorithm parameters. Applying an algorithm to an 
image is then conceptually a two-step process. First, an instance of an 
algorithm class is created, yielding an algorithm object. In a second step, the 
algorithm object is applied to an image. An algorithm can be applied in-place 
(using the :meth:`~ost.img.ImageHandle.ApplyIP` method), modifying the image, or 
out-of-place, (using :meth:`~ost.img.ImageHandle.Apply` ), leaving the original image 
untouched, and returning the result as a new image. 

Here is an example. All the algorithms used in the following are described in the :ref:`algorithms` section.

.. code-block:: python
 
  # creates an algorithm object
  rand_alg = img.alg.Randomize() 
  # applies algorithm object in place, overwriting the image
  im.ApplyIP( rand_alg )

Sometimes, there is no need to create a permanent instance of an algorithm 
object. A temporary object enough:

.. code-block:: python

  # applies temporary algorithm object in-place
  im.ApplyIP( img.alg.GaussianFilter(4.0) )

When used this way, the algorithm class will cease to exist as soon as the 
algorithm is applied. However, some algorithm are stateful and store 
information. One good example is the :class:`Stat` algorithm, which does not
modify the image when applied, but change its internal state to store 
information extracted from the image, which can be recovered later. For example:

.. code-block:: python

  # creates and applies an algorithm object
  stat=img.alg.Stat()
  im.ApplyIP(stat)
  # extracts information from the algorithm
  mean=stat.GetMean()

It is important to remember that when the algorithms ceases to exist, all information it stores is lost.

Fourier Transforming Images
----------------------------

An image can be Fourier-transformed using either the :class:`FFT` algorithm or 
the :class:`DFT` algorithm. The difference between the two is that the 
:class:`DFT` algorithm honors the :ref:`spatial-origin` of the image, and 
applies the corresponding phase shift in Fourier space. The :class:`FFT` does 
not follow this behavior. 

.. code-block:: python

  # create an instance of the Dft algorithm object
  dft=img.alg.DFT() 
  # do the actual Fourier transformation
  im_ft=im.Apply(dft) 
  # back-transform
  im2 = im_ft.Apply(dft) 

The :class:`FFT` and :class:`DFT` algorithms do not require a direction to be 
given (forward or back transform). This is implicitly determined by the current 
:ref:`data-domain` of the image being transformed. The following rules apply. 

* :obj:`SPATIAL` -> :obj:`HALF_FREQUENCY`
* :obj:`HALF_FREQUENCY` -> :obj:`SPATIAL`
* :obj:`FREQUENCY` -> :obj:`COMPLEX_SPATIAL`
* :obj:`COMPLEX_SPATIAL` -> :obj:`FREQUENCY`

.. _filters:

Filters
-------

OpenStructure makes several image filters available. Most of them are Fourier 
space filters, others are real space ones. However, since the 
:class:`~ost.img.ImagerHandle` class is aware of its own :ref:`data-domain`, 
the user does not need to convert the image to Fourier space or to real space. 
Irrespective of which domain the filter applies to, OpenStructure  will 
internally convert the image to the appropriate domain, apply the filter, and 
then return the image to its original conditions.

The following filters are available (their are described in the :ref:`algorithms` section below)

Fourier space filters:

* :class:`LowPassFilter`
* :class:`HighPassFilter`
* :class:`GaussianLowPassFilter`
* :class:`GaussianHighPassFilter`
* :class:`FermiLowPassFilter`
* :class:`FermiHighPassFilter`
* :class:`ButterworthLowPassFilter`
* :class:`ButterworthHighPassFilter`
* :class:`FermiLowPassFilter`

Real space filters:

* :class:`GaussianFilter`

.. _algorithms:

Selected Algorithms
--------------------------------------------------------------------------------

Many algorithms are available for image manipulation. What follows is a description of the 
most important ones.

.. class:: DFT()

   This algorithm performs a Fourier Transform of the image, honoring its 
   :ref:`spatial-origin`, thus applying the corresponding phase shift in Fourier 
   space.

.. class:: DiscreteShrink(block_size)

   The algorithm performs a scaling of the original image by merging adjacent 
   blocks of pixels. The block size is passed in the constructor in the form of 
   a :class:`~ost.img.Size` but can be changed later using the relevant method. 
   The :class:`~ost.img.Size` and the :class:`~ost.img.Extent` of the image are 
   changed when the algorithm is applied. The :ref:`pixel-sampling` of the image 
   is also adjusted according to the scaling, so that the size of the image in 
   the absolute reference system used by OpenStructure stays constant.
   
   :param block_size: Size of the blocks to be merged
   :type block_size: :class:`~ost.img.Size`

   .. method:: GetBlocksize()

     Returns the current size of the blocks to be merged

     :rtype: :class:`~ost.img.Size`

   .. method:: SetBlocksize(block size)

    Sets the size of the blocks to be shrunk to the specified value

    :param block_size:
    :type  block_size: :class:`~ost.img.Size`
 
.. class:: FFT()

    This algorithm performs a Fourier Transform of the image, without honoring 
    its :ref:`spatial-origin` (See :class:`DFT`)
	
.. class:: LowPassFilter(cutoff=1.0)

   This algorithm applies a Fourier low pass filter to the image. The filter cutoff frequency needs
   to be provided in sampling units (for example 8 Angstrom). Please notice that this filter features a sharp dropoff.

   :param cutoff: Frequency cutoff in sampling units
   :type  cutoff: float

   .. method:: GetLimit()
     
     Returns the current value of the filter cutoff frequency  (in sampling units).

     :rtype: float

   .. method:: SetLimit(cutoff)

     Sets the value of the filter cutoff frequency to the specified value (in sampling units).

     :param cutoff: Frequency cutoff in sampling units
     :type  cutoff: float

.. class:: HighPassFilter(cutoff=1.0)

   This algorithm applies a Fourier high pass filter to the image. The filter cutoff frequency needs
   to be provided in sampling units (for example 8 Angstrom). Please notice that this filter features a sharp dropoff.

   :param cutoff: Frequency cutoff in sampling units
   :type  cutoff: float

   .. method:: GetLimit()
     
     Returns the current value of the filter cutoff frequency  (in sampling units).

     :rtype: float

   .. method:: SetLimit(cutoff)

     Sets the value of the filter cutoff frequency to the specified value (in sampling units).

     :param cutoff: Frequency cutoff in sampling units
     :type  cutoff: float


.. class:: GaussianLowPassFilter(cutoff=1.0)

   This algorithm applies a Fourier `Gaussian low pass filter <http://en.wikipedia.org/wiki/Gaussian_filter>`_ to the
   image. The filter cutoff frequency needs to be provided in sampling units (for example 8 Angstrom). 

   :param cutoff: Frequency cutoff in sampling units
   :type  cutoff: float

   .. method:: GetLimit()
  
     Returns the current value of the filter cutoff frequency (in sampling units).

     :rtype: float

   .. method:: SetLimit(cutoff)

	 Sets the value of the filter cutoff frequency to the specified value (in sampling units).

	 :param cutoff: Frequency cutoff in sampling units
	 :type  cutoff: float

.. class:: GaussianHighPassFilter(cutoff=1.0)

   This algorithm applies a Fourier `Gaussian High pass filter <http://en.wikipedia.org/wiki/Gaussian_filter>`_ to the
   image. The filter cutoff frequency needs to be provided in sampling units (for example 8 Angstrom). 

   :param cutoff: Frequency cutoff in sampling units
   :type  cutoff: float

   .. method:: GetLimit()
  
     Returns the current value of the filter cutoff frequency (in sampling units).

     :rtype: float

   .. method:: SetLimit(cutoff)

	 Sets the value of the filter cutoff frequency to the specified value (in sampling units).

	 :param cutoff: Frequency cutoff in sampling units
	 :type  cutoff: float
	
.. class:: FermiLowPassFilter(cutoff=1.0,t=1.0)

   This algorithm applies a Fourier `Fermi low pass filter <http://en.wikipedia.org/wiki/Fermi_filter>`_ to the
   image. The filter cutoff frequency and the temperature parameter T need to be provided in sampling units 
   (for example 8 Angstrom). 

   :param cutoff: Frequency cutoff in sampling units
   :type  cutoff: float
   :param t: Temperature factor in sampling units
   :type  t: float

   .. method:: GetLimit()
  
     Returns the current value of the filter cutoff frequency in sampling units.

     :rtype: float

   .. method:: SetLimit(cutoff)

	 Sets the value of the filter cutoff frequency to the specified value (in sampling units).

	 :param cutoff: Frequency cutoff in sampling units
	 :type  cutoff: float
	
   .. method:: GetT()
  
     Returns the current value of the filter's T factor (in sampling units).

     :rtype: float

   .. method:: SetT(t_factor)

	 Sets the value of the filter's T factor to the specified value (in sampling units).

	 :param t_factor: Frequency cutoff in sampling units
	 :type  t_factor: float

.. class:: FermiHighPassFilter(cutoff=1.0,t=1.0)

   This algorithm applies a Fourier `Fermi high pass filter <http://en.wikipedia.org/wiki/Fermi_filter>`_ to the
   image. The filter cutoff frequency and the temperature parameter T need to be provided in sampling units 
   (for example 8 Angstrom). 

   :param cutoff: Frequency cutoff in sampling units
   :type  cutoff: float
   :param t: Temperature factor in sampling units
   :type  t: float

   .. method:: GetLimit()
  
     Returns the current value of the filter cutoff frequency in sampling units.

     :rtype: float

   .. method:: SetLimit(cutoff)

	 Sets the value of the filter cutoff frequency to the specified value (in sampling units).

	 :param cutoff: Frequency cutoff in sampling units
	 :type  cutoff: float
	
   .. method:: GetT()
  
     Returns the current value of the filter's T factor (in sampling units).

     :rtype: float

   .. method:: SetT(t_factor)

	 Sets the value of the filter's T factor to the specified value (in sampling units).

	 :param t_factor: Frequency cutoff in sampling units
	 :type  t_factor: float
	
.. class:: ButterworthLowPassFilter(passband=1.0,stopband=1.0)

   This algorithm applies a Fourier `Butterworth low pass filter <http://en.wikipedia.org/wiki/Butterworth_filter>`_ to
   the image. The filter passband and stopband frequencies need to be provided in sampling units (for example 8 Angstrom). 
   The default values of the Epsilon and Maximum Passband Gain parameters are set to 0.882 and 10.624 respectively.

   :param passband: Passband frequency in sampling units
   :type  passband: float
   :param stopband: Stopband frequency in sampling units
   :type  stopband: float

   .. method:: GetLimit()
  
     Returns the current value of the filter passband frequency in sampling units.

     :rtype: float

   .. method:: SetLimit(passband)

	 Sets the value of the filter passband frequency to the specified value (in sampling units).

	 :param passband: Frequency cutoff in sampling units
	 :type  passband: float
	
   .. method:: GetStop()
  
     Returns the current value of the filter's stopband frequency (in sampling units).

     :rtype: float

   .. method:: SetStop(stopband)

	 Sets the value of the filter's stopband frequency to the specified value (in sampling units).

	 :param stopband: Frequency cutoff in sampling units
	 :type  stopband: float	
	
   .. method:: GetEps()
  
     Returns the current value of the filter's Epsilon parameter.

     :rtype: float

   .. method:: SetEps(epsilon)

	 Sets the value of the filter's epsilon parameter to the specified value.

	 :param eps: Epsilon parameter
	 :type  eps: float
	
   .. method:: GetA()
  
     Returns the current value of the filter's Maximum Passband Gain parameter.

     :rtype: float

   .. method:: SetA(gain)

	 Sets the value of the filter's Maximum Passband Gain parameter to the specified value.

	 :param gain: Maximum Passband Gain parameter
	 :type  gain: float			
	
.. class:: ButterworthHighPassFilter(passband=1.0,stopband=1.0)

   This algorithm applies a Fourier `Butterworth high pass filter <http://en.wikipedia.org/wiki/Butterworth_filter>`_ 
   to the image. The filter passband and stopband frequencies need to be provided in sampling units (for example 8
   Angstrom). The default values of the Epsilon and Maximum Passband Gain parameters are set to 0.882 and 10.624
   respectively.

   :param passband: Passband frequency in sampling units
   :type  passband: float
   :param stopband: Stopband frequency in sampling units
   :type  stopband: float

   .. method:: GetLimit()
  
     Returns the current value of the filter passband frequency in sampling units.

     :rtype: float

   .. method:: SetLimit(passband)

	 Sets the value of the filter passband frequency to the specified value (in sampling units).

	 :param passband: Frequency cutoff in sampling units
	 :type  passband: float
	
   .. method:: GetStop()
  
     Returns the current value of the filter's stopband frequency (in sampling units).

     :rtype: float

   .. method:: SetStop(stopband)

	 Sets the value of the filter's stopband frequency to the specified value (in sampling units).

	 :param stopband: Frequency cutoff in sampling units
	 :type  stopband: float	
	
   .. method:: GetEps()
  
     Returns the current value of the filter's Epsilon parameter.

     :rtype: float

   .. method:: SetEps(epsilon)

	 Sets the value of the filter's epsilon parameter to the specified value.

	 :param eps: Epsilon parameter
	 :type  eps: float
	
   .. method:: GetA()
   
     Returns the current value of the filter's Maximum Passband Gain parameter.

     :rtype: float

   .. method:: SetA(gain)

	 Sets the value of the filter's Maximum Passband Gain parameter to the specified value.

	 :param gain: Maximum Passband Gain parameter
	 :type  gain: float	
	
.. class:: GaussianFilter(sigma=1.0)

	 This algorithm applies a real space Gaussian filter to the image, as defined in the following publication:
	
	 I.T.Young, L.J. van Vliet,"Recursive implementation of the Gaussian filter",Signal Processing, 44(1995), 139-151
	
	 :param sigma: Width of the Gaussian filter
	 :type  sigma: float

	 .. method:: GetSigma()

	   Returns the current value of the filter's width.

	   :rtype: float

	 .. method:: SetSigma(width)

	   Sets the value of the filter's width to the specified value.

	   :param sigma: Width of the Gaussian filter
	   :type  sigma: float			
	
	 .. method:: SetQ(q_param)

	   Sets the value of the filter's Q parameter (see publication) to the specified value.

	   :param q_param: Filter's Q parameter
	   :type  q_param: float			

.. class:: Histogram(bins, minimum, maximum)

   This algorithm performs an histogram analysis of the image. The minimum and 
   maximum pixel values of the histogram representation must be provided when 
   the algorithm object is created, as well as the number of bins in the 
   histogram. Bins are equally spaced and minimum and maximum values for each 
   bin are automatically computed.
   
   When the algorithm is applied to an image, the analysis is carried out. A 
   python 'list' object containing in sequence the pixel counts for all the bins 
   can the be recovered from the algorithm object.

   :param bins: Number of bins in the histogram
   :type  bins: int
   :param minimum: Minimum value in the histogram
   :type  minimum: float
   :param maximum: Maximum value in the histogram

   .. method:: GetBins()

     Returns the bins of the histogram representation

     :rtype: list of ints 

   :type  maximum: float

   .. method:: GetBins()

