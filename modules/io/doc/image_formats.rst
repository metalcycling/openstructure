Supported Image File Formats
================================================================================

CCP4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by the ccp4 software package.

* Recognized file extensions: *.ccp4, .map*  
* Recognized based on content: *yes*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.CCP4(normalize_on_save=false, endianess_on_save=OST_LOCAL_ENDIAN)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool
  :param endianess_on_save: Byte order for saving.
  :type  endianess_on_save: ost.io.Endianess


DAT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Simple binary format for square images.

* Recognized file extensions: *.dat, .img*  
* Recognized based on content: *no*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.DAT(normalize_on_save=false, format=OST_DEFAULT_FORMAT, signed=false, endianess=OST_LOCAL_ENDIAN)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool
  :param format: Data format for loading and saving.
  :type  format: ost.io.Format
  :param signed: Determines whether data is interpreted as signed or unsigned values during loading and saving.
  :type  signed: bool
  :param endianess: Byte order for loading and saving.
  :type  endianess: ost.io.Endianess

DF3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PovRay Density file format.

* Recognized file extensions: *.df3*  
* Recognized based on content: *no*
* File import: *no*
* File export: *yes*

.. class:: ost.io.DF3(normalize_on_save=false)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool


DM3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by Gatan Inc.'s Digital Micrograph software.

* Recognized file extensions: *.dm3*  
* Recognized based on content: *yes*
* File import: *yes*
* File export: *no*

.. class:: ost.io.DM3( )

DX
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by the OpenDX software package.

* Recognized file extensions: *.dx*  
* Recognized based on content: *no*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.DX(normalize_on_save=false)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool


IPL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Ditabis Micron Image Plate Scanner Format.

* Recognized file extensions: *.ipl*  
* Recognized based on content: *yes*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.IPL(normalize_on_save=false, format=OST_DEFAULT_FORMAT)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool
  :param format: Data format for loading and saving.
  :type  format: ost.io.Format

JPK
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by JPK Instruments AG's software (Customized Tiff format).

* Recognized file extensions: *.jpk*  
* Recognized based on content: *no (jpk files without file extension are recognized as tif)*
* File import: *yes*
* File export: *yes*



.. class:: ost.io.JPK(normalize_on_save=indeterminate, format=OST_DEFAULT_FORMAT, signed=false, phasecolor=false, subimage=-1)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: tribool
  :param format: Data format for loading and saving.
  :type  format: ost.io.Format
  :param signed: Determines whether data is interpreted as signed or unsigned values during loading and saving.
  :type  signed: bool
  :param phasecolor: Determines whether the phases of complex values are encoded as color information.
  :type  phasecolor: bool
  :param subimage: Determines which subimage is read from a multi-image jpk (-1 = default image).
  :type  subimage: int

MRC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by the MRC software package.

* Recognized file extensions: *.mrc*  
* Recognized based on content: *yes*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.MRC(normalize_on_save=false, subformat=MRC_NEW_FORMAT, endianess_on_save=OST_LOCAL_ENDIAN)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool
  :param subformat: Determines whether old style or new style mrc files should be saved.
  :type  subformat: ost.io.Subformat
  :param endianess_on_save: Byte order for saving.
  :type  endianess_on_save: ost.io.Endianess
  
.. py:class:: ost.io.Subformat

   .. py:attribute:: MRC_AUTO_FORMAT
   .. py:attribute:: MRC_OLD_FORMAT
   .. py:attribute:: MRC_NEW_FORMAT

Nanoscope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by software from Veeco.

* Recognized file extensions: *.mod*  
* Recognized based on content: *no*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.Nanoscope(subimage=-1)
 
  :param subimage: Determines which subimage is read from a multi-image nanoscope file (-1 = default image).
  :type  subimage: int


PNG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Portable Network Graphic image format.

* Recognized file extensions: *.png*  
* Recognized based on content: *yes*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.PNG(normalize_on_save=false)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool

Situs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by the Situs software package.

* Recognized file extensions: *.situs, .sit*  
* Recognized based on content: *no*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.Situs(normalize_on_save=false)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool

Spider
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Format used by the Spider software package.

* Recognized file extensions: *.spi*  
* Recognized based on content: *no*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.Spider(normalize_on_save=false, endianess_on_save=OST_LOCAL_ENDIAN)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: bool
  :param endianess_on_save: Byte order for saving.
  :type  endianess_on_save: ost.io.Endianess

TIF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Tagged Image File Format.

* Recognized file extensions: *.tif, .tiff*  
* Recognized based on content: *yes*
* File import: *yes*
* File export: *yes*

.. class:: ost.io.TIF(normalize_on_save=indeterminate, format=OST_DEFAULT_FORMAT, signed=false, phasecolor=false, subimage=-1)
 
  :param normalize_on_save: Determines whether to renormalize image upon saving.
  :type  normalize_on_save: tribool
  :param format: Data format for loading and saving.
  :type  format: ost.io.Format
  :param signed: Determines whether data is interpreted as signed or unsigned values during loading and saving.
  :type  signed: bool
  :param phasecolor: Determines whether the phases of complex values are encoded as color information.
  :type  phasecolor: bool
  :param subimage: Determines which subimage is read from a multi-image tif (-1 = default image).
  :type  subimage: int

Common parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:class:: ost.io.Endianess

   .. py:attribute:: OST_BIG_ENDIAN
   .. py:attribute:: OST_LITTLE_ENDIAN
   .. py:attribute:: OST_VAX_DATA
   .. py:attribute:: OST_LOCAL_ENDIAN

.. py:class:: ost.io.Format

   .. py:attribute:: OST_BIT8_FORMAT
   .. py:attribute:: OST_BIT16_FORMAT
   .. py:attribute:: OST_BIT32_FORMAT
   .. py:attribute:: OST_FLOAT_FORMAT
   .. py:attribute:: OST_DOUBLE_FORMAT
   .. py:attribute:: OST_DEFAULT_FORMAT

