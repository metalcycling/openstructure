The Sequence Viewer
================================================================================
.. currentmodule:: ost.gui

The sequence viewer of OpenStructure can be used to display sequences or
alignments.

The following example adds an alignment to a sequence viewer and shows it in a 
new sequence viewer instance:

.. code-block:: python
  
  aln=io.LoadAlignment('sh2.aln')
  v=gui.SequenceViewer()
  v.AddAlignment(aln)
  v.Show()

Core Features
--------------------------------------------------------------------------------

Display styles
^^^^^^^^^^^^^^
It is possible to display the data in different ways. To change the display 
style, click the tool icon of the menubar and select the 
display style from the list.

  .. image:: images/100615_select_display_style.png

For the moment there are four different display styles:



  **1. Highlight propteries**
    
    .. image:: images/100615_property_type.png
       :width: 500pt
    
    Color groups of amino acids
  
    =====  ======  ======  =====  ===  ====  ====
    grey   orange  yellow  green  red  blue  cyan
    =====  ======  ======  =====  ===  ====  ====
    G      F       C       S      K    D     P
    A      Y       M       T      R    E      
    V      W                      H    N      
    L                                  Q      
    I                                         
    =====  ======  ======  =====  ===  ====  ====

  
  **2. Secondary structure**
  
    .. image:: images/100615_secondary_structure.png
       :width: 500pt
      
    Display secondary structure information
    
    This mode is only available, when a :class:`mol.EntityView` is connected to 
    the sequence.
   
  **3. Highlight conservation 1**
  
    .. image:: images/100615_conservation_1.png
       :width: 500pt
  
    Display conservation
     
    This display style uses the :meth:`seq.alg.Conservation`-Algorithm which calculates
    the Conservation for each column. This Mode is only available for alignments.
  
  **4. Highlight conservation 2**
  
    .. image:: images/100615_conservation_2.png
       :width: 500pt
      
    Display conservation
    
    This implementation of conservation colors the columns grey, if all amino acids are 
    from the same type. If all amino acids are from the same group, it colors the column with 
    a light grey. Otherwise the column is colored white. This Mode is also only available for
    alignments.


Zoom
^^^^

With the sequence viewer of Openstructure it is possible to zoom in and out `<ctrl + mousewheel>`.

 .. image:: images/100616_zoom_in.png
 
 Zoom in, if you lost your glasses at home 
 
 .. image:: images/100616_zoom_out.png
 
 Zoom out, if you need a better overview of the sequence 
 
Copy and Paste
^^^^^^^^^^^^^^

It is possible to copy sequence text from the sequence viewer by pressing `<ctrl + c>`. If you paste 
the text (for example into a text editor), all the gaps will be filled with `-`.

 .. image:: images/100616_copy_paste.png

Search
^^^^^^

With `<ctrl + f>` you can show or hide the search bar which help you finding parts of a sequence.

 .. image:: images/100617_search_bar.png

.. class:: SequenceViewer
  
The SequenceViewer class represents a graphical sequence viewer. To see the sequences of a :class:`gfx.Entity`, you
can use the :meth:`AddEntity` method. If you want to load an alignment, the :meth:`AddAlignment` method can be used.

  .. method:: AddEntity(entity)
        
    Add :class:`gfx.Entity` to SequenceViewer.  
    
    :param entity: an entity instance 
    :type  entity: :class:`gfx.Entity`

  .. method:: RemoveEntity(entity)
  
    Remove entity from sequence viewer.
    
    :param entity: an entity instance
    :type  entity: :class:`gfx.Entity`

  .. method:: AddAlignment(alignment)
        
    Add :class:`seq.AlignmentHandle` to SequenceViewer.  
    
    :param alignment: an alignment instance
    :type  alignment: :class:`seq.AlignmentHandle`

  .. method:: RemoveAlignment(alignment)
  
    Remove alignment from sequence viewer.
    
    :param alignment: an alignment instance
    :type  alignment: :class:`seq.AlignmentHandle`
    :rtype: None

  .. method:: ChangeDisplayMode(style [, obj : None])
    
    Change the display style of all alignments and sequences to the given style. If the obj parameter is given, it 
    changes only the display style of the given object. 

    If the style is not known, nothing happens
    
    See also :meth:`GetDisplayModes`, :meth:`GetCurrentDisplayMode`
    
    :param style: The display style
    :type  style: :class:`string`
    :param obj: The object you would like to change
    :type obj: :class:`gfx.Entity`, :class:`seq.AlignmentHandle`
    

  .. method:: GetCurrentDisplayMode([obj : None])
    
    Returns the current display style. If there are objects with different display styles, " " will be returned 
    otherwise the string of the display style.

    :rtype: str
    :param obj: The object you would like to get the display style
    :type obj: :class:`gfx.Entity`, :class:`seq.AlignmentHandle`

  .. method:: GetDisplayModes([obj : None])

    Returns a list with all available display styles. If you pass an object, you will get only the display styles for 
    the given object.

    :rtype: StringList
    :param obj: The object you would like to get all display styles
    :type obj: :class:`gfx.Entity`, :class:`seq.AlignmentHandle`

  .. method:: GetQObject()
    
    Get the SIP-QObject (QWidget), learn more about :doc:`python_cpp`.
    
    :rtype: QWidget

  .. method:: Hide()
  
    Hide the Widget

  .. method:: Show()

    Shows the Widget

