Scene Window
================================================================================
.. currentmodule:: ost.gui

The scene window holds a list of all graphical objects currently registered to
the :class:`gfx.Scene`.

.. image:: images/100208_Scene_Win.png

Every node of the graphical scene tree is shown as a subnode of the scene. It 
is possible to select one or more objects in the scenewin (`ctrl + <left mouse>`
). The available actions of the context menu are applied to all currently 
selected objects.

.. class:: SceneWin

  .. method:: GetContextMenu()
  
    Returns the :class:`ContextMenu`-instance.
     
    :rtype: ContextMenu

  .. method:: GetQObject()

    Get the SIP-QObject (QWidget), learn more about :doc:`python_cpp`.
    
    :rtype: QWidget

  .. method:: Hide()
  
    Hide the Widget

  .. method:: Show()

    Shows the Widget
    
Context Menu
--------------------------------------------------------------------------------

The Context Menu of the Scene Window is context sensitive. So, dependent on what
is selected in the :class:`SceneWin` the context menu changes.

Context menu of an entity:
  .. image:: images/100614_entity_context_menu.png
 
Context menu of an entity view:
  .. image:: images/100614_context_menu_view.png
  
It is possible to extend the Context Menu from python, by creating a QAction and 
add it to the ContextMenu-Class. If you want to display the action only for 
certain objects / states, you can pass flags which marks when the action should 
be shown.

  .. code-block:: python
  
    from PyQt5 import QtWidgets
    cm = gui.GostyApp.Instance().scene_win.GetContextMenu()
    action = QtWidgets.QAction("Test single entity", cm.qobject)
    cm.AddAction(action, gui.ContextActionType.ENTITY | gui.ContextActionType.SINGLE)

.. class:: ContextMenu
  
  The ContextMenu is shown whenever the user presses `<right click>` on the :class:`SceneWin` 
  
  .. method:: AddAction(action, flags)
  
     Adds the given action to the context menu of the scene window. 

     :param action: The Action which should be displayed in the context menu
    
     :type  action: :class:`QAction`
     
     :param flags: Flags that indicates, when the context_menu should be active. The action 
                   will be shown, when all flags are true.
    
     :type  flags: :data:`ContextActionTypes`


  .. method:: GetQObject()

    Get the SIP-QObject (QObject), learn more about :doc:`python_cpp`.
    
    :rtype: QObject

    
.. data:: ContextActionType
 
  It is possible to specify flags for new :class:`ContextMenu` actions of the :class:`SceneWin`
  These flags are used, to specify when the action should be shown on the scene window. The 
  actions can be combined by `oring` them (|). When flags are combined with an `or`, the 
  resulting object is a :data:`ContextActionTypes`-object.
  The following flags are defined by default: 

  * NOT_SCENE
     All selected objects are not the root node  
  * GFX_OBJECT
     All selected objects are :class:`gfx.GfxObj`
  * ENTITY
      All selected objects are :class:`gfx.Entity`
  * MAP
      All selected objects are :class:`gfx.MapIso`
  * ENTITY_VIEW
      All selected objects are parts of the same Entity (either Queries or Views)
  * CUSTOM_VIEW
      All selected objects are named views (created by the user)  
  * NOT_HIDDEN
      All selected objects are visible (can not be combined with the ENTITY_VIEW flag)
  * NOT_VISIBLE
      All selected objects are hidden (can not be combined with the ENTITY_VIEW flag)
  * SINGLE
      There is only one selected object in the :class:`SceneWin`
  * MULTI
      There are multiple objects selected object in the :class:`SceneWin`
   
.. data:: ContextActionTypes

  When multiple :data:`ContextActionType` objects are being combined with an or (|) you will get an 
  :data:`ContextActionTypes`-Object.
  
  .. code-block:: python
    
    ctx_act_types = gui.ContextActionType.ENTITY | gui.ContextActionType.SINGLE

  
  
