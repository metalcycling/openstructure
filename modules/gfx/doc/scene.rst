The Scene
================================================================================

.. currentmodule:: ost.gfx

The scene is the central registry for graphical objects and manages rendering 
parameters. Among other parameters, it is used to set up the lighting, fog, 
background color and the position and viewing direction of the user. The scene
is a singleton, meaning that there is only one scene available. The instance can
be accessed via :func:`gfx.Scene`. Because the scene is so important and
commonly used, the scene is also available as the `scene` variable in the 
interactive Python shell as well as from scripts.

Manipulating the Camera
--------------------------------------------------------------------------------

The users' position and orientation in the scene is specified by the camera. The camera is represented by a center, a rotation and an offset from the center. These 3 attributes can be set independently. Manipulation of the camera is done by assigning a new :class:`~ost.mol.Transform` to the :attr:`Scene.transform` attribute. 

.. _camera-orbit:

Orbiting around a Point
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following code example will let the camera orbit around the center of the camera.

.. code-block:: python

  # tweakable parameters:
  #  - axis of rotation
  #  - angular step size. Making it bigger speeds up the rotation
  #  - number of frames
  axis = geom.Vec3(0, 1, 0)
  step_size = 0.01
  num_frames = 100
  
  angle = 0.0
  for i in range(num_frames):
    angle += step_size
    if angle > math.pi*2:
      angle -= math.pi*2
    camera = scene.transform
    rot = geom.AxisRotation(axis, angle)
    camera.SetRot(rot)
    camera.SetTrans(geom.Vec3(25.0, 0.0, -50))
    scene.transform = camera
    scene.Export('frame-%03d.png' % i, 500, 500)

It is interesting to note that the offset from center (`trans`) is given in rotated coordinates. Transforming along z shifts the camera along the viewing direction. Setting x and y to non-zero causes the center of the camera not to be projected onto the center of the screen any longer. For example, setting the value `trans` to `geom.Vec3(50, 0, 0)` gives a viewing direction perpendicular to the vector from the camera position to the center.

.. class:: Scene

  .. attribute:: background

    The background color of the scene. By default the background color is
    transparent black.
  
    :type: :class:`Color`
  
  .. attribute:: center
  
    The center of the scene. Read-write. Changing the center affects the 
    orientation of the camera, but not its position. To change the camera
    position, use :attr:`transform`. Also available as 
    :meth:`GetCenter`/:meth:`SetCenter`.
    
    See also :meth:`CenterOn`.
    
    :type: :class:`~ost.geom.Vec3`
  
  .. attribute:: transform
  
    The camera transformation. Read-write. For an example usage see
    :ref:`camera-orbit`. Also available as
    :meth:`GetTransform`/:meth:`SetTransform`.
        
    :type: :class:`~ost.mol.Transform`
    
    
  .. attribute:: fov
  
    The field of view angle (in degrees) in vertical direction. Read-write.
    Also available as :meth:`GetFOV`/:meth:`SetFOV`.
    
    :type: float. 

  .. method:: Add(obj[, replace_existing])

    Add an object to the root node of the scene. This means that
    
     * The object will be rendered on the screen.    
     * The object's bounding box affects the slabbing operations.     

    :param obj: The object to be added.
    :type  obj: :class:`GfxNode`
    :param replace_existing: If true, existing objects of the same name are
       silently replaced. If false (the default), trying to add an object with
       an existing name will raise a RuntimeError.
       
    :type replace_existing: bool


  .. method:: AttachObserver(obs)

    Attach a :class:`scene observer <SceneObserver>`. The new scene observer 
    will get notified when the scene's state changes or objects get updated. 
    
    :param obs:
    :type  obs: :class:`SceneObserver`


  .. method:: AutoAutoslab(flag)

    Enable/disable autoslabbing. If autoslabbing is enabled, the near and far 
    clipping planes will automatically be adjusted to contain all the objects in 
    the scene. This calculation is done before every redraw. Note that, 
    autoslabbing is not rotation invariant. See :meth:`AutoslabMax` for a
    rotation-invariant version.
    
    :param flag:
    :type  flag: bool


  .. method:: Autoslab(fast[, force])

    .. note::
      
      This method looks stale. Remove it?

    :param fast:
    :type  fast: bool
    :param force:
    :type  force: bool

  .. method:: AutoslabMax()

    Adjust the near and far clipping planes in such a way that they contain the 
    bounding sphere of all the objects in the scene. In constrast to AutoSlab() 
    the near/far clipping planes calculated with :meth:`AutoslabMax` are 
    invariant to rotation.

  .. method:: CenterOn(obj)
    
    Center the camera on the center of the obj. This does not update offset and 
    rotation of the camera. However, since the offset and rotation are applied 
    after the centering, the position and viewing direction are affected by the 
    change of center.
    
    :param obj: Object, or name of the object
    :type  obj: str, or :class:`GfxNode`


  .. method:: Export(filename, width, height[, transparent])
              Export(filename[, transparent])
              
    Renders (exports) the scene into a PNG image.
    
    :param filename: The output filename
    :type  filename: str
    :param width: The width of the image. Defaults to the width of the 
      OpenGL window.
    :type  width: int
    :param height: The height of the image. Defaults to the height of the OpenGL
        window.
    :type  height: int
    :param transparent: If true, and the background color of the scene is 
      transparent, will produce a transparent image. If false, the alpha 
      component will be set to opaque.
    :type  transparent: bool

  .. method:: ExportPov(filename[, working_dir])
    
    Export the scene to POV-Ray format. The export will generate two files, one 
    containing a general description of the scene, including camera position and 
    materials and a second containing the geometric description of the objects.
    The first file will be named `filename.pov`, the second `filename.inc`.
    
    The files will be placed in working directory
    
    .. note::
    
      This method is highly experimental and does not work for more complex
      objects. Stay tuned for updates.
      

    :param filename: The base filename without file extension. `.pov` and `.inc`
       will automatically be appended to te filename.
    :type  filename: str
    :param working_dir: The working directory. Defaults to the current
       directory.
    :type  working_dir: str


  .. method:: GetFOV()

    Get the field of view angle in the y direction (in degrees).
    
    :rtype: float

  .. method:: GetRTC()
    
    .. note::
    
      Looks stale. Remove it?
      
    :rtype: :class:`~ost.geom.Mat4`

  .. method:: GetTransform()

    See :attr:`transform`

  .. method:: Remove(obj)

    Remove the given object from the scene.
    
    :param obj:
    :type  obj: :class:`GfxNode`
  
  .. method:: RemoveAll()
  
    Remove all objects from the scene.
    
    
  .. method:: RenderGL()

    Renders the scene.

  .. method:: RequestRedraw()

    Request a redraw of the scene. 


  .. method:: Resize(width, height)
    
    Resize the OpenGL context.
    
    :param width: The new width
    :type  width: int
    :param height: The new height
    :type  height: int

  .. method:: SetBackground(color)

    See :attr:`background`

  .. method:: SetBlur(arg2)

    :param arg2:
    :type  arg2: int

  .. method:: SetFOV(angle)

    See :attr:`fov`.
    
  .. method:: SetFog(arg2)

    :param arg2:
    :type  arg2: bool

  .. method:: SetFogColor(arg2)

    :param arg2:
    :type  arg2: :class:`Color`

  .. method:: SetFogOffsets(arg2, arg3)

    :param arg2:
    :type  arg2: float
    :param arg3:
    :type  arg3: float

  .. method:: SetLightDir(arg2)

    :param arg2:
    :type  arg2: :class:`~ost.geom.Vec3`

  .. method:: SetLightProp(arg2, arg3, arg4)

    :param arg2:
    :type  arg2: :class:`Color`
    :param arg3:
    :type  arg3: :class:`Color`
    :param arg4:
    :type  arg4: :class:`Color`

  .. method:: SetNearFar(near, far)

    Manually sets the near and far clipping planes to the given values. Before 
    using this method, make sure :meth:`AutoAutoslab` is disabled.
    
    :param near:
    :type  near: float
    :param far:
    :type  far: float

  .. method:: UnProject(point[, ignore_vp])

    Apply the inverse of the current camera and perspective transformation. This 
    essentially gives the screen (or normalized view coordinates) of a point in
    global coordinates. The inverse operations is available as  :meth:`Project`.
    
    :param point: The point in global coordinates.
    :type  point: :class:`~ost.geom.Vec3`
    :param ignore_vp: If set to false (the default), also performs the viewport 
      transformation. Points inside the viewing frustum will have x and y 
      coordinates in the range [0, width)x[0, height]. If true, the returned 
      coordinates will be between [-1,1]x[-1,1].
    :type  ignore_vp: bool
    :rtype: :class:`~ost.geom.Vec3`

  .. method:: SetRTC(arg2)

    :param arg2:
    :type  arg2: :class:`~ost.geom.Mat4`

  .. method:: SetShadow(arg2)

    :param arg2:
    :type  arg2: bool

  .. method:: SetShadowQuality(arg2)

    :param arg2:
    :type  arg2: int

  .. method:: SetStereoEye(arg2)

    :param arg2:
    :type  arg2: int

  .. method:: SetStereoInverted(arg2)

    :param arg2:
    :type  arg2: bool

  .. method:: SetTransform(new_transform)

    See :attr:`transform`

  .. method:: Stereo(mode)

    :param mode:
    :type  mode: int

  
