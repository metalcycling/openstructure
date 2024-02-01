import math,traceback

from PyQt5 import QtCore

import ost
from ost import gfx, gui

class SpacenavControl(QtCore.QObject):
  def __init__(self, spnav,
               parent=None):
    QtCore.QObject.__init__(self, parent)
    QtCore.QObject.connect(spnav,QtCore.SIGNAL("deviceTransformed(int,int,int,int,int,int)"), self.Changed)
    QtCore.QObject.connect(spnav,QtCore.SIGNAL("deviceButtonPressed(int)"), self.Toggle)     
    
    self.trans = True
    self.rot = True
    self.speed = 480.0
    
  def Changed(self, tx,ty,tz,rx,ry,rz):
    scene=gfx.Scene()
    tf = scene.transform
    def d(r):
      if r==0.0:
        return 0.0
      rr=r/abs(r)*max(0.0,abs(r)-0.9)
      if abs(rr)>0:
        return rr/abs(rr)*(math.pow(1.01,abs(rr))-1.0)*40.0/self.speed
      return 0.0
    if(self.trans):
      tf.ApplyXAxisTranslation(d(tx))
      tf.ApplyYAxisTranslation(d(ty))
      # adjust translation speed to distance from viewpoint
      currz=tf.trans[2];
      delta=currz*math.pow(1.01,d(tz))-currz;
      tf.ApplyZAxisTranslation(delta);
      # adjust near and far clipping planes together with z trans
      scene.SetNearFar(scene.near-delta,scene.far-delta);
    if(self.rot):
      tf.ApplyXAxisRotation(d(rx))
      tf.ApplyYAxisRotation(d(ry))
      tf.ApplyZAxisRotation(d(rz))

    scene.transform=tf
    scene.RequestRedraw()

  def Toggle(self, button):
    if button == 0:
      self.trans = not self.trans
      ost.LogVerbose("SpaceNav: translation %s" % (self.strans and "enabled" or "disabled"))
    elif button == 1:
      self.rot = not self.rot
      #ost.LogVerbose("SpaceNav: rotation %s"%("enabled" if self.rot else "disabled"))
      ost.LogVerbose("SpaceNav: rotation %s"%(self.rot and "enabled" or "disabled"))
    elif button == 12:
      if self.speed > 20:
        self.speed *= 0.8
        ost.LogVerbose("SpaceNav: speed increased to "+str(self.speed))
    elif button == 13:
      self.speed /= 0.8
      ost.LogVerbose("SpaceNav: speed reduced to "+str(self.speed))
    else:
      ost.LogVerbose("SpaceNav: unmapped button press ["+str(button)+"]")

      
def _InitSpaceNav(app):
  try:
    spnav = gui.SpnavInput.GetQThread()
    if spnav:
      spnav.start()
      parent = app.gl_win.qobject
      SpacenavControl(spnav,parent)
      ost.LogInfo("SpaceNav: device found and connected")
    else:
      ost.LogInfo("SpaceNav: no device found, or could not connect to device socket")
  except:
    ost.LogInfo("SpaceNav: caught exception during initialization: %s"%(traceback.format_exc()))
