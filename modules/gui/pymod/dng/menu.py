from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from ost import *
from ost import gui
from ost.gui.init_splash import _InitSplash
from ost.gui.dng import termuse
from ost.gui.dng import superpositiondialog

import sys
class FileMenu(QMenu):
  def __init__(self, parent=None):
    QMenu.__init__(self, parent)
    self.setTitle('File')
    gui.AddMenuAction(self, 'Open', self._Open, shortcut='Ctrl+O')
    gui.AddMenuAction(self, 'Save As...', self._SaveAs, 
                  shortcut='Ctrl+Shift+S', 
                  enabled=gui.OneOf(gfx.Entity))
  def _Open(self):
    filename,_=QFileDialog.getOpenFileName(None, 'Open file','')
    if(QFileInfo(filename).isFile()):
      gui.FileLoader.LoadObject(str(filename))

  def _SaveAs(self):
    sel_node=gui.SceneSelection.Instance().GetActiveNode(0)
    filename=QFileDialog.getSaveFileName(None, 'Save As', 
                                         '%s.pdb' % sel_node.name)
    io.SavePDB(sel_node.view, str(filename))
    
class ClipWidget(QWidget):
  def __init__(self, width=500, height=500):
    QWidget.__init__(self)
    self.setWindowFlags(Qt.Tool)
    l = QGridLayout(self)
    l.addWidget(QLabel("Near"), 0, 0)
    l.addWidget(QLabel("Far"), 1, 0)
    bounds_near = QLineEdit(str(0))
    bounds_near.setValidator(QIntValidator(0, 9999, bounds_near))
    bounds_near.setMaxLength(4)
    bounds_near.setMaximumWidth(50)
    bounds_far = QLineEdit(str(200))
    bounds_far.setValidator(QIntValidator(0, 9999, bounds_far))
    bounds_far.setMaxLength(4)
    bounds_far.setMaximumWidth(50)
    l.addWidget(bounds_near, 0, 1, 2, 1)
    l.addWidget(bounds_far, 0, 3, 2, 1)
    self.near_ = QSlider(Qt.Horizontal)
    self.near_.setMinimum(int(bounds_near.text()))
    self.near_.setMaximum(int(bounds_far.text()))
    self.near_.setValue(int(gfx.Scene().near))
    self.far_ = QSlider(Qt.Horizontal)
    self.far_.setMinimum(int(bounds_near.text()))
    self.far_.setMaximum(int(bounds_far.text()))
    far = int(gfx.Scene().far)
    if far>sys.maxsize:
      far = sys.maxsize
    self.far_.setValue(far)
    self.auto_ = QCheckBox("Continuous Automatic Clipping")
    self.auto_.setChecked(gfx.Scene().GetAutoAutoslab())

    l.addWidget(self.near_, 0, 2)
    l.addWidget(self.far_, 1, 2)
    l.addWidget(self.auto_, 2, 0, 1, 4)
    self.near_.valueChanged.connect(self.SetNear)
    self.far_.valueChanged.connect(self.SetFar)
    self.auto_.stateChanged.connect(self.SetAuto)
    bounds_near.textEdited.connect(self.SetNearBounds)
    bounds_far.textEdited.connect(self.SetFarBounds)

  def SetNear(self, val):
    gfx.Scene().near = val

  def SetFar(self, val):
    gfx.Scene().far = val

  def SetAuto(self, val):
    gfx.Scene().AutoAutoslab(val)
    gfx.Scene().near = int(self.near_.value())
    gfx.Scene().far = int(self.far_.value())

  def SetNearBounds(self, text):
    if text!='':
      self.near_.setMinimum(int(text))
      self.far_.setMinimum(int(text))

  def SetFarBounds(self, text):
    if text!='':
      self.near_.setMaximum(int(text))
      self.far_.setMaximum(int(text))

class ExportSceneDialog(QDialog):
  def __init__(self, width=500, height=500):
    QDialog.__init__(self)
    l=QGridLayout(self)
    l.setColumnMinimumWidth(0, 100)
    l.addWidget(QLabel("Width (px)"), 0, 0)
    l.addWidget(QLabel("Height (px)"), 1, 0)
    self.width_=QLineEdit(str(width))
    self.height_=QLineEdit((str(height)))
    self.width_.setValidator(QIntValidator(50, 3000, self.width_))
    self.height_.setValidator(QIntValidator(50, 3000, self.height_))
    self.opaque_=QCheckBox("Force Opaque Background")
    l.addWidget(self.width_, 0, 1)
    l.addWidget(self.height_, 1, 1)
    l.addWidget(self.opaque_, 2, 1)
    hbox=QHBoxLayout()
    cancel=QPushButton("Cancel")
    cancel.clicked.connect(self.reject)
    hbox.addWidget(cancel)
    export=QPushButton("Export")
    hbox.addWidget(export)
    export.setDefault(True)
    l.addLayout(hbox, 3, 1, 2, 1)
    export.clicked.connect(self.accept)
  @property
  def transparent(self):
    return not self.opaque_.isChecked()

  @property
  def width(self):
    return int(self.width_.text())

  @property
  def height(self):
    return int(self.height_.text())

class SceneMenu(QMenu):
  def __init__(self, parent=None):
    QMenu.__init__(self, parent)
    self.setTitle('Scene')
    self.aboutToShow.connect(self._AboutToShow)
    gui.AddMenuAction(self, 'Background Color', self._SetSceneBackground)
    self.fog_action=gui.AddMenuAction(self, 'Depth Cueing', self._ToggleFog, 
                                  shortcut='Ctrl+Shift+F', checkable=True, 
                                  checked=gfx.Scene().fog)
    gui.AddMenuAction(self, 'Center', self._Center,
                  enabled=gui.ManyOf(gfx.GfxObj))
    gui.AddMenuAction(self, 'Fit To Screen', self._FitToScreen,
                  enabled=gui.OneOf(gfx.Entity))
    gui.AddMenuAction(self, 'Superpose', self._SuperposeDialog,
                      enabled=gui.TwoOf(gfx.Entity))
    gui.AddMenuAction(self, 'Save Snapshot', self._ExportScene)
    gui.AddMenuAction(self, 'Scene Clipping', self._ClipScene, shortcut='Ctrl+Shift+C')
    
  def _ExportScene(self):
    qd=ExportSceneDialog()
    if not qd.exec_():
      return

    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    filename, _ = QFileDialog.getSaveFileName(self,"Save Snapshot","snapshot.png", "", options=options)

    if filename:
      gfx.Scene().Export(str(filename), qd.width, qd.height, qd.transparent)

  def _ClipScene(self):
    self.cw = ClipWidget()
    self.cw.show()

  def _AboutToShow(self):
    self.fog_action.setChecked(gfx.Scene().fog)
    
  def _Center(self):
    sel=gui.SceneSelection.Instance()
    gfx.Scene().center=sel.GetActiveNode(0).center

  def _SetSceneBackground(self):
    new_color=gui.PickColor(gfx.Scene().bg)
    if new_color:
      gfx.Scene().bg=new_color

  def _ToggleFog(self):
    gfx.Scene().fog=not gfx.Scene().fog
    self.fog_action.setChecked(gfx.Scene().fog)

  def _FitToScreen(self):
    sel=gui.SceneSelection.Instance()
    gfx.FitToScreen(sel.GetActiveNode(0))

  def _SuperposeDialog(self):
    sel=gui.SceneSelection.Instance()
    act_count=sel.GetActiveNodeCount()
    # we now that there have to be 2 gfx.Entities, because of using TwoOf to
    # enable menu entry!
    i = 0;
    gfx_ent_1 = sel.GetActiveNode(i)
    while not isinstance(gfx_ent_1, gfx.Entity):
      i += 1
      gfx_ent_1 = sel.GetActiveNode(i)
    i += 1
    gfx_ent_2 = sel.GetActiveNode(i)
    while not isinstance(gfx_ent_2, gfx.Entity):
      i += 1
      gfx_ent_2 = sel.GetActiveNode(i)
    sd = superpositiondialog.SuperpositionDialog(gfx_ent_1, gfx_ent_2)
    if sd.rmsd_superposed_atoms != None:
      if sd.reference == 0:
        gfx_ent_2.UpdatePositions()
        gfx.Scene().CenterOn(gfx_ent_1)
      else:
        gfx_ent_1.UpdatePositions()
        gfx.Scene().CenterOn(gfx_ent_2)
      LogScript('RMSD: %.3f'%sd.rmsd)  
      LogScript('RMSD superposed atoms: %.3f'%sd.rmsd_superposed_atoms) 
      LogScript('fraction superposed: %.3f'%sd.fraction_superposed)   
    elif sd.rmsd != None:
      if sd.reference == 0:
        gfx_ent_2.UpdatePositions()
        gfx.Scene().CenterOn(gfx_ent_1)
      else:
        gfx_ent_1.UpdatePositions()
        gfx.Scene().CenterOn(gfx_ent_2)
      LogScript('RMSD: %.3f'%sd.rmsd)
    elif sd.superposition_error != None:
      LogScript('Superposition Failed: ' + sd.superposition_error)
    else:
      LogScript('Superposition Failed!')

class WindowMenu(QMenu):
  def __init__(self, parent=None):
    QMenu.__init__(self, parent)
    self.setTitle('Window')
    gosty=gui.GostyApp.Instance()
    self.aboutToShow.connect(self._AboutToShow)
    inspector_visible=gosty.GetWidget('InspectorDialog').isVisible()
    self._gl_win_action=gui.AddMenuAction(self, 'GL Window',  
                                      self._ToggleShowGLWindow, checkable=True,
                                      checked=gosty.gl_win.qobject.isVisible(),
                                      shortcut='Ctrl+G')
    self._inspector_action=gui.AddMenuAction(self, 'Inspector', 
                                         self._ToggleShowInspector,
                                         checkable=True,
                                         checked=inspector_visible,
                                         shortcut='Ctrl+I')
    self.addSeparator()
    self.addMenu(gosty.perspective.panels.menu)
    gui.AddMenuAction(self, 'Reset View', self._ResetView)
  def _AboutToShow(self):
    gosty=gui.GostyApp.Instance()
    self._gl_win_action.setChecked(gosty.gl_win.qobject.isVisible())
    inspector_visible=gosty.GetWidget('InspectorDialog').isVisible()
    self._inspector_action.setChecked(inspector_visible)
    
  def _ToggleShowGLWindow(self):
    gosty=gui.GostyApp.Instance()
    gl_win=gosty.GetGLWin()
    if gl_win and gl_win.qobject.isHidden():
      gl_win.Show()
    else:
      gl_win.Hide()

  def _ToggleShowInspector(self):
    gosty=gui.GostyApp.Instance()
    inspector=gosty.GetWidget('InspectorDialog')
    if inspector and inspector.isHidden():
      inspector.show()
    else:
      inspector.hide()
      
  def _ResetView(self):
    msg_box = QMessageBox()
    msg_box.setWindowTitle("Reset the Panels and Widget");
    msg_box.setIcon(QMessageBox.Question)
    msg_box.setText("Do you really want to reset the Panels and Widgets?");
    msg_box.setStandardButtons(QMessageBox.Yes |
                               QMessageBox.Cancel);
    msg_box.setDefaultButton(QMessageBox.Cancel);
    ret = msg_box.exec_();
    if(ret == QMessageBox.Yes):
      settings = QSettings()
      settings.setValue("restore_settings",QVariant(False))
      info_box = QMessageBox()
      info_box.setStandardButtons(QMessageBox.Ok)
      info_box.setIcon(QMessageBox.Information)
      info_box.setWindowTitle("Restart OpenStructure")
      info_box.setText("You must restart OpenStructure for the changes to take effect!");
      info_box.exec_();

class HelpMenu(QMenu):
  def __init__(self, parent=None):
    QMenu.__init__(self, parent)
    self.setTitle('Help')
    gui.AddMenuAction(self, 'Documentation', self._VisitDocs)
    gui.AddMenuAction(self, 'About', self._ShowAboutDialog)
    if sys.platform=='darwin':
      gui.AddMenuAction(self, 'Install Command Line Tool',
                    termuse.InstallTerminalPrograms)
  def _VisitDocs(self):
    QDesktopServices.openUrl(QUrl("http://www.openstructure.org/docs/"))
  
  def _ShowAboutDialog(self):
    _InitSplash()

def _InitMenu():
  _InitMenu.mbar=gui.GostyApp.Instance().perspective.GetMenuBar()
  file_menu=FileMenu(_InitMenu.mbar)
  scene_menu=SceneMenu(_InitMenu.mbar)
  win_menu=WindowMenu(_InitMenu.mbar)
  help_menu=HelpMenu(_InitMenu.mbar)
  _InitMenu.mbar.addMenu(file_menu)
  _InitMenu.mbar.addMenu(scene_menu)
  _InitMenu.mbar.addMenu(win_menu)
  _InitMenu.mbar.addMenu(help_menu)
