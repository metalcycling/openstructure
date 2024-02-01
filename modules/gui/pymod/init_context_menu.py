import platform

from PyQt5 import QtCore, QtWidgets

from ost import geom, gfx, gui, seq
from ost import settings
from ost import LogError, mol
from ost.bindings import WrappedTMAlign
from ost.bindings import msms
from ost.seq import alg
import ost
from ost.gui.scene.query_editor import QueryEditorWidget,QueryDialog

class SelectRefDialog(QtWidgets.QDialog):
  def __init__(self, ent_list, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.ent_list_ = ent_list
    vb = QtWidgets.QVBoxLayout()
    self.setLayout(vb)
    self.setWindowTitle("Select Reference Object")
    self.label = QtWidgets.QLabel("Please Select the Reference Object")
    self.list = QtWidgets.QTableWidget(self)
    self.list.horizontalHeader().setStretchLastSection(True)
    self.list.setColumnCount(2)
    self.list.verticalHeader().setVisible(False)
    self.list.horizontalHeader().setVisible(False)
    self.list.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    self.list.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    vb.addWidget(self.label)
    vb.addWidget(self.list)
    self.show_scores = QtWidgets.QCheckBox(self)
    self.show_scores.setText("Show Scores")
    self.show_scores.setChecked(True)
    vb.addWidget(self.show_scores)
    self.show_alignment = QtWidgets.QCheckBox(self)
    self.show_alignment.setText("Display Alignment")
    self.show_alignment.setChecked(False)
    vb.addWidget(self.show_alignment)
    hb = QtWidgets.QHBoxLayout()
    hb.setDirection(QtWidgets.QBoxLayout.LeftToRight)
    cancel_btn = QtWidgets.QPushButton("Cancel", self)
    load_btn = QtWidgets.QPushButton("Select", self)
    hb.addStretch(1)
    hb.addWidget(cancel_btn)
    hb.addWidget(load_btn)
    vb.addItem(hb)
    load_btn.setDefault(True)
    load_btn.clicked.connect(self.Select)
    cancel_btn.clicked.connect(self.reject)

    row = 0
    for ent in self.ent_list_:
      variant = QtCore.QVariant(ent)
      self.list.insertRow(row)
      new_item = QtWidgets.QTableWidgetItem("%i"%(row+1))
      new_item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
      new_item.setData(QtCore.Qt.UserRole,variant)
      self.list.setItem(row, 0, new_item)
      new_item = QtWidgets.QTableWidgetItem(ent.GetName())
      new_item.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
      self.list.setItem(row, 1, new_item)
      row += 1  

    self.list.resizeColumnsToContents()
    
  def Select(self):
    items = self.list.selectedItems()
    for item in items:
      if(item.column()==0):
        ent = item.data(QtCore.Qt.UserRole)
        self.ent_list_.remove(ent)
        self.ent_list_.insert(0,ent)
    self.accept()
     
  def GetShowScores(self):
    return self.show_scores.isChecked()
  
  def GetDisplayAlignment(self):
    return self.show_alignment.isChecked()
  
  def GetEntities(self):
    return self.ent_list_

class ShowResultDialog(QtWidgets.QDialog):
  def __init__(self, ent_list, res_list, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.ent_list_ = ent_list
    vb = QtWidgets.QVBoxLayout()
    self.setLayout(vb)
    self.setWindowTitle("Alignment result")
    self.label = QtWidgets.QLabel("Alignment results with %s as reference"%ent_list[0].GetName())
    self.list = QtWidgets.QTableWidget(self)
    self.list.horizontalHeader().setStretchLastSection(True)
    self.list.setColumnCount(4)
    self.list.verticalHeader().setVisible(False)
    self.list.setHorizontalHeaderLabels (["Name","RMSD","TMScore",""])
    self.list.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    vb.addWidget(self.label)
    vb.addWidget(self.list)
    hb = QtWidgets.QHBoxLayout()
    hb.setDirection(QtWidgets.QBoxLayout.LeftToRight)
    ok_btn = QtWidgets.QPushButton("OK", self)
    hb.addStretch(1)
    hb.addWidget(ok_btn)
    vb.addItem(hb)
    ok_btn.setDefault(True)
    ok_btn.clicked.connect(self.accept)

    for i in range(0, len(res_list)):
      self.list.insertRow(i)
      new_item = QtWidgets.QTableWidgetItem(ent_list[i+1].GetName())
      new_item.setFlags(QtCore.Qt.ItemIsEnabled)
      self.list.setItem(i, 0, new_item)
      new_item = QtWidgets.QTableWidgetItem("%.2f"%res_list[i].rmsd)
      new_item.setFlags(QtCore.Qt.ItemIsEnabled)
      self.list.setItem(i, 1, new_item)
      new_item = QtWidgets.QTableWidgetItem("%.2f"%res_list[i].tm_score)
      new_item.setFlags(QtCore.Qt.ItemIsEnabled)
      self.list.setItem(i, 2, new_item)
      new_item = QtWidgets.QTableWidgetItem()
      new_item.setFlags(QtCore.Qt.ItemIsEnabled)
      self.list.setItem(i, 3, new_item)
          
    self.list.resizeColumnsToContents()

class CalculateSurfaceSettingsDialog(QtWidgets.QDialog):
  def __init__(self, executable, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    vb = QtWidgets.QGridLayout()
    self.setLayout(vb)
    self.setWindowTitle("MSMS Surface Settings")
    msmsexe_label=QtWidgets.QLabel("executable")
    self.msmsexe_field=QtWidgets.QLineEdit()
    self.msmsexe_field.setText(executable)
    msmsexe_browsebutton=QtWidgets.QPushButton("Browse")
    vb.addWidget(msmsexe_label, 0, 0)
    vb.addWidget(self.msmsexe_field, 0, 1)
    vb.addWidget(msmsexe_browsebutton, 0, 2)
    surfname_label=QtWidgets.QLabel("surface name")
    self.surfname_field=QtWidgets.QLineEdit()
    self.surfname_field.setText("surface")
    vb.addWidget(surfname_label, 1, 0)
    vb.addWidget(self.surfname_field, 1, 1, 1, 2)
    density_label=QtWidgets.QLabel("density")
    self.density_spinbox=QtWidgets.QSpinBox()
    self.density_spinbox.setRange(1, 10)
    self.density_spinbox.setValue(4)
    vb.addWidget(density_label, 2, 0)
    vb.addWidget(self.density_spinbox, 2, 1, 1, 2)
    probe_label=QtWidgets.QLabel("probe radius")
    self.probe_spinbox=QtWidgets.QDoubleSpinBox()
    self.probe_spinbox.setDecimals(1)
    self.probe_spinbox.setSingleStep(0.1)
    self.probe_spinbox.setRange(0.3, 5.0)
    self.probe_spinbox.setValue(1.4)
    vb.addWidget(probe_label, 3, 0)
    vb.addWidget(self.probe_spinbox, 3, 1, 1, 2)
    selection_label=QtWidgets.QLabel("selection")
    self.selection_field=QtWidgets.QLineEdit()
    self.selection_field.setText("")
    vb.addWidget(selection_label, 4, 0)
    vb.addWidget(self.selection_field, 4, 1, 1, 2)
    self.noh_box=QtWidgets.QCheckBox("no hydrogens")
    vb.addWidget(self.noh_box, 5, 0)
    self.nohet_box=QtWidgets.QCheckBox("no hetatoms")
    vb.addWidget(self.nohet_box, 5, 1)
    self.nowat_box=QtWidgets.QCheckBox("no waters")
    vb.addWidget(self.nowat_box, 5, 2)
    
    cancel_btn = QtWidgets.QPushButton("Cancel", self)
    ok_btn = QtWidgets.QPushButton("OK", self)
    vb.addWidget(cancel_btn, 6, 1)
    vb.addWidget(ok_btn, 6, 2)
    
    msmsexe_browsebutton.clicked.connect(self.GetPath)
    ok_btn.clicked.connect(self.accept)
    cancel_btn.clicked.connect(self.reject)
    
  def GetPath(self):
    path, _ =QtWidgets.QFileDialog.getOpenFileName(self, "Choose MSMS Executable")
    if path!='':
      self.msmsexe_field.setText(path)

class SurfaceContextMenu(QtCore.QObject):
  def __init__(self, context_menu):
    try:
      settings_name="msms"
      self.executable=settings.Locate(settings_name)
    except settings.FileNotFound:
      self.executable=""
    QtCore.QObject.__init__(self, context_menu.qobject)
    self.action = QtWidgets.QAction("Calculate Surface", self)
    self.action.triggered.connect(self.CalculateSurface)
    context_menu.AddAction(self.action, gui.ContextActionType.ENTITY)
    
  def CalculateSurface(self):
    scene_selection = gui.SceneSelection.Instance()
    ent_list = list()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      ent_list.append(scene_selection.GetActiveNode(i))
    cssd = CalculateSurfaceSettingsDialog(self.executable)
    if(cssd.exec_()):
      self.__CalculateSurface(ent_list,
                              str(cssd.surfname_field.text()),
                              str(cssd.msmsexe_field.text()),
                              cssd.density_spinbox.value(),
                              cssd.probe_spinbox.value(),
                              str(cssd.selection_field.text()),
                              cssd.noh_box.isChecked(),
                              cssd.nohet_box.isChecked(),
                              cssd.nowat_box.isChecked())
      
  def __CalculateSurface(self,ent_list,name,msms_exe,density,
                         radius,selection,noh,nohet,nowat):
    for entity in ent_list:
      if isinstance(entity, gfx.Entity):
        try:
          s=msms.CalculateSurface(entity.view,
                                  msms_exe=msms_exe,
                                  density=density,
                                  radius=radius,
                                  selection=selection,
                                  no_hydrogens=noh,
                                  no_hetatoms=nohet,
                                  no_waters=nowat)[0]
          gfx.Scene().Add(gfx.Surface("%s_%s"%(entity.GetName(),name),s))
        except (RuntimeError, msms.MsmsProcessError):
          LogError("WARNING: Surface could not be calculated")
          return
        except UserWarning:
          LogError("WARNING: Entry with the same name already present in scene")
          return


class SecStructContextMenu(QtCore.QObject):
  def __init__(self, context_menu):
    QtCore.QObject.__init__(self, context_menu.qobject)
    self.action = QtWidgets.QAction("Assign Secondary Structure", self)
    self.action.triggered.connect(self.AssignSecondaryStructure)
    context_menu.AddAction(self.action, gui.ContextActionType.ENTITY)
    
  def AssignSecondaryStructure(self):
    scene_selection = gui.SceneSelection.Instance()
    ent_list = list()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      node = scene_selection.GetActiveNode(i)
      mol.alg.AssignSecStruct(node.view)
      node.UpdateView()

  
class AlignmentContextMenu(QtCore.QObject):

  def __init__(self, context_menu):
   
    QtCore.QObject.__init__(self, context_menu.qobject)
    self.action = QtWidgets.QAction("Align", self)
    self.action.triggered.connect(self.Align)
    context_menu.AddAction(self.action, gui.ContextActionType.ENTITY | gui.ContextActionType.MULTI)
    self.seq_viewer = None
    
  def Align(self):
    scene_selection = gui.SceneSelection.Instance()
    ent_list = list()
    for i in range(0,scene_selection.GetActiveNodeCount()):
      ent_list.append(scene_selection.GetActiveNode(i))
    sd = SelectRefDialog(ent_list)
    if(sd.exec_()):
      self.__Align(sd.GetEntities(),sd.GetShowScores(), sd.GetDisplayAlignment())
        
  def __Align(self, ent_list, show_scores=True, display_alignment=False):
    node = ent_list[0]
    res_list = list()
    if isinstance(node, gfx.Entity):
      ref = node.view.chains[0]
      for i in range(1,len(ent_list)):
        node = ent_list[i]
        if isinstance(node, gfx.Entity):
          tm_result = WrappedTMAlign(node.view.chains[0], ref)
          res_list.append(tm_result)
          node.view.handle.EditXCS().ApplyTransform(tm_result.transform)
          node.UpdatePositions()
    if show_scores:
      self.__ShowScore(ent_list, res_list)
    if display_alignment:
      self.__DisplayAlignment(ent_list, res_list)
    
  def __ShowScore(self, ent_list, res_list):
    if(len(res_list)==1):
      res = res_list[0]
      string = "RMSD: %.2f, TMScore: %.2f"%(res.rmsd, res.tm_score)
      gui.GostyApp.Instance().perspective.StatusMessage(string)
    elif(len(res_list)>1):
      ShowResultDialog(ent_list, res_list).exec_()
      
  def __DisplayAlignment(self, ent_list, res_list):
    if(len(res_list)>0):
      ref_seq = seq.CreateSequence("%s (ref)"%ent_list[0].GetName(),
                                   res_list[0].alignment.GetSequence(1).GetGaplessString())
      aln_list = seq.AlignmentList()
      if(ref_seq.IsValid()):
        for i in range(0, len(res_list)):
          # WrappedTMAlign returns an alignment with second sequence
          # being reference... let's swap...
          new_aln = seq.CreateAlignment()
          new_aln.AddSequence(res_list[i].alignment.GetSequence(1))
          new_aln.AddSequence(res_list[i].alignment.GetSequence(0))
          new_aln.SetSequenceName(1, ent_list[i+1].GetName())
          aln_list.append(new_aln)
        alignment = alg.MergePairwiseAlignments(aln_list, ref_seq)
        gosty = gui.GostyApp.Instance()
        main_area = gosty.perspective.GetMainArea()
        if self.seq_viewer:
          self.seq_viewer.qobject.close()
        self.seq_viewer = gui.SequenceViewer(True)
        self.seq_viewer.AddAlignment(alignment)
        self.seq_viewer.ChangeDisplayMode("Highlight conservation 1")
        self.seq_viewer.Show()

class SelectMenuPoints(QtCore.QObject):
  def __init__(self, context_menu):
    QtCore.QObject.__init__(self, context_menu.qobject)
    action=QtWidgets.QAction("Select...", self)
    action.triggered.connect(self._Select)
    context_menu.AddAction(action, gui.ENTITY)
    action=QtWidgets.QAction("Copy Selection...", self)
    action.triggered.connect(self._CopyViews)
    context_menu.AddAction(action, gui.ENTITY)
    action=QtWidgets.QAction('Select...', self)
    action.triggered.connect(self._SelectViewsSameEntity)
    context_menu.AddAction(action, gui.ENTITY_VIEW|gui.VIEWS_SAME_OBJECT)    
  def _Select(self):
    scene_selection=gui.SceneSelection.Instance()
    ent=scene_selection.GetActiveNode(0)
    dialog=QueryDialog('Select...')
    if dialog.exec_():
      q=mol.Query(dialog.query)
      if q.IsValid():
        ent.selection=ent.view.Select(dialog.query, dialog.query_flags)
      else:
        ost.LogError("invalid query: %s" % q.error)

  def _UniqueName(self, ent):
    """
    Returns a name based on ent that is unique within the scene
    """
    ent_name=ent.GetName()
    num=2
    while True:
      candidate_name='%s-%d' % (ent_name, num)
      if not gfx.Scene().HasNode(candidate_name):
        return candidate_name
      num+=1

  def _SelectViewsSameEntity(self):

    union=gui.SceneSelection.Instance().GetViewUnion()
    dialog=QueryDialog('Select...')
    if dialog.exec_():
      q=mol.Query(dialog.query)    
      if q.IsValid():
        ve=gui.SceneSelection.Instance().GetViewEntity()
        ve.selection=union.Select(q, dialog.query_flags)
      else:    
        ost.LogError("invalid query: %s" % q.error)    
  
  def _CopyViews(self):
    views_to_add=[]
    scene_selection=gui.SceneSelection.Instance()
    ent=scene_selection.GetActiveNode(0)
    dialog=QueryDialog('Select...')
    if dialog.exec_():
      q=mol.Query(dialog.query)
      if q.IsValid():
        for i in range(scene_selection.GetActiveNodeCount()):
          ent=scene_selection.GetActiveNode(i)
          selected=ent.view.Select(q, dialog.query_flags)
          gfx_ent=gfx.Entity(self._UniqueName(ent),selected)
          gfx.Scene().Add(gfx_ent)
      else:    
        ost.LogError("invalid query: %s" % q.error)

def _InitContextMenu(app):
  _InitContextMenu.cm=app.scene_win.GetContextMenu()
  _InitContextMenu.am=AlignmentContextMenu(_InitContextMenu.cm)
  _InitContextMenu.sec=SecStructContextMenu(_InitContextMenu.cm)
  _InitContextMenu.sc=SurfaceContextMenu(_InitContextMenu.cm)
  _InitContextMenu.sm=SelectMenuPoints(_InitContextMenu.cm)
