from ost import mol
from PyQt5 import QtCore, QtWidgets, QtGui

class QueryEditorWidget(QtWidgets.QWidget):
  def __init__(self, parent=None): 
    QtWidgets.QWidget.__init__(self, parent)
    self.default_font_=QtGui.QTextCharFormat()
    self.default_font_.setForeground(QtGui.QBrush(QtGui.QColor(0,0,0)))
    self.error_font_=QtGui.QTextCharFormat()
    self.error_font_.setForeground(QtGui.QBrush(QtGui.QColor(255,0,0))) 
    self.selection_edit_=QtWidgets.QTextEdit(self)
    self.selection_edit_.setFixedHeight(40)
    self.selection_edit_.updateGeometry()
    self.status_=QtWidgets.QLabel(" ",self);
    self.status_.setWordWrap(True)
    self.status_.setMargin(0)
    self.status_.setAlignment(QtCore.Qt.AlignRight)
    self.status_.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum,
                                                 QtWidgets.QSizePolicy.Expanding))
    vl=QtWidgets.QVBoxLayout()
    vl.addWidget(self.selection_edit_)
    self.no_bonds_=QtWidgets.QRadioButton('none')
    self.ex_bonds_=QtWidgets.QRadioButton('exclusive')
    self.in_bonds_=QtWidgets.QRadioButton('inclusive')
    self.in_bonds_.setChecked(True)
    self.match_res_=QtWidgets.QCheckBox('match residues')
    
    vl.setContentsMargins(0,0,0,0)
    vl.setSpacing(0)
    self.setLayout(vl)
    vl.addWidget(self.status_)
    hl=QtWidgets.QHBoxLayout()
    l=QtWidgets.QLabel("bonds:")

    hl.addWidget(l)
    hl.addSpacing(5)
    hl.addWidget(self.no_bonds_)
    hl.addWidget(self.ex_bonds_)
    hl.addWidget(self.in_bonds_)
    hl.addStretch(1)
    vl.addLayout(hl)
    vl.addWidget(self.match_res_)

    self.changing_text_=False;
    self.selection_edit_.textChanged.connect(self._StartTimer)
    self.timer_=QtCore.QTimer()
    self.timer_.timeout.connect(self._UpdateMessage)
    
  def GetQueryFlags(self):
    flags=0
    if self.no_bonds_.isChecked():
      flags|=mol.NO_BONDS
    if self.ex_bonds_.isChecked():
      flags|=mol.EXCLUSIVE_BONDS
    if self.match_res_.isChecked():
      flags|=mol.MATCH_RESIDUES
    return flags
  
  def GetQuery(self):
    return mol.Query(str(self.selection_edit_.toPlainText()))

  def GetQueryText(self):
    return str(self.selection_edit_.toPlainText())
  
  def SetQueryFlags(self,flags):
    self.no_bonds_.setChecked(flags & mol.NO_BONDS)
    self.ex_bonds_.setChecked(flags & mol.EXCLUSIVE_BONDS)
    self.in_bonds_.setChecked(not (flags & mol.NO_BONDS|mol.EXCLUSIVE_BONDS))
    self.match_res_.setChecked(flags & mol.MATCH_RESIDUES)
  
  def SetQuery(self,query):
    self.selection_edit_.setText(query)
    
  def _StartTimer(self):
    self.timer_.stop()
    self.timer_.start(500)

  def _UpdateMessage(self):
    if self.changing_text_:
      return
    self.changing_text_ = True
    query=self.GetQuery()
    
    cursor=self.selection_edit_.textCursor()
    cursor.select(QtGui.QTextCursor.Document)
    cursor.setCharFormat(self.default_font_)  
    
    if query.IsValid():
      self.status_.setText("")
    else:      
      d=query.GetErrorDescription()
      self.status_.setText("<font color='red'>%s</font>"%d.msg.strip())
      #self.status_.setFixedSize(self.width(),self.status_.height())
      cursor.movePosition(QtGui.QTextCursor.Start)
      if d.range.Loc<len(query.string):
        
        cursor.movePosition(QtGui.QTextCursor.NextCharacter, 
                            QtGui.QTextCursor.MoveAnchor, d.range.Loc)
        cursor.movePosition(QtGui.QTextCursor.NextCharacter, 
                            QtGui.QTextCursor.KeepAnchor, d.range.Length)
        cursor.setCharFormat(self.error_font_)
    self.changing_text_=False
    
class QueryDialog(QtWidgets.QDialog):
  def __init__(self, title, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    l=QtWidgets.QVBoxLayout(self)
    self.setWindowTitle(title)
    self.editor=QueryEditorWidget(self)
    l.addWidget(self.editor)
    l.addSpacing(10)
    l3=QtWidgets.QHBoxLayout()
    ab=QtWidgets.QPushButton('OK')
    ab.setDefault(True)
    cb=QtWidgets.QPushButton('Cancel')
    l3.addStretch(1)
    l3.addWidget(cb, 0)
    l3.addWidget(ab, 0)
    l.addLayout(l3)
    cb.clicked.connect(self.reject)
    ab.clicked.connect(self.accept)
    self.editor.selection_edit_.textChanged.connect(self._CheckNewline)

  @property 
  def query_flags(self):
    return self.editor.GetQueryFlags()
  @property
  def query(self):
    return self.editor.GetQueryText()
  
  def event(self, e):
    if e.type() == QtCore.QEvent.KeyPress and (e.key () == QtCore.Qt.Key_Return or e.key () == QtCore.Qt.Key_Enter):
      self.accept()
      return True
    else:
      return QtWidgets.QDialog.event(self, e)
  
  def _CheckNewline(self):
    if self.editor.GetQueryText().endswith("\n"):
      self.accept()
