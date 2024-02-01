from PyQt5 import QtGui, QtCore, QtWidgets
from ost.gui import AdminRights
import ost
import os
#from ost.gui import AdminRights

usage='''The DNG application bundle contains two shell commands ('ost' and 'dng') that lets you use the OpenStructure command-line interpreter and dng from the terminal. If you want to use these commands, it is recommended that a symbolic link is created. 
'''
class TerminalUsageDialog(QtWidgets.QDialog):
  def __init__(self, parent=None):
    QtWidgets.QDialog.__init__(self, parent)
    self.setWindowTitle('Enhanced Terminal Usage')
    self.setFixedSize(QSize(480, 300))
    l=QVBoxLayout(self)
    title=QLabel('Enhanced Terminal Usage')
    font=title.font()
    font.setPointSize(20)
    title.setFont(font)
    l.addWidget(title)
    text=QLabel(usage)
    l.addWidget(text)
    l2=QHBoxLayout()
    l2.addWidget(QLabel('If you proceed, the link will be created in: '))
    self.path_combo=QComboBox()
    self.path_combo.setFixedWidth(150)
    for path in os.getenv('PATH').split(':'):
      exp_path=os.path.expanduser(path)
      if os.path.exists(exp_path) and exp_path.find('DNG.app')==-1:
        self.path_combo.addItem(path)
    l2.addWidget(self.path_combo)
    l.addLayout(l2)
    l3=QHBoxLayout()
    ab=QPushButton('Create Link')
    ab.setDefault(True)
    cb=QPushButton('Don\'t Create')
    l3.addStretch(1)
    l3.addWidget(cb, 0)
    l3.addWidget(ab, 0)
    l.addLayout(l3)
    text.setWordWrap(True)
    QObject.connect(cb, SIGNAL('clicked()'), self.reject)
    QObject.connect(ab, SIGNAL('clicked()'), self.accept)    
  def GetSelectedPath(self):
    return str(self.path_combo.currentText())

def _CreateLinks(bin_dir, sel_dir):
  for bin in ('ost', 'dng','lddt', 'chemdict_tool'):
    if os.path.exists(os.path.join(sel_dir, bin)):
      os.unlink(os.path.join(sel_dir, bin))
    os.system('ln -s "%s" "%s"' % (os.path.join(bin_dir, bin), 
                                   os.path.join(sel_dir, bin)))
def InstallTerminalPrograms():
  """
  Installs symlinks to the 'ost' and 'dng' command line programs into a 
  user-specified directory in the path.
  """  
  term_use=TerminalUsageDialog()
  if term_use.exec_():
    prefix=ost.GetPrefixPath()
    bin_dir=os.path.join(prefix, 'bin')
    sel_path=term_use.GetSelectedPath()    
    if not os.access(sel_path, os.W_OK):
      admin_rights=AdminRights()
      if admin_rights.Acquire():
        for bin in ('ost', 'dng', 'lddt', 'chemdict_tool'):
          admin_rights.CreateLink(os.path.join(bin_dir, bin), 
                                  os.path.join(sel_path, bin))
      admin_rights.Release()
    else:
      _CreateLinks(bin_dir, sel_path)
