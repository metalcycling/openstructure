from PyQt5 import QtWidgets
from ost import gui
class DocWidget(QtWidgets.QWidget):
  def __init__(self, parent=None):
    QtWidgets.QWidget.__init__(self, parent)
    l=QtWidgets.QVBoxLayout(self)
    l.setContentsMargins(0,0,0,0)
    self.setWindowTitle('OpenStructure Help')
    self.searchbox=QtWidgets.QLineEdit(self)
    l.addWidget(self.searchbox)
    self.doctext=QtWidgets.QTextEdit('', self)
    l.addWidget(self.doctext)
    self.setLayout(l)
    self.searchbox.returnPressed.connect(self.OnReturnPressed)
  def OnReturnPressed(self):
    help(str(self.searchbox.text()))
doc_widget=DocWidget()
doc_widget_for_panel=gui.Widget(doc_widget)
panels=gui.GostyApp.Instance().perspective.panels
panels.AddWidgetToPool("OpenStructure Help",doc_widget_for_panel)


class Help:
  def __call__(self, what=None):
    import pydoc
    d=pydoc.HTMLDoc()
    if what==None:
      doc_widget.doctext.setHtml('')
      doc_widget.searchbox.setText('')
    else:
      try:
        s=d.document(pydoc.resolve(what)[0])
      except ImportError as detail:
        s=str(detail)
      doc_widget.doctext.setHtml(s)
      doc_widget.searchbox.setText('')
    if not doc_widget.isVisible():
      panels.AddWidget(gui.PanelPosition.RIGHT_PANEL, doc_widget_for_panel, False)
  def __repr__(self):
    return 'Type help() for interactive help, or help(object) for help about object.'

help=Help()
