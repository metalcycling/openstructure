import __main__
import sys
import os.path
import optparse
from ost import io, mol, seq, geom, conop, gui, settings, gfx

import ost
try: 
  from ost import img
  import ost.img.alg
  _img_present=True
except ImportError:
  _img_present=False
  pass

import http.client

from PyQt5 import QtGui, QtWidgets, QtCore
from ost.gui.scene.init_inspector import _InitInspector
from ost.gui.init_menubar import _InitMenuBar
from ost.gui.init_spacenav import _InitSpaceNav
from ost.gui.init_context_menu import _InitContextMenu
from ost.gui.init_splash import _InitSplash
from ost.gui.dng import termuse
from ost.gui.dng import superpositiondialog
from ost.gui.scene.remote import RemoteLoader

import ost.gui.dng.menu

def _my_exit(code):
  QtWidgets.QApplication.instance().quit()
  gui.GostyApp.Instance().ProcessEvents()
  sys._exit(code)

sys._exit=sys.exit
sys.exit=_my_exit

def _InitRuleBasedProcessor():
  compound_lib_path=os.path.join(ost.GetSharedDataPath(), 'compounds.chemlib')
  if os.path.exists(compound_lib_path):
    conop_inst=conop.Conopology.Instance()
    compound_lib=conop.CompoundLib.Load(compound_lib_path)
    conop_inst.SetDefaultLib(compound_lib)
    io.profiles['DEFAULT'].processor = conop.RuleBasedProcessor(compound_lib)

# switch to rule-based processor for high fidelity if compounds.chemlib is 
# available
_InitRuleBasedProcessor()

def _CheckRestore():
  settings = QtCore.QSettings()
  restore = settings.value("restore_settings",QtCore.QVariant(False))
  if not restore:
    settings.clear()
  settings.setValue("restore_settings",QtCore.QVariant(True))

def _InitPanels(app):
  panels = app.perspective.panels
  panels.AddWidgetToPool('ost.gui.FileBrowser', -1)
  panels.AddWidgetToPool('ost.gui.PythonShell', 1)
  panels.AddWidgetToPool('ost.gui.SceneWin', 1)
  panels.AddWidgetToPool('ost.gui.SequenceViewer', 1)
  panels.AddWidgetToPool('ost.gui.MessageWidget', 1)
  if not panels.Restore("ui/perspective/panels"):
    panels.AddWidget(gui.PanelPosition.LEFT_PANEL, app.scene_win)
    panels.AddWidgetByName(gui.PanelPosition.LEFT_PANEL, 
                           'ost.gui.FileBrowser', False)    
    panels.AddWidget(gui.PanelPosition.BOTTOM_PANEL, app.seq_viewer)
    panels.AddWidget(gui.PanelPosition.BOTTOM_PANEL, app.py_shell)
    panels.AddWidget(gui.PanelPosition.RIGHT_PANEL, app.message_widget)
    return False
  return True


def _InitFrontEnd(try_stereo):
  _CheckRestore()
  app=gui.GostyApp.Instance()
  app.SetAppTitle("DNG")
  app.TryStereo(try_stereo)
  main_area=app.perspective.main_area
  _InitSpaceNav(app)
  _InitContextMenu(app)
  main_area.AddPersistentWidget("3D Scene", "gl_win" , app.gl_win,
                                int(QtCore.Qt.WindowMaximized))
  _InitInspector(app)
  ost.gui.dng.menu._InitMenu()
  app.perspective.Restore()
  additional_modules=getattr(__main__, 'ADDITIONAL_GUI_MODULES', [])
  for module_name in additional_modules:
    __import__(module_name)
  app.ProcessEvents()
  


  if not _InitPanels(app):
    _InitSplash()
  
def _load_files():
  input_files=[_SplitIDSel(arg) for arg in loading_list]
  for f in input_files:
    selection=_get_selection_query(f[1])
    gui.FileLoader.LoadObject(f[0],selection)

def _get_selection_query(sel):
  if len(options.query)>0:
    if len(sel)>0:
      return '(%s) and (%s)' % (options.query, sel)
    else:
      return options.query
  elif len(sel)>0:
    return sel
  return ""

def _execute_script():
  script=script_argv[0]
  sys_argv_backup=sys.argv
  sys.argv=script_argv
  try:
    exec(compile(open(script).read(), script, 'exec'), __main__.__dict__)
  finally:
    sys.argv=sys_argv_backup     

def show_help(option, opt, value, parser):
  parser.print_help()
  QtWidgets.QApplication.instance().exit()
  sys.exit(-1)

def parse_script_option(option, opt, value, parser):
  script_argv.append(value)
  for arg in parser.rargs:
     script_argv.append(arg)
  del parser.rargs[0:len(parser.rargs)]

def _SplitIDSel(name):
  pos=name.find('[')
  if pos>-1:
    return name[:pos], name[pos+1:-1]
  return name, ''  

def stop():
  gui.GostyApp.Instance().StopScript()


loading_list=[]
script_argv=[]
#images=[]
#viewers=[]
usage = 'usage: dng [options] [files to load]'
class OstOptionParser(optparse.OptionParser):
  def __init__(self, **kwargs):
    optparse.OptionParser.__init__(self, **kwargs)
  def exit(self, status_code, error_message):
    print(error_message, end=' ')
    QtWidgets.QApplication.instance().exit()
    sys.exit(-1)

parser=OstOptionParser(usage=usage,conflict_handler="resolve")
parser.add_option("-h", "--help", action="callback", callback=show_help, help="show this help message and exit")
parser.add_option("-v", "--verbosity_level", action="store", type="int", dest="vlevel", default=2, 
                  help="sets the verbosity level [default: %default]")
parser.add_option("-s", "--script", action="callback", default=[], dest="script", type="string", callback=parse_script_option, help="executes a script (syntax: -s SCRIPT [options] [args]) Anything that follows this option is passed to the script")
parser.add_option("-p", "--processor", dest="processor", default="HEURISTIC", help="Type of processor used by the progam (either RULE_BASED or HEURISTIC) [default: %default]")
parser.add_option("-c", "--compound_library", dest="complib", default="compounds.chemlib", help="Compound library for the RULE_BASED processor (only used if --processor option is set to RULE_BASED, otherwise ignored [default: %default]")
parser.add_option("-q", "--query", dest="query", default="", help="Selection query to be highlighted automatically upon loading (only used when a PDB file is loaded, otherwise ignored [default: None]")
parser.add_option("-S","--stereo", dest="try_stereo", default=False, action="store_true",help="try to get a quad-buffer stereo visual")
parser.disable_interspersed_args()
(options, args) = parser.parse_args()

if len(parser.rargs)!=0:
  for rargs_string in parser.rargs:
    if not rargs_string.endswith('.py'):  
      loading_list.append(rargs_string)
    else:
      print('Error:  one of the files to load is a Python script, use -s flag to execute it\n')
      QtWidgets.QApplication.instance().exit()
      sys.exit(-1)    

if len(options.script)!=0:
  script_argv=options.script

if options.processor=="RULE_BASED":
  from ost import conop
  if os.path.exists(options.complib):
    compound_lib=conop.CompoundLib.Load(options.complib)
    conop.SetDefaultLib(compound_lib)
    io.profiles['DEFAULT'].processor = conop.RuleBasedProcessor(compound_lib)

home = os.getenv('HOME') or os.getenv('USERPROFILE')
_ostrc=os.path.join(home, '.ostrc')
if os.path.exists(_ostrc):
  try:
    exec(compile(open(_ostrc).read(), _ostrc, 'exec'))
  except Exception as e:
    print(e)
else:
  rcfile=open(_ostrc,"w")
  print('# This python file is parsed by ost and dng at startup', file=rcfile)
  print('# Its content is made available in the global namespace', file=rcfile)
  print('# It can be used to define custom variables and functions', file=rcfile)
  print('# For example:', file=rcfile)
  print('# IMPORTANT_DIR="path/to/important/dir"', file=rcfile)
  rcfile.close()

ost.gui.PushVerbosityLevel(options.vlevel)
working_dir=settings.GetValue("DNG_WORKING_DIR",None)

if working_dir != None and os.path.isdir(working_dir):
  os.chdir(working_dir)

_InitFrontEnd(options.try_stereo)

if len(loading_list)!=0:
  _load_files()
  gfx.Scene().Autoslab()
if len(script_argv)!=0:
  _execute_script()
