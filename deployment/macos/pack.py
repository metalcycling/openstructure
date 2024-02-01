"""
deploy.py helps you package your code into a standalone application.

for now, macOS and Linux are supported...
"""
import os
import subprocess
import shutil
import sys
import glob

def _lib_name(component):
  return 'lib%s.dylib' % component

def _deps_for_lib(lib, pool, recursive=True):
  if lib in pool:
    return
  otool=subprocess.Popen(['otool', '-L', lib], stdout=subprocess.PIPE)
  output=otool.communicate()[0]
  lines=output.split('\n')[1:]
  for line in lines:
    d=line.split(' ')[0].strip()
    if len(d)>0:
      if d==lib:
        continue
      if d not in pool:
        if d.startswith('/System') or d.startswith('/usr/lib'):
          continue
        if recursive:
          _deps_for_lib(d, pool)
        pool.add(d)        
  return

def collect_deps(stage_dir, components, binaries, libexec_binaries,
                 site_packages, site_packages_dir, libexec_path='openstructure'):
  """
  Collect the dependencies for the given components and returns a list of 
  frameworks/libraries that the component depends on.
  """
  pool=set()
  for component in components:
    lib_name=os.path.abspath(os.path.join(stage_dir, 'lib', 
                                          _lib_name(component)))  
    if not os.path.exists(lib_name):
      print('WARNING:', lib_name, 'does not exist')
    if lib_name not in pool:
      _deps_for_lib(lib_name, pool)
      pool.add(lib_name)    
  for bin in binaries:  
    bin_name=os.path.abspath(os.path.join(stage_dir, 'bin', 
                                          bin))  
    if not os.path.exists(bin_name):
      print('WARNING:', bin_name, 'does not exist')
      continue
    if bin_name not in pool:
      _deps_for_lib(bin_name, pool)
  for bin in libexec_binaries:
    bin_name=os.path.abspath(os.path.join(stage_dir, 'libexec', libexec_path,
                                          bin))
    if not os.path.exists(bin_name):
      print('WARNING:', bin_name, 'does not exist')
      continue
    if bin_name not in pool:
      _deps_for_lib(bin_name, pool)
  for site_package in site_packages:
    full_path=get_python_module_path(site_package)
    print(full_path)
    if not os.path.exists(full_path):
      print('WARNING:', site_package, 'does not exists')
      continue
    if os.path.isdir(full_path):
      for so_file in glob.glob(os.path.join(full_path, '*.so')):
        _deps_for_lib(so_file, pool)
  return pool

LIBEXEC_SCRIPTS=['ost_config']
LIBEXEC_BINARIES=[]
GUI_LIBEXEC_BINARIES=['gosty']
BINARIES=['ldt', 'chemdict_tool', 'tmalign', 'tmscore']
GUI_BINARIES=[]
GUI_COMPONENTS=['gfx', 'gui', 'info']
COMPONENTS=['mol', 'geom', 'conop', 'seq_alg', 'seq',
            'img', 'img_alg', 'io', 'db', 'base']
GUI_SCRIPTS=['dng']
SCRIPTS=['ost']
CHANGE_ID_RPATH='install_name_tool -id @rpath/%s %s'   
CHANGE_ID='install_name_tool -id @rpath/%s %s'
CHANGE_LOAD_CMD_RPATH='install_name_tool -change %s @rpath/%s %s'
CHANGE_LOAD_CMD='install_name_tool -change %s @executable_path/%s %s'
ADD_RPATH='install_name_tool -add_rpath %s %s 2> /dev/null'
SITE_PACKAGES=[]
GUI_SITE_PACKAGES=['sip.so', 'sipconfig.py', 'sipdistutils.py', 'PyQt4']
REMOVE_HEADERS='rm -rf `find %s/lib -type d -name Headers`'
REMOVE_CURRENT='rm -rf `find %s/lib -type d -name Current`'
# collect libs of non-standard libraries/frameworks we depend on

def copy_binaries(stage_dir, outdir, binary_names, scripts, bin_dir,
                  append_bin=True):

  exe_path=os.path.abspath(os.path.join(outdir, bin_dir))
  for binary_name in binary_names:
    if append_bin:
      bin_name=os.path.join(stage_dir, bin_dir, binary_name)
    else:
      bin_name=os.path.join(stage_dir, binary_name)
    if not os.path.exists(bin_name):
      print('WARNING:', binary_name, 'does not exist')
      continue
    dst_name=os.path.join(outdir, bin_dir, os.path.basename(bin_name))
    shutil.copy(bin_name, dst_name)
    update_load_commands(dst_name, True, exe_path)
  for script in scripts:
    shutil.copy(os.path.join(stage_dir, bin_dir, script),
                os.path.join(outdir,bin_dir, script))

def split_framework_components(abs_path):
    """
    Splits the path pointing to a dynamic library within a framework
    
    '/System/Frameworks/A.framework/Versions/4/A' =>
    ['/System/Frameworks/A.framework', 'Versions/4/A']
    """
    parts=abs_path.split('/')
    for i, s in enumerate(parts):
      if s.endswith('.framework'):
        lead=os.path.join('/', *parts[:i+1])
        trail=os.path.join(*parts[i+1:])
        return lead, trail

def change_id(id, lib):
  os.chmod(lib, 0o666)
  os.system(CHANGE_ID_RPATH % (id,lib))
  os.chmod(lib, 0o444)

def update_load_commands(lib, exe, exe_path):
  direct_deps=set()
  _deps_for_lib(lib, direct_deps, recursive=False)
  os.chmod(lib, 0o666)
  for direct_dep in direct_deps:
    if direct_dep.endswith('.dylib'):
      new_name=os.path.basename(direct_dep)
      os.system(CHANGE_LOAD_CMD_RPATH % (direct_dep, new_name, lib))
    else:
      assert direct_dep.find('.framework/')>=0
      framework_path, rel_path=split_framework_components(direct_dep)
      framework_name=os.path.basename(framework_path)
      new_name=os.path.join(framework_name, rel_path)
      os.system(CHANGE_LOAD_CMD_RPATH % (direct_dep, new_name, lib))
  if exe:
    os.chmod(lib, 0o555)
  else:
    os.chmod(lib, 0o444)

def copy_deps(dependencies, outdir):
  exe_path=os.path.join(outdir, 'bin')
  for dep in dependencies:
    if dep.endswith('.dylib'):
      dst_name=os.path.join(outdir, 'lib', os.path.basename(dep))
      if not os.path.exists(dep):
        continue
      if os.path.exists(dst_name):
        continue
      shutil.copy(dep, dst_name)
      change_id(os.path.basename(dep), dst_name)
      update_load_commands(dst_name, False, exe_path)
    else:
      assert dep.find('.framework/')>=0
      framework_path, rel_path=split_framework_components(dep)
      framework_name=os.path.basename(framework_path)
      dst_name=os.path.join(outdir, 'lib', framework_name)
      shutil.copytree(framework_path, dst_name)
      change_id(os.path.join(dst_name, rel_path),
                os.path.join(dst_name, rel_path))
      os.unlink(os.path.join(dst_name, os.path.splitext(framework_name)[0]))
      update_load_commands(os.path.join(dst_name, rel_path), False, 
                           exe_path)

def update_pymod_shared_objects(lib_path, path, files):
  exe_path=os.path.abspath(os.path.join(lib_path, '../bin'))
  for f in files:
    if not os.path.exists(os.path.join(path, f)):
      continue
    base, ext=os.path.splitext(f)
    if  ext=='.so':
      path_to_lib_path=os.path.relpath(lib_path, path)
      abs_name=os.path.join(path, f)
      os.system(ADD_RPATH % (path_to_lib_path, abs_name))
      update_load_commands(abs_name, False, exe_path)
    elif ext in ('.pyc', '.pyo'):
      os.unlink(os.path.join(path, f))

def merge_tree(src, dst):
  """
  Similar to shutil.copytree, but does not complain when the destination
  directory already exists.
  """
  names = os.listdir(src)
  if not os.path.exists(dst):
    os.makedirs(dst)
  errors = []
  for name in names:
    srcname = os.path.join(src, name)
    dstname = os.path.join(dst, name)
    try:
        if os.path.islink(srcname):
            linkto = os.readlink(srcname)
            os.symlink(linkto, dstname)
        elif os.path.isdir(srcname):
            merge_tree(srcname, dstname)
        else:
            shutil.copy2(srcname, dstname)
    except (IOError, os.error) as why:
        errors.append((srcname, dstname, str(why)))
    except shutil.Error as err:
        errors.extend(err.args[0])
  try:
      shutil.copystat(src, dst)
  except OSError as why:
      if WindowsError is not None and isinstance(why, WindowsError):
          # Copying file access times may fail on Windows
          pass
      else:
          errors.extend((src, dst, str(why)))
  if errors:
      raise shutil.Error(errors)

def get_site_package_dir():
  """
  Get site-package directory of this python installation. This assumes 
  that ost was linked against the same version of Python (which is a very
  reasonable thing to do, as this script is most likely run with ost).
  """
  for p in sys.path:
    pattern='/site-packages'
    index=p.find(pattern)
    if index>=0:
      return p[:index+len(pattern)]
  raise RuntimeError("Couldn't determine site-packages location")

def get_python_module_path(module):
  for path in sys.path:
    full_path=os.path.join(path, module)
    if os.path.exists(full_path):
      return full_path
  return None
  
  
def get_python_home():
  """
  Derive Python home by looking at the location of the os module
  """
  return os.path.dirname(sys.modules['os'].__file__)


class Package(object):
  def __init__(self, name, root_dir, binaries=[], scripts=[], 
               modules=[], libraries=[], libexec_dir=None, libexec_scripts=[]):
    self.root_dir=root_dir
    self.name=name
    self.binaries=binaries
    self.scripts=scripts
    self.libraries=libraries
    self.libexec_dir=libexec_dir
    self.libexec_scripts=libexec_scripts
    self.pymod_dir=os.path.join('lib', 'python%d.%d' % sys.version_info[0:2],
                                'site-packages')
    self.modules=modules
    self.libexec_binaries=[]
    self.site_packages=[]
    self.site_packages_dir=''
  def status(self, message):
    print('%s: %s' % (self.name, message))

  def _prepare_output_dir(self, output_dir):
    """
    Prepares the output directory structure, including lib, bin and an optional
    libexec directory.
    """
    #if os.path.exists(output_dir):
    #  shutil.rmtree(output_dir)
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    if not os.path.exists(os.path.join(output_dir, 'bin')):
      os.makedirs(os.path.join(output_dir, 'bin'))
    if not os.path.exists(os.path.join(output_dir, 'lib')):
      os.makedirs(os.path.join(output_dir, 'lib'))
    if self.libexec_dir:
      out_exec_dir=os.path.join(output_dir, 'libexec', self.libexec_dir)
      if not os.path.exists(out_exec_dir):
        print('making...', out_exec_dir)
        os.makedirs(out_exec_dir)
  def _copy_site_packages(self, output_dir):
    for sp in SITE_PACKAGES:
      src=get_python_module_path(sp)
      if os.path.isdir(src):
        merge_tree(src, os.path.joini(output_dir, self.pymod_dir, sp))
      else:
        shutil.copy(src, os.path.join(output_dir, self.pymod_dir, sp))
    print('updating link commands of python shared objects')
    os.path.walk(os.path.join(output_dir, 'lib'), 
                 update_pymod_shared_objects, 
                 os.path.join(output_dir, 'lib'))

  def ship(self, output_dir):
    self._prepare_output_dir(output_dir)
    if os.path.exists(os.path.join(self.root_dir, 'share')):
      self.status('copying shared data files')
      merge_tree(os.path.join(self.root_dir, 'share'), 
                 os.path.join(output_dir, 'share'))
    self.status('collecting dependencies')
    deps=collect_deps(self.root_dir, self.libraries, self.binaries, 
                      self.libexec_binaries, self.site_packages, 
                      self.site_packages_dir)
    # when running in non-gui mode, we are most likely missing the boost
    # python library. Let's add it to the list of dependencies by
    # inspecting "_ost_base.so".
    pymod_dir='lib/python%d.%d/site-packages' % sys.version_info[0:2]
    _deps_for_lib(os.path.join(self.root_dir, pymod_dir, 'ost/_ost_base.so'),
                  deps, recursive=False)
    self.status('copying dependencies')
    copy_deps(deps, output_dir)
    if self.libexec_dir:
      self.status('copying libexec binaries')
      copy_binaries(self.root_dir, output_dir, self.libexec_binaries, 
                    self.libexec_scripts,
                    os.path.join('libexec', self.libexec_dir))
    self.status('copying binaries')
    copy_binaries(self.root_dir, output_dir, self.binaries, 
                self.scripts, 'bin')
    self.status('copying pymod')
    merge_tree(os.path.join(self.root_dir,self.pymod_dir), 
               os.path.join(output_dir, self.pymod_dir))
    self._copy_site_packages(output_dir)

class OpenStructure(Package):
  def __init__(self, stage_dir, minimal=True):
    libs=['ost_mol', 'ost_geom', 'ost_conop', 'ost_seq_alg',
          'ost_io', 'ost_db', 'ost_base', 'ost_seq', 'ost_mol_alg']
    super(OpenStructure, self).__init__('OpenStructure', stage_dir, 
                                        binaries=['ldt', 'chemdict_tool'],
                                        libexec_scripts=['ost_config'],
                                        scripts=['ost'],
                                        libraries=libs,
                                        libexec_dir='openstructure')

class Qmean(Package):
  def __init__(self, stage_dir):
    super(Qmean, self).__init__('Qmean', stage_dir,
                                scripts=['qmean'],
                                libexec_scripts=['qmean_script.py'],
                                libexec_dir='qmean',
                                libraries=['qmean'])
ost=OpenStructure('../../stage')
qmean=Qmean('../../../../../qmean/stage')
ost.ship('qmean')
qmean.ship('qmean')

