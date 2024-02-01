import os
import subprocess
import shutil
import sys
import glob

def _lib_name(component):
  return 'libost_%s.dylib' % component

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
    bin_name=os.path.abspath(os.path.join(stage_dir, 'libexec/openstructure',
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
BINARIES=['lddt', 'chemdict_tool', 'tmalign', 'tmscore']
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

def get_site_package_dir():
  """
  Get site-package directory of this python installation. This assumes 
  that ost was linked against the same version of Python
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

def make_standalone(stage_dir, outdir, no_includes, force_no_rpath=False,
                    macports_workaround=False, no_gui=False):
  site_packages_dir=get_site_package_dir()

  if os.path.exists(outdir):
    shutil.rmtree(outdir)
  os.system('mkdir -p "%s"' % outdir)
  os.system('mkdir -p "%s/lib"' % outdir)
  os.system('mkdir -p "%s/bin"' % outdir)
  os.system('mkdir -p "%s/libexec/openstructure"' % outdir)
  print('copying shared datafiles')
  shutil.copytree(os.path.join(stage_dir, 'share'), 
                  os.path.join(outdir, 'share'))
  scripts=SCRIPTS
  binaries=BINARIES
  components=COMPONENTS
  libexec_scripts=LIBEXEC_SCRIPTS
  site_packages=SITE_PACKAGES
  libexec_binaries=LIBEXEC_BINARIES
  if not no_gui:
    scripts+=GUI_SCRIPTS
    binaries+=GUI_BINARIES
    components+=GUI_COMPONENTS
    libexec_binaries+=GUI_LIBEXEC_BINARIES
    site_packages+=GUI_SITE_PACKAGES
  print('collecting dependencies')
  deps=collect_deps(stage_dir, components, binaries, libexec_binaries,
                    site_packages, site_packages_dir)
  # when running in non-gui mode, we are most likely missing the boost
  # python library. Let's add it to the list of dependencies by
  # inspecting "_ost_base.so".
  pymod_dir='lib/python%d.%d/site-packages' % sys.version_info[0:2]
  _deps_for_lib(os.path.join(stage_dir, pymod_dir, 'ost/_ost_base.so'),
                deps, recursive=False)
  print('copying dependencies')
  copy_deps(deps, outdir)
  print('copying libexec binaries')
  copy_binaries(stage_dir, outdir, libexec_binaries, libexec_scripts,
                'libexec/openstructure')
  print('copying binaries')
  copy_binaries(stage_dir, outdir, binaries, scripts, 'bin')
  print('copying pymod')
  shutil.copytree(os.path.join(stage_dir,pymod_dir),
                  os.path.join(outdir, pymod_dir))
  copied_py_framework=False
  non_std_python=False
  if os.path.exists(os.path.join(outdir, 'lib/Python.framework')):
    framework_path=os.path.join(outdir, 'lib/Python.framework')
    nonstd_python=True
    copied_py_framework=True
  if len(glob.glob(os.path.join(outdir, 'lib', 'libpython*')))>0:
    non_std_python=True
  if non_std_python:
    print('looks like we are using a non-standard python.')
    python_home=get_python_home()    
    if not copied_py_framework:
      print('also copying python modules from %s' % python_home)
      modules_dst=os.path.join(outdir, 'lib', os.path.basename(python_home))
      shutil.copytree(python_home, modules_dst)
      if os.path.exists(os.path.join(modules_dst, 'site-packages')):
        shutil.rmtree(os.path.join(modules_dst, 'site-packages'))
      copy_binaries(os.path.join(python_home, '../..'), outdir, 
                    ['python'], [], 'bin')
      python_bin=os.path.abspath(os.path.join(python_home, '../../bin/python'))
    else:
      # For MacPorts it's even more involved. Python is not executed directly 
      # but rather uses a wrapper executable that calls the actual python exe.
      # We have to include that one into the bundle.
      if macports_workaround:
        path_to_app='../../Resources/Python.app/Contents/MacOS/'
        exe_path=os.path.join(python_home, path_to_app)
        copy_binaries(exe_path, outdir, ['python'], [], 
                      append_bin=False)
        python_bin=os.path.join('/opt/local/bin/python')
      else:
        copy_binaries(os.path.join(python_home, '../..'), outdir, 
                      ['python'], [])
        python_bin=os.path.abspath(os.path.join(python_home, '../../bin/python'))
      # remove all versions but the one we are using
      version_string=sys.version[0:3]
      prefix, postfix=split_framework_components(python_home)
      site_packages_dir=os.path.join(outdir, 'lib', 'Python.framework', 
                                 postfix, 'site-packages')
      if os.path.exists(site_packages_dir):
        shutil.rmtree(site_packages_dir)
      for directory in glob.glob(os.path.join(framework_path, 'Versions/*')):
        if os.path.basename(directory)!=version_string:
          shutil.rmtree(directory)
    # replace the python executable
    ost_script=os.path.join(outdir, 'libexec', 'openstructure', 'ost_config')
    os.chmod(ost_script, 0o666)
    script=''.join(open(ost_script, 'r').readlines())
    script=script.replace(python_bin, '$BIN_DIR/python')
    open(ost_script, 'w').write(script)
    os.chmod(ost_script, 0o555)

  if no_includes:
    os.system(REMOVE_HEADERS % outdir)
    os.system(REMOVE_CURRENT % outdir)  
  print('copying site-packages')
  for sp in SITE_PACKAGES:
    src=get_python_module_path(sp)
    if os.path.isdir(src):
      shutil.copytree(src, os.path.join(outdir, pymod_dir, sp))
    else:
      shutil.copy(src, os.path.join(outdir, pymod_dir, sp))
  print('updating link commands of python shared objects')
  os.path.walk(os.path.join(outdir, 'lib'), 
               update_pymod_shared_objects, 
               os.path.join(outdir, 'lib'))
