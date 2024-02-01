import os, platform
import __main__
 
def GetPlatform():
  """
  Returns platform.system(). If the system call is interrupted, it is repeated.
  This function is a workaround for the buggy system call handling in Python.
  """
  system=None
  while not system:
    try:
      system=platform.system()
    except IOError:
      pass
  return system

def GetValue(val_key,val_default=None,prefix='OST'):
  """
  Returns the value of the variable val_key if defined, otherwise returns the 
  default value provided by the user (if provided). Search order: 
  
   * environment variable called $prefix_$val_key 
   * variable called val_key in .ostrc file
  """
  if prefix:
    env_var_name='%s_%s' % (prefix, val_key)
  else:
    env_var_name=val_key
  env_value=os.getenv(env_var_name)
  if env_value:
    return env_value
  else:
    main_attr=getattr(__main__, val_key, None)
    if main_attr:
      return main_attr
    else:
      return val_default

class FileNotFound(RuntimeError):
  """
  Raised when :func:`Locate` is unable to locate a file. The exception contains
  detailed information on what was tried to locate the file, i.e. search paths, 
  environment variables and also provides useful hints on how to let Locate know
  where to find the file.
  """
  def __init__(self, name, reason):
    self.name=name
    self.reason=reason
  def __str__(self):
    return 'Could not find "%s": %s' % (self.name, self.reason)

def Locate(file_name, explicit_file_name=None, search_paths=[],
           env_name=None, search_system_paths=True):
  """
  Helper function to locate files. To get the full name of an executable, let's 
  say qmake, use
  
  .. code-block:: python

    abs_qmake_path=Locate('qmake', env_name='QMAKE_EXECUTABLE')

  First the function checks if an environment variable with the name 
  QMAKE_EXECUTABLE is set. If so, the value of this variable is returned. Next, 
  each directory listed in search_paths is searched. If the executable could 
  still not be found and search_system_paths is set to True, the binary search 
  paths are searched.

  If the file could not be located, a :exc:`~ost.settings.FileNotFound` 
  exception will be raised containing a detail description why Locate failed. The 
  error message is formatted in such a way that it can directly be presented to 
  the user.
  """
  def _is_executable(filename):
    return os.path.exists(filename) and os.access(filename, os.X_OK)
  if type(file_name) is str:
    file_names=[file_name]
  else:
    file_names=file_name
  env_var_inexistent='env variable %s points to inexistent file %s'
  epxl_inexistent='explicitly set file "%s" does not exist'
  set_env_var='set the environment variable %s to the absolute path to %s or '
  if explicit_file_name:
    if _is_executable(explicit_file_name):
      return explicit_file_name
    else:
      raise FileNotFound(file_name, epxl_inexistent % explicit_file_name)
  if env_name:
    file_env_name=os.getenv(env_name, None)
    if file_env_name:
      if _is_executable(file_env_name):
        return file_env_name
      else:
        raise FileNotFound(file_name, 
                           env_var_inexistent % (env_name, file_env_name))
  searched=list(search_paths)
  for search_path in search_paths:
    for file_name in file_names:
      full_file_name=os.path.join(search_path, file_name)
      if _is_executable(full_file_name):
        return full_file_name

  if search_system_paths:
    paths=os.getenv('PATH')
    if GetPlatform() == "Windows":
      searched+=paths.split(';')
    else:
      searched+=paths.split(':')
    for path in searched:
      for file_name in file_names:
        full_file_name=os.path.join(path, file_name)
        if _is_executable(full_file_name):
          return full_file_name
  msg=''        
  if len(searched)>0:
    msg='searched in \n%s\n' % ( '\n'.join([' - %s' % s for s in searched]))
  if env_name:
    msg+=set_env_var % (env_name, ', ' % file_names)
  msg+='put %s into one of the search paths' % ', '.join(file_names)
  raise FileNotFound(file_name, msg)
