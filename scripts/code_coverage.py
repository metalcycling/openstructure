"""
Code coverage report generator for Ost
"""

import optparse
import glob
import os, os.path

def Which(program_name):
  """
  emulate unit which program
  """
  for path in os.environ['PATH'].split(os.pathsep):
    path=os.path.expanduser(os.path.expandvars(path))
    full=os.path.join(path, program_name)
    if os.path.exists(full) and os.access(full, os.X_OK):
      return full
  return None

def SubDirsOf(directory):
  directories=[]
  for d in glob.glob(os.path.join(directory, '*')):
      if os.path.isdir(d):
        directories.append(os.path.basename(d))
  return directories
  
def ListModules():
  return ','.join(SubDirsOf('modules'))

def ParseOptions():
  parser=optparse.OptionParser()
  parser.add_option('-m', '--modules', dest='modules',
                    help='modules to test, separated by \',\'')
  parser.add_option('-o', '--output-directory', dest='output_dir',
                    type='string', help='output directory')

  parser.set_defaults(modules=ListModules())
  parser.set_defaults(output_dir='LCOV')
  (options, args)=parser.parse_args()
  options.modules=options.modules.split(',')
  return options

def CheckLCOV():
  if not Which('lcov'):
    print('please install lcov to run this script')
    sys.exit(-1)

CheckLCOV()

opts=ParseOptions()

def Print(message):
  print('~~~~~~')
  print('|  *  |', message)
  print('~~~~~~')  
  
def RemoveFiles(directory, extension):
  glob_pattern='*.%s' % extension
  for f in glob.glob(os.path.join(directory, glob_pattern)):
    os.unlink(f)

def Cleanup(modules):
  for module in modules:
    sub_dirs=SubDirsOf(os.path.join('modules', module))
    for sub_dir in sub_dirs:
      RemoveFiles(os.path.join('modules',module, sub_dir), 'gcda')

def RunLCOV(module):
  Print('Running Tests for %s' % module)  
  sub_dirs=SubDirsOf(os.path.join('modules', module))
  for sub_dir in sub_dirs:
    if sub_dir=='pymod':
      continue
    full_dir=os.path.join('modules', module, sub_dir)
    out='lcov-tmp/%s_%s.info' % (module, sub_dir)
    cmd='lcov --capture --b . --directory "%s" -o "%s"' % (full_dir, out)
    os.system(cmd+' > /dev/null')
  files_to_merge=glob.glob('lcov-tmp/%s_*.info' % module)
  output_name='lcov-tmp/%s_merge.info' % module  

  merge_cmd='lcov %s -o %s' % (' '.join(['-a %s' % f for f in files_to_merge]),
                               output_name)
  os.system(merge_cmd+' > /dev/null')
  
  libs_to_remove=['"*/boost/*"', '"/usr/*"', '"/opt*"']
  remove_cmd='lcov -r %s %s -o %s ' % (output_name, ' '.join(libs_to_remove),
                                       output_name)
  os.system(remove_cmd+' > /dev/null')
  

def MergeAll(modules):
  filtered_modules=[]
  for f in modules:
    if os.access('lcov-tmp/%s_merge.info' % f, os.F_OK):
      filtered_modules.append(f)
  
  f=' '.join(['-a lcov-tmp/%s_merge.info' % f for f in filtered_modules])
  merge_cmd='lcov %s -o %s' % (f, 'lcov-tmp/merge.info')
  os.system(merge_cmd+' > /dev/null')

def GenReport(output_dir):
  os.system('rm -rf LCOV')
  cmd='genhtml -o %s lcov-tmp/merge.info' % (output_dir)
  os.system(cmd+'> /dev/null')

def RunTests(modules):
  try:
    os.mkdir('lcov-tmp')
  except:
    pass
  for module in modules:
    test_dir=os.path.join('modules', module, 'tests')
    if os.path.exists(test_dir):
      test_binary=os.path.join(os.path.join(test_dir, 'tests'))
      if os.path.exists(test_binary) and os.access(test_binary, os.X_OK):
        print('running tests for module', module)        
        Cleanup(modules)
        cmd='cd "%s"; ./tests' % test_dir
        os.system(cmd)
        RunLCOV(module)
  MergeAll(modules)
  
Print('recompiling project with --test-coverage, if neccessary')
os.system('scons --test-coverage unit_tests > /dev/null')  
RunTests(opts.modules)
Print('Generating report in HTML format')
GenReport(opts.output_dir)
