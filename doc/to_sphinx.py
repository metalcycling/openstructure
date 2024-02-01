import re
to_document=sys.argv[1]
parts=to_document.split('.')
module='.'.join(parts[:-1])
class_name=parts[-1]
__import__(module, globals(), locals())
mod=sys.modules[module]
the_class=getattr(mod, class_name)
print('.. currentmodule:: %s' % module)
print('')
print('.. class:: %s' % class_name)

def _arg_type(arg_type):
  if arg_type not in ('str', 'int', 'bool', 'float', 'None'):
    return ':class:`%s`' % arg_type
  return arg_type
class TypedArgument:
  def __init__(self, name, arg_type):
    self.name=name
    self.type=_arg_type(arg_type)

class Method:
  def __init__(self, name, rtype, args, optional):
    self.name=name
    self.rtype=rtype
    self.args=args
    self.opt=optional
  def to_sphinx(self):
    if len(self.opt)>0:
      signature='%s(%s%s)' % (self.name, ', '.join([a.name for a in self.args]), 
                            '[, %s]' % ', '.join([a.name for a in self.opt]))
    else:
      signature='%s(%s)' % (self.name, ', '.join([a.name for a in self.args]))
    val='  .. method:: %s\n\n' % signature
    for arg in self.args+self.opt:
      val+='    :param %s:\n' % arg.name
      val+='    :type  %s: %s\n' % (arg.name, arg.type)
    val+='    :rtype: %s\n' % self.rtype
    return val

def parse_signature(signature):
  arg=re.compile(r'\s*(\((?P<type>\w+)\)(?P<name>\w+),?\s*)*')
  opt_args=re.compile(r'(?P<args>.*)(:?\[,(?P<kwargs>.*)\])')
  method=re.compile(r'(?P<name>\w+)\((?P<args>.*)\) -> (?P<rtype>\w+)')
  method_match=method.match(signature)
  if method_match:
    optional=opt_args.match(method_match.groupdict()['args'])
    opt_args=[]
    if optional:
      for opt_arg in optional.groupdict()['kwargs'].split(','):
        a=arg.match(opt_arg)
        assert a
        opt_args.append(TypedArgument(a.groupdict()['name'], a.groupdict()['type']))
      arg_string=optional.groupdict()['args']
    else:
      arg_string=method_match.groupdict()['args']
    args=[]
    for aarg in arg_string.split(','):
      a=arg.match(aarg)
      assert a
      args.append(TypedArgument(a.groupdict()['name'], a.groupdict()['type']))  
    
    return Method(method_match.groupdict()['name'], 
                  _arg_type(method_match.groupdict()['rtype']), 
                  args[1:], opt_args)
  print(signature, 'not matched')

if '--no-derived' in sys.argv:
  members=list(the_class.__dict__.keys())
else:
  members=dir(the_class)
for m in members:
  if m.startswith('__'):
    continue
  member_doc=getattr(the_class, m).__doc__
  if member_doc:
    method=parse_signature(member_doc.split('\n')[1])
    print(method.to_sphinx())