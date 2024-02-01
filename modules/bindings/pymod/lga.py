"""
Interface to LGA

Author: Marco Biasini
"""
import tempfile
import os
from ost import io
import subprocess
from ost import geom
from ost import settings

def _FindLGABinary(lga_bin):
  return settings.Locate('lga', explicit_file_name=lga_bin, 
                         env_name='LGA_BINARY')

def _PrepareInput(pdb1, pdb2, output_pdb):
  """
  Deal with the brain-dead input requirements of LGA.
  """
  mol1=os.path.basename(os.path.splitext(pdb1)[0])
  mol2=os.path.basename(os.path.splitext(pdb2)[0])  
  os.system('echo MOLECULE %s > %s' % (mol1, output_pdb))
  os.system('cat %s >> %s' % (pdb1, output_pdb))
  os.system('echo "\n" >> %s' % (output_pdb))  
  os.system('echo MOLECULE %s >> %s' % (mol2, output_pdb))
  os.system('cat %s >> %s' % (pdb2, output_pdb))  

class GDTResult:
  def __init__(self, rotation, shift, gdt_ts, gdt_ha):
    self.rotation=rotation
    self.shift=shift
    self.gdt_ts=gdt_ts
    self.gdt_ha=gdt_ha
  def GetTransform(self):
    transform=geom.Mat4()
    transform.PasteTranslation(self.shift)
    return transform*geom.Mat4(self.rotation)

def _ParseRotationAndShift(lines):
  t=[l.split() for l in lines]
  rot=geom.Mat3()
  shift=geom.Vec3()
  for i, x in enumerate(t):
    rot[(i, 0)]=+float(x[2])
    rot[(i, 1)]=+float(x[6])
    rot[(i, 2)]=+float(x[10])
    shift[i]=float(x[14])
  return rot, shift

def _ParseGDTSection(section, residue_count):
  cutoffs=[float(e) for e in section[0].split()[2:]]
  num_ca=[int(e) for e in section[1].split()[2:]]
  gdt_ts=[float(e) for e in section[2].split()[2:]]
  scores=dict(list(zip(cutoffs, gdt_ts)))
  numbers=dict(list(zip(cutoffs, num_ca)))
  factor=(1.0/(4*residue_count))*100
  ts_cutoffs=(1.0, 2.0, 4.0, 8.0)
  ha_cutoffs=(0.5, 1.0, 2.0, 4.0)  
  gdt_ts=(sum([numbers[c] for c in ts_cutoffs]))*factor
  gdt_ha=(sum([numbers[c] for c in ha_cutoffs]))*factor  
  return gdt_ts, gdt_ha 

def _ParseLGAOutput(output, residue_count):
  result=GDTResult(geom.Mat3(), geom.Vec3(), 0.0, 0.0)
  found_gdt_section=False
  found_transform_section=False
  for index, line in enumerate(output):
    if line.startswith('GLOBAL_DISTANCE_TEST'): 
      next_lines=output[index+1:index+5]
      result.gdt_ts, result.gdt_ha=_ParseGDTSection(next_lines, residue_count)
      found_gdt_section=True      
    if line.startswith('Unitary ROTATION matrix'):
      next_lines=output[index+1:index+4]
      result.rotation, result.shift=_ParseRotationAndShift(next_lines)
      found_transform_section=True
      break
  assert found_transform_section and found_gdt_section
  return result
  
def GDT(pdb1, pdb2, chain1='', chain2='', reference_length=None, lga_bin=None):
  """
  Calculate GDT value between pdb1 and pdb2. It is assumed that the
  corresponding residues in pdb1 and pdb2 have the same residue numbers.
  """
  lga_bin=_FindLGABinary(lga_bin)  
  temp_d=tempfile.mkdtemp(prefix='lga_gdt_ts_')
  for d in ('MOL2', 'TMP', 'RESULTS'):
    os.mkdir(os.path.join(temp_d, d))

  if not chain1:
    chain1=pdb1.chains[0].name
  if not chain2:
    chain2=pdb2.chains[0].name
  pdb_one_name=os.path.join(temp_d, 'MOL2', 'one.pdb')
  pdb_two_name=os.path.join(temp_d, 'MOL2', 'two.pdb')
  io.SaveEntity(pdb1, pdb_one_name)
  io.SaveEntity(pdb2, pdb_two_name)
  _PrepareInput(pdb_one_name, pdb_two_name, 
                os.path.join(temp_d, 'MOL2', 'input.pdb'))
  output_file=os.path.join(temp_d, 'out.txt')
  ch1, ch2=('', '')
  if len(chain1.strip()): ch1=' -ch1:%s ' % chain1
  if len(chain2.strip()): ch2=' -ch2:%s ' % chain2
    
  params=(temp_d, lga_bin, 'input.pdb', ch1, ch2)  
  command='cd %s; %s %s %s%s-ie -3 -d:4 -sda'

  expanded_cmd=command % params
  lga_proc=subprocess.Popen(expanded_cmd, shell=True, 
                            stdout=subprocess.PIPE)
  stdout, _ = lga_proc.communicate()

  length=reference_length or max(pdb1.residue_count, pdb2.residue_count)
  result=_ParseLGAOutput(stdout.decode().splitlines(), reference_length)
  os.system('rm -r %s' % temp_d)
  return result
