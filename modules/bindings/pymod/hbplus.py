"""
Interface to HB+ command line program

Author: Marco Biasini
"""
import tempfile
from ost import settings
import subprocess
import re
import os
import shutil
from ost import io
from ost import mol

def _LocateHBPlus(hbplus_bin):
  return settings.Locate('hbplus', explicit_file_name=hbplus_bin,
                         env_name="HBPLUS")

class HBond:
  def __init__(self, donor, acceptor):
    self.donor=donor
    self.acceptor=acceptor

def _ParseOutput(ent, output):
  hbonds=[]
  # skip header
  lines=list(output)
  for index, line in enumerate(lines):
    if line.startswith('<---DONOR---> <-ACCEPTOR-->'):
      lines=lines[index+4:]
      break
  for line in lines:
    if len(line.strip())==0:
      continue
    don_chain=line[0]
    don_rnum=int(line[1:5])
    don_ins_c=line[5]
    don_atom=line[9:13]
    if don_chain=='-':
      don_chain=' '
    if don_ins_c=='-':
      don_ins_c='\0'
    acc_chain=line[14]
    acc_rnum=int(line[15:19])
    acc_ins_c=line[19]
    acc_atom=line[24:28]
    if acc_chain=='-':
      acc_chain=' '
    if acc_ins_c=='-':
      acc_ins_c='\0'
    donor=ent.FindAtom(don_chain, mol.ResNum(don_rnum, don_ins_c),
                       don_atom.strip())
    acc=ent.FindAtom(acc_chain, mol.ResNum(acc_rnum, acc_ins_c),
                     acc_atom.strip())                       
    assert donor.IsValid()
    assert acc.IsValid()
    hbonds.append(HBond(donor, acc))
  return hbonds

def HBondList(ent, hbplus_bin=None):
  """
  returns a list of HBonds found in the given entity (handle or view)
  """
  full_bin=_LocateHBPlus(hbplus_bin)
  temp_d=tempfile.mkdtemp(prefix='hbplus_')
  file_name=os.path.join(temp_d, 'ent.pdb')
  io.SaveEntity(ent, file_name)
  hb_proc=subprocess.Popen(full_bin, shell=True, stdout=subprocess.PIPE,
                           stdin=subprocess.PIPE)
  hb_proc.stdin.write(('%s\n' % temp_d).encode())
  hb_proc.stdin.write(('%s\n\n' % file_name).encode())
  stdout,_ = hb_proc.communicate()

  for line in stdout.decode().splitlines():
    match=re.match(r'Configured for (\d+) atoms and\s+(\d+) residues\.', line)
    if match:
      assert ent.atom_count<int(match.group(1))
      assert ent.residue_count<int(match.group(2))
  hb_out=open(os.path.join(temp_d, 'ent.hb2'), 'r')
  hbonds=_ParseOutput(ent, hb_out)
  hb_out.close()
  shutil.rmtree(temp_d)
  return hbonds

def HBondScore(ent1, ent2, hbplus_bin=None):
  """
  returns percentage of hydrogen bonds in ent1 that are also present in ent2 in
  the range of 0 to 1.
  
  this function is slow as hell, and could be improved drastically by first
  sorting the hydrogen bond lists by donor residue number.
  """

  hbonds_a=HBondList(ent1, hbplus_bin=hbplus_bin)
  hbonds_b=HBondList(ent2, hbplus_bin=hbplus_bin)
  matching=0
  for a in hbonds_a:
    found=False
    for b in hbonds_b:
      acc_names_match=b.acceptor.name==a.acceptor.name
      don_names_match=b.donor.name==a.donor.name
      don_num=b.donor.residue.number==a.donor.residue.number
      acc_num=b.acceptor.residue.number==a.acceptor.residue.number
      if acc_names_match and don_names_match and acc_num and don_num:
        found=True
        break
    if found:
      matching+=1
  if len(hbonds_a)>0:
    return float(matching)/len(hbonds_a)
  else:
    return 0.0

if __name__=='__main__':
  print('HBond Score:', HBondScore(io.LoadPDB(sys.argv[1]), 
                                   io.LoadPDB(sys.argv[2])))