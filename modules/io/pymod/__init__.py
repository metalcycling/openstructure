#------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------
import os, tempfile, ftplib, http.client

from ._ost_io import *
from ost import mol, geom, conop, seq

profiles=None

class IOProfiles:
  def __init__(self):
     self._dict={}

  def __getitem__(self, key):
    return IOProfileRegistry.Instance().Get(key)

  def __setitem__(self, key, value):
    if isinstance(value, str):
      value=self[value].Copy()
    IOProfileRegistry.Instance().Set(key, value)
    self._dict[key]=value

  def __len__(self):
    return len(self._dict)

  def __iter__(self):
    return self._dict.__iter__()

if not profiles:
  profiles=IOProfiles()
  if conop.GetDefaultLib():
    processor = conop.RuleBasedProcessor(conop.GetDefaultLib())
  else:
    processor = conop.HeuristicProcessor()
  profiles['STRICT']=IOProfile(dialect='PDB', fault_tolerant=False,
                               quack_mode=False, processor=processor.Copy())
  profiles['SLOPPY']=IOProfile(dialect='PDB', fault_tolerant=True,
                               quack_mode=True, processor=processor.Copy())
  profiles['CHARMM']=IOProfile(dialect='CHARMM', fault_tolerant=True,
                               quack_mode=False, processor=processor.Copy())
  profiles['DEFAULT']='STRICT'

def _override(val1, val2):
  if val2!=None:
    return val2
  else:
    return val1

def LoadPDB(filename, restrict_chains="", no_hetatms=None,
            fault_tolerant=None, load_multi=False, quack_mode=None,
            join_spread_atom_records=None, calpha_only=None,
            profile='DEFAULT', remote=False, remote_repo='pdb',
            dialect=None, seqres=False, bond_feasibility_check=None):
  """
  Load PDB file from disk and return one or more entities. Several options 
  allow to customize the exact behaviour of the PDB import. For more information 
  on these options, see :doc:`profile`.
  
  Residues are flagged as ligand if they are mentioned in a HET record.

  :param restrict_chains: If not an empty string, only chains listed in the
     string will be imported.

  :param fault_tolerant: Enable/disable fault-tolerant import. If set, overrides 
     the value of :attr:`IOProfile.fault_tolerant`.

  :param no_hetatms: If set to True, HETATM records will be ignored. Overrides 
      the value of :attr:`IOProfile.no_hetatms`

  :param load_multi: If set to True, a list of entities will be returned instead
      of only the first. This is useful when dealing with multi-PDB files.

  :param join_spread_atom_records: If set, overrides the value of 
      :attr:`IOProfile.join_spread_atom_records`.
  
  :param remote: If set to True, the method tries to load the pdb from the 
     remote repository given as *remote_repo*. The filename is then 
     interpreted as the entry id as further specified for the *remote_repo*
     parameter.

  :param remote_repo: Remote repository to fetch structure if *remote* is True. 
                      Must be one in ['pdb', 'smtl', 'pdb_redo']. In case of 
                      'pdb' and 'pdb_redo', the entry must be given as lower 
                      case pdb id, which loads the deposited assymetric unit 
                      (e.g. '1ake'). In case of 'smtl', the entry must also 
                      specify the desired biounit (e.g. '1ake.1').
     
  :rtype: :class:`~ost.mol.EntityHandle` or a list thereof if `load_multi` is 
      True.

  :param dialect: Specifies the particular dialect to use. If set, overrides 
    the value of :attr:`IOProfile.dialect`

  :param seqres: Whether to read SEQRES records. If set to True, the loaded 
    entity and seqres entry will be returned as a tuple.

  :type dialect: :class:`str`
  

  :raises: :exc:`~ost.io.IOException` if the import fails due to an erroneous or 
      inexistent file
  """
  def _override(val1, val2):
    if val2!=None:
      return val2
    else:
      return val1
  if isinstance(profile, str):
    prof=profiles[profile].Copy()
  elif isinstance(profile, IOProfile):
    prof=profile.Copy()
  else:
    raise TypeError('profile must be of type string or IOProfile, '+\
                    'instead of %s'%type(profile))
  if dialect not in (None, 'PDB', 'CHARMM',):
    raise ValueError('dialect must be PDB or CHARMM')
  prof.calpha_only=_override(prof.calpha_only, calpha_only)
  prof.no_hetatms=_override(prof.no_hetatms, no_hetatms)
  prof.dialect=_override(prof.dialect, dialect)
  prof.quack_mode=_override(prof.quack_mode, quack_mode)
  if prof.processor:
    prof.processor.check_bond_feasibility=_override(prof.processor.check_bond_feasibility, 
                                                    bond_feasibility_check)
  prof.fault_tolerant=_override(prof.fault_tolerant, fault_tolerant)
  prof.join_spread_atom_records=_override(prof.join_spread_atom_records,
                                          join_spread_atom_records)

  tmp_file = None # avoid getting out of scope
  if remote:
    if remote_repo not in ['pdb', 'smtl', 'pdb_redo']:
      raise IOError("remote_repo must be in ['pdb', 'smtl', 'pdb_redo']")
    from ost.io.remote import RemoteGet
    tmp_file =RemoteGet(filename, from_repo=remote_repo)
    filename = tmp_file.name
  
  conop_inst=conop.Conopology.Instance()
  if prof.processor:
    if prof.dialect=='PDB':
      prof.processor.dialect=conop.PDB_DIALECT
    elif prof.dialect=='CHARMM':
      prof.processor.dialect=conop.CHARMM_DIALECT
  reader=PDBReader(filename, prof)
  reader.read_seqres=seqres
  try:
    if load_multi:
      ent_list=[]
      while reader.HasNext():
        ent=mol.CreateEntity()
        reader.Import(ent, restrict_chains)
        if prof.processor:
          prof.processor.Process(ent)
        ent_list.append(ent)
      if len(ent_list)==0:
        raise IOError("File '%s' doesn't contain any entities" % filename)
      return ent_list
    else:
      ent=mol.CreateEntity()
      if reader.HasNext():
        reader.Import(ent, restrict_chains)
        if prof.processor:
          prof.processor.Process(ent)
      else:
        raise IOError("File '%s' doesn't contain any entities" % filename)
      if seqres:
        return ent, reader.seqres
      return ent
  except:
    raise

def SavePDB(models, filename, dialect=None,  pqr=False, profile='DEFAULT'):
  """
  Save entity or list of entities to disk. If a list of entities is supplied
  the PDB file will be saved as a multi PDB file. Each of the entities is
  wrapped into a MODEL/ENDMDL pair.

  If the atom number exceeds 99999, '*****' is used.

  :param models: The entity or list of entities (handles or views) to be saved
  :param filename: The filename
  :type  filename: string
  :raises: IOException if the restrictions of the PDB format are not satisfied
           (with the exception of atom numbers, see above):

             * Chain names with more than one character
             * Atom positions with coordinates outside range [-999.99, 9999.99]
             * Residue names longer than three characters
             * Atom names longer than four characters
             * Numeric part of :class:`ost.mol.ResNum` outside range [-999, 9999] 
             * Alternative atom indicators longer than one character
  """
  if not getattr(models, '__len__', None):
    models=[models]
  if isinstance(profile, str):
    profile=profiles[profile].Copy()
  elif isinstance(profile, IOProfile):
    profile.Copy()
  else:
    raise TypeError('profile must be of type string or IOProfile, '+\
                    'instead of %s'%type(profile))
  profile.dialect=_override(profile.dialect, dialect)
  writer=PDBWriter(filename, profile)
  writer.SetIsPQR(pqr)
  if len(models)>1:
    writer.write_multi_model=True
  for model in models:
    writer.Write(model)

try:
  from ost import img
  LoadMap = LoadImage
  SaveMap = SaveImage
except ImportError:
  pass

 ## loads several images and puts them in an ImageList
 # \sa \example fft_li.py "View Fourier Transform Example"
def LoadImageList (files):
  image_list=img.ImageList()
  for file in files:
    image=LoadImage(file)
    image_list.append(image)
  return image_list

LoadMapList=LoadImageList

def LoadCHARMMTraj(crd, dcd_file=None, profile='CHARMM',
                   lazy_load=False, stride=1, 
                   dialect=None, detect_swap=True,swap_bytes=False):
  """
  Load CHARMM trajectory file.
  
  :param crd: EntityHandle or filename of the (PDB) file containing the
      structure. The structure must have the same number of atoms as the 
      trajectory
  :param dcd_file: The filename of the DCD file. If not set, and crd is a 
      string, the filename is set to the <crd>.dcd
  :param layz_load: Whether the trajectory should be loaded on demand. Instead 
      of loading the complete trajectory into memory, the trajectory frames are 
      loaded from disk when requested.
  :param stride: The spacing of the frames to load. When set to 2, for example, 
      every second frame is loaded from the trajectory. By default, every frame 
      is loaded.
  :param dialect: The dialect for the PDB file to use. See :func:`LoadPDB`. If 
      set, overrides the value of the profile
  :param profile: The IO profile to use for loading the PDB file. See 
      :doc:`profile`.
  :param detect_swap: if True (the default), then automatic detection of endianess
      is attempted, otherwise the swap_bytes parameter is used
  :param swap_bytes: is detect_swap is False, this flag determines whether bytes
      are swapped upon loading or not
  """
  if not isinstance(crd, mol.EntityHandle):
    if dcd_file==None:
      dcd_file='%s.dcd' % os.path.splitext(crd)[0]    
    crd=LoadPDB(crd, profile=profile, dialect=dialect)

  else:
    if not dcd_file:
      raise ValueError("No DCD filename given")
  return LoadCHARMMTraj_(crd, dcd_file, stride, lazy_load, detect_swap, swap_bytes)

def LoadMMCIF(filename, fault_tolerant=None, calpha_only=None, profile='DEFAULT', remote=False, seqres=False, info=False):
  """
  Load a mmCIF file and return one or more entities. Several options allow to
  customize the exact behaviour of the mmCIF import. For more information on
  these options, see :doc:`profile`.
  
  Residues are flagged as ligand if they are mentioned in a HET record.

  :param fault_tolerant: Enable/disable fault-tolerant import. If set, overrides
     the value of :attr:`IOProfile.fault_tolerant`.

  :param remote: If set to True, the method tries to load the pdb from the 
     remote pdb repository www.pdb.org. The filename is then interpreted as the 
     pdb id.

  :rtype: :class:`~ost.mol.EntityHandle` (or tuple if *seqres* or *info* are
          True).

  :param seqres: Whether to read SEQRES records. If True, a
                 :class:`~ost.seq.SequenceList` object is returned as the second
                 item. The sequences in the list are named according to the
                 mmCIF chain name.
                 This feature requires a default
                 :class:`compound library <ost.conop.CompoundLib>`
                 to be defined and accessible via
                 :func:`~ost.conop.GetDefaultLib` or an empty list is returned.

  :param info: Whether to return an info container with the other output.
               If True, a :class:`MMCifInfo` object is returned as last item.

  :raises: :exc:`~ost.io.IOException` if the import fails due to an erroneous
           or non-existent file.
  """
  def _override(val1, val2):
    if val2!=None:
      return val2
    else:
      return val1
  if isinstance(profile, str):
    prof = profiles[profile].Copy()
  else:
    prof = profile.Copy()

  prof.calpha_only=_override(prof.calpha_only, calpha_only)
  prof.fault_tolerant=_override(prof.fault_tolerant, fault_tolerant)

  if remote:
    from ost.io.remote import RemoteGet
    tmp_file = RemoteGet(filename, from_repo='cif')
    filename = tmp_file.name
  
  try:
    ent = mol.CreateEntity()
    reader = MMCifReader(filename, ent, prof)
    reader.read_seqres = seqres
    
    # NOTE: to speed up things, we could introduce a restrict_chains parameter
    #       similar to the one in LoadPDB. Here, it would have to be a list/set
    #       of chain-name-strings.

    #if reader.HasNext():
    reader.Parse()
    if prof.processor:
      prof.processor.Process(ent)
      reader.info.ConnectBranchLinks()
    #else:
    #  raise IOError("File doesn't contain any entities")
    if seqres and info:
      return ent, reader.seqres, reader.info
    if seqres:
      return ent, reader.seqres
    if info:
      return ent, reader.info
    return ent
  except:
    raise

# this function uses a dirty trick: should be a member of MMCifInfoBioUnit
# which is totally C++, but we want the method in Python... so we define it
# here (__init__) and add it as a member to the class. With this, the first
# arguement is the usual 'self'.
# documentation for this function was moved to mmcif.rst,
# MMCifInfoBioUnit.PDBize, since this function is not included in SPHINX.
def _PDBize(biounit, asu, seqres=None, min_polymer_size=None,
            transformation=False, peptide_min_size=10, nucleicacid_min_size=10,
            saccharide_min_size=10):
  if min_polymer_size is not None:
    pdbizer = mol.alg.PDBize(min_polymer_size=min_polymer_size)
  else:
    pdbizer = mol.alg.PDBize(peptide_min_size=peptide_min_size,
                             nucleicacid_min_size=nucleicacid_min_size,
                             saccharide_min_size=saccharide_min_size)

  chains = biounit.GetChainList()
  c_intvls = biounit.GetChainIntervalList()
  o_intvls = biounit.GetOperationsIntervalList()
  ss = seqres
  if not ss:
    ss = seq.CreateSequenceList()
  # create list of operations
  # for cartesian products, operations are stored in a list, multiplied with
  # the next list of operations and re-stored... until all lists of operations
  # are multiplied in an all-against-all manner.
  operations = biounit.GetOperations()
  for i in range(0,len(c_intvls)):
    trans_matrices = geom.Mat4List()
    l_operations = operations[o_intvls[i][0]:o_intvls[i][1]]
    if len(l_operations) > 0:
      for op in l_operations[0]:
        rot = geom.Mat4()
        rot.PasteRotation(op.rotation)
        trans = geom.Mat4()
        trans.PasteTranslation(op.translation)
        tr = geom.Mat4()
        tr = trans * rot
        trans_matrices.append(tr)
      for op_n in range(1, len(l_operations)):
        tmp_ops = geom.Mat4List()
        for o in l_operations[op_n]:
          rot = geom.Mat4()
          rot.PasteRotation(o.rotation)
          trans = geom.Mat4()
          trans.PasteTranslation(o.translation)
          tr = geom.Mat4()
          tr = trans * rot
          for t_o in trans_matrices:
            tp = t_o * tr
            tmp_ops.append(tp)
        trans_matrices = tmp_ops
    # select chains into a view as basis for each transformation
    assu = asu.Select('cname='+','.join(mol.QueryQuoteName(name) \
                                        for name in \
                                        chains[c_intvls[i][0]:c_intvls[i][1]]))
    pdbizer.Add(assu, trans_matrices, ss)
  pdb_bu = pdbizer.Finish(transformation)
  if transformation:
    return pdb_bu, pdb_bu.GetTransform().GetMatrix()
  return pdb_bu

MMCifInfoBioUnit.PDBize = _PDBize
