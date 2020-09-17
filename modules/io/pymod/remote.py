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

import urllib.request, urllib.error, urllib.parse
import tempfile

from ost.io import LoadPDB, LoadMMCIF

class RemoteRepository:
  """
  A remote repository represents a structural database accessible through the 
  internet, e.g. the PDB or SWISS-MODEL template library.

  :param name:        Name of the repository
  :param url_pattern: URL pattern for repository. Required format is described
                      in :func:`URLForID`
  :param type:        Data format to expect at resolved URL must be in 
                      ('pdb', 'cif')
  :param id_transform: Transformation to apply to ID before resolving URL
                       in :func:`URLForID`. Must be in ('lower', 'upper')

  :type name:         :class:`str`
  :type url_pattern:  :class:`str`
  :type type:         :class:`str`
  :type id_transform: :class:`str`
  """
  def __init__(self, name, url_pattern, type, id_transform='upper'):
    self.name = name
    self.url_pattern = url_pattern
    self.type = type
    if type not in ('cif', 'pdb'):
      raise ValueError('only cif and pdb types are supported')
    self.id_transform = id_transform

  def URLForID(self, id):
    """
    Resolves URL given *url_pattern* and *id_transform* provided at object 
    initialization.
    The *url_pattern* must contain substring '$ID'. Given *id*, the URL to 
    the structure gets constructed by applying *id_transform* and inserting it
    at the location of '$ID'. e.g. 'https://files.rcsb.org/view/$ID.pdb' given 
    1ake as *id* and 'upper' as *id_transform* resolves to:
    'https://files.rcsb.org/view/1AKE.pdb'
    """
    if self.id_transform == 'upper':
      id = id.upper()
    if self.id_transform == 'lower':
      id = id.lower()
    return self.url_pattern.replace('$ID', id)

  def Get(self, id):
    """
    Resolves URL with :func:`URLForID`, dumps the content in a temporary file 
    and returns its path.
    
    :param id:          ID to resolve
    :type id:           :class:`str`
    """
    remote_url = self.URLForID(id)
    tmp_file_suffix = '.%s' % self.type
    if remote_url.endswith('.gz'):
      tmp_file_suffix+='.gz'

    try:
      connection = urllib.request.urlopen(remote_url)
      if hasattr(connection, 'code'):
        status = connection.code
      else:
        status = connection.getcode()
    except urllib.error.HTTPError as e:
      status = e.code
    if status != 200:
      raise IOError('Could not load %s from %s (status code %d, url %s)' \
                    % (id, self.name, status, remote_url))
    tmp_file = tempfile.NamedTemporaryFile(suffix=tmp_file_suffix)
    tmp_file.write(connection.read())
    tmp_file.flush()
    return tmp_file

  def Load(self, id):
    """
    Resolves URL with :func:`URLForID` and directly loads/returns the according 
    :class:`ost.mol.EntityHandle`. Loading invokes the 
    :func:`ost.io.LoadPDB`/:func:`ost.io.LoadMMCIF` with default parameterization. If you need
    custom settings, you might want to consider to call :func:`Get` and do the
    loading manually.
    
    :param id:          ID to resolve
    :type id:           :class:`str`
    """
    tmp_file = self.Get(id)
    if self.type == 'pdb':
      return LoadPDB(tmp_file.name)
    if self.type == 'cif':
      return LoadMMCIF(tmp_file.name)

REMOTE_REPOSITORIES = {
    'pdb' : RemoteRepository('pdb.org (PDB)', 'https://www.pdb.org/pdb/files/$ID.pdb.gz',
                   type='pdb', id_transform='upper'),
    'smtl' : RemoteRepository('SMTL', 'https://swissmodel.expasy.org/templates/$ID.pdb',
                   type='pdb', id_transform='lower'),
    'cif' : RemoteRepository('pdb.org (mmCIF)', 'https://www.pdb.org/pdb/files/$ID.cif.gz',
                   type='cif', id_transform='lower'),
    'pdb_redo' : RemoteRepository('pdbredo', 'https://pdb-redo.eu/db/$ID/$ID_besttls.pdb.gz',
                   type='pdb', id_transform='lower'),
}

def RemoteGet(id, from_repo='pdb'):
  """
  Invokes :func:`RemoteRepository.Get` on predefined repositories 
  ('pdb', 'smtl', 'cif', 'pdb_redo')

  :param from_repo:    One of the predefined repositories
  :type from_repo:     :class:`str`
  """
  remote_repo = REMOTE_REPOSITORIES.get(from_repo, None) 
  if not remote_repo:
    raise ValueError('%s is not a valid repository' % from_repo)
  return remote_repo.Get(id)

def RemoteLoad(id, from_repo='pdb'):
  """
  Invokes :func:`RemoteRepository.Load` on predefined repositories 
  ('pdb', 'smtl', 'cif', 'pdb_redo')

  :param from_repo:    One of the predefined repositories
  :type from_repo:     :class:`str`
  """
  remote_repo = REMOTE_REPOSITORIES.get(from_repo, None) 
  if not remote_repo:
    raise ValueError('%s is not a valid repository' % from_repo)
  return remote_repo.Load(id)
