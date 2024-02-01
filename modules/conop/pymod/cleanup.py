from ost import conop, mol
  
def Cleanup(entity, strip_water=True, canonicalize=True, remove_ligands=True):
  """
  This function returns a cleaned-up (simplified) version of the protein
  structure. Different parameters affect the behaviour of the function.

  :param strip_water: Whether to remove water from the structure
  :param canonicalize: Whether to strip off modifications of amino acids and map
     them back to their parent standard amino acid, e.g. selenium methionine to
     methionine.For more complex amino acids, where the relation between the
     modified and the standard parent amino acid is not known, sidechain atoms
     are removed. D-peptide-linking residues are completely removed as well.
  :param remove_ligands: Whether to remove ligands from the structure

  :return: a cleaned version of the entity
  """
  #setup
  lib = conop.GetDefaultLib()
  if not lib:
    raise RuntimeError("Cleanup requires a compound library.")
  clean_entity = entity.Copy()
  ed = clean_entity.EditXCS()
  #remove water residues
  if strip_water:
    _StripWater(clean_entity, ed)
  #replace modified residues before removing ligands to avoid removing MSE and others
  if canonicalize:
    _CanonicalizeResidues(clean_entity, ed, lib)
  #remove all hetatoms that are not water
  if remove_ligands:
    _RemoveLigands(clean_entity, ed)
  return clean_entity


def _StripWater(clean_entity, ed) :
  """
  This function removes water residues from the structure
  """
  for res in clean_entity.residues:
    if res.IsValid():
      if res.chem_class == mol.WATER:
        ed.DeleteResidue(res.handle)
  ed.UpdateICS()
  return
  
def _RemoveLigands(clean_entity, ed) :
  """
  This function removes ligands from the structure
  """
  for res in clean_entity.residues:
    if res.IsValid():
      #WHEN mmCIF WILL BE USED, CHANGE IsPeptideLinking() TO IsProtein()
      if not res.IsPeptideLinking() and res.atoms[0].is_hetatom and res.chem_class != mol.WATER:
        ed.DeleteResidue(res.handle)
  ed.UpdateICS()
  return
  
def _CanonicalizeResidues(clean_entity, ed, compound_lib) :
  """
  This function strips off modifications of amino acids and maps
  them back to their parent standard amino acid, e.g. selenium methionine to
  methionine.For more complex amino acids, where the relation between the 
  modified and the standard parent amino acid is not known, sidechain atoms 
  are removed. D-peptide-linking residues are completely removed as well.
  """

  for res in clean_entity.residues:
    if res.IsValid() and res.IsPeptideLinking() :
      parent_olc = res.one_letter_code
      if parent_olc == "X" :
        _DeleteSidechain(res, ed)
        for atom in res.atoms:
          atom.is_hetatom = False
      else:
        parent_tlc = conop.OneLetterCodeToResidueName(parent_olc)
        parent_res = compound_lib.FindCompound(parent_tlc)
        if not parent_res:
          _DeleteSidechain(res, ed)
          for atom in res.atoms:
            atom.is_hetatom = False
          print("Removing sidechain of %s, beacuse it has not been found in the compound library"% parent_tlc)
        else:
          #collect atom's names
          modif_atom_names = set([atom.name for atom in res.atoms
                                                  if atom.element != "H" and atom.element != "D" ])
          #if the res is the first or last take all the atoms from the parent res
          if res.FindAtom("OXT").IsValid() :
            parent_atom_names = set([atom.name for atom in parent_res.atom_specs
                                                  if atom.element != "H" and atom.element != "D" ])
          else:
            parent_atom_names = set([atom.name for atom in parent_res.atom_specs
                                                    if atom.element != "H" and atom.element != "D" and not atom.is_leaving ])
          additional_parent_atoms = parent_atom_names - modif_atom_names
          additional_modif_atoms = modif_atom_names - parent_atom_names
          #WHEN mmCIF WILL BE USED, CHANGE IsPeptideLinking() TO IsProtein(), TO EXCLUDE LIGANDS FROM CANONICALISATION
          if res.atoms[0].is_hetatom :
            old_name = res.name
            ed.RenameResidue(res, parent_tlc)
            if additional_parent_atoms:
              if additional_modif_atoms:
                #replacement
                _Replacement(res, ed, old_name)
              else:
                #deletion
                _Deletion(res, ed)
            elif additional_modif_atoms:
              #addition
              _Addition(res, ed, additional_modif_atoms)
            else:
              #unchanged, later check stereochemistry or H atoms
              _Unchanged(res, ed)
          #the res is a peptide but not a ligand (is a protein res)
          else:
            if additional_parent_atoms:# if the sidechain is incomplete
              _DeleteSidechain(res, ed)
  ed.UpdateICS()
  return
  
def _Replacement(res, ed, old_name) :
  #TEMP ONLY MSE
  if old_name == "MSE" :
    for atom in res.atoms:
      atom.is_hetatom = False
    sel = res.FindAtom("SE")
    if sel.IsValid() :
      ed.InsertAtom( res, "SD", sel.pos, "S", sel.occupancy, sel.b_factor ) #S radius=~1;SE=~1.2
      ed.DeleteAtom( sel )
    else:
      _DeleteSidechain(res, ed)
  else:
    _DeleteSidechain(res, ed)
  return
  
def _Deletion(res, ed) :
  _DeleteSidechain(res, ed)
  for atom in res.atoms :
    atom.is_hetatom = False
  return
  
def _Addition(res, ed, additional_modif_atoms) :
  for add_atom_name in additional_modif_atoms:
    add_atom = res.FindAtom( add_atom_name )
    if add_atom.IsValid() :
      ed.DeleteAtom( add_atom )
  for atom in res.atoms:
    atom.is_hetatom = False
  return
  
def _Unchanged(res, ed) :
  if res.chem_class == mol.D_PEPTIDE_LINKING:
    ed.DeleteResidue(res)
  else:
    _DeleteSidechain(res, ed)
    for atom in res.atoms :
      atom.is_hetatom = False
  return
  
def _DeleteSidechain(res, ed) :
  for atom in res.atoms:
    if not atom.name in ['CA','CB','C','N','O']:
      ed.DeleteAtom(atom)
  return
  
#visible functions
__all__ = [Cleanup]
