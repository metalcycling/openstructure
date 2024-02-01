model_file='sdh.pdb'
e=io.LoadEntity(model_file)


def StripSidechain(edi, residue):
  atoms = residue.GetAtomList()
  for atom in atoms:
    if atom.name not in ['CA', 'C', 'N', 'O']:
      edi.DeleteAtom(atom)

def StripSidechains(entity):
  edi=entity.EditXCS(mol.EditMode.BUFFERED_EDIT)
  for c in e.GetChainList():
    residues=c.GetResidueList()
    for residue in residues:
      StripSidechain(edi, residue)

StripSidechains(e)
g=gfx.Entity('bare backbone', e)
scene.Add(g)
scene.CenterOn(g)


