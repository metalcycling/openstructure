"""
A bunch of algorithms operating on two views. 

Authors: Marco Biasini, Pascal Benkert
"""

from ost import io,mol

def _GetChain(view):
  if view.chain_count!=1:
    raise RuntimeError("chain count of view must be one")
  return view.chains[0]
  
def _GetChains(view_a, view_b):
  return _GetChain(view_a), _GetChain(view_b)


def _EmptyView(view):
  if isinstance(view, mol.EntityHandle):
    return view.CreateEmptyView()
  return view.handle.CreateEmptyView()


def PairResiduesByNum(view_a, view_b, 
                      view_add_flags=mol.ViewAddFlag.INCLUDE_ATOMS):
  """
  Pair residues by residue number.
  """
  not_supported="PairResiduesByNum has no support for unsorted chains"
  chain_a, chain_b=_GetChains(view_a, view_b)
  if not chain_a.InSequence() or not chain_b.InSequence():
    raise RuntimeError(not_supported)
  residues_a=chain_a.residues
  residues_b=chain_b.residues
  residues_a=iter(chain_a.residues)
  residues_b=iter(chain_b.residues)
  result_a=_EmptyView(view_a)
  result_b=_EmptyView(view_b)
  try:
    while True:
      r1=next(residues_a)
      r2=next(residues_b)
      while r1.number<r2.number:
        r1=next(residues_a)
      while r2.number<r1.number:
        r2=next(residues_b)
      assert r1.number==r2.number
      result_a.AddResidue(r1, view_add_flags)
      result_b.AddResidue(r2, view_add_flags)
  except StopIteration:
    pass
  return result_a, result_b
  
def RepresentativeAtoms(ent, chain=None, alpha_and_beta=False):
  """
  Returns a view with one or two representative atom per amino acid residue. 

  There are two basic modes, controlled by the alpha_and_beta parameter:

  When the parameter is false, for residues with a sidechain, the C-beta atom is 
  used, for residues without  sidechain, the C-alpha atom is used. Note that this 
  is different from using the selection

    (aname=CA and rname=GLY) or (aname=CB and rname!=GLY)

  When the alpha_and_beta parameter is true, both C-alpha and C-beta (if
  available) are added to the view.

  If chain is not equal to None, only atoms of the chain with that chain name
  will be added to the view.
  """
  e_view=ent.CreateEmptyView()
  if ent.IsValid():
    # all chains:
    if chain==None:
      for res in ent.residues:
        res=res.handle
        if res.IsPeptideLinking():
          atom =  res.FindAtom('CB')
          if atom.IsValid():
            e_view.AddAtom(atom)   

            #also add C-alpha for hybrid-potential:
            if alpha_and_beta == 1:
              atom =  res.FindAtom('CA')
              if atom.IsValid():
                e_view.AddAtom(atom)

          else: # Cbeta does not exist, Calpha?
            atom =  res.FindAtom('CA')
            if atom.IsValid():
              e_view.AddAtom(atom)

              # count C_alphas twice for hybrid-potential:
              if alpha_and_beta == 1:
                e_view.AddAtom(atom)

    elif chain != "" and ent.FindChain(chain).IsValid():
      for res in ent.FindChain(chain).GetResidueList():
        res=res.handle
        if res.IsPeptideLinking():
          atom =  res.FindAtom('CB')
          if atom.IsValid():
            e_view.AddAtom(atom)

            if alpha_and_beta == 1:
              atom =  res.FindAtom('CA')
              if atom.IsValid():
                e_view.AddAtom(atom)

          else: # Cbeta does not exist, Calpha?
            atom =  res.FindAtom('CA')
            if atom.IsValid():
              e_view.AddAtom(atom)

              # count C_alphas twice for glycine:
              if alpha_and_beta == 1:
                e_view.AddAtom(atom)
  return e_view
