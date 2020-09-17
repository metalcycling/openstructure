from ._ost_seq_alg import *
from ost.seq.alg.mat import *

def ValidateSEQRESAlignment(aln, chain=None):
  """
  Checks if sequence in alignment has same connectivity as residues in chain.
  This looks for connected stretches in both the sequence and the chain and
  returns False if they don't match. This uses the connectivity of the protein
  backbone.

  :param aln: Alignment of two sequences with second one expected to map to
              residues in *chain*.
  :type aln: :class:`~ost.seq.AlignmentHandle`
  :param chain: Source of the sequence.
  :type chain: :class:`~ost.mol.ChainHandle`

  :returns: True if all residues (beside gapped ones) are connected, False
            otherwise.
  """
  from ost import LogWarning
  from ost import seq
  from ost import mol
  if aln.GetCount() != 2:
    raise ValueError('Alignment contains more than 2 sequences!')
  sequence = aln.GetSequence(1)
  if len(sequence) == 0:
    return True
  if chain == None:
    if sequence.HasAttachedView() == False:
      raise ValueError("Alignment is missing an attached chain view.")
    chain = sequence.GetAttachedView()
  residues = chain.residues
  # eat up all beginning gaps
  j = 1
  for s in sequence:
    if s != '-':
      break
    j += 1;
  l = sequence[j-1]
  i = 0
  # run over sequence & alignment
  for s in sequence[j:]:
    if s != '-':
      i += 1
      r1 = residues[i-1]
      r2 = residues[i]
      if r1.one_letter_code=='?' or r2.one_letter_code=='?':
        continue
      if l != '-':
        if not mol.InSequence(r1.handle, r2.handle):
          LogWarning('%s and %s are not connected by peptide bond' % (str(r1), str(r2)))
          return False
      else:
        if mol.InSequence(r1.handle, r2.handle):
          LogWarning('%s and %s are connected by peptide bond' % (str(r1), str(r2)))
          return False
    l = s
  return True

def AlignToSEQRES(chain, seqres, try_resnum_first=False, validate=True):
  """
  Aligns the residues of chain to the SEQRES sequence, inserting gaps where 
  needed. The function uses the connectivity of the protein backbone to find 
  consecutive peptide fragments. These fragments are then aligned to the SEQRES 
  sequence.
  
  All the non-ligand, peptide-linking residues of the chain must be listed in 
  SEQRES. If there are any additional residues in the chain, the function
  raises a ValueError.

  :param chain: Source of the sequence
  :type chain: :class:`~ost.mol.ChainHandle`
  :param seqres: SEQRES sequence
  :type seqres: :class:`str`
  :param try_resnum_first: If set to True, this first builds an alignment using
                           residue numbers and checks if the one-letter-codes
                           match. If they all match, this alignment is used
                           (and possibly validated). Otherwise, it displays a
                           warning and falls back to the connectivity-based
                           alignment.
  :type try_resnum_first: :class:`bool`
  :param validate: If set to True, the alignment is additionally checked by
                   :func:`~ost.seq.alg.ValidateSEQRESAlignment` and raises
                   a ValueError if the validation failed.
  :type validate: :class:`bool`

  :returns: The alignment of the residues in the chain and the SEQRES entries.
  :rtype: :class:`~ost.seq.AlignmentHandle`
  """

  def IsEqual(olc1, olc2):
    return olc1 in ('X', '?') or olc2 in ('X', '?') or olc1 == olc2

  from ost import seq
  from ost import mol
  from ost import LogWarning
  view=chain
  residues=view.residues
  if len(residues)==0:
    return seq.CreateAlignment()
  if try_resnum_first:
    aln_seq = seq.CreateSequence('atoms', '-'*len(seqres))
    aligned_resnums = set()
    for r1 in residues:
      if r1.number.num in aligned_resnums:
        LogWarning('Residue numbers must be unique. Already observed %i, ' \
                   'cannot align %s anymore.'%(r1.number.num, r1.qualified_name))
        try_resnum_first = False
        break
      if r1.number.num <= len(seqres) and r1.number.num > 0:
        if IsEqual(seqres[r1.number.num - 1], r1.one_letter_code):
          aln_seq[r1.number.num - 1] = r1.one_letter_code
          aligned_resnums.add(r1.number.num)
        else:
          LogWarning('Sequence mismatch: chain has "' + r1.one_letter_code +
                     '", while SEQRES is "' + seqres[r1.number.num - 1] +
                     '" at the corresponding position.')
          try_resnum_first = False
          break
      else:
        warning = 'Residue with number %i is outside of the range covered by '\
                  'SEQRES [1, %i]'%(r1.number.num, len(seqres)) 
        LogWarning(warning)
        try_resnum_first = False
        break
  if not try_resnum_first:
    fragments=[residues[0].one_letter_code]
    for r1, r2 in zip(residues[:-1], residues[1:]):
      if not mol.InSequence(r1.handle, r2.handle):
        fragments.append('')
      fragments[-1]+=r2.one_letter_code
    ss=str(seqres)
    pos=0
    aln_seq=''
    for frag in fragments:
      new_pos=ss.find(frag, pos)
      if new_pos==-1:
        raise ValueError('"%s" is not a substring of "%s"' % (frag, ss))
      aln_seq+='-'*(new_pos-pos)+frag
      pos=new_pos+len(frag)
    aln_seq = seq.CreateSequence('atoms',
                                 aln_seq+('-'*(len(seqres)-len(aln_seq))))
  alignment = seq.CreateAlignment(seq.CreateSequence('SEQRES', str(seqres)), 
                                  aln_seq)
  if validate and not ValidateSEQRESAlignment(alignment, view):
    raise ValueError("SEQRES cannot be aligned with its corresponding chain.")
  return alignment


def AlignmentFromChainView(chain, handle_seq_name='handle', 
                           view_seq_name='view'):
  """
  Creates and returns the sequence alignment of the given chain view to the 
  chain handle. The alignment contains two sequences, the first containing all 
  non-ligand peptide-linking residues, the second containing all non-ligand 
  peptide-linking residues that are part of the view. 
  
  :param chain: A valid chain
  :type chain: :class:`~ost.mol.ChainView`
  
  :param handle_seq_name: Name of the handle sequence in the output alignment
  :param view_seq_name: Name of the view sequence in the output alignment
  :returns: The alignment
  :rtype: :class:`~ost.seq.AlignmentHandle`
  
  """
  from ost import seq
  v0=chain.handle.Select('ligand=false and peptide=true')
  v1=chain.Select('ligand=false and peptide=true')
  s0=seq.CreateSequence(handle_seq_name, '')
  s1=seq.CreateSequence(view_seq_name, '')
  s0.AttachView(v0)
  s1.AttachView(v1)
  res0=v0.residues
  res1=v1.residues
  idx0, idx1=(0, 0)
  while idx0<len(res0):
    s0.Append(res0[idx0].one_letter_code)
    if idx1<len(res1) and res1[idx1].handle==res0[idx0].handle:
      s1.Append(res1[idx1].one_letter_code)
      idx1+=1
    else:
      s1.Append('-')
    idx0+=1
  return seq.CreateAlignment(s0, s1)

def PredictContacts(ali):
  """
  Predicts contacts from a multiple sequence alignment using a combination
  of Mutual Information (*MI*) and the Contact Substitution Score (*CoEvoSc*).
  MI is calculated with the APC and small number corrections as well as with a 
  transformation into Z-scores. The *CoEvoSc* is calculated using the default 
  PairSubstWeightMatrix (see seq.alg.LoadDefaultPairSubstWeightMatrix).
  The final score for a pair of columns *(i,j)* of **ali** is obtained from:
  
  Sc(i,j)=MI(i,j)exp(CoEvoSc(i,j))      if *(i,j)* >=0
  
  Sc(i,j)=MI(i,j)exp(1-CoEvoSc(i,j))    if *(i,j)* <0

  :param ali: The multiple sequence alignment
  :type ali: :class:`~ost.seq.AlignmentHandle`
  """
  import math
  from ost import seq
  if not type(ali)==type(seq.AlignmentHandle()):
    print("Parameter should be an AlignmentHandle")
    return
  mi=CalculateMutualInformation(ali)
  CoEvoSc=CalculateContactSubstitutionScore(ali)
  ncol=ali.GetLength()
  for i in range(ncol):
    for j in range(ncol):
      if mi.matrix[i][j]>=0:
        mi.matrix[i][j]=mi.matrix[i][j]*(math.exp(CoEvoSc.matrix[i][j]))
      else:
        mi.matrix[i][j]=mi.matrix[i][j]*(math.exp(1.-CoEvoSc.matrix[i][j]))
  mi.RefreshSortedIndices()
  return mi

def CalculateContactProbability(cpred_res,method):
  """
  Calculate the probability of a predicted contact to be correct.
  This simply transforms the score associated with a prediction into a probability.

  :param cpred_res: A contact prediction result
  :param method: The method which was used for contact prediction. Should be one
   of "MI", "MIp", "MIpz", "cevoMI", "cevo"

  :type cpred_res: :class:`ost.seq.alg.ContactPredictionScoreResult`
  :type method: :class:`str`
  """
  import math
  def _growth_function(x,K,x0,B,nu,Q):
    p=K/(1+Q*math.exp(-B*(x-x0)))**(1/nu)
    p=min(1,max(0,p))
    return p
  
  def _decay_function(x,A,k):
    p=A*math.exp(-k*x)
    p=min(1,max(0,p))
    return p

  prob_params={}
  prob_params["MI"]=[0.05858455345122933, 0.8930350957023122]
  prob_params["MIp"]=[0.10019621004607637, 0.9065429261332942]
  prob_params["MIpz"]=[0.7368147563063437, 4.638820470770171, 0.6383539191372934, 1.1627300209863376, 19.63060874042289]
  prob_params["cevoMI"]=[0.779979757231944, 1.642357937131157, 0.3267847173036033, 0.3742848849873358, 5.1816922372446]
  prob_params["cevo"]=[0.7564532665623617, 0.5947472274271304, 3.8548775389879166, 0.46203017320927053, 1.6198602780123705]

  if not method in prob_params:
    raise ValueError("method should be one of MI, MIp, MIpz, cevoMI, cevo")
  params=prob_params[method]
  cpred_res.RefreshSortedIndices()
  nres=len(cpred_res.matrix)
  probabilities=[[0 for i in range(nres)] for j in range(nres)]

  if len(params)==2:func=_decay_function
  else:func=_growth_function

  nres=float(nres)
  if len(params)==2:
    for k,(i,j) in enumerate(cpred_res.sorted_indices):
      p=_decay_function(math.log10(0.01+k/nres),*params)
      probabilities[i][j]=p
      probabilities[j][i]=p
  else:
    for k,(i,j) in enumerate(cpred_res.sorted_indices):
      p=_growth_function(cpred_res.GetScore(i,j),*params)
      probabilities[i][j]=p
      probabilities[j][i]=p
  cpred_res.probabilities=probabilities
  return cpred_res




