from ost import seq, mol

def _RenumberSeq(seq_handle, old_number_label=None):
  if not seq_handle.HasAttachedView():
    raise RuntimeError("Sequence Handle has no attached view")
  ev = seq_handle.attached_view.CreateEmptyView()
  new_numbers = mol.ResNumList()
  for pos in range(len(seq_handle)):
    if seq_handle[pos] != '-':
      r = seq_handle.GetResidue(pos)
      if r.IsValid():
          ev.AddResidue(r, mol.INCLUDE_ALL)
          new_numbers.append(pos+1)
          if old_number_label is not None:
            r.SetIntProp(old_number_label, r.number.GetNum())
      else:
        raise RuntimeError('Error: renumbering failed at position %s' % pos)
  return ev, new_numbers

def _RenumberAln(aln, seq_index, old_number_label=None):
  if not aln.sequences[seq_index].HasAttachedView():
    raise RuntimeError("Sequence Handle has no attached view")
  counter=0
  ev = aln.sequences[seq_index].attached_view.CreateEmptyView()
  new_numbers = mol.ResNumList()
  for col in aln:
    if col[0] != '-' and col[seq_index] != '-':
      if col[0] != col[seq_index]:
        raise RuntimeError("residue mismatch at position %d (%s vs %s) "\
                           "(renumbering failed)" % (counter, col[0],
                                                     col[seq_index]))
      rnum = aln.GetSequence(seq_index).GetResidueIndex(counter)
      r = aln.GetSequence(seq_index).GetResidue(counter)
      if not r.IsValid():
        raise RuntimeError("invalid residue at postion %s (renumbering failed)"\
                           % (counter))
      ev.AddResidue(r, mol.INCLUDE_ALL)
      new_numbers.append(counter+1)
      if old_number_label is not None:
        r.SetIntProp(old_number_label, r.number.GetNum())
    counter += 1
  return ev, new_numbers


def Renumber(seq_handle, sequence_number_with_attached_view=1, old_number_label=None):
  """
  Function to renumber an entity according to an alignment between the model
  sequence and the full-length target sequence. The aligned model sequence or
  the alignment itself with an attached view needs to be provided. Upon
  succcess, the renumbered entity is returned.
  If an alignment is given, the first sequence of the alignment is considered
  the full-length sequence and it must match the model sequence wherever it is
  aligned (i.e. excluding gaps).

  .. code-block:: python

    from ost.seq.alg import renumber
    from ost.bindings.clustalw import *
    ent = io.LoadPDB("path_to_model")
    s = io.LoadSequence("path_to_full_length_fasta_seqeunce")
    pdb_seq = seq.SequenceFromChain("model", ent.chains[0])
    aln = ClustalW(s, pdb_seq)
    aln.AttachView(1, ent.chains[0].Select(""))
    e = Renumber(aln.sequences[1])
    io.SavePDB(e, "renum.pdb")
   
  :param seq_handle: Sequence or alignment handle with attached view.
  :type seq_handle:  :class:`SequenceHandle` / :class:`AlignmentHandle`
  :param sequence_number_with_attached_view: Sequence number for the aln. handle
                                             (not used if seq. handle given)
  :type sequence_number_with_attached_view:  :class:`int`
  :raises: :exc:`RuntimeError` if unknown type of *seq_handle* or if attached
           view is missing or if the given alignment sequence is inconsistent.
  """
  if isinstance(seq_handle, seq.SequenceHandle) \
     or isinstance(seq_handle, seq.ConstSequenceHandle):
    ev, new_numbers = _RenumberSeq(seq_handle, old_number_label)
  elif isinstance(seq_handle, seq.AlignmentHandle):
    ev, new_numbers = _RenumberAln(seq_handle, sequence_number_with_attached_view, old_number_label)
  else:
    raise RuntimeError("Unknown input type " + str(type(seq_handle)))

  ev.AddAllInclusiveBonds()
  new_ent = mol.CreateEntityFromView(ev, False);
  new_ent.EditXCS().RenumberChain(new_ent.chains[0], new_numbers)
  return new_ent

# choose visible interface
__all__ = ('Renumber', )