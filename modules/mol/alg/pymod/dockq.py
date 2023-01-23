from ost import geom
from ost import mol
from ost import seq

def _PreprocessStructures(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2,
                          ch1_aln = None, ch2_aln = None):
    """ Preprocesses *mdl* and *ref*

    Sets int properties to each residue in mdl_ch1, mdl_ch2 as well as
    the respective reference chains.
    dockq_mapped: 1 if a residues in mdl_ch1 is mapped to a residue in ref_ch1 
                  and vice versa, 0 otherwise. Same is done for mdl_ch2 and
                  ref_ch2.
    dockq_idx: If a pair of residue is mapped, the same index will be set
               to both residues. The index is unique otherwise.

    By default, mapping happens with residue numbers but you can enforce a
    mapping with an alignment. In the example of ch1_aln, the first sequence
    corresponds to the ATOMSEQ of ref_ch1 in ref and the second sequence to
    the ATOMSEQ of mdl_ch1 in mdl.
    """

    # set default values for dockq_mapped and dockq_idx properties
    # => makes sure that we have a clean slate if stuff has been set in
    # previous runs
    for cname in [ref_ch1, ref_ch2]:
        ch = ref.FindChain(cname)
        for r in ch.residues:
            r.SetIntProp("dockq_mapped", 0)
            r.SetIntProp("dockq_idx", -1)
    for cname in [mdl_ch1, mdl_ch2]:
        ch = mdl.FindChain(cname)
        for r in ch.residues:
            r.SetIntProp("dockq_mapped", 0)
            r.SetIntProp("dockq_idx", -1)

    dockq_idx = 0
    if ch1_aln is not None and ch2_aln is not None:
        # there are potentially already views attached to the alns but
        # we go for *mdl* and *ref* here
        if ch1_aln.GetCount() != 2 or ch2_aln.GetCount() != 2:
            raise RuntimeError("Expect exactly two sequences in alns provided "
                               "to DockQ!")
        tmp = ch1_aln.GetSequence(0)
        ref_s1 = seq.CreateSequence(tmp.GetName(), str(tmp))
        ref_s1.SetOffset(tmp.GetOffset())
        ref_s1.AttachView(ref.Select(f"cname={ref_ch1}"))
        tmp = ch1_aln.GetSequence(1)
        mdl_s1 = seq.CreateSequence(tmp.GetName(), str(tmp))
        mdl_s1.SetOffset(tmp.GetOffset())
        mdl_s1.AttachView(mdl.Select(f"cname={mdl_ch1}"))
        new_ch1_aln = seq.CreateAlignment(ref_s1, mdl_s1)
        for col in new_ch1_aln:
            if col[0] != '-' and col[1] != '-':
                ref_r = col.GetResidue(0)
                mdl_r = col.GetResidue(1)
                if not (ref_r.IsValid() and ref_r.one_letter_code == col[0]):
                    raise RuntimeError("DockQ: mismatch between provided "
                                       "alignments and ATOMSEQ in structures")
                if not (mdl_r.IsValid() and mdl_r.one_letter_code == col[1]):
                    raise RuntimeError("DockQ: mismatch between provided "
                                       "alignments and ATOMSEQ in structures")
                ref_r.SetIntProp("dockq_idx", dockq_idx)
                mdl_r.SetIntProp("dockq_idx", dockq_idx)
                ref_r.SetIntProp("dockq_mapped", 1)
                mdl_r.SetIntProp("dockq_mapped", 1)
                dockq_idx += 1

        tmp = ch2_aln.GetSequence(0)
        ref_s2 = seq.CreateSequence(tmp.GetName(), str(tmp))
        ref_s2.SetOffset(tmp.GetOffset())
        ref_s2.AttachView(ref.Select(f"cname={ref_ch2}"))
        tmp = ch2_aln.GetSequence(1)
        mdl_s2 = seq.CreateSequence(tmp.GetName(), str(tmp))
        mdl_s2.SetOffset(tmp.GetOffset())
        mdl_s2.AttachView(mdl.Select(f"cname={mdl_ch2}"))
        new_ch2_aln = seq.CreateAlignment(ref_s2, mdl_s2)
        for col in new_ch2_aln:
            if col[0] != '-' and col[1] != '-':
                ref_r = col.GetResidue(0)
                mdl_r = col.GetResidue(1)
                if not (ref_r.IsValid() and ref_r.one_letter_code == col[0]):
                    raise RuntimeError("DockQ: mismatch between provided "
                                       "alignments and ATOMSEQ in structures")
                if not (mdl_r.IsValid() and mdl_r.one_letter_code == col[1]):
                    raise RuntimeError("DockQ: mismatch between provided "
                                       "alignments and ATOMSEQ in structures")
                ref_r.SetIntProp("dockq_idx", dockq_idx)
                mdl_r.SetIntProp("dockq_idx", dockq_idx)
                ref_r.SetIntProp("dockq_mapped", 1)
                mdl_r.SetIntProp("dockq_mapped", 1)
                dockq_idx += 1
    else:
        # go by residue numbers
        for mdl_r in mdl.Select(f"cname={mdl_ch1}").residues:
            ref_r = ref.FindResidue(ref_ch1, mdl_r.GetNumber())
            if ref_r.IsValid():
                ref_r.SetIntProp("dockq_idx", dockq_idx)
                mdl_r.SetIntProp("dockq_idx", dockq_idx)
                ref_r.SetIntProp("dockq_mapped", 1)
                mdl_r.SetIntProp("dockq_mapped", 1)
                dockq_idx += 1
        for mdl_r in mdl.Select(f"cname={mdl_ch2}").residues:
            ref_r = ref.FindResidue(ref_ch2, mdl_r.GetNumber())
            if ref_r.IsValid():
                ref_r.SetIntProp("dockq_idx", dockq_idx)
                mdl_r.SetIntProp("dockq_idx", dockq_idx)
                ref_r.SetIntProp("dockq_mapped", 1)
                mdl_r.SetIntProp("dockq_mapped", 1)
                dockq_idx += 1

    # set unique dockq_idx property for all residues that are still -1
    for cname in [ref_ch1, ref_ch2]:
        ch = ref.FindChain(cname)
        for r in ch.residues:
            if r.GetIntProp("dockq_idx") == -1:
                r.SetIntProp("dockq_idx", dockq_idx)
                dockq_idx += 1
    for cname in [mdl_ch1, mdl_ch2]:
        ch = mdl.FindChain(cname)
        for r in ch.residues:
            if r.GetIntProp("dockq_idx") == -1:
                r.SetIntProp("dockq_idx", dockq_idx)
                dockq_idx += 1

def _GetContacts(ent, ch1, ch2, dist_thresh):
    int1 = ent.Select(f"cname={ch1} and {dist_thresh} <> [cname={ch2}]")
    int2 = ent.Select(f"cname={ch2} and {dist_thresh} <> [cname={ch1}]")
    contacts = set()
    int1_p = [geom.Vec3List([a.pos for a in r.atoms]) for r in int1.residues]
    int2_p = [geom.Vec3List([a.pos for a in r.atoms]) for r in int2.residues]
    for r1, p1 in zip(int1.residues, int1_p):
        for r2, p2 in zip(int2.residues, int2_p):
            if p1.IsWithin(p2, dist_thresh):
                contacts.add((r1.GetIntProp("dockq_idx"),
                              r2.GetIntProp("dockq_idx")))
    return contacts

def _ContactScores(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2, dist_thresh=5.0):
    ref_contacts = _GetContacts(ref, ref_ch1, ref_ch2, dist_thresh)
    mdl_contacts = _GetContacts(mdl, mdl_ch1, mdl_ch2, dist_thresh)

    nnat = len(ref_contacts)
    nmdl = len(mdl_contacts)

    fnat = len(ref_contacts.intersection(mdl_contacts))
    if nnat > 0:
        fnat /= nnat

    fnonnat = len(mdl_contacts.difference(ref_contacts))
    if len(mdl_contacts) > 0:
        fnonnat /= len(mdl_contacts)

    return (nnat, nmdl, fnat, fnonnat)

def _RMSDScores(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2, dist_thresh=10.0):

    # make mapped residues accessible by the dockq_idx property
    mapped_mdl = mdl.Select(f"cname={mdl_ch1},{mdl_ch2} and grdockq_mapped=1")
    mapped_ref = ref.Select(f"cname={ref_ch1},{ref_ch2} and grdockq_mapped=1")
    mdl_ch1_residues = mapped_mdl.FindChain(mdl_ch1).residues
    mdl_ch1_residues = {r.GetIntProp("dockq_idx"): r for r in mdl_ch1_residues}
    mdl_ch2_residues = mapped_mdl.FindChain(mdl_ch2).residues
    mdl_ch2_residues = {r.GetIntProp("dockq_idx"): r for r in mdl_ch2_residues}
    ref_ch1_residues = mapped_ref.FindChain(ref_ch1).residues
    ref_ch1_residues = {r.GetIntProp("dockq_idx"): r for r in ref_ch1_residues}
    ref_ch2_residues = mapped_ref.FindChain(ref_ch2).residues
    ref_ch2_residues = {r.GetIntProp("dockq_idx"): r for r in ref_ch2_residues}

    # iRMSD
    #######
    int1 = ref.Select(f"cname={ref_ch1} and {dist_thresh} <> "
                      f"[cname={ref_ch2}]")
    int2 = ref.Select(f"cname={ref_ch2} and {dist_thresh} <> "
                      f"[cname={ref_ch1}]")

    int1_indices = [r.GetIntProp("dockq_idx") for r in int1.residues]
    int2_indices = [r.GetIntProp("dockq_idx") for r in int2.residues]
    ref_pos = geom.Vec3List()
    mdl_pos = geom.Vec3List()
    atom_names = ['CA','C','N','O']
    for idx in int1_indices:
        if idx in ref_ch1_residues and idx in mdl_ch1_residues:
            ref_r = ref_ch1_residues[idx]
            mdl_r = mdl_ch1_residues[idx]
            for aname in atom_names:
                ref_a = ref_r.FindAtom(aname)
                mdl_a = mdl_r.FindAtom(aname)
                if ref_a.IsValid() and mdl_a.IsValid():
                    ref_pos.append(ref_a.pos)
                    mdl_pos.append(mdl_a.pos)
    for idx in int2_indices:
        if idx in ref_ch2_residues and idx in mdl_ch2_residues:
            ref_r = ref_ch2_residues[idx]
            mdl_r = mdl_ch2_residues[idx]
            for aname in atom_names:
                ref_a = ref_r.FindAtom(aname)
                mdl_a = mdl_r.FindAtom(aname)
                if ref_a.IsValid() and mdl_a.IsValid():
                    ref_pos.append(ref_a.pos)
                    mdl_pos.append(mdl_a.pos)

    if len(mdl_pos) >= 3:
        sup_result = mol.alg.SuperposeSVD(mdl_pos, ref_pos)
        irmsd = sup_result.rmsd
    else:
        irmsd = 0.0

    # lRMSD
    #######
    # receptor is by definition the larger chain in ref
    n_ch1 = len(ref.FindChain(ref_ch1).residues)
    n_ch2 = len(ref.FindChain(ref_ch2).residues)
    if n_ch1 > n_ch2:
        ref_receptor_residues = ref_ch1_residues.values()
        ref_ligand_residues = ref_ch2_residues.values()
        mdl_receptor_residues = \
        [mdl_ch1_residues[idx] for idx in ref_ch1_residues.keys()]
        mdl_ligand_residues = \
        [mdl_ch2_residues[idx] for idx in ref_ch2_residues.keys()] 
    else:
        ref_receptor_residues = ref_ch2_residues.values()
        ref_ligand_residues = ref_ch1_residues.values()
        mdl_receptor_residues = \
        [mdl_ch2_residues[idx] for idx in ref_ch2_residues.keys()]
        mdl_ligand_residues = \
        [mdl_ch1_residues[idx] for idx in ref_ch1_residues.keys()]
    ref_receptor_positions = geom.Vec3List()
    mdl_receptor_positions = geom.Vec3List()
    ref_ligand_positions = geom.Vec3List()
    mdl_ligand_positions = geom.Vec3List()
    for ref_r, mdl_r in zip(ref_receptor_residues, mdl_receptor_residues):
        for aname in atom_names:
            ref_a = ref_r.FindAtom(aname)
            mdl_a = mdl_r.FindAtom(aname)
            if ref_a.IsValid() and mdl_a.IsValid():
                ref_receptor_positions.append(ref_a.pos)
                mdl_receptor_positions.append(mdl_a.pos)
    for ref_r, mdl_r in zip(ref_ligand_residues, mdl_ligand_residues):
        for aname in atom_names:
            ref_a = ref_r.FindAtom(aname)
            mdl_a = mdl_r.FindAtom(aname)
            if ref_a.IsValid() and mdl_a.IsValid():
                ref_ligand_positions.append(ref_a.pos)
                mdl_ligand_positions.append(mdl_a.pos)

    if len(mdl_receptor_positions) >= 3:
        sup_result = mol.alg.SuperposeSVD(mdl_receptor_positions,
                                          ref_receptor_positions)
        mdl_ligand_positions.ApplyTransform(sup_result.transformation)
        lrmsd = mdl_ligand_positions.GetRMSD(ref_ligand_positions)
    else:
        lrmsd = 0.0

    return (irmsd, lrmsd)

def _ScaleRMSD(rmsd, d):
    return 1.0/(1+(rmsd/d)**2)

def _DockQ(fnat, lrmsd, irmsd, d1, d2):
    """ The final number chrunching as described in the DockQ manuscript
    """
    return (fnat + _ScaleRMSD(lrmsd, d1) + _ScaleRMSD(irmsd, d2))/3

def DockQ(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2,
          ch1_aln=None, ch2_aln=None):
    """ Computes DockQ for specified interface

    DockQ is described in: Sankar Basu and Bjoern Wallner (2016), "DockQ: A
    Quality Measure for Protein-Protein Docking Models", PLOS one 

    Residues are mapped based on residue numbers by default. If you provide
    *ch1_aln* and *ch2_aln* you can enforce an arbitrary mapping.

    :param mdl: Model structure
    :type mdl: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param ref: Reference structure, i.e. native structure
    :type ref: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param mdl_ch1: Specifies chain in model constituting first part of
                    interface
    :type mdl_ch1: :class:`str`
    :param mdl_ch2: Specifies chain in model constituting second part of
                    interface
    :type mdl_ch2: :class:`str`
    :param ref_ch1: ref equivalent of mdl_ch1
    :type ref_ch1: :class:`str`
    :param ref_ch2: ref equivalent of mdl_ch2
    :type ref_ch2: :class:`str`
    :param ch1_aln: Alignment with two sequences to map *ref_ch1* and *mdl_ch1*.
                    The first sequence must match the sequence in *ref_ch1* and
                    the second to *mdl_ch1*.
    :type ch1_aln: :class:`ost.seq.AlignmentHandle`
    :param ch2_aln: Alignment with two sequences to map *ref_ch2* and *mdl_ch2*.
                    The first sequence must match the sequence in *ref_ch2* and
                    the second to *mdl_ch2*.
    :type ch2_aln: :class:`ost.seq.AlignmentHandle`
    :returns: :class:`dict` with keys nnat, nmdl, fnat, fnonnat, irmsd, lrmsd,
              DockQ which corresponds to the equivalent values in the original
              DockQ implementation.
    """
    _PreprocessStructures(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2,
                          ch1_aln = ch1_aln, ch2_aln = ch2_aln)
    nnat, nmdl, fnat, fnonnat = _ContactScores(mdl, ref, mdl_ch1, mdl_ch2,
                                               ref_ch1, ref_ch2)
    irmsd, lrmsd = _RMSDScores(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2)
    return {"nnat": nnat,
            "nmdl": nmdl,
            "fnat": fnat,
            "fnonnat": fnonnat,
            "irmsd": round(irmsd, 3),
            "lrmsd": round(lrmsd, 3),
            "DockQ": round(_DockQ(fnat, lrmsd, irmsd, 8.5, 1.5), 3)}
