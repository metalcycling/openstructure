from ost import geom
from ost import mol
from ost import seq

def _PreprocessStructures(mdl, ref, mdl_ch1, mdl_ch2, ref_ch1, ref_ch2,
                          ch1_aln = None, ch2_aln = None):
    """ Preprocesses *mdl* and *ref*

    Returns two entity views with the exact same number of residues. I.e. the
    residues correspond to a one-to-one mapping. Additionally, each residue gets
    the int property "dockq_map" assigned, which corresponds to the residue
    index in the respective chain of the processed structures.
    """
    mdl_residues_1 = list()
    mdl_residues_2 = list()
    ref_residues_1 = list()
    ref_residues_2 = list()

    if ch1_aln is not None and ch2_aln is not None:
        # there are potentially already views attached to the alns but
        # we go for *mdl* and *ref* here
        if ch1_aln.GetCount() != 2 or ch2_aln.GetCount() != 2:
            raise RuntimeError("Expect exactly two sequences in provided alns!")

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
                ref_residues_1.append(ref_r)
                mdl_residues_1.append(mdl_r)

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
                ref_residues_2.append(ref_r)
                mdl_residues_2.append(mdl_r)
    else:
        # go by residue numbers
        for mdl_r in mdl.Select(f"cname={mdl_ch1}").residues:
            ref_r = ref.FindResidue(ref_ch1, mdl_r.GetNumber())
            if ref_r.IsValid():
                mdl_residues_1.append(mdl_r)
                ref_residues_1.append(ref_r)
        for mdl_r in mdl.Select(f"cname={mdl_ch2}").residues:
            ref_r = ref.FindResidue(ref_ch2, mdl_r.GetNumber())
            if ref_r.IsValid():
                mdl_residues_2.append(mdl_r)
                ref_residues_2.append(ref_r)

    new_mdl = mdl.handle.CreateEmptyView()
    new_ref = ref.handle.CreateEmptyView()
    for r in mdl_residues_1:
        new_mdl.AddResidue(r.handle, mol.INCLUDE_ALL)
    for r in mdl_residues_2:
        new_mdl.AddResidue(r.handle, mol.INCLUDE_ALL)
    for r in ref_residues_1:
        new_ref.AddResidue(r.handle, mol.INCLUDE_ALL)
    for r in ref_residues_2:
        new_ref.AddResidue(r.handle, mol.INCLUDE_ALL)

    # set dockq_map property
    ch = new_mdl.FindChain(mdl_ch1)
    for r_idx, r in enumerate(ch.residues):
        r.SetIntProp("dockq_map", r_idx)
    ch = new_mdl.FindChain(mdl_ch2)
    for r_idx, r in enumerate(ch.residues):
        r.SetIntProp("dockq_map", r_idx)
    ch = new_ref.FindChain(ref_ch1)
    for r_idx, r in enumerate(ch.residues):
        r.SetIntProp("dockq_map", r_idx)
    ch = new_ref.FindChain(ref_ch2)
    for r_idx, r in enumerate(ch.residues):
        r.SetIntProp("dockq_map", r_idx)

    return (new_mdl, new_ref)

def _GetContacts(ent, ch1, ch2, dist_thresh):
    int1 = ent.Select(f"cname={ch1} and {dist_thresh} <> [cname={ch2}]")
    int2 = ent.Select(f"cname={ch2} and {dist_thresh} <> [cname={ch1}]")
    contacts = set()
    int1_p = [geom.Vec3List([a.pos for a in r.atoms]) for r in int1.residues]
    int2_p = [geom.Vec3List([a.pos for a in r.atoms]) for r in int2.residues]
    for r1, p1 in zip(int1.residues, int1_p):
        for r2, p2 in zip(int2.residues, int2_p):
            if p1.IsWithin(p2, dist_thresh):
                contacts.add((r1.GetIntProp("dockq_map"), r2.GetIntProp("dockq_map")))
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

    mdl_ch1_residues = mdl.FindChain(mdl_ch1).residues
    mdl_ch2_residues = mdl.FindChain(mdl_ch2).residues
    ref_ch1_residues = ref.FindChain(ref_ch1).residues
    ref_ch2_residues = ref.FindChain(ref_ch2).residues

    # iRMSD
    #######
    int1 = ref.Select(f"cname={ref_ch1} and {dist_thresh} <> [cname={ref_ch2}]")
    int2 = ref.Select(f"cname={ref_ch2} and {dist_thresh} <> [cname={ref_ch1}]")
    int1_indices = [r.GetIntProp("dockq_map") for r in int1.residues]
    int2_indices = [r.GetIntProp("dockq_map") for r in int2.residues]
    ref_pos = geom.Vec3List()
    mdl_pos = geom.Vec3List()
    atom_names = ['CA','C','N','O']
    for idx in int1_indices:
        ref_r = ref_ch1_residues[idx]
        mdl_r = mdl_ch1_residues[idx]
        for aname in atom_names:
            ref_a = ref_r.FindAtom(aname)
            mdl_a = mdl_r.FindAtom(aname)
            if ref_a.IsValid() and mdl_a.IsValid():
                ref_pos.append(ref_a.pos)
                mdl_pos.append(mdl_a.pos)

    for idx in int2_indices:
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
    # receptor is by definition the larger chain
    if len(ref_ch1_residues) > len(ref_ch2_residues):
        ref_receptor_residues = ref_ch1_residues
        ref_ligand_residues = ref_ch2_residues
        mdl_receptor_residues = mdl_ch1_residues
        mdl_ligand_residues = mdl_ch2_residues
    else:
        ref_receptor_residues = ref_ch2_residues
        ref_ligand_residues = ref_ch1_residues
        mdl_receptor_residues = mdl_ch2_residues
        mdl_ligand_residues = mdl_ch1_residues

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
    mapped_model, mapped_ref = _PreprocessStructures(mdl, ref, mdl_ch1, mdl_ch2,
                                                     ref_ch1, ref_ch2,
                                                     ch1_aln = ch1_aln,
                                                     ch2_aln = ch2_aln)
    nnat, nmdl, fnat, fnonnat = _ContactScores(mapped_model, mapped_ref,
                                         mdl_ch1, mdl_ch2, ref_ch1, ref_ch2)
    irmsd, lrmsd = _RMSDScores(mapped_model, mapped_ref,
                               mdl_ch1, mdl_ch2, ref_ch1, ref_ch2)
    return {"nnat": nnat,
            "nmdl": nmdl,
            "fnat": fnat,
            "fnonnat": fnonnat,
            "irmsd": round(irmsd, 3),
            "lrmsd": round(lrmsd, 3),
            "DockQ": round(_DockQ(fnat, lrmsd, irmsd, 8.5, 1.5), 3)}
