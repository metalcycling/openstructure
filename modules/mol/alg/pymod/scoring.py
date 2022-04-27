from ost import mol
from ost import seq
from ost import io
from ost.mol.alg import lddt
from ost.mol.alg import qsscoring


class lDDTBSScorer:
    """Scorer specific for a reference/model pair

    Computes lDDT score on residues that constitute a binding site and can deal
    with oligos using a chain mapping derived from
    :class:`ost.mol.alg.qsscoring.QSscorer.chain_mapping`

    There are two options to initialize :class:`lDDTBSScorer`

    * provide a *reference* and *model* structure, will be used to internally
      setup a :class:`ost.mol.alg.qsscoring.QSscorer` from which a chain mapping
      is derived.
    * provide a :class:`ost.mol.alg.qsscoring.QSscorer` directly to make use of a
      potentially cached chain mapping.

    In both cases, the actually evaluated structures are derived from
    :class:`ost.mol.alg.qsscoring.QSscorer.qs_ent_1` (reference) and
    :class:`ost.mol.alg.qsscoring.QSscorer.qs_ent_2` (model). That means they
    are cleaned as described in
    :class:`ost.mol.alg.qsscoring.QSscoreEntity.ent`.


    :param reference: Reference structure
    :type reference: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param model: Model structure
    :type model: :class:`ost.mol.EntityView`/:class:`ost.mol.EntityHandle`
    :param qs_scorer: QSscorer object where you potentially already have a 
                      chain mapping cached - *model* and *reference* will be
                      neglected if this is provided.
    :type qs_scorer: :class:`ost.mol.alg.qsscoring.QSscorer`
    :param residue_number_alignment: Passed to QSscorer constructor if it needs
                                     to be initialized with *reference* and
                                     *model*.
    :type residue_number_alignment: :class:`bool`
    :raises: :class:`RuntimeError` if you don't provide *qs_scorer* or
             *reference* and *model*,
             :class:`ost.mol.alg.qsscoring.QSscoreError` if QSscorer
             constructor raises.
    """
    def __init__(self, reference=None, model=None,
                 qs_scorer=None, residue_number_alignment=False):
        if qs_scorer is not None:
            self.qs_scorer = qs_scorer
        elif model is not None and reference is not None:
            self.qs_scorer = qsscoring.QSscorer(reference, model,
                                                residue_number_alignment)
        else:
            raise RuntimeError("Must either provide qs_scorer or reference and "
                               "model")
        self.ref = self.qs_scorer.qs_ent_1.ent.Select("peptide=true")
        self.mdl = self.qs_scorer.qs_ent_2.ent.Select("peptide=true")

    def ScoreBS(self, ligand, radius = 4.0, lddt_radius=10.0,
                return_mapping=False):
        """Computes binding site lDDT score given *ligand*

        :param ligand: Defines the scored binding site, i.e. provides positions
                       to perform proximity search
        :type ligand: r'((Residue)|(Chain)|(Entity))((View)|(Handle))'
        :param radius: Reference residues with any atom position within *radius*
                       of *ligand* consitute the scored binding site
        :type radius: :class:`float`
        :param lddt_radius: Passed as *inclusion_radius* to
                            :class:`ost.mol.alg.lddt.lDDTScorer`
        :type lddt_radius: :class:`float`
        :param return_mapping: If true, returns binding site mapping information
                               in addition to the raw lDDTBS score, i.e. returns
                               a tuple with 1: lDDTBS score 2: list of qualified
                               residue names in reference 3: same for model
        :type return_mapping: :class:`bool`

        :returns: lDDTBS score or tuple if *return_mapping* is True
        """

        # create view of reference binding site
        ref_residues_hashes = set() # helper to keep track of added residues
        for ligand_at in ligand.atoms:
            close_atoms = self.ref.FindWithin(ligand_at.GetPos(), radius)
            for close_at in close_atoms:
                ref_res = close_at.GetResidue()
                h = ref_res.handle.GetHashCode()
                if h not in ref_residues_hashes:
                    ref_residues_hashes.add(h)

        # reason for doing that separately is to guarantee same ordering of
        # residues as in underlying entity. (Reorder by ResNum seems only
        # available on ChainHandles)
        ref_bs = self.ref.CreateEmptyView()
        for ch in self.ref.chains:
            for r in ch.residues:
                if r.handle.GetHashCode() in ref_residues_hashes:
                    ref_bs.AddResidue(r, mol.ViewAddFlag.INCLUDE_ALL)

        # create view of model binding site using residue mapping from qs_scorer
        # build up ref to mdl alignments for each chain on the go (alns)
        mdl_bs = self.mdl.CreateEmptyView()
        alns = dict()
        rmapping = self.qs_scorer.mapped_residues
        chain_mapping = self.qs_scorer.chain_mapping
        for ref_chain in ref_bs.chains:
            ref_cname = ref_chain.GetName()
            ref_olcs = list()
            mdl_olcs = list()
            for ref_r in ref_chain.residues:
                ref_rnum = ref_r.GetNumber().GetNum()
                ref_olcs.append(ref_r.one_letter_code)
                mdl_res_found = False
                if ref_cname in rmapping and ref_rnum in rmapping[ref_cname]:
                    mdl_cname = chain_mapping[ref_cname]
                    mdl_rnum = rmapping[ref_cname][ref_rnum]
                    mdl_r = self.mdl.FindResidue(mdl_cname, mol.ResNum(mdl_rnum))
                    if mdl_r.IsValid():
                        mdl_res_found = True
                        mdl_bs.AddResidue(mdl_r, mol.ViewAddFlag.INCLUDE_ALL)
                        mdl_olcs.append(mdl_r.one_letter_code)
                if not mdl_res_found:
                    mdl_olcs.append('-')
            if list(set(mdl_olcs)) != ['-']:
                mdl_cname = chain_mapping[ref_cname]
                a = seq.CreateAlignment()
                a.AddSequence(seq.CreateSequence(ref_cname, ''.join(ref_olcs)))
                a.AddSequence(seq.CreateSequence(mdl_cname, ''.join(mdl_olcs)))
                alns[mdl_cname] = a

        scorer = lddt.lDDTScorer(ref_bs, inclusion_radius = lddt_radius)
        # lddt wants model chains mapped on target chain => invert
        inv_chain_mapping = {v:k for k,v in chain_mapping.items()}
        # additionally, lddt only wants chains in that mapping that
        # actually exist in the provided structures
        lddt_chain_mapping = {k: inv_chain_mapping[k] for k in alns.keys()}
        score, _ = scorer.lDDT(mdl_bs, chain_mapping = lddt_chain_mapping,
                               residue_mapping = alns)

        if return_mapping:
            trg_residues = [str(r) for r in ref_bs.residues]
            mdl_residues = [str(r) for r in mdl_bs.residues]
            return (score, trg_residues, mdl_residues)
        else:
            return score
